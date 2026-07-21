import tools.bio
import argparse
import collections
import json
import logging
import os
import re
import sys
import time
import warnings
import pandas as pd
import multiprocessing as mp
from functools import partial
from pathlib import Path
import tools.fileOps
import tools.intervals
import tools.mathOps
import tools.misc
import tools.nameConversions
import tools.procOps
import tools.sqlInterface
import tools.transcripts
from tools.defaultOrderedDict import DefaultOrderedDict
from sqlalchemy import inspect

# Suppress warnings for better performance
warnings.filterwarnings('ignore')
warnings.filterwarnings('ignore', category=FutureWarning, module='pandas')
warnings.filterwarnings('ignore', category=DeprecationWarning, module='pandas')

DENOVO_PREFIXES = ('augPB-', 'strg-')

def _is_denovo(aln_id):
    """Check if an alignment ID belongs to a de novo mode (augPB or strg)."""
    return isinstance(aln_id, str) and aln_id.startswith(DENOVO_PREFIXES)

logger = logging.getLogger(__name__)
ID_TEMPLATE = '{genome:.10}_{tag_type}{unique_id:07d}'


def _import_consensus_parallel():
    """Import shared consensus logic from this repo (not ~/.local site-packages cat)."""
    import importlib

    repo_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    if repo_root in sys.path:
        sys.path.remove(repo_root)
    sys.path.insert(0, repo_root)
    for mod in list(sys.modules):
        if mod == 'cat' or mod.startswith('cat.'):
            del sys.modules[mod]

    cp = importlib.import_module('cat.consensus_parallel')
    cp_file = getattr(cp, '__file__', '')
    if repo_root not in os.path.realpath(cp_file or ''):
        raise ImportError(
            f"Expected consensus_parallel from {repo_root}, got {cp_file}. "
            "Remove or upgrade the pip-installed `cat` package if it shadows this repo."
        )

    names = [
        'process_chromosome',
        'track_multi_locus_mappings',
        'select_consensus_with_cnv',
        'create_transcript_attributes',
        'classify_augpb_transcript',
        'check_novel_splices',
        'score_transcripts',
        'initialize_metrics',
        'load_transmap_evals',
        'load_metrics_from_db',
        'load_evaluations_from_db',
        'deduplicate_consensus',
        'finalize_consensus_after_source_gene_resolution',
        'resolve_overlapping_cds_intervals',
        'calculate_completeness',
        'write_consensus_gps',
        'write_consensus_gff3',
        'write_consensus_fastas',
        'apply_reference_gene_biotype_policy',
        'rescue_missing_reference_pc_genes',
        'rescue_missing_reference_noncoding_genes',
        'rescue_missing_reference_transcripts',
        'rescue_alternative_source_isoforms',
        'build_alignment_coverage_map',
        'norm_ensg',
    ]
    out = {n: getattr(cp, n) for n in names}
    missing = [n for n, fn in out.items() if fn is None]
    if missing:
        raise ImportError(
            f"cat.consensus_parallel from {cp_file} is missing: {missing}. "
            "Use the cat2 repo checkout, not an older pip install."
        )
    return out


_cp = _import_consensus_parallel()
process_chromosome = _cp['process_chromosome']
track_multi_locus_mappings = _cp['track_multi_locus_mappings']
select_consensus_with_cnv = _cp['select_consensus_with_cnv']
create_transcript_attributes = _cp['create_transcript_attributes']
classify_augpb_transcript = _cp['classify_augpb_transcript']
check_novel_splices = _cp['check_novel_splices']
score_transcripts = _cp['score_transcripts']
initialize_metrics = _cp['initialize_metrics']
load_transmap_evals = _cp['load_transmap_evals']
load_metrics_from_db = _cp['load_metrics_from_db']
load_evaluations_from_db = _cp['load_evaluations_from_db']
deduplicate_consensus = _cp['deduplicate_consensus']
finalize_consensus_after_source_gene_resolution = _cp['finalize_consensus_after_source_gene_resolution']
resolve_overlapping_cds_intervals = _cp['resolve_overlapping_cds_intervals']
calculate_completeness = _cp['calculate_completeness']
write_consensus_gps = _cp['write_consensus_gps']
write_consensus_gff3 = _cp['write_consensus_gff3']
write_consensus_fastas = _cp['write_consensus_fastas']
apply_reference_gene_biotype_policy = _cp['apply_reference_gene_biotype_policy']
rescue_missing_reference_pc_genes = _cp['rescue_missing_reference_pc_genes']
rescue_missing_reference_noncoding_genes = _cp['rescue_missing_reference_noncoding_genes']
rescue_missing_reference_transcripts = _cp['rescue_missing_reference_transcripts']
rescue_alternative_source_isoforms = _cp['rescue_alternative_source_isoforms']
norm_ensg = _cp['norm_ensg']


def main():
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    
    logger.info("="*80)
    logger.info(f"Starting consensus generation for {args.genome} (LOCAL VERSION)")
    logger.info("="*80)
    
    # Run the main consensus logic
    metrics = generate_consensus(args)

    # Write the final metrics file (exclude internal keys)
    logger.info(f"Writing metrics to {args.metrics_json}")
    output_metrics = {k: v for k, v in metrics.items() if not k.startswith('_')}
    with open(args.metrics_json, 'w') as outf:
        json.dump(output_metrics, outf, indent=4)
    logger.info(f"Successfully generated consensus gene set and metrics for {args.genome}.")
    logger.info("="*80)


def add_arguments(parser):
    """Adds all command-line arguments to the argparse parser."""
    # --- Input Files ---
    parser.add_argument("--gp-list", nargs='+', required=True, help="Space-separated list of all genePred files to consider.")
    parser.add_argument("--db-path", required=True, help="Path to the genome's primary database.")
    parser.add_argument("--ref-db-path", required=True, help="Path to the reference genome's database.")
    parser.add_argument("--fasta", required=True, help="Path to the target genome FASTA file.")
    parser.add_argument("--genome", required=True, help="Name of the target genome.")
    
    # --- Output Files ---
    parser.add_argument("--consensus-gp", required=True, help="Output path for consensus genePred.")
    parser.add_argument("--consensus-gp-info", required=True, help="Output path for consensus gp_info TSV.")
    parser.add_argument("--consensus-gff3", required=True, help="Output path for consensus GFF3.")
    parser.add_argument("--consensus-fasta", required=True, help="Output path for consensus transcript FASTA.")
    parser.add_argument("--protein-fasta", required=True, help="Output path for consensus protein FASTA.")
    parser.add_argument("--metrics-json", required=True, help="Output path for consensus metrics JSON.")
    
    # --- Control and Filtering Parameters ---
    parser.add_argument("--intron-rnaseq-support", type=float, default=0.0, help="Percent of introns that must be supported by RNA-seq.")
    parser.add_argument("--exon-rnaseq-support", type=float, default=0.0, help="Percent of exons supported by RNA-seq.")
    parser.add_argument("--intron-annot-support", type=float, default=0.0, help="Percent of introns supported by reference annotation.")
    parser.add_argument("--exon-annot-support", type=float, default=0.0, help="Percent of exons supported by reference annotation.")
    parser.add_argument("--original-intron-support", type=float, default=0.0, help="Percent of original introns that must be preserved.")
    parser.add_argument("--in-species-rna-support-only", action="store_true", help="Use in-species RNA-seq support only, ignoring cross-species evidence.")
    parser.add_argument("--filter-overlapping-genes", action="store_true", help="Filter out overlapping CDS intervals from different genes.")
    parser.add_argument("--overlapping-ignore-bases", type=int, default=0, help="Number of bases to ignore when clustering for overlap.")
    parser.add_argument("--cnv-score-similarity", type=float, default=0.80, help="Keep multi-locus transcripts if scores are within this fraction of max (default: 0.80 = 80%%)")
    parser.add_argument("--fragment-max-coverage", type=float, default=30.0, help="A non-denovo transcript with coverage AND identity below their respective thresholds is reclassified as 'fragment'. Coverage threshold (percent); set to 0 to disable fragment reclassification (keeps them as regular orthologs).")
    parser.add_argument("--fragment-max-identity", type=float, default=30.0, help="Identity (percent) companion to --fragment-max-coverage; set to 0 to disable fragment reclassification.")
    parser.add_argument(
        "--keep-protein-only-novel",
        action="store_true",
        help="Retain (and cleanly label) lineage-specific genes found ONLY in protein/augMP "
             "evidence, i.e. augMP models whose protein is not a reference transcript. When set, "
             "such orphan models are relabeled 'putative_novel' protein_coding with MP-NOVEL-<genome>-N "
             "gene ids, redundant per-locus copies are collapsed, and (unless "
             "--protein-novel-keep-overlapping) only intergenic loci are kept. Off by default; "
             "enabled automatically under high_recall.",
    )
    parser.add_argument(
        "--protein-novel-min-coverage",
        type=float,
        default=0.0,
        help="Minimum miniprot protein-space coverage (percent) for a protein-only novel gene "
             "to be kept (0 = no extra filter beyond miniprot's own).",
    )
    parser.add_argument(
        "--protein-novel-min-identity",
        type=float,
        default=0.0,
        help="Minimum miniprot protein-space identity (percent) for a protein-only novel gene "
             "to be kept (0 = no extra filter).",
    )
    parser.add_argument(
        "--protein-novel-keep-overlapping",
        action="store_true",
        help="Also keep protein-only novel genes whose exons overlap an existing reference-projected "
             "feature (any biotype, same strand). Default: keep only intergenic loci, which is safer "
             "for precision (drops projected pseudogenes/retrocopies mistaken for novel genes).",
    )
    parser.add_argument(
        "--rescue-expressed-noncoding-to-pc",
        action="store_true",
        help="Promote an expressed, ORF-intact model to protein_coding when it only carries a "
             "non-coding gene_biotype because it was inherited from a non-coding reference ortholog "
             "(lncRNA/pseudogene). Requires a valid ORF, >= --rescue-expressed-min-cds-aa residues, "
             "multi-exon (unless --rescue-expressed-allow-single-exon), and IsoSeq/RNA support. "
             "Recovers real PC genes without adding intergenic false positives.",
    )
    parser.add_argument(
        "--rescue-expressed-min-cds-aa",
        type=int,
        default=100,
        help="Minimum CDS length (amino acids) for the expressed-noncoding->PC rescue [100].",
    )
    parser.add_argument(
        "--rescue-expressed-allow-single-exon",
        action="store_true",
        help="Allow single-exon models to be rescued to protein_coding (default: require multi-exon "
             "to exclude processed pseudogenes/retrocopies).",
    )
    parser.add_argument(
        "--rescue-noncoding-require-protein-evidence",
        dest="rescue_noncoding_require_protein_evidence",
        action="store_true",
        default=True,
        help="When rescuing expressed models whose reference biotype is confidently non-coding "
             "(lncRNA / pseudogene), additionally require independent protein evidence: a same-strand "
             "augMP (miniprot-hinted Augustus) CDS overlapping the model's CDS. This keeps ~85-95%% of "
             "genuine PC recoveries (RefSeq-confirmed) while dropping ~90%% of expressed lncRNAs that "
             "merely carry an incidental ORF. Denovo/unknown and IG/TR biotypes are never gated. "
             "Default: on.",
    )
    parser.add_argument(
        "--no-rescue-noncoding-require-protein-evidence",
        dest="rescue_noncoding_require_protein_evidence",
        action="store_false",
        help="Disable the protein-evidence gate for lncRNA/pseudogene rescue (more permissive, "
             "recovers ~40 more real genes per genome but adds ~800 incidental-ORF lncRNAs).",
    )
    parser.add_argument(
        "--rescue-dropped-augMP",
        dest="rescue_dropped_augMP",
        action="store_true",
        help="Recover high-quality augMP models that consensus selection dropped, at loci where "
             "the consensus is otherwise EMPTY on that strand. augMP is frequently the only mode "
             "that finds lineage-specific genes (it uses species-specific proteins), so this closes "
             "a real recall gap. Requires multi-exon + ORF >= --rescue-augMP-min-cds-aa. Off by default.",
    )
    parser.add_argument(
        "--rescue-augMP-min-exons",
        type=int,
        default=2,
        help="Minimum exon count for augMP empty-locus recovery [2].",
    )
    parser.add_argument(
        "--rescue-augMP-min-cds-aa",
        type=int,
        default=100,
        help="Minimum CDS length (aa) for augMP empty-locus recovery [100].",
    )
    parser.add_argument(
        "--rescue-augMP-single-exon-min-cds-aa",
        type=int,
        default=300,
        help="Single-exon augMP/augPB models are recovered only if their ORF is at least this "
             "many residues (keeps intronless genes like SMEK1, excludes short retrocopies) [300].",
    )
    parser.add_argument(
        "--rescue-augMP-min-coverage",
        type=float,
        default=0.0,
        help="Minimum miniprot coverage (%%) for augMP empty-locus recovery, if metrics exist [0].",
    )
    parser.add_argument(
        "--rescue-augMP-min-identity",
        type=float,
        default=0.0,
        help="Minimum miniprot identity (%%) for augMP empty-locus recovery, if metrics exist [0].",
    )
    parser.add_argument(
        "--protein-novel-min-exons",
        type=int,
        default=2,
        help="Minimum exon count for a protein-only novel gene to be kept. Default 2 removes "
             "single-exon models (dominated by retrocopies/spurious hits). Set to 1 to keep "
             "single-exon novel genes (high_recall does this).",
    )
    parser.add_argument(
        "--protein-novel-min-cds-aa",
        type=int,
        default=100,
        help="Minimum CDS length (amino acids) for a protein-only novel gene to be kept. "
             "Default 100 removes very short ORFs unlikely to be real genes. Set to 0 to disable "
             "(high_recall does this).",
    )
    parser.add_argument(
        "--txTM-min-coverage",
        type=float,
        default=0.0,
        help="Drop txTM transcripts with AlnCoverage_mRNA <= this value (percent). "
             "NaN coverage is kept unless --txTM-strict-metrics is set.",
    )
    parser.add_argument(
        "--txTM-strict-metrics",
        action="store_true",
        help="Require txTM transcripts to have AlnCoverage_mRNA in the metrics DB "
             "(do not treat missing coverage as 100%%).",
    )
    parser.add_argument(
        "--txTM-transmap-anchor-overlap",
        type=float,
        default=0.0,
        help="When transMap exists for a transcript, drop txTM copies at other loci "
             "whose exons overlap transMap by less than this fraction (0 disables). "
             "Paralog loci without transMap are unaffected.",
    )
    parser.add_argument(
        "--txTM-min-coverage-no-transmap",
        type=float,
        default=None,
        help="Min AlnCoverage (%%) for txTM when this transcript accession has NO transMap "
             "anchor (default: 80). Uses >= threshold.",
    )
    parser.add_argument(
        "--txTM-strict-metrics-no-transmap",
        action="store_true",
        default=None,
        help="Require metrics for txTM-only transcripts; drop NaN coverage (default: true).",
    )
    parser.add_argument(
        "--no-txTM-strict-metrics-no-transmap",
        dest="txTM_strict_metrics_no_transmap",
        action="store_false",
        help="Explicitly allow missing txTM metrics for txTM-only transcripts.",
    )
    parser.add_argument(
        "--min-pc-len-ratio-txTM-only-rescue",
        type=float,
        default=None,
        help="Min span ratio vs reference for rescuing PC genes with only txTM candidates "
             "(default: same as --min-pc-len-ratio-vs-reference).",
    )
    parser.add_argument(
        "--augMP-min-coverage-no-anchor",
        type=float,
        default=80.0,
        help="Min AlnCoverage (%%) for augMP when no transMap anchor exists for that "
             "transcript (default: 80). Uses >= threshold; set 0 to disable.",
    )
    parser.add_argument(
        "--augMP-strict-metrics-no-anchor",
        action="store_true",
        default=None,
        help="Require metrics for augMP-only transcripts; drop NaN coverage (default: true).",
    )
    parser.add_argument(
        "--no-augMP-strict-metrics-no-anchor",
        dest="augMP_strict_metrics_no_anchor",
        action="store_false",
        help="Allow missing augMP coverage metrics for augMP-only transcripts.",
    )
    parser.add_argument(
        "--rescue-reference-isoforms",
        action="store_true",
        default=True,
        help="After consensus selection, re-add missing reference PC isoforms from "
             "transMap/txTM when they pass rescue coverage rules (default: true).",
    )
    parser.add_argument(
        "--no-rescue-reference-isoforms",
        dest="rescue_reference_isoforms",
        action="store_false",
        help="Disable reference isoform rescue pass.",
    )
    parser.add_argument(
        "--rescue-min-txTM-coverage",
        type=float,
        default=80.0,
        help="Min txTM AlnCoverage for rescued reference isoforms only (default: 80). "
             "transMap rescues are not coverage-filtered here.",
    )
    parser.add_argument(
        "--rescue-reference-noncoding-genes",
        action="store_true",
        default=True,
        help="After consensus, re-add missing ref pseudogene/lncRNA genes from transMap/txTM "
             "(default: true).",
    )
    parser.add_argument(
        "--no-rescue-reference-noncoding-genes",
        dest="rescue_reference_noncoding_genes",
        action="store_false",
        help="Disable reference non-coding gene rescue pass.",
    )
    parser.add_argument(
        "--rescue-min-txTM-coverage-noncoding",
        type=float,
        default=50.0,
        help="Min txTM AlnCoverage for rescued ref pseudogene/lncRNA genes (default: 50).",
    )
    parser.add_argument(
        "--min-nc-len-ratio-txTM-only-rescue",
        type=float,
        default=0.25,
        help="Min span ratio vs reference for rescuing non-coding genes with only txTM "
             "candidates (default: 0.25).",
    )
    parser.add_argument(
        "--min-nc-len-ratio-vs-reference",
        type=float,
        default=0.0,
        help="Min span ratio vs reference for non-coding rescue when transMap exists "
             "(default: 0 = disabled).",
    )
    parser.add_argument(
        "--txTM-min-coverage-noncoding",
        type=float,
        default=None,
        help="txTM min coverage for non-coding pool only (default: same as --txTM-min-coverage).",
    )
    parser.add_argument(
        "--txTM-min-coverage-no-transmap-noncoding",
        type=float,
        default=None,
        help="txTM-only min coverage for non-coding pool (default: same as coding orphan threshold).",
    )
    parser.add_argument(
        "--txTM-strict-metrics-no-transmap-noncoding",
        action="store_true",
        default=None,
        help="Require metrics for txTM-only non-coding transcripts (default: false when unset).",
    )
    parser.add_argument(
        "--no-txTM-strict-metrics-no-transmap-noncoding",
        dest="txTM_strict_metrics_no_transmap_noncoding",
        action="store_false",
        help="Allow missing txTM metrics for txTM-only non-coding transcripts.",
    )
    parser.add_argument(
        "--rescue-alternative-isoforms",
        action="store_true",
        default=True,
        help="Promote reference isoforms listed only in alternative_source_transcripts "
             "to full transcript records (default: true).",
    )
    parser.add_argument(
        "--no-rescue-alternative-isoforms",
        dest="rescue_alternative_isoforms",
        action="store_false",
        help="Disable alternative-source isoform promotion.",
    )

    # --- Parallelization Parameters (LOCAL ONLY) ---
    parser.add_argument("--num-workers", type=int, default=None, help="Number of parallel workers (default: auto-detect, max 23)")
    
    # --- De Novo Parameters ---
    parser.add_argument("--denovo-tx-modes", nargs='*', default=[], help="List of de novo modes to consider (e.g., augPB).")
    parser.add_argument("--denovo-num-introns", type=int, default=1, help="De novo isoforms (augPB and strg) must have at least this many introns (default 1 → ≥2 exons).")
    parser.add_argument("--denovo-splice-support", type=float, default=1.0, help="Percent of de novo splices that must be RNA-seq supported.")
    parser.add_argument("--denovo-exon-support", type=float, default=1.0, help="Percent of de novo exons that must be RNA-seq supported.")
    parser.add_argument("--denovo-ignore-novel-genes", action="store_true", help="If set, only incorporate de novo transcripts as novel isoforms, not novel genes.")
    parser.add_argument("--denovo-only-novel-genes", action="store_true", help="If set, only incorporate de novo transcripts if they are novel genes.")
    parser.add_argument("--denovo-allow-unsupported", action="store_true", help="Allow de novo transcripts with novel splices that lack RNA-seq support.")
    parser.add_argument("--denovo-allow-bad-annot-or-tm", action="store_true", help="Allow de novo models flagged as 'badAnnotOrTm'.")
    parser.add_argument("--denovo-allow-novel-ends", action="store_true", help="Allow de novo models with novel 5' or 3' ends.")
    parser.add_argument("--denovo-novel-end-distance", type=int, default=0)
    parser.add_argument("--strg-min-single-exon-len", type=int, default=500, help="Deprecated; strg now uses --denovo-num-introns like augPB.")

    # --- PacBio Parameters ---
    parser.add_argument("--require-pacbio-support", action="store_true", help="If set, remove any consensus transcript not validated by Iso-Seq data.")
    parser.add_argument("--hints-db-has-rnaseq", action="store_true", help="Flag if the hints DB contains RNA-seq, for tagging purposes.")
    
    # --- RNA-seq / Annotation Support ---
    parser.add_argument("--bam-files", nargs='*', default=[], help="BAM files for computing real RNA-seq splice/exon support.")
    parser.add_argument("--isoseq-bam-files", nargs='*', default=[], help="IsoSeq BAM files for computing splice/exon support.")
    parser.add_argument("--ref-gp", default=None, help="Reference annotation genePred for computing annotation support.")

    # --- Postprocess (split runaway genes + drop weak paralog copies) ---
    parser.add_argument(
        "--no-consensus-postprocess",
        dest="consensus_postprocess",
        action="store_false",
        default=True,
        help=(
            "Disable the post-consensus cleanup that (a) splits pc gene records whose "
            "transcripts have ballooned to span multiple distinct loci into per-cluster "
            "records, and (b) reclassifies single-exon copies of multi-exon reference "
            "genes as processed_pseudogene and drops weak duplicate copies when a "
            "stronger ortholog exists at another locus. Default: ON (requires --ref-gp "
            "and an adjacent <ref_gp>_attrs file)."
        ),
    )
    parser.add_argument("--postprocess-min-introns-low-support", type=int, default=3,
                        help="Min intron count before the low-intron-support drop rule can fire (default 3).")
    parser.add_argument("--postprocess-low-support-fraction", type=float, default=0.3,
                        help="Supported/total intron ratio below which a low-support duplicate may be dropped (default 0.3). Lower to drop fewer.")
    parser.add_argument("--postprocess-augpb-chimera-exon-ratio", type=float, default=1.5,
                        help="augPB duplicate is treated as a chimera (droppable) when its exon count exceeds the reference by this ratio (default 1.5).")
    parser.add_argument("--postprocess-allow-drop-strong-modes", dest="postprocess_protect_strong_modes",
                        action="store_false", default=True,
                        help="Allow the low-intron-support drop rule to also drop loci backed by a strong mode (transMap/augTM/augPB/strg). Default: protect them (do not drop).")

    # --- Diagnostics / auditing ---
    parser.add_argument(
        "--pc-audit-tsv",
        type=Path,
        default=None,
        help=(
            "Optional: write a TSV auditing protein-coding gene inclusion (per mode: input vs final consensus). "
            "Useful for spotting which PC genes were provided to consensus but not retained."
        ),
    )
    parser.add_argument(
        "--pc-audit-chrom",
        default=None,
        help=(
            "Optional: restrict the protein-coding audit to a single chromosome/contig (e.g. chr20). "
            "This does not change the consensus results; it only filters the audit report."
        ),
    )
    parser.add_argument(
        "--only-chrom",
        default=None,
        help=(
            "Optional: run consensus computation on only one chromosome/contig (e.g. chr20). "
            "This is intended for quick testing/debugging and will produce a partial consensus gene set."
        ),
    )
    parser.add_argument(
        "--disregard-long-mode-ratio",
        type=float,
        default=2.0,
        help=(
            "When computing gene spans for overlap/conflict resolution, ignore any alignment_mode whose "
            "span for a gene is >= ratio * median(other modes' spans) AND exceeds --disregard-long-mode-min-bp "
            "in absolute extra length. Default: 2.0."
        ),
    )
    parser.add_argument(
        "--disregard-long-mode-min-bp",
        type=int,
        default=50_000,
        help=(
            "Minimum extra length (bp) beyond the median span required to disregard an outlier mode for a gene "
            "during overlap/conflict resolution. Default: 50000."
        ),
    )
    parser.add_argument(
        "--spurious-overlap-min-assembly-overlap-bp",
        type=int,
        default=2000,
        help=(
            "For the 'spurious 2-way overlap (not in reference)' filter, require at least this many bp of "
            "actual overlap between the two gene intervals in the target assembly before removing one. "
            "Smaller overlaps are treated as boundary noise and both genes are kept. Default: 2000."
        ),
    )
    parser.add_argument(
        "--spurious-overlap-min-reciprocal",
        type=float,
        default=0.02,
        help=(
            "For the 'spurious 2-way overlap (not in reference)' filter, require at least this reciprocal "
            "overlap fraction (min of overlap/gene_len) before removing one. Default: 0.02 (2%%)."
        ),
    )
    parser.add_argument(
        "--min-pc-len-ratio-vs-reference",
        type=float,
        default=0.4,
        help=(
            "Remove transcripts when their span is below this fraction of the reference source_gene span. "
            "Applies to all gene biotypes. Set to 0 to disable. Default: 0.4."
        ),
    )
    parser.add_argument(
        "--filter-spurious-pc-overlaps-not-in-reference",
        action="store_true",
        help=(
            "If set, remove overlapping protein-coding gene loci when the overlap is not supported by the "
            "reference gene coordinates (including different reference chromosomes). This is a strict cleanup "
            "pass to reduce spurious loci."
        ),
    )

def map_mode_to_db_table(mode):
    """
    Map display mode names to their database table names.
    Each mode now has its own table, so just return the mode as-is.
    """
    return mode


def normalize_alignment_id(aln_id, mode):
    """
    Normalize alignment IDs according to mode:
    - transMap/transMap_pairwise: keep as-is
    - txTM: strip underscore and numbers after (e.g., ENST00000123_1 -> ENST00000123)
    - augTM/augTMR/augMP (and pairwise variants): strip prefix (e.g., augTM-ENST00000123 -> ENST00000123)
    """
    if mode in ['transMap', 'transMap_pairwise']:
        return aln_id
    elif mode == 'txTM':
        # Strip cross-mode _cp suffix first, then txTM CNV copy _N (e.g. XM_….1_14_cp9 → XM_….1).
        import re
        base = re.sub(r'_cp\d*$', '', aln_id)
        base = re.sub(r'_\d+$', '', base)
        return base
    elif mode in ['augTM', 'augTM_pairwise', 'augTMR', 'augTMR_pairwise', 'augMP']:
        # Strip the prefix
        if aln_id.startswith('augTM-'):
            return aln_id[6:]  # len('augTM-') = 6
        elif aln_id.startswith('augTMR-'):
            return aln_id[7:]  # len('augTMR-') = 7
        elif aln_id.startswith('augMP-'):
            return aln_id[6:]  # len('augMP-') = 6
        return aln_id
    elif mode == 'augPB':
        # Keep augPB IDs as-is for now
        return aln_id
    return aln_id


def normalize_gene_id(gene_id, mode):
    """
    Normalize gene IDs according to mode (same as transcript ID normalization).
    This is important for txTM CNV copies which append _N to both transcript and gene IDs.
    
    - transMap/transMap_pairwise: keep as-is
    - txTM: strip ONLY if pattern is _N_M (two consecutive number suffixes)
    - augTM/augTMR/augMP: keep as-is
    - augPB: keep as-is
    
    Examples for txTM:
    - hg002_chrY_paternal_691_1 → hg002_chrY_paternal_691 (CNV copy, strip _1)
    - hg002_chrY_paternal_691 → hg002_chrY_paternal_691 (original, keep as-is)
    - hg002_chrY_paternal_157_10 → hg002_chrY_paternal_157 (CNV copy 10, strip _10)
    - hg002_chrY_paternal_166_cp2 → hg002_chrY_paternal_166_cp2 (paralog, PRESERVE _cp2)
    """
    if mode == 'txTM':
        # Only strip the LAST _N if there are TWO consecutive _N patterns
        # Pattern: ends with _digits_digits (e.g., _691_1 or _157_10)
        import re
        match = re.search(r'_\d+_(\d+)$', gene_id)
        if match:
            # Has double _N_M pattern, strip the last _M
            base = re.sub(r'_\d+$', '', gene_id)
            return base
    # For all other modes or single _N pattern, keep gene ID as-is
    return gene_id


def identify_mode(aln_id, gp_file=None):
    """Identify the alignment mode from the alignment ID or file path"""
    if aln_id.startswith('augPB-'):
        return 'augPB'
    elif aln_id.startswith('strg-'):
        return 'strg'
    elif aln_id.startswith('augTMR-'):
        # Check if it's from pairwise file
        if gp_file and 'pairwise' in gp_file:
            return 'augTMR_pairwise'
        return 'augTMR'
    elif aln_id.startswith('augTM-'):
        # Check if it's from pairwise file
        if gp_file and 'pairwise' in gp_file:
            return 'augTM_pairwise'
        return 'augTM'
    elif gp_file and 'txTM' in gp_file:
        return 'txTM'
    elif gp_file and 'transMap_pairwise' in gp_file:
        return 'transMap_pairwise'
    elif gp_file and 'transMap' in gp_file:
        return 'transMap'
    elif gp_file and 'augMP' in gp_file:
        return 'augMP'
    # Default fallback
    return 'transMap'


def _is_txtm_cnv_copy_id(aln_id):
    """True when aln_id ends with a txTM CNV copy suffix (_N), not a version dot."""
    return bool(re.search(r'_\d+$', str(aln_id)))


def _strip_cross_mode_cp_suffix(aln_id):
    """Strip consensus _cp/_cp2 suffix from cross-mode duplicate renaming."""
    return re.sub(r'_cp\d*$', '', str(aln_id))


def _build_txtm_metrics_match_keys(valid_aln_ids, alignment_source_map):
    """
    AlignmentIds usable to find txTM rows in the metrics DB.

    tx_dict may use cross-mode _cp suffixes (XM_….1_cp9) while the DB stores XM_….1.
    """
    keys = set()
    for aid in valid_aln_ids:
        if alignment_source_map.get(aid) != 'txTM':
            continue
        keys.add(str(aid))
        base_no_cp = _strip_cross_mode_cp_suffix(aid)
        if base_no_cp != aid:
            keys.add(base_no_cp)
    return keys


def _txtm_db_id_matches(db_id, match_keys, valid_aln_ids):
    """True if a txTM metrics DB AlignmentId corresponds to a tx_dict entry."""
    db_id = str(db_id)
    if db_id in match_keys:
        return True
    for n in range(1, 100):
        if f"{db_id}_{n}" in match_keys:
            return True
    return False


def _txtm_allowed_db_ids(match_keys):
    """Vectorized-filter equivalent of ``_txtm_db_id_matches`` over a whole column.

    A DB AlignmentId matches iff it is in ``match_keys`` OR one of ``<id>_1`` ..
    ``<id>_99`` is in ``match_keys`` (a txTM CNV copy). Rather than test that per
    row with an up-to-99-iteration Python loop (~200M string ops over millions of
    metric rows), precompute the full set of acceptable DB ids once: every
    ``match_keys`` entry, plus the ``<prefix>`` of every ``<prefix>_<1..99>`` key.
    Callers then use a single vectorized ``Series.isin(...)``.
    """
    allowed = set(match_keys)
    for k in match_keys:
        idx = k.rfind("_")
        if idx <= 0 or idx == len(k) - 1:
            continue
        suffix = k[idx + 1:]
        # Only numeric CNV suffixes 1..99 count (matches the old range(1, 100)).
        if suffix.isdigit() and 1 <= int(suffix) <= 99:
            allowed.add(k[:idx])
    return allowed


def remap_metrics_to_txtm_cp_aliases(metrics_df, valid_aln_ids, alignment_source_map):
    """
    Attach metrics to txTM AlignmentIds renamed with _cp for cross-mode conflicts.

    GenePred loading renames e.g. XM_….1 → XM_….1_cp9 when another mode used XM_….1 first.
    Metrics remain keyed by the DB id; duplicate rows under the tx_dict id before merge.
    """
    if metrics_df is None or len(metrics_df) == 0 or 'Mode' not in metrics_df.columns:
        return metrics_df
    if not (metrics_df['Mode'] == 'txTM').any():
        return metrics_df

    existing = set(metrics_df['AlignmentId'].astype(str))
    alias_targets = {}
    for aid in valid_aln_ids:
        if alignment_source_map.get(aid) != 'txTM':
            continue
        base_no_cp = _strip_cross_mode_cp_suffix(aid)
        if base_no_cp == aid:
            continue
        alias_targets.setdefault(base_no_cp, []).append(str(aid))

    # Index the base rows we need ONCE instead of rescanning metrics_df per db_id
    # (the old `metrics_df[metrics_df['AlignmentId'] == db_id]` inside the loop was
    # O(n_alias_ids * n_rows); with ~500k cross-mode _cp duplicates that stalls just
    # like the CNV expansion did). Filter to the needed db_ids first, then bucket.
    new_rows = []
    needed_ids = set(alias_targets.keys())
    rows_by_id = {}
    if needed_ids:
        sub = metrics_df[metrics_df['AlignmentId'].astype(str).isin(needed_ids)]
        for rec in sub.to_dict('records'):
            rows_by_id.setdefault(str(rec['AlignmentId']), []).append(rec)
    for db_id, tx_ids in alias_targets.items():
        base_recs = rows_by_id.get(db_id)
        if not base_recs:
            continue
        for tx_id in tx_ids:
            if tx_id in existing:
                continue
            for rec in base_recs:
                r = dict(rec)
                r['AlignmentId'] = tx_id
                r['Mode'] = 'txTM'
                new_rows.append(r)
            existing.add(tx_id)

    if not new_rows:
        return metrics_df

    logger.info(
        f"    Remapped txTM metrics onto {len(new_rows)} cross-mode _cp alias rows "
        f"({len(alias_targets)} DB ids with aliases)"
    )
    out = pd.concat([metrics_df, pd.DataFrame(new_rows)], ignore_index=True)
    return out.drop_duplicates(subset=['AlignmentId'], keep='first')


def _txtm_cnv_copy_ids(base_id, valid_aln_ids, max_copy=99):
    """Numbered CNV copies in genePred, including cross-mode _cp-renamed IDs."""
    base_id = str(base_id)
    root = _strip_cross_mode_cp_suffix(base_id)
    copies = []
    seen = set()
    for n in range(1, max_copy + 1):
        for candidate in (f"{root}_{n}", f"{base_id}_{n}"):
            if candidate in valid_aln_ids and candidate not in seen:
                copies.append(candidate)
                seen.add(candidate)
    for aid in valid_aln_ids:
        aid = str(aid)
        if aid in seen or aid == base_id:
            continue
        no_cp = _strip_cross_mode_cp_suffix(aid)
        if no_cp == root:
            continue
        if no_cp.startswith(f"{root}_") and re.search(r'_\d+$', no_cp):
            if aid not in seen:
                copies.append(aid)
                seen.add(aid)
    return copies


def expand_txtm_cnv_metrics(metrics_df, valid_aln_ids):
    """
    Propagate txTM base alignment metrics onto numbered CNV copy IDs in genePred.

    On Y (and similar loci), genePred often has both a bare accession (XM_… .1) and
    numbered copies (XM_….1_14). The DB stores metrics on every ID. Older logic
    deleted bare base rows whenever any _N copy existed, which left bare genePred
    IDs without metrics (NaN coverage after merge).

    Rules:
    - Keep a base metrics row when that base AlignmentId is in valid_aln_ids.
    - Remove a base row only when it is not in genePred but numbered copies are
      (metrics live only under the base accession in the DB).
    - Add propagated rows only for copy IDs that lack metrics rows already.
    """
    if metrics_df is None or len(metrics_df) == 0 or 'Mode' not in metrics_df.columns:
        return metrics_df
    if not (metrics_df['Mode'] == 'txTM').any():
        return metrics_df

    non_txtm = metrics_df[metrics_df['Mode'] != 'txTM']
    txtm = metrics_df[metrics_df['Mode'] == 'txTM'].copy()
    existing_ids = set(txtm['AlignmentId'].astype(str))

    # Precompute genePred CNV copies per txTM base id in a SINGLE pass over
    # valid_aln_ids. The old code called _txtm_cnv_copy_ids() once per base, and
    # each call scanned ALL of valid_aln_ids -> O(n_bases * n_valid), which at
    # panprimate scale (hundreds of thousands of txTM bases x ~2.2M ids) is ~10^11+
    # iterations and effectively never finishes. This builds the same mapping in a
    # single O(n_valid) pass: `c` is a CNV copy of base `b` iff strip_cp(c) starts
    # with strip_cp(b)+"_" and ends in _<digits> (exactly the old Part-B condition,
    # which also subsumes the numbered Part-A candidates). Bucket each qualifying
    # `c` under every underscore-boundary prefix that is a known txTM base root.
    txtm_roots = {_strip_cross_mode_cp_suffix(str(a))
                  for a in txtm['AlignmentId'].unique()
                  if not _is_txtm_cnv_copy_id(a)}
    copies_by_root = collections.defaultdict(list)
    if txtm_roots:
        _cnv_suffix_re = re.compile(r'_\d+$')
        for c in valid_aln_ids:
            c = str(c)
            s = _strip_cross_mode_cp_suffix(c)
            if not _cnv_suffix_re.search(s):
                continue
            idx = s.find('_')
            while idx != -1:
                prefix = s[:idx]
                if prefix and prefix in txtm_roots:
                    copies_by_root[prefix].append(c)
                idx = s.find('_', idx + 1)

    expanded_rows = []
    bases_with_copies = []
    for aid in txtm['AlignmentId'].unique():
        aid = str(aid)
        if _is_txtm_cnv_copy_id(aid):
            continue
        copies = copies_by_root.get(_strip_cross_mode_cp_suffix(aid), [])
        if not copies:
            continue
        bases_with_copies.append(aid)
        base_metrics = txtm[txtm['AlignmentId'] == aid]
        if len(base_metrics) == 0:
            continue
        for cnv_id in copies:
            if cnv_id in existing_ids:
                continue
            for _, row in base_metrics.iterrows():
                cnv_row = row.copy()
                cnv_row['AlignmentId'] = cnv_id
                cnv_row['Mode'] = 'txTM'
                expanded_rows.append(cnv_row)
            existing_ids.add(cnv_id)

    bases_to_drop = [b for b in bases_with_copies if b not in valid_aln_ids]
    if bases_to_drop:
        txtm = txtm[~txtm['AlignmentId'].isin(bases_to_drop)]

    n_added = 0
    if expanded_rows:
        txtm = pd.concat([txtm, pd.DataFrame(expanded_rows)], ignore_index=True)
        n_added = len(expanded_rows)

    if len(txtm) > 0:
        txtm = txtm.drop_duplicates(subset=['AlignmentId'], keep='first')

    if n_added > 0 or bases_to_drop:
        logger.info(
            f"    txTM CNV metrics expansion: {len(bases_with_copies)} bases with copies, "
            f"+{n_added} propagated rows, {len(bases_to_drop)} unused base-only rows removed"
        )

    return pd.concat([non_txtm, txtm], ignore_index=True)


def backfill_cnv_metrics(mrna_metrics_df, cds_metrics_df, valid_aln_ids, alignment_source_map, db_path):
    """
    Backfill metrics for CNV/paralog copies that don't have metrics in the database.
    This handles:
    1. TxTM CNV copies with _N suffix (e.g., ENST00000123.1_1)
    2. Internal paralog copies with _cp suffix (e.g., ENST00000123.1_cp2)
    
    We duplicate the base metrics for all copies.
    """
    if len(mrna_metrics_df) == 0 and len(cds_metrics_df) == 0:
        return mrna_metrics_df, cds_metrics_df
    
    # Find transcripts with suffix that are missing metrics
    existing_mrna_ids = set(mrna_metrics_df['AlignmentId'].values) if len(mrna_metrics_df) > 0 else set()
    existing_cds_ids = set(cds_metrics_df['AlignmentId'].values) if len(cds_metrics_df) > 0 else set()
    
    missing_mrna = []
    missing_cds = []
    
    import re
    
    for aln_id in valid_aln_ids:
        mode = alignment_source_map.get(aln_id, 'unknown')
        base_id = None
        
        # Case 1: TxTM CNV copy (_N suffix)
        if mode == 'txTM' and '_' in aln_id.split('.')[-1]:
            # This is a txTM CNV copy (_N suffix)
            base_id = normalize_alignment_id(aln_id, mode)  # Strips _N
        
        # Case 2: Internal paralog copy (_cp, _cp2, _cp3, etc. suffix)
        elif re.search(r'_cp\d*$', aln_id):
            # Strip _cp suffix to get base ID
            base_id = re.sub(r'_cp\d*$', '', aln_id)
        
        if base_id and base_id != aln_id:
            # Check if this copy is missing metrics but base has them
            if aln_id not in existing_mrna_ids and base_id in existing_mrna_ids:
                missing_mrna.append((aln_id, base_id))
            
            if aln_id not in existing_cds_ids and base_id in existing_cds_ids:
                missing_cds.append((aln_id, base_id))
    
    if missing_mrna or missing_cds:
        logger.info(f"  Backfilling metrics for {len(set([x[0] for x in missing_mrna + missing_cds]))} CNV/paralog copies...")
    
    # OPTIMIZED: Duplicate metrics for missing CNV copies using pandas merge
    if missing_mrna and len(mrna_metrics_df) > 0:
        # Create mapping dataframe: base_id -> cnv_id
        mapping_df = pd.DataFrame(missing_mrna, columns=['cnv_id', 'base_id'])
        
        # Merge with metrics to duplicate rows
        # For each base_id in metrics, create copies for all corresponding cnv_ids
        merged = pd.merge(mapping_df, mrna_metrics_df, left_on='base_id', right_on='AlignmentId', how='inner')
        
        if len(merged) > 0:
            # Replace AlignmentId with cnv_id
            merged['AlignmentId'] = merged['cnv_id']
            merged = merged.drop(['cnv_id', 'base_id'], axis=1)
            
            mrna_metrics_df = pd.concat([mrna_metrics_df, merged], ignore_index=True)
            logger.info(f"    ✓ Added {len(merged)} mRNA metric records for CNV/paralog copies")
    
    if missing_cds and len(cds_metrics_df) > 0:
        # Create mapping dataframe: base_id -> cnv_id
        mapping_df = pd.DataFrame(missing_cds, columns=['cnv_id', 'base_id'])
        
        # Merge with metrics to duplicate rows
        merged = pd.merge(mapping_df, cds_metrics_df, left_on='base_id', right_on='AlignmentId', how='inner')
        
        if len(merged) > 0:
            # Replace AlignmentId with cnv_id
            merged['AlignmentId'] = merged['cnv_id']
            merged = merged.drop(['cnv_id', 'base_id'], axis=1)
            
            cds_metrics_df = pd.concat([cds_metrics_df, merged], ignore_index=True)
            logger.info(f"    ✓ Added {len(merged)} CDS metric records for CNV/paralog copies")
    
    return mrna_metrics_df, cds_metrics_df


def reclassify_protein_only_novel(final_consensus, tx_dict, ref_df, metrics, genome,
                                  mrna_metrics_df=None, min_coverage=0.0,
                                  min_identity=0.0, intergenic_only=True,
                                  min_exons=2, min_cds_aa=100):
    """Promote orphan augMP models (protein-only, no reference transcript) to a
    clean protein-only novel-gene class so lineage-specific genes found ONLY in
    protein/miniprot evidence survive and are reported as novel.

    An augMP alignment whose (stripped) query name is not a reference transcript
    is a gene found solely from protein homology -- typically a lineage-specific
    gene with no reference ortholog. In the normal path these reach consensus but
    are labeled 'ortholog' under a synthetic ``UNKNOWN_GENE_*`` id, so downstream
    reports mistake them for reference loci instead of counting them as novel.

    Operating on the finalized ``final_consensus`` list (so the reference CDS
    footprint already includes rescued reference genes) this:

      * ignores poor models already downgraded to ``fragment`` /
        ``processed_pseudogene`` (only genuinely good models are promoted),
      * drops orphans below ``min_coverage`` / ``min_identity`` (miniprot mRNA
        metrics; both default to 0 = no extra filter, since augMP coverage is
        already gated upstream by ``augMP_min_coverage_no_anchor``),
      * collapses redundant per-locus copies (several species' orthologous
        proteins hitting one locus) to one representative (longest ORF),
      * when ``intergenic_only`` (default) drops any locus overlapping a
        reference-anchored CDS, keeping only genuinely novel intergenic loci,
      * relabels survivors as ``putative_novel`` protein_coding with stable
        ``MP-NOVEL-<genome>-N`` source ids (numbered by genomic position) and a
        blank ``source_gene_biotype`` so they are counted as novel genes.

    Returns the (possibly shorter) consensus list.
    """
    nc = tools.nameConversions

    def _versionless(x):
        return re.sub(r"\.[0-9]+$", "", str(x))

    ref_tx_ids = set()
    if ref_df is not None and 'TranscriptId' in getattr(ref_df, 'columns', []):
        ref_tx_ids = {_versionless(t) for t in ref_df['TranscriptId'].astype(str)}

    cov_map = {}
    if mrna_metrics_df is not None and 'AlignmentId' in getattr(mrna_metrics_df, 'columns', []):
        cov_col = 'AlnCoverage' if 'AlnCoverage' in mrna_metrics_df.columns else None
        id_col = 'AlnIdentity' if 'AlnIdentity' in mrna_metrics_df.columns else None
        if cov_col or id_col:
            for _, r in mrna_metrics_df.iterrows():
                cov_map[r['AlignmentId']] = (r.get(cov_col) if cov_col else None,
                                             r.get(id_col) if id_col else None)

    orphan_idx = []
    for i, (aln_id, attrs) in enumerate(final_consensus):
        if attrs.get('alignment_mode') != 'augMP' and not nc.aln_id_is_augustus_mp(str(aln_id)):
            continue
        if attrs.get('transcript_class') in ('fragment', 'processed_pseudogene'):
            continue
        base = _versionless(nc.strip_alignment_numbers(str(aln_id)))
        if base in ref_tx_ids:
            continue  # augMP refining a reference transcript -> genuine ortholog
        tx_obj = tx_dict.get(aln_id)
        if tx_obj is None or getattr(tx_obj, 'cds_size', 0) <= 0:
            continue
        orphan_idx.append(i)

    if not orphan_idx:
        return final_consensus

    orphan_pos = set(orphan_idx)

    # Exonic footprint of every reference-anchored (non-orphan) consensus model,
    # merged per (chrom, strand), for the intergenic test. A protein-only novel
    # model overlapping ANY projected exon -- coding OR non-coding (lncRNA,
    # pseudogene, retrocopy) -- on the same strand is almost always that projected
    # feature, not a genuinely new gene, so it is dropped. Merged intervals + bisect
    # keep this O(N log N) instead of O(orphans x reference).
    import bisect as _bisect
    fp_raw = collections.defaultdict(list)
    for i, (aln_id, attrs) in enumerate(final_consensus):
        if i in orphan_pos:
            continue
        tx_obj = tx_dict.get(aln_id)
        if tx_obj is None:
            continue
        strand = getattr(tx_obj, 'strand', None)
        for ex in getattr(tx_obj, 'exon_intervals', []) or []:
            fp_raw[(ex.chromosome, strand)].append((ex.start, ex.stop))
    footprint = {}
    for key, ivs in fp_raw.items():
        ivs.sort()
        starts, stops = [], []
        cs, ce = ivs[0]
        for s, e in ivs[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                starts.append(cs); stops.append(ce); cs, ce = s, e
        starts.append(cs); stops.append(ce)
        footprint[key] = (starts, stops)

    def _overlaps_reference(ci):
        f = footprint.get((ci.chromosome, ci.strand))
        if not f:
            return False
        starts, stops = f
        j = _bisect.bisect_left(starts, ci.stop) - 1
        return j >= 0 and stops[j] > ci.start

    drop_counts = collections.Counter()
    drop_pos = set()
    kept = []
    for i in orphan_idx:
        aln_id, attrs = final_consensus[i]
        tx_obj = tx_dict.get(aln_id)
        # Structural quality gates: single-exon models are dominated by retrocopies
        # / spurious hits, and very short ORFs are rarely real novel genes.
        n_ex = len(getattr(tx_obj, 'exon_intervals', []) or [])
        if min_exons and n_ex < int(min_exons):
            drop_pos.add(i); drop_counts['single_exon'] += 1; continue
        if min_cds_aa and getattr(tx_obj, 'cds_size', 0) < int(min_cds_aa) * 3:
            drop_pos.add(i); drop_counts['short_cds'] += 1; continue
        cov, ident = cov_map.get(aln_id, (None, None))
        if min_coverage and cov is not None and pd.notna(cov) and float(cov) < float(min_coverage):
            drop_pos.add(i); drop_counts['low_coverage'] += 1; continue
        if min_identity and ident is not None and pd.notna(ident) and float(ident) < float(min_identity):
            drop_pos.add(i); drop_counts['low_identity'] += 1; continue
        kept.append(i)

    clusters = []
    for i in kept:
        aln_id, attrs = final_consensus[i]
        ci = tx_dict[aln_id].coding_interval
        placed = False
        for cl in clusters:
            if cl['chrom'] == ci.chromosome and cl['strand'] == ci.strand and \
               ci.start < cl['stop'] and ci.stop > cl['start']:
                cl['members'].append(i)
                cl['start'] = min(cl['start'], ci.start)
                cl['stop'] = max(cl['stop'], ci.stop)
                placed = True
                break
        if not placed:
            clusters.append({'chrom': ci.chromosome, 'strand': ci.strand,
                             'start': ci.start, 'stop': ci.stop, 'members': [i]})

    clusters.sort(key=lambda c: (str(c['chrom']), c['start'], c['stop']))

    n_kept = 0
    novel_idx = 0
    for cl in clusters:
        members = cl['members']
        if intergenic_only:
            union = tools.intervals.ChromosomeInterval(cl['chrom'], cl['start'], cl['stop'], cl['strand'])
            if _overlaps_reference(union):
                drop_pos.update(members); drop_counts['overlaps_reference'] += len(members); continue
        rep = max(members, key=lambda i: getattr(tx_dict[final_consensus[i][0]], 'cds_size', 0))
        for i in members:
            if i != rep:
                drop_pos.add(i)
        novel_idx += 1
        gene_id = "MP-NOVEL-{}-{}".format(genome, novel_idx)
        attrs = final_consensus[rep][1]
        attrs['transcript_class'] = 'putative_novel'
        attrs['gene_biotype'] = 'protein_coding'
        attrs['transcript_biotype'] = 'protein_coding'
        attrs['source_gene'] = gene_id
        attrs['source_gene_biotype'] = 'N/A'
        attrs['source_gene_common_name'] = None
        attrs['protein_only_novel'] = 'True'
        n_kept += 1

    new_consensus = [entry for i, entry in enumerate(final_consensus) if i not in drop_pos]
    metrics.setdefault('Protein-only novel', {})
    metrics['Protein-only novel']['kept'] = n_kept
    metrics['Protein-only novel']['dropped'] = len(drop_pos)
    metrics['Protein-only novel']['dropped_single_exon'] = drop_counts['single_exon']
    metrics['Protein-only novel']['dropped_short_cds'] = drop_counts['short_cds']
    metrics['Protein-only novel']['dropped_low_coverage'] = drop_counts['low_coverage']
    metrics['Protein-only novel']['dropped_low_identity'] = drop_counts['low_identity']
    metrics['Protein-only novel']['dropped_overlaps_reference'] = drop_counts['overlaps_reference']
    return new_consensus


def rescue_augMP_at_empty_loci(final_consensus, tx_dict, ref_df, metrics, genome,
                               mrna_metrics_df=None, min_exons=2, min_cds_aa=100,
                               min_coverage=0.0, min_identity=0.0,
                               single_exon_min_cds_aa=300):
    """Recover high-quality augMP CDS models that consensus dropped, at loci that
    lack a same-strand protein-coding CDS in the consensus.

    Rationale (validated vs RefSeq): most RefSeq PC genes we fail to recover are
    still found by a raw mode -- for genes with no reference ortholog, augMP is
    frequently the ONLY mode that finds them, because it uses species/lineage-
    specific proteins from the expanded protein DB. Two failure modes are covered:

      1. the locus is completely empty (a truly novel lineage-specific gene), and
      2. the locus is occupied only by a NON-coding model (a pseudogene/lncRNA
         projection with no ORF) that beat the CDS-bearing augMP model -- confirmed
         on real genes (WASHC2A, MOXD2, KIR3DX1, PRSS45, ...). We treat a locus as
         available whenever it has no same-strand protein-coding CDS, so the augMP
         CDS model is added as a protein_coding gene alongside the non-coding call.

    Any augMP model in ``tx_dict`` that was NOT selected, has an ORF >=
    ``min_cds_aa`` residues, passes optional miniprot coverage/identity, is either
    multi-exon (>= ``min_exons``) OR single-exon with a long ORF (>=
    ``single_exon_min_cds_aa`` residues, to keep intronless genes like SMEK1 while
    excluding short retrocopies), and lands on a locus with no same-strand PC CDS is
    added back. Overlap with an existing same-strand PC gene is permitted only for
    retrocopies (single-exon models, e.g. retrogenes nested in a host gene's
    intron); multi-exon models overlapping an annotated PC gene are excluded as
    likely fragments of that gene. Per-locus duplicates collapse to the longest ORF. Models whose
    protein is a reference transcript are labelled recovered orthologs (source gene
    = that reference gene, counting toward reference recall / paralogs); the rest
    are labelled protein-only ``putative_novel``.

    Returns the (possibly longer) consensus list.
    """
    nc = tools.nameConversions
    import bisect as _bis

    def _versionless(x):
        return re.sub(r"\.[0-9]+$", "", str(x))

    # reference transcript -> (GeneId, biotype) for labelling recovered orthologs
    ref_tx_ids = set()
    tx_to_gene = {}
    if ref_df is not None and 'TranscriptId' in getattr(ref_df, 'columns', []):
        base = ref_df['TranscriptId'].astype(str).map(_versionless)
        ref_tx_ids = set(base)
        gid = ref_df['GeneId'] if 'GeneId' in ref_df.columns else base
        gbt = ref_df['GeneBiotype'] if 'GeneBiotype' in ref_df.columns else None
        for b, g, bt in zip(base, gid, (gbt if gbt is not None else [None] * len(base))):
            tx_to_gene[b] = (g, bt)

    cov_map = {}
    if mrna_metrics_df is not None and 'AlignmentId' in getattr(mrna_metrics_df, 'columns', []):
        cov_col = 'AlnCoverage' if 'AlnCoverage' in mrna_metrics_df.columns else None
        id_col = 'AlnIdentity' if 'AlnIdentity' in mrna_metrics_df.columns else None
        if cov_col or id_col:
            for _, r in mrna_metrics_df.iterrows():
                cov_map[r['AlignmentId']] = (r.get(cov_col) if cov_col else None,
                                             r.get(id_col) if id_col else None)

    # "Occupied" footprint: only same-strand PROTEIN-CODING CDS of consensus models.
    # A locus counts as available if it has no PC CDS here -- so pseudogene / lncRNA
    # only loci (no ORF) become eligible for a CDS-bearing augMP/augPB model, while
    # loci already annotated protein_coding are left untouched.
    selected_ids = set()
    fp_raw = collections.defaultdict(list)
    for aln_id, attrs in final_consensus:
        selected_ids.add(aln_id)
        if attrs.get('gene_biotype') != 'protein_coding':
            continue
        tx_obj = tx_dict.get(aln_id)
        if tx_obj is None:
            continue
        ci0 = getattr(tx_obj, 'coding_interval', None)
        if ci0 is None or ci0.start is None or ci0.stop is None or ci0.stop <= ci0.start:
            continue
        strand = getattr(tx_obj, 'strand', None)
        for ex in getattr(tx_obj, 'exon_intervals', []) or []:
            a = max(ex.start, ci0.start); b = min(ex.stop, ci0.stop)
            if b > a:
                fp_raw[(ex.chromosome, strand)].append((a, b))
    occupied = {}
    for key, ivs in fp_raw.items():
        ivs.sort(); starts, stops = [], []
        cs, ce = ivs[0]
        for s, e in ivs[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                starts.append(cs); stops.append(ce); cs, ce = s, e
        starts.append(cs); stops.append(ce)
        occupied[key] = (starts, stops)

    def _occupied(ci):
        f = occupied.get((ci.chromosome, ci.strand))
        if not f:
            return False
        starts, stops = f
        j = _bis.bisect_left(starts, ci.stop) - 1
        return j >= 0 and stops[j] > ci.start

    drop_counts = collections.Counter()
    candidates = []
    for aln_id, tx_obj in tx_dict.items():
        if aln_id in selected_ids:
            continue
        # augMP (protein evidence) or augPB (IsoSeq denovo): both carry a real ORF.
        if not (nc.aln_id_is_augustus_mp(str(aln_id)) or nc.aln_id_is_pb(str(aln_id))):
            continue
        if tx_obj is None or getattr(tx_obj, 'cds_size', 0) <= 0:
            continue
        n_ex = len(getattr(tx_obj, 'exon_intervals', []) or [])
        cds_aa = getattr(tx_obj, 'cds_size', 0) // 3
        # single-exon allowed only with a long ORF (keeps SMEK1-type intronless
        # genes, excludes short retrocopies); multi-exon uses the normal floor.
        if n_ex < int(min_exons):
            if n_ex >= 1 and single_exon_min_cds_aa and cds_aa >= int(single_exon_min_cds_aa):
                pass
            else:
                drop_counts['single_exon'] += 1; continue
        if min_cds_aa and getattr(tx_obj, 'cds_size', 0) < int(min_cds_aa) * 3:
            drop_counts['short_cds'] += 1; continue
        cov, ident = cov_map.get(aln_id, (None, None))
        if min_coverage and cov is not None and pd.notna(cov) and float(cov) < float(min_coverage):
            drop_counts['low_coverage'] += 1; continue
        if min_identity and ident is not None and pd.notna(ident) and float(ident) < float(min_identity):
            drop_counts['low_identity'] += 1; continue
        ci = getattr(tx_obj, 'coding_interval', None)
        if ci is None or ci.start is None or ci.stop is None:
            continue
        # Overlap with an existing same-strand PC gene is allowed ONLY for
        # retrocopies (single-exon, intronless models -- e.g. a retrogene inserted
        # into a host gene's intron/exon). Multi-exon models overlapping an already
        # annotated PC gene are almost always fragments of that same gene, so keep
        # excluding them.
        if _occupied(ci) and n_ex > 1:
            drop_counts['locus_occupied'] += 1; continue
        candidates.append((aln_id, tx_obj, ci))

    if not candidates:
        metrics.setdefault('augMP empty-locus recovery', {})
        metrics['augMP empty-locus recovery'].update(
            {'recovered_genes': 0, 'recovered_orthologs': 0, 'recovered_novel': 0,
             'dropped_single_exon': drop_counts['single_exon'],
             'dropped_short_cds': drop_counts['short_cds'],
             'dropped_locus_occupied': drop_counts['locus_occupied']})
        return final_consensus

    # Cluster candidates by CDS overlap on the same strand; keep the longest ORF.
    candidates.sort(key=lambda c: (str(c[2].chromosome), c[2].strand or '.', c[2].start, c[2].stop))
    clusters = []
    for aln_id, tx_obj, ci in candidates:
        if clusters:
            cl = clusters[-1]
            if cl['chrom'] == ci.chromosome and cl['strand'] == ci.strand and ci.start < cl['stop']:
                cl['members'].append((aln_id, tx_obj, ci))
                cl['stop'] = max(cl['stop'], ci.stop)
                continue
        clusters.append({'chrom': ci.chromosome, 'strand': ci.strand,
                         'start': ci.start, 'stop': ci.stop,
                         'members': [(aln_id, tx_obj, ci)]})

    n_orth = 0
    n_novel = 0
    novel_idx = 0
    for cl in clusters:
        aln_id, tx_obj, ci = max(cl['members'], key=lambda m: getattr(m[1], 'cds_size', 0))
        base = _versionless(nc.strip_alignment_numbers(str(aln_id)))
        attrs = {
            'alignment_id': aln_id,
            'alignment_mode': 'augMP',
            'gene_biotype': 'protein_coding',
            'transcript_biotype': 'protein_coding',
            'score': 0,
            'augMP_recovered': 'True',
        }
        if base in ref_tx_ids:
            g, bt = tx_to_gene.get(base, (base, 'protein_coding'))
            attrs['source_gene'] = g
            attrs['source_gene_biotype'] = bt or 'protein_coding'
            attrs['source_gene_common_name'] = None
            attrs['transcript_class'] = 'ortholog'
            n_orth += 1
        else:
            novel_idx += 1
            gene_id = "MP-RECOVERED-{}-{}".format(genome, novel_idx)
            attrs['source_gene'] = gene_id
            attrs['source_gene_biotype'] = 'N/A'
            attrs['source_gene_common_name'] = None
            attrs['transcript_class'] = 'putative_novel'
            attrs['protein_only_novel'] = 'True'
            n_novel += 1
        final_consensus.append((aln_id, attrs))

    metrics.setdefault('augMP empty-locus recovery', {})
    metrics['augMP empty-locus recovery'].update(
        {'recovered_genes': n_orth + n_novel,
         'recovered_orthologs': n_orth,
         'recovered_novel': n_novel,
         'dropped_single_exon': drop_counts['single_exon'],
         'dropped_short_cds': drop_counts['short_cds'],
         'dropped_locus_occupied': drop_counts['locus_occupied']})
    return final_consensus


def rescue_expressed_noncoding_to_pc(final_consensus, tx_dict, metrics,
                                     min_cds_aa=100, require_multiexon=True,
                                     require_protein_evidence_for_noncoding=True):
    """Promote an expressed, ORF-intact model to protein_coding when it carries a
    NON-coding gene_biotype only because it was inherited from a non-coding
    reference ortholog (lncRNA / pseudogene / etc.).

    This recovers real protein-coding genes whose reference (human) ortholog is
    annotated non-coding but which are protein-coding in this species -- a genuine
    recall gap confirmed against RefSeq. It is deliberately strict to avoid the
    intergenic false positives seen with protein-only calls: a model is promoted
    only if it has an intact ORF (valid start/stop, proper frame), a CDS of at
    least ``min_cds_aa`` residues, is multi-exon (unless disabled), AND has direct
    transcription evidence (IsoSeq isoform support or RNA intron/exon support).
    Models with no expression evidence (e.g. IsoSeq-less genomes) are never
    promoted, keeping this precise. Promotion is applied at the gene level (all
    isoforms of a qualifying gene become protein_coding); each isoform's own
    transcript_biotype is set to protein_coding only if it is itself ORF-intact.

    Returns the consensus list (modified in place).
    """
    coding = 'protein_coding'

    def has_expression(attrs):
        if attrs.get('pacbio_isoform_supported') is True:
            return True
        for k in ('exon_rna_support', 'intron_rna_support'):
            v = attrs.get(k, '')
            if isinstance(v, str) and '1' in v.split(','):
                return True
        return False

    def orf_intact(attrs, aln_id):
        if not (attrs.get('proper_orf') and attrs.get('valid_start') and attrs.get('valid_stop')):
            return False
        tx = tx_dict.get(aln_id)
        return tx is not None and getattr(tx, 'cds_size', 0) >= int(min_cds_aa) * 3

    # Protein-evidence gate for confidently non-coding reference biotypes.
    # A reference lncRNA/pseudogene that merely carries an incidental ORF should NOT
    # become protein_coding on expression alone: validated against RefSeq, only ~9%
    # of such expressed-ORF lncRNAs are real PC, vs ~84% of those that ALSO have a
    # same-strand augMP (miniprot-hinted Augustus) CDS. Requiring augMP support keeps
    # the genuine recoveries and drops the incidental-ORF ones. Denovo/unknown and
    # IG/TR biotypes are protein-coding-competent by construction and are never gated.
    def _needs_protein_gate(attrs):
        src = attrs.get('source_gene_biotype') or ''
        return src == 'lncRNA' or 'pseudogene' in src or src in ('misc_RNA', 'TEC')

    # Merged augMP CDS-span footprint per (chrom, strand) from ALL augMP inputs.
    nc = tools.nameConversions
    augmp_fp = {}
    if require_protein_evidence_for_noncoding:
        import bisect as _bis
        raw = collections.defaultdict(list)
        for aid, tx in tx_dict.items():
            if not nc.aln_id_is_augustus_mp(str(aid)):
                continue
            ci = getattr(tx, 'coding_interval', None)
            if ci is None or ci.start is None or ci.stop is None or ci.stop <= ci.start:
                continue
            raw[(ci.chromosome, ci.strand)].append((ci.start, ci.stop))
        for key, ivs in raw.items():
            ivs.sort()
            starts, stops = [], []
            cs, ce = ivs[0]
            for s, e in ivs[1:]:
                if s <= ce:
                    ce = max(ce, e)
                else:
                    starts.append(cs); stops.append(ce); cs, ce = s, e
            starts.append(cs); stops.append(ce)
            augmp_fp[key] = (starts, stops)

        def _has_augmp_support(aln_id):
            tx = tx_dict.get(aln_id)
            ci = getattr(tx, 'coding_interval', None) if tx is not None else None
            if ci is None or ci.start is None or ci.stop is None:
                return False
            f = augmp_fp.get((ci.chromosome, ci.strand))
            if not f:
                return False
            starts, stops = f
            j = _bis.bisect_left(starts, ci.stop) - 1
            return j >= 0 and stops[j] > ci.start
    else:
        def _has_augmp_support(aln_id):
            return True

    n_gated = 0

    # Pass 1: find qualifying (coding-competent, expressed) non-PC transcripts.
    promote_genes = set()
    qualifying_idx = set()
    for i, (aln_id, attrs) in enumerate(final_consensus):
        if attrs.get('gene_biotype') == coding:
            continue
        if attrs.get('transcript_class') in ('fragment', 'processed_pseudogene'):
            continue
        if not orf_intact(attrs, aln_id):
            continue
        if require_multiexon:
            tx = tx_dict.get(aln_id)
            if len(getattr(tx, 'exon_intervals', []) or []) < 2:
                continue
        if not has_expression(attrs):
            continue
        if _needs_protein_gate(attrs) and not _has_augmp_support(aln_id):
            n_gated += 1
            continue
        qualifying_idx.add(i)
        g = attrs.get('source_gene')
        if g and g != 'N/A':
            promote_genes.add(g)

    # Pass 2: promote qualifying genes (gene-level) + orphan qualifying transcripts.
    genes_promoted = set()
    n_tx = 0
    for i, (aln_id, attrs) in enumerate(final_consensus):
        g = attrs.get('source_gene')
        gene_hit = (g in promote_genes) if (g and g != 'N/A') else (i in qualifying_idx)
        if not gene_hit or attrs.get('gene_biotype') == coding:
            continue
        attrs['gene_biotype'] = coding
        if orf_intact(attrs, aln_id):
            attrs['transcript_biotype'] = coding
        attrs['expressed_pc_rescue'] = 'True'
        genes_promoted.add(g if (g and g != 'N/A') else aln_id)
        n_tx += 1

    metrics.setdefault('Expressed noncoding->PC rescue', {})
    metrics['Expressed noncoding->PC rescue']['genes'] = len(genes_promoted)
    metrics['Expressed noncoding->PC rescue']['transcripts'] = n_tx
    metrics['Expressed noncoding->PC rescue']['gated_no_protein_evidence'] = n_gated
    return final_consensus


def generate_consensus(args):
    """
    Main consensus finding logic with chromosome-based processing.
    Uses local multiprocessing only (no SLURM).
    """
    start_time = time.time()
    
    logger.info("\n" + "="*80)
    logger.info("STEP 1: Loading Input Data")
    logger.info("="*80)
    logger.info(f"Input genePred files: {len(args.gp_list)}")
    for gp_file in args.gp_list:
        logger.info(f"  - {gp_file}")
    
    logger.info("\nBuilding alignment source map...")
    
    # Create a mapping from alignment IDs to their source files
    alignment_source_map = {}
    tx_by_mode = collections.defaultdict(dict)
    
    # Load genePreds and track their sources
    for gp_idx, gp_file in enumerate(args.gp_list, 1):
        logger.info(f"  Loading {gp_idx}/{len(args.gp_list)}: {gp_file}...")
        tx_count = 0
        
        # Determine mode from filename
        # Check more specific patterns first to avoid false matches
        if 'transMap_pairwise' in gp_file:
            mode = 'transMap_pairwise'  # Keep as separate mode
        elif 'transMap' in gp_file:
            mode = 'transMap'
        elif 'txTM' in gp_file:
            mode = 'txTM'
        elif 'augTMR_pairwise' in gp_file:
            mode = 'augTMR_pairwise'  # Keep as separate mode
        elif 'augTMR' in gp_file:
            mode = 'augTMR'
        elif 'augTM_pairwise' in gp_file:
            mode = 'augTM_pairwise'  # Keep as separate mode
        elif 'augTM' in gp_file:
            mode = 'augTM'
        elif 'augMP' in gp_file:
            mode = 'augMP'  # Keep as separate mode
        elif 'augPB' in gp_file:
            mode = 'augPB'
        else:
            mode = 'unknown'
        
        # Load transcripts - handle duplicates within the same mode (e.g., paralogs, CNV copies)
        for t in tools.transcripts.gene_pred_iterator(gp_file):
            tx_count += 1
            
            # VALIDATE: Check if transcript has valid coordinates
            # Skip transcripts with suspiciously large spans (likely chimeric from gff3ToGenePred merging
            # transcripts with duplicate IDs on different chromosomes)
            transcript_span = t.stop - t.start
            if transcript_span > 10_000_000:  # 10 Mb threshold
                logger.warning(f"  Skipping transcript {t.name} on {t.chromosome}: span too large ({transcript_span:,} bp). "
                             f"Coordinates: {t.start}-{t.stop}. Likely chimeric from duplicate GFF3 IDs on different chromosomes.")
                continue
            
            # Use pre-determined mode or infer from ID if unknown
            if mode == 'unknown':
                inferred_mode = identify_mode(t.name, gp_file)
                actual_mode = inferred_mode
            else:
                actual_mode = mode
            
            # Normalize txTM IDs BEFORE duplicate detection to catch CNV copies at different loci
            # TxTM uses _N suffix (e.g., gene_1, gene_2) for CNV copies that should be preserved
            if actual_mode == 'txTM':
                normalized_name = normalize_alignment_id(t.name, actual_mode)
                if normalized_name != t.name:
                    #logger.debug(f"    Normalized txTM ID: {t.name} → {normalized_name}")
                    t.name = normalized_name
                # Also normalize the gene ID (name2 field) to strip txTM's _N suffix
                if t.name2:
                    normalized_gene = normalize_gene_id(t.name2, actual_mode)
                    if normalized_gene != t.name2:
                        #logger.debug(f"    Normalized txTM gene ID: {t.name2} → {normalized_gene}")
                        t.name2 = normalized_gene
            
            # Handle duplicates within the same mode (isoforms vs paralogs at different loci)
            original_name = t.name
            original_gene_name = t.name2  # Save original gene ID
            if original_name in tx_by_mode[actual_mode]:
                existing = tx_by_mode[actual_mode][original_name]
                # Check if transcripts are identical (same hash AND same location)
                same_location = (existing.chromosome == t.chromosome and 
                               existing.start == t.start and 
                               existing.stop == t.stop)
                if hash(existing) == hash(t) and same_location:
                    # Truly identical transcript at same location - skip
                    continue
                
                # Check if transcripts overlap (isoforms) or are at different loci (paralogs)
                overlaps = (existing.chromosome == t.chromosome and
                           t.start <= existing.stop and t.stop >= existing.start)
                
                # Different transcript - add with unique suffix to transcript ID
                suffix_num = 2
                new_name = f"{original_name}_cp{suffix_num}"
                while new_name in tx_by_mode[actual_mode]:
                    suffix_num = suffix_num + 1
                    new_name = f"{original_name}_cp{suffix_num}"
                t.name = new_name
                
                # CRITICAL: Only update gene ID (name2) if transcripts DON'T overlap
                # - Overlapping transcripts → isoforms of same gene → keep same gene ID
                # - Non-overlapping transcripts → paralogs → different gene IDs
                if not overlaps and t.name2:
                    t.name2 = f"{original_gene_name}_cp{suffix_num}"
                    #logger.info(f"    Found paralog '{original_name}' at non-overlapping locus - renamed to '{t.name}' (gene: {t.name2})")
                #else:
                    # Keep same gene ID for overlapping isoforms
                    #logger.info(f"    Found isoform '{original_name}' at overlapping locus - renamed to '{t.name}' (same gene: {t.name2})")
            
            # CROSS-MODE duplicate handling: If this transcript ID exists in alignment_source_map 
            # (meaning another mode already has it), add a suffix to keep both versions
            # This ensures ALL transcripts from ALL modes are included in consensus generation
            if t.name in alignment_source_map and alignment_source_map[t.name] != actual_mode:
                # Transcript ID exists in another mode - add suffix to distinguish
                original_name = t.name
                suffix_num = 2
                new_name = f"{original_name}_cp{suffix_num}"
                # Check across ALL modes to find an unused name
                while any(new_name in mode_txs for mode_txs in tx_by_mode.values()) or new_name in alignment_source_map:
                    suffix_num += 1
                    new_name = f"{original_name}_cp{suffix_num}"
                logger.debug(f"    Cross-mode duplicate: {t.name} exists in {alignment_source_map[t.name]}, renaming {actual_mode} version to {new_name}")
                t.name = new_name
            
            # Store the transcript
            tx_by_mode[actual_mode][t.name] = t
            
            # Track in alignment_source_map (now each transcript has a unique ID)
            alignment_source_map[t.name] = actual_mode
        
        logger.info(f"    ✓ Loaded {tx_count} transcripts (mode: {actual_mode})")
    
    # Report cross-mode duplicate handling
    cross_mode_duplicates = sum(1 for tx_id in alignment_source_map.keys() if re.search(r'_cp\d+$', tx_id))
    if cross_mode_duplicates > 0:
        logger.info(f"\n  ℹ️  Found {cross_mode_duplicates} cross-mode duplicate transcript IDs")
        logger.info(f"     (renamed with _cp suffix to include all versions in consensus)")
    
    # Set the global mapping for use by alignment_type function
    tools.nameConversions.set_alignment_source_map(alignment_source_map)
    
    # Build tx_dict from tx_by_mode (which has all the _cp suffix handling)
    logger.info("\nBuilding tx_dict from loaded transcripts...")
    tx_dict = {}
    for mode, txs in tx_by_mode.items():
        tx_dict.update(txs)
    
    logger.info(f"\n✓ Loaded {len(tx_dict)} total transcripts")
    
    # Load .gp_attrs files for biotype information (especially for txTM)
    logger.info("\nLoading genePred attribute files (.gp_attrs) for biotype information...")
    gp_attrs_transcript_biotypes = {}  # transcript_id -> transcript biotype
    gp_attrs_gene_biotypes = {}  # transcript_id -> gene biotype
    for gp_file in args.gp_list:
        attrs_file = gp_file + '_attrs'
        if os.path.exists(attrs_file):
            logger.info(f"  Loading {attrs_file}...")
            with open(attrs_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        transcript_id, attr_name, attr_value = parts[0], parts[1], parts[2]
                        # Store transcript-level and gene-level biotypes separately
                        if attr_name in ['transcript_biotype', 'transcript_type']:
                            gp_attrs_transcript_biotypes[transcript_id] = attr_value
                        elif attr_name in ['gene_biotype', 'gene_type']:
                            gp_attrs_gene_biotypes[transcript_id] = attr_value
    
    logger.info(f"✓ Loaded biotype info for {len(gp_attrs_transcript_biotypes)} transcripts from .gp_attrs files")
    logger.info(f"  Gene-level biotypes: {len(gp_attrs_gene_biotypes)} transcripts")
    logger.info(f"  Modes found: {list(tx_by_mode.keys())}")
    for mode, txs in tx_by_mode.items():
        logger.info(f"    {mode}: {len(txs)} transcripts")
    
    # Update alignment_source_map for transcripts that were renamed during loading
    # (e.g., cross-mode conflicts where load_gps added "_cp", "_cp2", "_cp3" suffix)
    # Strategy: The _cp version was created when load_gps found a duplicate ID
    # We need to find which mode had the transcript that got renamed
    # by comparing hash and location with all modes' transcripts
    logger.info("  Updating alignment_source_map for renamed transcripts...")
    renamed_count = 0
    
    # Build a faster lookup: map (chrom, start, stop) -> mode for original IDs
    # This avoids nested loops
    location_to_mode = {}
    for mode, mode_txs in tx_by_mode.items():
        for tx_id, tx_obj in mode_txs.items():
            # Only index non-_cp versions (originals)
            if not re.search(r'_cp\d*$', tx_id):
                location_key = (tx_obj.chromosome, tx_obj.start, tx_obj.stop)
                if location_key not in location_to_mode:
                    location_to_mode[location_key] = mode
    
    # Now process only _cp transcripts
    _cp_tx_ids = [tx_id for tx_id in tx_dict.keys() if re.search(r'_cp\d*$', tx_id) and tx_id not in alignment_source_map]
    #logger.info(f"    Found {len(_cp_tx_ids)} _cp transcripts to process...")
    
    for tx_id in _cp_tx_ids:
        cp_tx = tx_dict[tx_id]
        location_key = (cp_tx.chromosome, cp_tx.start, cp_tx.stop)
        
        # Fast lookup
        if location_key in location_to_mode:
            found_mode = location_to_mode[location_key]
            alignment_source_map[tx_id] = found_mode
            renamed_count += 1
    
    if renamed_count > 0:
        logger.info(f"  ✓ Updated alignment_source_map for {renamed_count} renamed transcripts (cross-mode conflicts)")
    else:
        logger.info(f"  ✓ No cross-mode conflicts found")
    
    # Track unique genes from each source
    logger.info("\n  Gene counts by source:")
    genes_by_source = {}
    for mode, txs in tx_by_mode.items():
        logger.info(f"    Processing {mode} ({len(txs)} transcripts)...")
        mode_genes = set()
        for tx_id, tx_obj in txs.items():
            # Extract gene from transcript ID
            # For transMap/txTM, the transcript ID itself can be used to get gene from ref_df later
            # For now, just track transcript IDs as proxy
            mode_genes.add(tx_id)
        genes_by_source[mode] = mode_genes
        logger.info(f"    {mode}: {len(mode_genes)} transcript IDs")
    
    # Load reference annotation information
    logger.info("\nLoading reference annotation...")
    ref_df = tools.sqlInterface.load_annotation(args.ref_db_path)
    logger.info(f"✓ Loaded {len(ref_df)} reference annotations")
    ref_biotype_counts = collections.Counter(ref_df.TranscriptBiotype)
    coding_count = ref_biotype_counts['protein_coding']
    non_coding_count = sum(y for x, y in ref_biotype_counts.items() if x != 'protein_coding')
    
    # Create readthrough gene map from ExtraTags (check both readthrough_gene and readthrough_transcript)
    readthrough_gene_set = set()
    if 'ExtraTags' in ref_df.columns:
        # Some annotations have readthrough_gene tag, some have readthrough_transcript tag
        readthrough_mask = (ref_df['ExtraTags'].str.contains('readthrough_gene', na=False) | 
                           ref_df['ExtraTags'].str.contains('readthrough_transcript', na=False))
        readthrough_gene_set = set(ref_df[readthrough_mask]['GeneId'].unique())
        logger.info(f"  Found {len(readthrough_gene_set)} readthrough genes in reference")
    
    # Build reference gene coordinate map for checking overlaps
    # This will be used to check if overlapping genes in the target also overlap in the reference
    logger.info("  Building reference gene coordinate map...")
    ref_gene_coords = {}
    genes_with_overlaps_in_ref = set()  # Track which genes overlap with others in reference
    
    if 'GeneId' in ref_df.columns and 'ExtraTags' in ref_df.columns:
        # Extract coordinates from ExtraTags using vectorized string operations
        # ExtraTags format: "...;Seqid=chr1;...;Start=12345;...;End=67890;..."
        ref_df_with_coords = ref_df[ref_df['ExtraTags'].notna()].copy()
        
        # Extract Seqid, Start, End using regex
        ref_df_with_coords['Seqid'] = ref_df_with_coords['ExtraTags'].str.extract(r'Seqid=([^;]+)', expand=False)
        ref_df_with_coords['Start'] = ref_df_with_coords['ExtraTags'].str.extract(r'Start=(\d+)', expand=False)
        ref_df_with_coords['End'] = ref_df_with_coords['ExtraTags'].str.extract(r'End=(\d+)', expand=False)
        
        # Filter rows where all three fields are present
        ref_df_with_coords = ref_df_with_coords[
            ref_df_with_coords['Seqid'].notna() & 
            ref_df_with_coords['Start'].notna() & 
            ref_df_with_coords['End'].notna()
        ]
        
        if len(ref_df_with_coords) > 0:
            # Convert Start and End to integers
            ref_df_with_coords['Start'] = ref_df_with_coords['Start'].astype(int)
            ref_df_with_coords['End'] = ref_df_with_coords['End'].astype(int)
            
            # Group by GeneId and get min Start, max End
            gene_coords_grouped = ref_df_with_coords.groupby('GeneId').agg({
                'Seqid': 'first',  # Chromosome should be same for all transcripts of a gene
                'Start': 'min',
                'End': 'max'
            })
            
            # Convert to dictionary (vectorized - avoid iterrows)
            ref_gene_coords = {
                gene_id: (row.Seqid, row.Start, row.End)
                for gene_id, row in zip(gene_coords_grouped.index, gene_coords_grouped.itertuples(index=False))
            }
            
            # Find which genes overlap with others in the reference
            # This helps identify gene families that naturally overlap
            logger.info("  Identifying genes that overlap in reference...")
            ref_genes_by_chrom = collections.defaultdict(list)
            for gene_id, (chrom, start, end) in ref_gene_coords.items():
                ref_genes_by_chrom[chrom].append((start, end, gene_id))
            
            for chrom, gene_list in ref_genes_by_chrom.items():
                gene_list.sort()  # Sort by start position
                for i in range(len(gene_list)):
                    start_i, end_i, gene_i = gene_list[i]
                    for j in range(i + 1, len(gene_list)):
                        start_j, end_j, gene_j = gene_list[j]
                        if start_j >= end_i:  # No more overlaps possible
                            break
                        # Check if they overlap
                        if not (end_i <= start_j or end_j <= start_i):
                            genes_with_overlaps_in_ref.add(gene_i)
                            genes_with_overlaps_in_ref.add(gene_j)
            
            logger.info(f"  Found {len(genes_with_overlaps_in_ref)} genes that overlap with others in reference")
            
        logger.info(f"  Built reference gene coordinates for {len(ref_gene_coords)} genes")
    
    gene_biotype_map = tools.sqlInterface.get_gene_biotype_map(args.ref_db_path)
    transcript_biotype_map = tools.sqlInterface.get_transcript_biotype_map(args.ref_db_path)
    logger.info(f"  Reference has {len(gene_biotype_map)} genes, {len(transcript_biotype_map)} transcripts")
    logger.info(f"  Biotype breakdown: {dict(ref_biotype_counts)}")
    
    # Map transcript IDs to gene IDs for tracking (vectorized for speed)
    ref_df['TranscriptId_base'] = ref_df['TranscriptId'].astype(str).str.replace(r"\.[0-9]+$", "", regex=True)
    tx_to_gene_map = dict(zip(ref_df['TranscriptId_base'], ref_df['GeneId']))
    
    # Now track actual gene IDs from each source
    logger.info("\n  Mapping to actual gene IDs:")
    gene_ids_by_source = {}
    pc_gene_ids_by_source = {}
    for mode, tx_ids in genes_by_source.items():
        mode_gene_ids = set()
        mode_pc_gene_ids = set()
        tx_found = 0
        tx_not_found = 0
        for tx_id in tx_ids:
            # First try to get gene ID directly from the transcript object (tx.name2)
            # This is much more reliable than mapping through reference database
            if tx_id in tx_dict:
                tx_obj = tx_dict[tx_id]
                gene_id = tx_obj.name2  # Gene ID from genePred name2 field
                tx_found += 1
            else:
                # Fallback: try to map through reference database (for backwards compatibility)
                normalized_id = normalize_alignment_id(tx_id, mode)
                normalized_id_base = re.sub(r'\.[0-9]+$', '', normalized_id)
                gene_id = tx_to_gene_map.get(normalized_id_base)
                tx_not_found += 1
            
            if gene_id:
                # Strip _cp suffix from gene_id for normalization (used for grouping unique genes)
                # The suffix is internal bookkeeping and doesn't represent a different gene
                gene_id_normalized = re.sub(r'_cp\d+$', '', gene_id)
                mode_gene_ids.add(gene_id_normalized)
                
                # Check if protein coding
                # Strip _cp from transcript ID and gene ID for biotype lookup
                tx_id_for_lookup = re.sub(r'_cp\d+$', '', tx_id)
                gene_id_for_lookup = gene_id_normalized
                
                # First check gene_biotype from .gp_attrs (most reliable for TxTM)
                is_pc = gp_attrs_gene_biotypes.get(tx_id_for_lookup) == 'protein_coding'
                # Fallback to gene_biotype_map from reference database
                if not is_pc:
                    is_pc = gene_biotype_map.get(gene_id_for_lookup) == 'protein_coding'
                if is_pc:
                    mode_pc_gene_ids.add(gene_id_normalized)
        gene_ids_by_source[mode] = mode_gene_ids
        pc_gene_ids_by_source[mode] = mode_pc_gene_ids
        logger.info(f"    {mode}: {len(mode_gene_ids)} genes ({len(mode_pc_gene_ids)} protein-coding) - {tx_found} tx found in dict, {tx_not_found} not found")
    
    # Load transMap evaluation data
    logger.info("\nLoading transMap evaluation data...")
    tm_eval_df = load_transmap_evals(args.db_path)
    logger.info(f"✓ Loaded {len(tm_eval_df)} transMap evaluations")
    
    # Determine which modes are available
    tx_modes = list(tx_by_mode.keys())
    tx_modes_with_metrics = [x for x in tx_modes if x not in ('augPB', 'strg')]
    
    logger.info("\n" + "="*80)
    logger.info("STEP 2: Loading Alignment Metrics")
    logger.info("="*80)
    # Load metrics (filter to actual transcripts for speed)
    valid_aln_ids = set(tx_dict.keys())
    
    if tx_modes_with_metrics:
        logger.info(f"Loading metrics for modes: {tx_modes_with_metrics}")
        txtm_metrics_match_keys = _build_txtm_metrics_match_keys(valid_aln_ids, alignment_source_map)
        # Precompute the acceptable txTM DB ids once so the per-mode filters below
        # can use a single vectorized isin() instead of a row-wise Python loop.
        txtm_allowed_db_ids = _txtm_allowed_db_ids(txtm_metrics_match_keys)
        
        mrna_dfs = []
        for tx_mode in tx_modes_with_metrics:
            # Map display mode to database table name
            db_mode = map_mode_to_db_table(tx_mode)
            df = load_metrics_from_db(args.db_path, db_mode, 'mRNA')
            # Add mode column to track source
            df['Mode'] = tx_mode
            # Filter early to reduce memory and speed up merges
            # For txTM, database has base IDs (ENST00000002596.6) but tx_dict has CNV copies (ENST00000002596.6_1, etc.)
            # So we need to check if the base ID or any CNV copy is in valid_aln_ids
            if tx_mode == 'txTM':
                df = df[df['AlignmentId'].astype(str).isin(txtm_allowed_db_ids)]
            else:
                df = df[df['AlignmentId'].isin(valid_aln_ids)]
            mrna_dfs.append(df)
        mrna_metrics_df = pd.concat(mrna_dfs) if mrna_dfs else pd.DataFrame()
        if len(mrna_metrics_df) > 0 and 'txTM' in tx_modes_with_metrics:
            mrna_metrics_df = expand_txtm_cnv_metrics(mrna_metrics_df, valid_aln_ids)
            mrna_metrics_df = remap_metrics_to_txtm_cp_aliases(
                mrna_metrics_df, valid_aln_ids, alignment_source_map
            )
        logger.info(f"✓ Loaded {len(mrna_metrics_df)} mRNA metrics (filtered to actual transcripts)")
        
        cds_dfs = []
        for tx_mode in tx_modes_with_metrics:
            # Map display mode to database table name
            db_mode = map_mode_to_db_table(tx_mode)
            df = load_metrics_from_db(args.db_path, db_mode, 'CDS')
            # Add mode column to track source
            df['Mode'] = tx_mode
            # Same logic as mRNA: for txTM, match base IDs to CNV copies
            if tx_mode == 'txTM':
                df = df[df['AlignmentId'].astype(str).isin(txtm_allowed_db_ids)]
            else:
                df = df[df['AlignmentId'].isin(valid_aln_ids)]
            cds_dfs.append(df)
        cds_metrics_df = pd.concat(cds_dfs) if cds_dfs else pd.DataFrame()
        if len(cds_metrics_df) > 0 and 'txTM' in tx_modes_with_metrics:
            cds_metrics_df = expand_txtm_cnv_metrics(cds_metrics_df, valid_aln_ids)
            cds_metrics_df = remap_metrics_to_txtm_cp_aliases(
                cds_metrics_df, valid_aln_ids, alignment_source_map
            )
        logger.info(f"✓ Loaded {len(cds_metrics_df)} CDS metrics (filtered to actual transcripts)")
        
        # Backfill metrics for txTM CNV copies (_N suffix) that don't have metrics
        mrna_metrics_df, cds_metrics_df = backfill_cnv_metrics(
            mrna_metrics_df, cds_metrics_df, valid_aln_ids, alignment_source_map, args.db_path
        )
        
        eval_dfs = []
        for tx_mode in tx_modes_with_metrics:
            # Map display mode to database table name
            db_mode = map_mode_to_db_table(tx_mode)
            df = load_evaluations_from_db(args.db_path, db_mode)
            df = df[df.index.isin(valid_aln_ids)]
            eval_dfs.append(df)
        eval_df = pd.concat(eval_dfs).reset_index() if eval_dfs else pd.DataFrame()
        logger.info(f"✓ Loaded {len(eval_df)} evaluation entries (filtered to actual transcripts)")
    else:
        logger.warning("No modes with metrics found")
        mrna_metrics_df = pd.DataFrame()
        cds_metrics_df = pd.DataFrame()
        eval_df = pd.DataFrame()
    
    # Create support dataframe
    logger.info("\nCreating support dataframe...")
    logger.info(f"  args.denovo_tx_modes = {args.denovo_tx_modes} (type: {type(args.denovo_tx_modes)})")
    support_df = create_support_dataframe(
        tx_dict, args.db_path, ref_df, alignment_source_map,
        denovo_tx_modes=args.denovo_tx_modes,
        bam_files=args.bam_files,
        isoseq_bam_files=args.isoseq_bam_files,
        ref_gp_path=args.ref_gp
    )
    
    # Process by chromosome
    logger.info("\n" + "="*80)
    logger.info("STEP 3: Grouping Transcripts by Chromosome")
    logger.info("="*80)
    tx_by_chrom = collections.defaultdict(list)
    for aln_id, tx_obj in tx_dict.items():
        tx_by_chrom[tx_obj.chromosome].append(aln_id)

    if args.only_chrom:
        only = str(args.only_chrom)
        if only not in {str(c) for c in tx_by_chrom.keys()}:
            available = ", ".join(sorted({str(c) for c in tx_by_chrom.keys()}))
            raise ValueError(f"--only-chrom {only} not found in inputs. Available: {available}")
        # Normalize keys by string compare, then filter.
        tx_by_chrom = collections.defaultdict(list, {c: v for c, v in tx_by_chrom.items() if str(c) == only})
    
    logger.info(f"✓ Found {len(tx_by_chrom)} chromosomes")
    for chrom in sorted(tx_by_chrom.keys()):
        logger.info(f"  {chrom}: {len(tx_by_chrom[chrom])} transcripts")
    
    logger.info("\n" + "="*80)
    logger.info("STEP 4: Processing Chromosomes (Local Multiprocessing)")
    logger.info("="*80)
    
    # Use multiprocessing mode (LOCAL ONLY)
    # OPTIMIZATION: Match workers to chromosomes to avoid extra serialization overhead
    # With Python multiprocessing, each worker requires pickling all shared data
    max_workers = args.num_workers if args.num_workers else min(mp.cpu_count(), len(tx_by_chrom))
    # Cap at chromosome count to avoid wasted serialization
    max_workers = min(max_workers, len(tx_by_chrom))
    logger.info(f"Using {max_workers} parallel workers for {len(tx_by_chrom)} chromosomes")
    
    all_consensus_transcripts = []
    metrics = initialize_metrics()
    
    # Prepare chromosome tasks
    chrom_tasks = []
    for chrom_num, chrom in enumerate(sorted(tx_by_chrom.keys()), 1):
        chrom_tx_ids = tx_by_chrom[chrom]
        chrom_tx_set = set(chrom_tx_ids)
        
        # Filter dataframes to this chromosome
        chrom_support_df = support_df[support_df['AlignmentId'].isin(chrom_tx_set)].copy()
        chrom_mrna_df = mrna_metrics_df[mrna_metrics_df['AlignmentId'].isin(chrom_tx_set)].copy() if len(mrna_metrics_df) > 0 else pd.DataFrame()
        chrom_cds_df = cds_metrics_df[cds_metrics_df['AlignmentId'].isin(chrom_tx_set)].copy() if len(cds_metrics_df) > 0 else pd.DataFrame()
        chrom_eval_df = eval_df[eval_df['AlignmentId'].isin(chrom_tx_set)].copy() if len(eval_df) > 0 else pd.DataFrame()
        
        # OPTIMIZATION: Filter large dataframes to reduce serialization overhead
        # tm_eval_df and ref_df are passed to ALL workers - filter them per chromosome
        chrom_tm_eval_df = tm_eval_df  # Keep full for now (small enough)
        # ref_df is large - but needed for all gene lookups, so keep full
        # Alternative: could filter by GeneId if we knew which genes are on this chromosome
        
        # OPTIMIZATION: Filter tx_dict to this chromosome only (HUGE savings!)
        chrom_tx_dict = {tx_id: tx_dict[tx_id] for tx_id in chrom_tx_ids if tx_id in tx_dict}
        
        # OPTIMIZATION: Filter alignment_source_map to this chromosome
        chrom_alignment_source_map = {tx_id: alignment_source_map[tx_id] for tx_id in chrom_tx_ids if tx_id in alignment_source_map}
        
        # OPTIMIZATION: Filter gp_attrs biotypes to this chromosome
        chrom_gp_attrs_transcript_biotypes = {tx_id: gp_attrs_transcript_biotypes[tx_id] for tx_id in chrom_tx_ids if tx_id in gp_attrs_transcript_biotypes} if gp_attrs_transcript_biotypes else {}
        chrom_gp_attrs_gene_biotypes = {tx_id: gp_attrs_gene_biotypes[tx_id] for tx_id in chrom_tx_ids if tx_id in gp_attrs_gene_biotypes} if gp_attrs_gene_biotypes else {}
        
        chrom_tasks.append((
            chrom, chrom_num, len(tx_by_chrom), chrom_tx_ids, chrom_tx_dict, chrom_support_df,
            chrom_mrna_df, chrom_cds_df, chrom_eval_df, chrom_tm_eval_df,
            ref_df, chrom_alignment_source_map, args, readthrough_gene_set, ref_gene_coords, genes_with_overlaps_in_ref, chrom_gp_attrs_transcript_biotypes, chrom_gp_attrs_gene_biotypes
        ))
    
    # Process chromosomes in parallel
    if max_workers > 1:
        logger.info("Processing chromosomes in parallel...")
        with mp.Pool(processes=max_workers) as pool:
            results = pool.starmap(process_chromosome_wrapper, chrom_tasks)
    else:
        logger.info("Processing chromosomes sequentially...")
        results = [process_chromosome_wrapper(*task) for task in chrom_tasks]
    
    # Aggregate results
    for chrom_consensus, chrom_metrics in results:
        all_consensus_transcripts.extend(chrom_consensus)
        merge_metrics(metrics, chrom_metrics)
    
    logger.info(f"\n✓ Total consensus transcripts across all chromosomes: {len(all_consensus_transcripts)}")
    elapsed = time.time() - start_time
    logger.info(f"  Chromosome processing took {elapsed:.1f} seconds")
    
    # Convert to consensus dict format
    logger.info("\n" + "="*80)
    logger.info("STEP 5: Final Filtering and Cleanup")
    logger.info("="*80)
    logger.info("Converting to consensus dictionary format...")
    consensus_dict = {}
    for aln_id, attrs in all_consensus_transcripts:
        consensus_dict[aln_id] = attrs
    logger.info(f"✓ Consensus dict has {len(consensus_dict)} entries")
    
    # Deduplication is now done per-chromosome (see process_chromosome Step 9)
    # This prevents genes on different chromosomes from being treated as duplicates
    logger.info("\nSkipping global deduplication (already done per-chromosome)...")
    deduplicated_consensus = consensus_dict
    if metrics['Duplicate transcripts']:
        logger.info(f"  Total duplicates removed across all chromosomes: {dict(metrics['Duplicate transcripts'])}")
    
    args.ref_pc_ensg = {norm_ensg(g) for g, b in gene_biotype_map.items() if b == 'protein_coding'}
    final_consensus = finalize_consensus_after_source_gene_resolution(
        consensus_dict,
        tx_dict,
        metrics,
        args,
        readthrough_gene_set,
        ref_gene_coords,
        genes_with_overlaps_in_ref,
        run_resolve_overlapping_different_genes=False,
    )
    final_consensus = apply_reference_gene_biotype_policy(
        final_consensus, gene_biotype_map, metrics=metrics
    )
    final_consensus = rescue_missing_reference_pc_genes(
        final_consensus,
        tx_dict,
        alignment_source_map,
        ref_gene_coords,
        gene_biotype_map,
        ref_df,
        args,
        mrna_metrics_df=mrna_metrics_df,
        metrics=metrics,
    )
    final_consensus = rescue_missing_reference_noncoding_genes(
        final_consensus,
        tx_dict,
        alignment_source_map,
        ref_gene_coords,
        gene_biotype_map,
        ref_df,
        args,
        mrna_metrics_df=mrna_metrics_df,
        metrics=metrics,
    )
    final_consensus = rescue_missing_reference_transcripts(
        final_consensus,
        tx_dict,
        alignment_source_map,
        ref_df,
        gene_biotype_map,
        args,
        mrna_metrics_df=mrna_metrics_df,
        metrics=metrics,
    )
    final_consensus = rescue_alternative_source_isoforms(
        final_consensus,
        tx_dict,
        alignment_source_map,
        ref_df,
        args,
        mrna_metrics_df=mrna_metrics_df,
        metrics=metrics,
    )

    # Guarantee retention + clean labeling of lineage-specific genes found ONLY in
    # protein/augMP evidence (no reference transcript). Off by default; enabled via
    # --keep-protein-only-novel (which the high_recall preset turns on). Runs after
    # the reference rescues so the intergenic test sees the full reference footprint,
    # and before the gene/PC counting + completeness below so survivors are counted
    # as novel protein-coding genes in the final stats.
    if getattr(args, 'keep_protein_only_novel', False):
        before = len(final_consensus)
        final_consensus = reclassify_protein_only_novel(
            final_consensus, tx_dict, ref_df, metrics, args.genome,
            mrna_metrics_df=mrna_metrics_df,
            min_coverage=getattr(args, 'protein_novel_min_coverage', 0.0),
            min_identity=getattr(args, 'protein_novel_min_identity', 0.0),
            intergenic_only=not getattr(args, 'protein_novel_keep_overlapping', False),
            min_exons=getattr(args, 'protein_novel_min_exons', 2),
            min_cds_aa=getattr(args, 'protein_novel_min_cds_aa', 100),
        )
        pn = metrics.get('Protein-only novel', {})
        logger.info(
            f"  Protein-only novel genes: kept {pn.get('kept', 0)} as putative_novel, "
            f"dropped {pn.get('dropped', 0)} (single_exon={pn.get('dropped_single_exon', 0)}, "
            f"short_cds={pn.get('dropped_short_cds', 0)}, low_cov={pn.get('dropped_low_coverage', 0)}, "
            f"low_id={pn.get('dropped_low_identity', 0)}, overlaps_ref={pn.get('dropped_overlaps_reference', 0)}); "
            f"consensus {before} -> {len(final_consensus)}"
        )

    # Recover high-quality augMP models that consensus selection dropped, at loci
    # where the consensus is otherwise empty. augMP is often the only mode that
    # finds lineage-specific genes (species-specific proteins), so this closes a
    # measured recall gap without touching occupied loci. Runs after the novel
    # path so its occupied footprint already includes the rescued novel genes.
    if getattr(args, 'rescue_dropped_augMP', False):
        before = len(final_consensus)
        final_consensus = rescue_augMP_at_empty_loci(
            final_consensus, tx_dict, ref_df, metrics, args.genome,
            mrna_metrics_df=mrna_metrics_df,
            min_exons=getattr(args, 'rescue_augMP_min_exons', 2),
            min_cds_aa=getattr(args, 'rescue_augMP_min_cds_aa', 100),
            min_coverage=getattr(args, 'rescue_augMP_min_coverage', 0.0),
            min_identity=getattr(args, 'rescue_augMP_min_identity', 0.0),
            single_exon_min_cds_aa=getattr(args, 'rescue_augMP_single_exon_min_cds_aa', 300),
        )
        am = metrics.get('augMP empty-locus recovery', {})
        logger.info(
            f"  augMP/augPB CDS recovery (PC-empty loci): added {am.get('recovered_genes', 0)} genes "
            f"({am.get('recovered_orthologs', 0)} recovered orthologs, "
            f"{am.get('recovered_novel', 0)} protein-only novel); dropped "
            f"single_exon={am.get('dropped_single_exon', 0)}, short_cds={am.get('dropped_short_cds', 0)}, "
            f"locus_has_PC_CDS={am.get('dropped_locus_occupied', 0)}; consensus {before} -> {len(final_consensus)}"
        )

    # Recover real protein-coding genes that inherited a non-coding biotype from a
    # non-coding reference ortholog but are expressed and ORF-intact in this
    # species (recall gap confirmed vs RefSeq). Strict + expression-gated, so it
    # does not re-introduce intergenic false positives. Runs before counting so
    # rescued genes are reflected in the PC stats.
    if getattr(args, 'rescue_expressed_noncoding_to_pc', False):
        final_consensus = rescue_expressed_noncoding_to_pc(
            final_consensus, tx_dict, metrics,
            min_cds_aa=getattr(args, 'rescue_expressed_min_cds_aa', 100),
            require_multiexon=not getattr(args, 'rescue_expressed_allow_single_exon', False),
            require_protein_evidence_for_noncoding=getattr(
                args, 'rescue_noncoding_require_protein_evidence', True),
        )
        r = metrics.get('Expressed noncoding->PC rescue', {})
        logger.info(
            f"  Expressed non-coding -> protein_coding rescue: promoted {r.get('genes', 0)} genes "
            f"({r.get('transcripts', 0)} transcripts) with intact ORF + IsoSeq/RNA support "
            f"(gated for lacking protein evidence: {r.get('gated_no_protein_evidence', 0)})"
        )

    # Count genes and protein coding genes
    unique_genes = set()
    unique_pc_genes = set()
    unique_pc_genes_by_source = collections.defaultdict(set)
    
    for aln_id, attrs in final_consensus:
        gene = attrs.get('source_gene')
        mode = alignment_source_map.get(aln_id, 'unknown')
        if gene and gene != 'N/A':
            unique_genes.add(gene)
            if attrs.get('gene_biotype') == 'protein_coding':
                unique_pc_genes.add(gene)
                unique_pc_genes_by_source[mode].add(gene)
    
    logger.info(f"  Total unique genes: {len(unique_genes)}")
    logger.info(f"  Protein-coding genes: {len(unique_pc_genes)}")
    logger.info(f"  Protein-coding genes by source: {dict((k, len(v)) for k, v in unique_pc_genes_by_source.items())}")
    
    # Compare input to output for each source
    logger.info("\n  Protein-coding gene comparison (input vs output):")
    for mode in pc_gene_ids_by_source.keys():
        input_genes = pc_gene_ids_by_source[mode]
        output_genes = unique_pc_genes_by_source.get(mode, set())
        missing_genes = input_genes - output_genes
        logger.info(f"    {mode}: {len(input_genes)} input → {len(output_genes)} output")
        if missing_genes:
            logger.info(f"      Missing {len(missing_genes)} genes: {sorted(list(missing_genes))[:10]}{'...' if len(missing_genes) > 10 else ''}")

    def _tx_on_audit_chrom(tx_id: str) -> bool:
        """Whether a transcript should be included in the optional PC audit."""
        if not args.pc_audit_chrom:
            return True
        tx = tx_dict.get(tx_id)
        if tx is None:
            return False
        return str(tx.chromosome) == str(args.pc_audit_chrom)

    # Optional: write a TSV to make "what was provided vs kept" easy to inspect.
    if args.pc_audit_tsv is not None:
        # Build output gene sets (any biotype and protein_coding), by mode, with optional chrom filter.
        out_any_by_mode = collections.defaultdict(set)
        out_pc_by_mode = collections.defaultdict(set)
        out_any_tx_by_mode = collections.defaultdict(list)
        out_pc_tx_by_mode = collections.defaultdict(list)

        for tx_id, attrs in final_consensus:
            if not _tx_on_audit_chrom(tx_id):
                continue
            mode = attrs.get("source") or alignment_source_map.get(tx_id, "unknown")
            gene = attrs.get("source_gene")
            if not gene and tx_id in tx_dict:
                gene = tx_dict[tx_id].name2
            if not gene:
                continue
            gene = re.sub(r"_cp\d+$", "", str(gene))
            out_any_by_mode[mode].add(gene)
            out_any_tx_by_mode[mode].append(tx_id)
            if attrs.get("gene_biotype") == "protein_coding":
                out_pc_by_mode[mode].add(gene)
                out_pc_tx_by_mode[mode].append(tx_id)

        # Recompute input genes by mode with optional chrom filter (using tx.name2 + biotype lookups).
        in_any_by_mode = {}
        in_pc_by_mode = {}
        in_any_tx_by_mode = {}
        in_pc_tx_by_mode = {}

        for mode, tx_ids in genes_by_source.items():
            mode_any_genes = set()
            mode_pc_genes = set()
            mode_any_txs = []
            mode_pc_txs = []
            for tx_id in tx_ids:
                if not _tx_on_audit_chrom(tx_id):
                    continue
                tx_obj = tx_dict.get(tx_id)
                gene_id = tx_obj.name2 if tx_obj else None
                if not gene_id:
                    continue
                gene_id_norm = re.sub(r"_cp\d+$", "", str(gene_id))
                mode_any_genes.add(gene_id_norm)
                mode_any_txs.append(tx_id)

                # Use gp_attrs biotype if present (best for txTM); else fall back to ref gene_biotype_map.
                tx_id_for_lookup = re.sub(r"_cp\d+$", "", tx_id)
                is_pc = gp_attrs_gene_biotypes.get(tx_id_for_lookup) == "protein_coding"
                if not is_pc:
                    is_pc = gene_biotype_map.get(gene_id_norm) == "protein_coding"
                if is_pc:
                    mode_pc_genes.add(gene_id_norm)
                    mode_pc_txs.append(tx_id)

            in_any_by_mode[mode] = mode_any_genes
            in_pc_by_mode[mode] = mode_pc_genes
            in_any_tx_by_mode[mode] = mode_any_txs
            in_pc_tx_by_mode[mode] = mode_pc_txs

        audit_rows = []
        modes = sorted(set(list(genes_by_source.keys()) + list(out_any_by_mode.keys())))
        for mode in modes:
            all_genes = sorted(in_any_by_mode.get(mode, set()) | out_any_by_mode.get(mode, set()))
            for gene in all_genes:
                audit_rows.append(
                    {
                        "chrom_filter": args.pc_audit_chrom or "",
                        "mode": mode,
                        "gene_name": gene,
                        "in_input_any": "yes" if gene in in_any_by_mode.get(mode, set()) else "no",
                        "in_input_pc": "yes" if gene in in_pc_by_mode.get(mode, set()) else "no",
                        "in_output_any": "yes" if gene in out_any_by_mode.get(mode, set()) else "no",
                        "in_output_pc": "yes" if gene in out_pc_by_mode.get(mode, set()) else "no",
                        "n_input_tx_any_mode": len(in_any_tx_by_mode.get(mode, [])),
                        "n_input_tx_pc_mode": len(in_pc_tx_by_mode.get(mode, [])),
                        "n_output_tx_any_mode": len(out_any_tx_by_mode.get(mode, [])),
                        "n_output_tx_pc_mode": len(out_pc_tx_by_mode.get(mode, [])),
                    }
                )

        args.pc_audit_tsv.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(audit_rows).to_csv(args.pc_audit_tsv, sep="\t", index=False)
        logger.info(f"\n✓ Wrote protein-coding audit TSV to {args.pc_audit_tsv}")
    
    # Calculate final metrics
    logger.info("\n" + "="*80)
    logger.info("STEP 6: Calculating Metrics")
    logger.info("="*80)
    calculate_completeness(final_consensus, metrics, gene_biotype_map, transcript_biotype_map)
    
    # Log summary metrics
    if 'Multi-locus mappings' in metrics and metrics['Multi-locus mappings']:
        logger.info("\nMulti-locus mapping summary:")
        for mode, count in metrics['Multi-locus mappings'].items():
            kept = metrics['Multi-locus kept'].get(mode, 0)
            logger.info(f"  {mode}: {count} transcripts map to multiple loci, kept {kept} total copies")
    
    if 'AugPB Classes' in metrics and metrics['AugPB Classes']:
        logger.info("\nAugPB classification summary:")
        for cls, count in metrics['AugPB Classes'].items():
            logger.info(f"  {cls}: {count}")
    
    # Write outputs
    logger.info("\n" + "="*80)
    logger.info("STEP 7: Writing Output Files")
    logger.info("="*80)
    logger.info(f"Writing consensus genePred to {args.consensus_gp}")
    consensus_gene_dict = write_consensus_gps(args.consensus_gp, args.consensus_gp_info,
                                              final_consensus, tx_dict, args.genome)
    logger.info(f"✓ Wrote {len(final_consensus)} transcripts to genePred")

    if getattr(args, 'consensus_postprocess', True) and getattr(args, 'ref_gp', None):
        logger.info("\nApplying consensus postprocess (split runaway pc genes + "
                    "reclassify/drop weak duplicate copies)...")
        try:
            from cat.consensus_postprocess import apply_postprocess
            report_path = args.consensus_gp_info + '.postprocess_report.tsv'
            stats = apply_postprocess(
                consensus_gene_dict, args.ref_gp,
                consensus_gp=args.consensus_gp,
                consensus_gp_info=args.consensus_gp_info,
                genome=args.genome,
                report_path=report_path,
                min_introns_for_low_support=getattr(args, 'postprocess_min_introns_low_support', 3),
                augpb_chimera_exon_ratio=getattr(args, 'postprocess_augpb_chimera_exon_ratio', 1.5),
                low_support_fraction=getattr(args, 'postprocess_low_support_fraction', 0.3),
                protect_strong_modes=getattr(args, 'postprocess_protect_strong_modes', True),
            )
            logger.info(f"✓ Postprocess: split={stats['split']}, "
                        f"reclassify={stats['reclassify']}, drop={stats['drop']}")
        except Exception as exc:
            logger.warning(f"Postprocess failed (continuing with raw consensus): {exc}")

    logger.info(f"\nWriting GFF3 to {args.consensus_gff3}")
    write_consensus_gff3(consensus_gene_dict, args.consensus_gff3)
    logger.info("✓ Wrote GFF3 file")
    
    logger.info(f"\nWriting FASTA files...")
    logger.info(f"  Transcript FASTA: {args.consensus_fasta}")
    logger.info(f"  Protein FASTA: {args.protein_fasta}")
    write_consensus_fastas(consensus_gene_dict, args.consensus_fasta, args.protein_fasta, args.fasta)
    logger.info("✓ Wrote FASTA files")
    
    total_time = time.time() - start_time
    logger.info(f"\n✓ Total consensus generation time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    
    return metrics


def load_alt_gene_assignments(db_path, denovo_tx_modes):
    """Load alternative gene assignments (AssignedGeneId) for augPB transcripts"""
    session = tools.sqlInterface.start_session(db_path)
    r = []
    for tx_mode in denovo_tx_modes:
        try:
            table = tools.sqlInterface.tables['alt_names'][tx_mode]
            inspector = inspect(session.bind)
            if inspector.has_table(table.__tablename__):
                alt_data = tools.sqlInterface.load_alternatives(table, session)
                if not alt_data.empty:
                    r.append(alt_data)
        except Exception as e:
            logger.warning(f"Error loading alternative names for {tx_mode}: {e}")
    
    session.close()
    
    if not r:
        # Return empty DataFrame with correct columns
        empty_df = pd.DataFrame(columns=['TranscriptId', 'AssignedGeneId', 'AlternativeGeneIds', 'ResolutionMethod'])
        empty_df.columns = [x if x != 'TranscriptId' else 'AlignmentId' for x in empty_df.columns]
        return empty_df
    
    df = pd.concat(r, ignore_index=True)
    # Rename TranscriptId to AlignmentId for merging
    df.columns = [x if x != 'TranscriptId' else 'AlignmentId' for x in df.columns]
    return df


def _extract_junctions_one_bam(bam_path, chromosomes=None):
    """Extract splice junctions from a single BAM file (worker for parallel extraction).

    Returns a set of (chrom, junction_start, junction_end) tuples, or an empty
    set if the BAM is missing/unindexed/unreadable.
    """
    import pysam

    junctions = set()
    if not os.path.exists(bam_path):
        logger.warning(f"  BAM file not found: {bam_path}")
        return junctions

    if not os.path.exists(bam_path + '.bai') and not os.path.exists(bam_path + '.csi'):
        logger.warning(f"  BAM index not found for {bam_path}, skipping")
        return junctions

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        logger.warning(f"  Could not open BAM {bam_path}: {e}")
        return junctions

    try:
        bam_chroms = set(bam.references)
        fetch_chroms = chromosomes if chromosomes else bam_chroms

        for chrom in fetch_chroms:
            if chrom not in bam_chroms:
                continue
            for read in bam.fetch(chrom):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.cigartuples is None:
                    continue

                ref_pos = read.reference_start
                for op, length in read.cigartuples:
                    if op == 0:    # M (match/mismatch)
                        ref_pos += length
                    elif op == 2:  # D (deletion)
                        ref_pos += length
                    elif op == 3:  # N (splice junction)
                        junctions.add((chrom, ref_pos, ref_pos + length))
                        ref_pos += length
                    elif op == 1:  # I (insertion)
                        pass
                    elif op == 4:  # S (soft clip)
                        pass
                    elif op == 5:  # H (hard clip)
                        pass
                    elif op == 7:  # = (sequence match)
                        ref_pos += length
                    elif op == 8:  # X (sequence mismatch)
                        ref_pos += length
    finally:
        bam.close()

    return junctions


def extract_splice_junctions_from_bams(bam_paths, chromosomes=None, max_workers=None):
    """Extract all splice junctions from BAM files.

    Returns a set of (chrom, junction_start, junction_end) tuples where
    junction_start is the last exonic base before the intron and
    junction_end is the first exonic base after the intron (0-based half-open).

    Work is parallelized across (BAM, chromosome) pairs — not just per BAM — so a
    handful of large IsoSeq BAMs can still saturate all available cores instead of
    leaving most idle. Each (BAM, chromosome) fetch uses the BAM index, so opening
    the same file concurrently is cheap. Set max_workers to bound the pool
    (default: min(#tasks, available CPUs)).
    """
    import multiprocessing as mp

    valid_bams = [p for p in bam_paths if p]
    if not valid_bams:
        return set()

    # Build fine-grained (bam, [chrom]) tasks when we know the chromosome list, so
    # #tasks = #bams * #chroms can use far more than #bams cores. Without a
    # chromosome list, fall back to one whole-file task per BAM.
    if chromosomes:
        chrom_list = sorted(chromosomes)
        tasks = [(bam_path, [c]) for bam_path in valid_bams for c in chrom_list]
    else:
        tasks = [(bam_path, None) for bam_path in valid_bams]

    # Serial fast-path for a single task (avoids pool overhead).
    if len(tasks) == 1:
        bam_path, chl = tasks[0]
        return _extract_junctions_one_bam(bam_path, chl)

    if max_workers is None:
        try:
            avail = len(os.sched_getaffinity(0))
        except AttributeError:
            avail = mp.cpu_count()
        max_workers = min(len(tasks), max(1, avail))

    junctions = set()
    if max_workers <= 1:
        for bam_path, chl in tasks:
            junctions |= _extract_junctions_one_bam(bam_path, chl)
        return junctions

    # Use 'fork' (the default). 'spawn' is NOT usable on these compute nodes:
    # the workers fail to re-open POSIX semaphores (SemLock._rebuild ->
    # FileNotFoundError) because of the restricted /dev/shm there. fork inherits
    # the semaphore fds and works. The former fork memory blow-up came from the
    # per-worker read buffers, which the streaming workers no longer allocate, so
    # each worker now touches almost nothing of the parent heap. gc.freeze() moves
    # the large parent objects out of the garbage collector's reach so a GC pass
    # in a forked worker can't dirty shared copy-on-write pages.
    import gc
    gc.freeze()
    with mp.Pool(processes=max_workers) as pool:
        for result in pool.starmap(_extract_junctions_one_bam, tasks):
            junctions |= result
    return junctions


def _coverage_intervals_one_bam_chrom(bam_path, chrom):
    """Worker: merged covered intervals (accounting for splicing) for one
    (BAM, chromosome). Returns (chrom, merged_starts, merged_ends).

    Only the (small) merged intervals are returned to the parent — never the raw
    reads — so pickling stays cheap even for deep IsoSeq BAMs.
    """
    import pysam

    if not os.path.exists(bam_path):
        return (chrom, [], [])
    if not os.path.exists(bam_path + '.bai') and not os.path.exists(bam_path + '.csi'):
        return (chrom, [], [])
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception:
        return (chrom, [], [])

    try:
        if chrom not in set(bam.references):
            return (chrom, [], [])

        # Memory-bounded streaming merge. The previous version buffered EVERY
        # aligned block of the whole chromosome (tens of millions of tuples for
        # deep IsoSeq) before sorting; with several workers running concurrently
        # that pinned the job at its memory ceiling and spent all its time in
        # cgroup reclaim. Instead we exploit that fetch() yields reads sorted by
        # reference_start: once we reach a read starting at frontier P, no future
        # block can fall below P, so any merged interval ending before P is final
        # and can be emitted and dropped. The 'active' set therefore stays tiny
        # (≈1 interval for contiguous coverage), independent of read depth.
        merged_starts = []
        merged_ends = []
        active = []  # sorted-by-start, disjoint [s, e]; every e >= current frontier

        def _merge_block(b_s, b_e):
            n = len(active)
            lo = 0
            while lo < n and active[lo][1] < b_s:  # strictly-before intervals: keep
                lo += 1
            hi = lo
            while hi < n and active[hi][0] <= b_e:  # overlapping/adjacent: absorb
                hi += 1
            if lo == hi:
                active.insert(lo, [b_s, b_e])
            else:
                ns = active[lo][0] if active[lo][0] < b_s else b_s
                ne = active[hi - 1][1] if active[hi - 1][1] > b_e else b_e
                active[lo:hi] = [[ns, ne]]

        for read in bam.fetch(chrom):
            if read.is_unmapped:
                continue
            frontier = read.reference_start
            k = 0
            while k < len(active) and active[k][1] < frontier:
                merged_starts.append(active[k][0])
                merged_ends.append(active[k][1])
                k += 1
            if k:
                del active[:k]
            for b_s, b_e in read.get_blocks():
                _merge_block(b_s, b_e)

        for s, e in active:
            merged_starts.append(s)
            merged_ends.append(e)

        if not merged_starts:
            return (chrom, [], [])
        return (chrom, merged_starts, merged_ends)
    finally:
        bam.close()


def compute_exon_coverage_from_bams(bam_paths, exon_intervals, min_reads=1, max_workers=None):
    """Check which exons have RNA-seq/IsoSeq read coverage.

    Parallelizes across (BAM, chromosome) pairs: each worker streams one
    chromosome of one BAM and returns its merged covered intervals. The parent
    unions the intervals per chromosome and binary-searches the exons against
    them. This replaces the old single-core serial pass, which was the dominant
    cost when many large IsoSeq BAMs are present.

    Args:
        bam_paths: list of BAM file paths
        exon_intervals: list of (chrom, start, end) tuples (0-based half-open)
        min_reads: kept for API compatibility; the streaming approach supports
                   min_reads=1 (any coverage marks the exon supported)
        max_workers: cap on the process pool (default: min(#tasks, CPUs))

    Returns:
        list of 0/1 values, one per exon
    """
    import bisect

    support = [0] * len(exon_intervals)
    if not exon_intervals:
        return support

    # Group exons by chromosome for one merged coverage pass per chromosome.
    chrom_exons = collections.defaultdict(list)
    for i, (chrom, start, end) in enumerate(exon_intervals):
        chrom_exons[chrom].append((start, end, i))
    for chrom in chrom_exons:
        chrom_exons[chrom].sort()

    valid_bams = [p for p in bam_paths if p and os.path.exists(p)]
    if not valid_bams:
        return support

    tasks = [(bam_path, chrom) for bam_path in valid_bams for chrom in chrom_exons.keys()]
    if not tasks:
        return support

    if max_workers is None:
        try:
            avail = len(os.sched_getaffinity(0))
        except AttributeError:
            import multiprocessing as _mp
            avail = _mp.cpu_count()
        # The streaming worker holds only its (tiny) running merged intervals, so
        # memory is no longer the constraint; use all available cores like the
        # junction pass, with a sane upper bound on concurrent BAM readers.
        max_workers = max(1, min(len(tasks), avail, 32))

    # chrom -> list of (starts, ends) contributed by each (bam, chrom) scan
    per_chrom = collections.defaultdict(list)
    if max_workers <= 1 or len(tasks) == 1:
        for bam_path, chrom in tasks:
            c, s, e = _coverage_intervals_one_bam_chrom(bam_path, chrom)
            if s:
                per_chrom[c].append((s, e))
    else:
        import multiprocessing as mp
        import gc
        # fork (default): 'spawn' can't start workers on these nodes (semaphore
        # re-open fails under the restricted /dev/shm). The streaming worker keeps
        # only its running merged intervals, so fork's copy-on-write stays cheap.
        # gc.freeze() keeps the collector from dirtying the parent's shared pages
        # inside the forked workers.
        gc.freeze()
        with mp.Pool(processes=max_workers) as pool:
            for c, s, e in pool.starmap(_coverage_intervals_one_bam_chrom, tasks):
                if s:
                    per_chrom[c].append((s, e))

    # Merge intervals across BAMs per chromosome, then binary-search exons.
    for chrom, exons in chrom_exons.items():
        interval_lists = per_chrom.get(chrom)
        if not interval_lists:
            continue

        all_intervals = []
        for starts, ends in interval_lists:
            all_intervals.extend(zip(starts, ends))
        all_intervals.sort()

        merged_starts = []
        merged_ends = []
        cur_s, cur_e = all_intervals[0]
        for s, e in all_intervals[1:]:
            if s <= cur_e:
                if e > cur_e:
                    cur_e = e
            else:
                merged_starts.append(cur_s)
                merged_ends.append(cur_e)
                cur_s, cur_e = s, e
        merged_starts.append(cur_s)
        merged_ends.append(cur_e)

        for ex_start, ex_end, orig_idx in exons:
            # Rightmost interval that starts strictly before ex_end.
            pos = bisect.bisect_left(merged_starts, ex_end) - 1
            if pos >= 0 and merged_ends[pos] > ex_start:
                support[orig_idx] = 1

    return support


def build_reference_intron_exon_sets(ref_gp_path):
    """Build sets of intron and exon intervals from the reference annotation genePred.
    
    Returns:
        ref_introns: set of (chrom, start, end) for each intron in reference
        ref_exons: dict of chrom -> sorted list of (start, end) for interval tree lookups
    """
    ref_introns = set()
    ref_exons = collections.defaultdict(list)
    
    with open(ref_gp_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 10:
                continue
            
            try:
                tx = tools.transcripts.GenePredTranscript(parts)
            except Exception:
                continue
            
            chrom = tx.chromosome
            for intron in tx.intron_intervals:
                ref_introns.add((chrom, intron.start, intron.stop))
            
            for exon in tx.exon_intervals:
                ref_exons[chrom].append((exon.start, exon.stop))
    
    # Sort exon lists for binary search
    for chrom in ref_exons:
        ref_exons[chrom].sort()
    
    return ref_introns, dict(ref_exons)


def exon_overlaps_reference(chrom, start, end, ref_exons, min_overlap_frac=0.1):
    """Check if an exon overlaps any reference annotation exon."""
    import bisect
    
    chrom_exons = ref_exons.get(chrom, [])
    if not chrom_exons:
        return False
    
    # Binary search for candidate overlapping exons
    idx = bisect.bisect_left(chrom_exons, (start,))
    # Check a window around the insertion point
    for i in range(max(0, idx - 1), min(len(chrom_exons), idx + 50)):
        ref_start, ref_end = chrom_exons[i]
        if ref_start >= end:
            break
        if ref_end <= start:
            continue
        # Overlap found
        overlap = min(end, ref_end) - max(start, ref_start)
        exon_len = end - start
        if exon_len > 0 and overlap / exon_len >= min_overlap_frac:
            return True
    
    return False


def compute_real_support(tx_dict, support_df, bam_files, isoseq_bam_files, ref_gp_path, alignment_source_map):
    """Compute real RNA-seq and annotation support from BAM files and reference.
    
    For annotation support, we use the non-denovo transcripts already in tx_dict
    (transMap, txTM, augTM, etc.) as the projected reference annotation on the
    target genome. The ref_gp_path (source genome coords) would not match target
    genome chromosomes.
    
    Updates support_df in-place with computed support vectors and percentages.
    """
    all_bams = list(bam_files or []) + list(isoseq_bam_files or [])
    has_bams = any(os.path.exists(b) for b in all_bams)
    
    if not has_bams:
        logger.info("  No BAM files provided — keeping default support vectors")
        return
    
    # Step 1: Extract splice junctions from BAMs
    splice_junctions = set()
    logger.info(f"  Extracting splice junctions from {len(all_bams)} BAM file(s)...")
    extract_start = time.time()
    
    chromosomes = set()
    for aln_id in support_df['AlignmentId'].values:
        tx = tx_dict.get(aln_id)
        if tx:
            chromosomes.add(tx.chromosome)
    
    splice_junctions = extract_splice_junctions_from_bams(all_bams, chromosomes)
    logger.info(f"    ✓ Found {len(splice_junctions)} unique splice junctions in {time.time() - extract_start:.1f}s")
    
    # Step 2: Build annotation intron/exon sets from non-denovo transcripts in tx_dict
    # These are transMap/txTM/augTM projections of the reference onto the target genome
    ref_introns = set()
    ref_exons = collections.defaultdict(list)
    annot_start = time.time()
    non_denovo_count = 0
    
    for aln_id, tx in tx_dict.items():
        if _is_denovo(aln_id):
            continue
        non_denovo_count += 1
        chrom = tx.chromosome
        for intron in tx.intron_intervals:
            ref_introns.add((chrom, intron.start, intron.stop))
        for exon in tx.exon_intervals:
            ref_exons[chrom].append((exon.start, exon.stop))
    
    # Sort and deduplicate exon lists for binary search
    for chrom in ref_exons:
        ref_exons[chrom] = sorted(set(ref_exons[chrom]))
    ref_exons = dict(ref_exons)
    
    logger.info(f"    ✓ Built annotation index from {non_denovo_count} non-denovo transcripts: "
                 f"{len(ref_introns)} unique introns, "
                 f"{sum(len(v) for v in ref_exons.values())} unique exons in {time.time() - annot_start:.1f}s")
    
    # Step 3: Compute per-transcript support
    logger.info("  Computing per-transcript support vectors...")
    compute_start = time.time()
    
    aln_ids = support_df['AlignmentId'].values
    new_intron_rna = []
    new_exon_rna = []
    new_intron_annot = []
    new_exon_annot = []
    new_cds_annot = []
    
    # Collect exon intervals for batch coverage check
    exon_intervals_for_coverage = []
    exon_interval_map = []  # (transcript_idx, exon_idx_within_transcript)
    
    for tx_idx, aln_id in enumerate(aln_ids):
        tx = tx_dict.get(aln_id)
        if tx is None:
            new_intron_rna.append([])
            new_exon_rna.append([])
            new_intron_annot.append([])
            new_exon_annot.append([])
            new_cds_annot.append([])
            continue
        
        chrom = tx.chromosome
        is_denovo = _is_denovo(aln_id)
        
        # Intron RNA support: check splice junctions (all transcripts, including denovo)
        intron_rna_vec = []
        for intron in tx.intron_intervals:
            if splice_junctions and (chrom, intron.start, intron.stop) in splice_junctions:
                intron_rna_vec.append(1)
            else:
                intron_rna_vec.append(0)
        new_intron_rna.append(intron_rna_vec)
        
        # Intron annotation support: check reference introns
        intron_annot_vec = []
        for intron in tx.intron_intervals:
            if ref_introns and (chrom, intron.start, intron.stop) in ref_introns:
                intron_annot_vec.append(1)
            else:
                intron_annot_vec.append(0)
        new_intron_annot.append(intron_annot_vec)
        
        # Exon annotation support: check reference exon overlap
        exon_annot_vec = []
        for exon in tx.exon_intervals:
            if ref_exons and exon_overlaps_reference(chrom, exon.start, exon.stop, ref_exons):
                exon_annot_vec.append(1)
            else:
                exon_annot_vec.append(0)
        new_exon_annot.append(exon_annot_vec)
        
        # CDS annotation support (same as exon for now)
        new_cds_annot.append(list(exon_annot_vec))
        
        # Exon RNA support: collect intervals for batch coverage query (all transcripts)
        exon_rna_placeholder = []
        for exon_idx, exon in enumerate(tx.exon_intervals):
            exon_intervals_for_coverage.append((chrom, exon.start, exon.stop))
            exon_interval_map.append((tx_idx, exon_idx))
            exon_rna_placeholder.append(0)
        new_exon_rna.append(exon_rna_placeholder)
    
    # Step 4: Batch compute exon coverage from BAMs
    if has_bams and exon_intervals_for_coverage:
        logger.info(f"    Computing exon coverage for {len(exon_intervals_for_coverage)} exons...")
        coverage_results = compute_exon_coverage_from_bams(all_bams, exon_intervals_for_coverage)
        
        for result_idx, (tx_idx, exon_idx) in enumerate(exon_interval_map):
            new_exon_rna[tx_idx][exon_idx] = coverage_results[result_idx]
    
    # Step 5: Update support_df
    support_df['IntronRnaSupport'] = new_intron_rna
    support_df['AllSpeciesIntronRnaSupport'] = new_intron_rna
    support_df['ExonRnaSupport'] = new_exon_rna
    support_df['AllSpeciesExonRnaSupport'] = new_exon_rna
    support_df['IntronAnnotSupport'] = new_intron_annot
    support_df['ExonAnnotSupport'] = new_exon_annot
    support_df['CdsAnnotSupport'] = new_cds_annot
    
    # Recalculate percentages
    def calc_pct(vectors):
        return [100.0 * sum(1 for x in vec if x > 0) / len(vec) if len(vec) > 0 else 0.0
                for vec in vectors]
    
    support_df['IntronRnaSupportPercent'] = calc_pct(new_intron_rna)
    support_df['AllSpeciesIntronRnaSupportPercent'] = calc_pct(new_intron_rna)
    support_df['ExonRnaSupportPercent'] = calc_pct(new_exon_rna)
    support_df['AllSpeciesExonRnaSupportPercent'] = calc_pct(new_exon_rna)
    support_df['IntronAnnotSupportPercent'] = calc_pct(new_intron_annot)
    support_df['ExonAnnotSupportPercent'] = calc_pct(new_exon_annot)
    support_df['CdsAnnotSupportPercent'] = calc_pct(new_cds_annot)
    
    elapsed = time.time() - compute_start
    
    # Log summary statistics
    non_denovo_mask = ~support_df['IsAugPB']
    if non_denovo_mask.any():
        intron_pcts = support_df.loc[non_denovo_mask, 'IntronRnaSupportPercent']
        exon_pcts = support_df.loc[non_denovo_mask, 'ExonRnaSupportPercent']
        intron_annot_pcts = support_df.loc[non_denovo_mask, 'IntronAnnotSupportPercent']
        exon_annot_pcts = support_df.loc[non_denovo_mask, 'ExonAnnotSupportPercent']
        
        logger.info(f"    Non-denovo support stats ({non_denovo_mask.sum()} transcripts):")
        logger.info(f"      Intron RNA:  mean={intron_pcts.mean():.1f}%, "
                     f"non-zero={(intron_pcts > 0).sum()}/{len(intron_pcts)}")
        logger.info(f"      Exon RNA:    mean={exon_pcts.mean():.1f}%, "
                     f"non-zero={(exon_pcts > 0).sum()}/{len(exon_pcts)}")
        logger.info(f"      Intron Annot: mean={intron_annot_pcts.mean():.1f}%, "
                     f"non-zero={(intron_annot_pcts > 0).sum()}/{len(intron_annot_pcts)}")
        logger.info(f"      Exon Annot:  mean={exon_annot_pcts.mean():.1f}%, "
                     f"non-zero={(exon_annot_pcts > 0).sum()}/{len(exon_annot_pcts)}")
    
    logger.info(f"  ✓ Computed real support vectors in {elapsed:.1f}s")


def create_support_dataframe(tx_dict, db_path, ref_df, alignment_source_map, denovo_tx_modes=None,
                             bam_files=None, isoseq_bam_files=None, ref_gp_path=None):
    """Create support dataframe with RNA-seq and annotation support"""
    logger.info("  Loading alignment evaluation data...")
    
    # Load evaluation data
    tm_eval = tools.sqlInterface.load_alignment_evaluation(db_path)
    tm_filter_eval = tools.sqlInterface.load_filter_evaluation(db_path)
    
    # Build a mapping of original IDs to renamed IDs for cross-mode conflicts
    # (e.g., when load_gps renamed transcripts by adding "_cp", "_cp2", "_cp3", etc.)
    import re as re_module
    id_rename_map = {}
    for tx_id in tx_dict.keys():
        # Check if this is a renamed transcript with _cp, _cp2, _cp3, etc. suffix
        cp_match = re_module.search(r'_cp\d*$', tx_id)
        if cp_match:
            original_id = tx_id[:cp_match.start()]
            # Check if the original ID exists in tx_dict and is from a different source
            if original_id in tx_dict:
                # We have a conflict: original_id exists in tx_dict (likely from a different mode)
                # and tx_id is the renamed version. The database likely has data under original_id
                # that should be associated with the renamed version.
                # We need to determine which mode each belongs to
                original_mode = alignment_source_map.get(original_id)
                renamed_mode = alignment_source_map.get(tx_id)
                
                # If they're from different modes, the database entry for original_id from renamed_mode
                # should be updated to use tx_id
                if original_mode != renamed_mode:
                    # The database has entries under original_id for renamed_mode's data
                    # We need to map: (original_id, renamed_mode) -> tx_id
                    id_rename_map[(original_id, renamed_mode)] = tx_id
    
    # Update AlignmentIds in database DataFrames to match tx_dict's naming
    if id_rename_map:
        logger.info(f"  Updating AlignmentIds for {len(id_rename_map)} cross-mode conflicts...")
        
        # For tm_eval and tm_filter_eval, we need to map based on the mode
        # But we don't have a mode column in these DataFrames, so we need to infer it
        # from the AlignmentId itself using alignment_source_map
        
        def update_alignment_ids(df):
            """Update AlignmentIds in a DataFrame based on the rename map (vectorized)"""
            if 'AlignmentId' not in df.columns or len(df) == 0:
                return df
            
            # Vectorized approach: check all at once
            aln_ids = df['AlignmentId'].values
            valid_ids = set(tx_dict.keys())
            
            # Create a mapping for IDs that need _cp suffix (_cp, _cp2, _cp3, etc.)
            def fix_id(aln_id):
                if aln_id not in valid_ids:
                    # Try _cp, _cp2, _cp3, etc.
                    for suffix in ['_cp', '_cp2', '_cp3', '_cp4', '_cp5', '_cp6', '_cp7', '_cp8', '_cp9', '_cp10']:
                        candidate = f"{aln_id}{suffix}"
                        if candidate in valid_ids:
                            return candidate
                return aln_id
            
            new_ids = [fix_id(aln_id) for aln_id in aln_ids]
            updated_count = sum(1 for old, new in zip(aln_ids, new_ids) if old != new)
            
            df['AlignmentId'] = new_ids
            if updated_count > 0:
                logger.info(f"    Updated {updated_count} AlignmentIds to match tx_dict naming")
            return df
        
        tm_eval = update_alignment_ids(tm_eval)
        tm_filter_eval = update_alignment_ids(tm_filter_eval)
    
    # CRITICAL: Duplicate database rows for _cp copies
    # If tx_dict has ENST123_cp2 and ENST123_cp3, but database only has ENST123,
    # we need to duplicate the ENST123 row for each _cp copy
    def duplicate_for_cp_copies(df):
        """Duplicate rows for transcripts that have multiple _cp copies in tx_dict"""
        if 'AlignmentId' not in df.columns or len(df) == 0:
            return df
        
        rows_to_add = []
        valid_ids = set(tx_dict.keys())
        
        for _, row in df.iterrows():
            base_id = row['AlignmentId']
            # Check if there are _cp2, _cp3, etc. versions in tx_dict
            for suffix in ['_cp2', '_cp3', '_cp4', '_cp5']:
                cp_id = f"{base_id}{suffix}"
                if cp_id in valid_ids:
                    # Create a copy of this row with the new AlignmentId
                    new_row = row.copy()
                    new_row['AlignmentId'] = cp_id
                    rows_to_add.append(new_row)
        
        if rows_to_add:
            logger.info(f"    Duplicating {len(rows_to_add)} database rows for _cp copies")
            df = pd.concat([df, pd.DataFrame(rows_to_add)], ignore_index=True)
        
        return df
    
    tm_eval = duplicate_for_cp_copies(tm_eval)
    tm_filter_eval = duplicate_for_cp_copies(tm_filter_eval)
    
    # Filter to only transcripts we have before merging (faster)
    valid_aln_ids = set(tx_dict.keys())
    tm_eval = tm_eval[tm_eval['AlignmentId'].isin(valid_aln_ids)]
    tm_filter_eval = tm_filter_eval[tm_filter_eval['AlignmentId'].isin(valid_aln_ids)]
    
    tm_eval_df = pd.merge(tm_eval, tm_filter_eval, on=['TranscriptId', 'AlignmentId'])
    logger.info(f"    Loaded {len(tm_eval_df)} evaluation records (filtered to actual transcripts)")
    
    # Create base support_df
    support_df = tm_eval_df[['GeneId', 'TranscriptId', 'AlignmentId']].copy()
    
    # CRITICAL: Update GeneId from tx_dict to get _cp suffixed gene IDs for paralogs
    # This must happen BEFORE normalization, as tx_dict has the corrected gene IDs
    def get_gene_id_from_tx_dict(aln_id):
        tx_obj = tx_dict.get(aln_id)
        return tx_obj.name2 if tx_obj and tx_obj.name2 else None
    
    support_df['GeneId_from_tx'] = support_df['AlignmentId'].apply(get_gene_id_from_tx_dict)
    # Update GeneId where we have a value from tx_dict (this includes _cp suffixed paralogs)
    has_tx_gene = support_df['GeneId_from_tx'].notna()
    if has_tx_gene.any():
        support_df.loc[has_tx_gene, 'GeneId'] = support_df.loc[has_tx_gene, 'GeneId_from_tx']
        logger.info(f"    Updated GeneId from tx_dict for {has_tx_gene.sum()} transcripts (includes _cp paralogs)")
    support_df.drop(columns=['GeneId_from_tx'], inplace=True)
    
    # Normalize GeneId for txTM transcripts to strip _N suffix (but NOT _cpN suffix!)
    # This ensures consistency with the normalized gene IDs in tx_dict
    # Note: _cp suffixes are preserved because they don't match the _\d+$ pattern
    txTM_mask = support_df['AlignmentId'].apply(lambda x: alignment_source_map.get(x, '') == 'txTM')
    if txTM_mask.any():
        support_df.loc[txTM_mask, 'GeneId'] = support_df.loc[txTM_mask, 'GeneId'].apply(
            lambda x: normalize_gene_id(x, 'txTM') if pd.notna(x) else x
        )
        logger.info(f"    Normalized GeneId for {txTM_mask.sum()} txTM transcripts in support_df")
    
    # Add missing transcripts (e.g., augPB, txTM without metrics)
    logger.info("  Adding missing transcripts to support dataframe...")
    existing_aln_ids = set(support_df['AlignmentId'].values)
    missing_aln_ids = set(tx_dict.keys()) - existing_aln_ids

    if missing_aln_ids:
        logger.info(f"    Found {len(missing_aln_ids)} transcripts not in evaluation data (likely augPB/txTM)")
        
        # Build a lookup dict for faster reference gene finding (avoid repeated DataFrame queries)
        ref_gene_lookup = dict(zip(ref_df['TranscriptId'], ref_df['GeneId']))
        
        missing_transcripts = []
        txTM_examples = []
        for tx_id in missing_aln_ids:
            mode = alignment_source_map.get(tx_id, 'unknown')
            
            # Get gene ID directly from the genePred object (preserves _N suffixes for txTM)
            tx_obj = tx_dict.get(tx_id)
            if tx_obj and tx_obj.name2:
                gene_id = tx_obj.name2
            else:
                gene_id = f'UNKNOWN_GENE_{tx_id}'
            
            # Normalize based on mode for TranscriptId lookup
            if mode == 'txTM':
                # Strip txTM's _N suffix first, then strip_alignment_numbers handles the rest
                normalized_id = normalize_alignment_id(tx_id, mode)
                base_tx_id = tools.nameConversions.strip_alignment_numbers(normalized_id)
                if len(txTM_examples) < 3:
                    txTM_examples.append(f"{tx_id} → gene_id={gene_id}, base_tx={base_tx_id}")
            else:
                base_tx_id = tools.nameConversions.strip_alignment_numbers(tx_id)
            
            missing_transcripts.append({
                'GeneId': gene_id,
                'TranscriptId': base_tx_id,
                'AlignmentId': tx_id
            })
        
        missing_df = pd.DataFrame(missing_transcripts)
        support_df = pd.concat([support_df, missing_df], ignore_index=True)
        logger.info(f"    ✓ Added {len(missing_transcripts)} missing transcripts")

        if txTM_examples:
            logger.info(f"    TxTM ID normalization examples: {txTM_examples}")
    
    # Debug: count txTM transcripts in final support_df
    txTM_in_support = support_df['AlignmentId'].apply(lambda x: alignment_source_map.get(x, '') == 'txTM').sum()
    
    import time
    vec_start = time.time()
    
    support_df['IsAugPB'] = support_df['AlignmentId'].apply(_is_denovo)
    
    logger.info(f"  Creating support vectors for {len(support_df)} transcripts...")
    
    # Build per-exon/intron support vectors for all transcripts
    aln_ids = support_df['AlignmentId'].values
    intron_supports = []
    exon_supports = []
    intron_annot_supports = []
    exon_annot_supports = []
    cds_annot_supports = []
    
    for aln_id in aln_ids:
        tx = tx_dict[aln_id]
        num_exons = len(tx.exon_frames)
        num_introns = num_exons - 1
        
        # Default: all zeros; real values computed from BAMs by compute_real_support()
        intron_supports.append([0] * num_introns)
        exon_supports.append([0] * num_exons)
        intron_annot_supports.append([0] * num_introns)
        exon_annot_supports.append([0] * num_exons)
        cds_annot_supports.append([0] * num_exons)
    
    support_df['AllSpeciesIntronRnaSupport'] = intron_supports
    support_df['AllSpeciesExonRnaSupport'] = exon_supports
    support_df['IntronRnaSupport'] = intron_supports
    support_df['ExonRnaSupport'] = exon_supports
    support_df['IntronAnnotSupport'] = intron_annot_supports
    support_df['CdsAnnotSupport'] = cds_annot_supports
    support_df['ExonAnnotSupport'] = exon_annot_supports
    
    # Calculate default percentages
    def calc_percent_batch(support_vectors):
        return [100.0 * sum(1 for x in vec if x > 0) / len(vec) if len(vec) > 0 else 0.0
                for vec in support_vectors]
    
    support_df['IntronAnnotSupportPercent'] = calc_percent_batch(intron_annot_supports)
    support_df['ExonAnnotSupportPercent'] = calc_percent_batch(exon_annot_supports)
    support_df['CdsAnnotSupportPercent'] = calc_percent_batch(cds_annot_supports)
    support_df['ExonRnaSupportPercent'] = calc_percent_batch(exon_supports)
    support_df['IntronRnaSupportPercent'] = calc_percent_batch(intron_supports)
    support_df['AllSpeciesExonRnaSupportPercent'] = calc_percent_batch(exon_supports)
    support_df['AllSpeciesIntronRnaSupportPercent'] = calc_percent_batch(intron_supports)
    
    vec_elapsed = time.time() - vec_start
    logger.info(f"  ✓ Created default support vectors in {vec_elapsed:.1f}s")
    
    # Load alternative gene assignments for augPB transcripts (if denovo_tx_modes provided)
    if denovo_tx_modes:
        logger.info(f"  Loading alternative gene assignments for augPB transcripts (denovo_tx_modes={denovo_tx_modes})...")
        alt_assignments_df = load_alt_gene_assignments(db_path, denovo_tx_modes)
        
        if len(alt_assignments_df) > 0:
            logger.info(f"    Loaded {len(alt_assignments_df)} alternative gene assignments")
            # Count non-empty AssignedGeneIds
            non_empty = alt_assignments_df[alt_assignments_df['AssignedGeneId'].notna() & (alt_assignments_df['AssignedGeneId'] != '')].shape[0]
            logger.info(f"    {non_empty} have non-empty AssignedGeneId")
            
            # Merge into support_df
            support_df = pd.merge(support_df, alt_assignments_df, on='AlignmentId', how='left')
            logger.info(f"    ✓ Merged alternative gene assignments into support dataframe")
            
            # Debug: check how many rows have AssignedGeneId after merge
            merged_with_assigned = support_df[(support_df['AssignedGeneId'].notna()) & (support_df['AssignedGeneId'] != '')].shape[0]
            logger.info(f"    After merge: {merged_with_assigned} rows have AssignedGeneId")
            
            # Look up source gene biotype from reference annotation
            logger.info(f"    Looking up source gene biotypes from reference annotation...")
            # Create a mapping of GeneId -> GeneBiotype from ref_df
            gene_biotype_map = dict(zip(ref_df['GeneId'], ref_df['GeneBiotype']))
            # Map AssignedGeneId to SourceGeneBiotype
            support_df['SourceGeneBiotype'] = support_df['AssignedGeneId'].map(gene_biotype_map)
            source_biotypes_found = support_df['SourceGeneBiotype'].notna().sum()
            logger.info(f"    Found biotypes for {source_biotypes_found} source genes")
        else:
            logger.info("    No alternative gene assignments found")
            support_df['AssignedGeneId'] = None
            support_df['AlternativeGeneIds'] = None
            support_df['ResolutionMethod'] = None
            support_df['SourceGeneBiotype'] = None
    else:
        logger.info("  No denovo_tx_modes specified - skipping alternative gene assignments")
        support_df['AssignedGeneId'] = None
        support_df['AlternativeGeneIds'] = None
        support_df['ResolutionMethod'] = None
        support_df['SourceGeneBiotype'] = None
    
    logger.info(f"Created support dataframe with {len(support_df)} entries")
    
    # Compute real RNA-seq and annotation support from BAM files
    if bam_files or isoseq_bam_files or ref_gp_path:
        logger.info("\nComputing real support from BAM files and reference annotation...")
        compute_real_support(tx_dict, support_df, bam_files, isoseq_bam_files, ref_gp_path, alignment_source_map)
    
    return support_df


def merge_metrics(target_metrics, source_metrics):
    """Merge metrics from one chromosome into the global metrics"""
    for key in source_metrics:
        if isinstance(source_metrics[key], collections.Counter):
            target_metrics[key].update(source_metrics[key])
        elif isinstance(source_metrics[key], collections.defaultdict):
            for subkey in source_metrics[key]:
                if isinstance(source_metrics[key][subkey], (list, collections.defaultdict)):
                    if isinstance(source_metrics[key][subkey], list):
                        target_metrics[key][subkey].extend(source_metrics[key][subkey])
                    else:
                        for subsubkey in source_metrics[key][subkey]:
                            if isinstance(source_metrics[key][subkey][subsubkey], list):
                                target_metrics[key][subkey][subsubkey].extend(source_metrics[key][subkey][subsubkey])
                            else:
                                target_metrics[key][subkey][subsubkey] += source_metrics[key][subkey][subsubkey]
                else:
                    target_metrics[key][subkey] += source_metrics[key][subkey]
        elif isinstance(source_metrics[key], dict) and not isinstance(source_metrics[key], collections.Counter):
            for subkey in source_metrics[key]:
                if subkey not in target_metrics.get(key, {}):
                    continue
                if isinstance(source_metrics[key][subkey], dict):
                    for subsubkey in source_metrics[key][subkey]:
                        target_metrics[key][subkey][subsubkey] += source_metrics[key][subkey][subsubkey]
                elif isinstance(source_metrics[key][subkey], (int, float)):
                    target_metrics[key][subkey] += source_metrics[key][subkey]
        elif isinstance(source_metrics[key], int):
            if key not in target_metrics:
                target_metrics[key] = 0
            target_metrics[key] += source_metrics[key]
        elif isinstance(source_metrics[key], list):
            if key not in target_metrics:
                target_metrics[key] = []
            target_metrics[key].extend(source_metrics[key])


def process_chromosome_wrapper(chrom, chrom_num, total_chroms, chrom_tx_ids, tx_dict, support_df, 
                              mrna_metrics_df, cds_metrics_df, eval_df, tm_eval_df, ref_df, 
                              alignment_source_map, args, readthrough_gene_set, ref_gene_coords,
                              genes_with_overlaps_in_ref, gp_attrs_transcript_biotypes, gp_attrs_gene_biotypes):
    """
    Wrapper for process_chromosome that handles logging and metrics initialization.
    Designed for parallel execution.
    """
    # Set up logging for this process
    logger = logging.getLogger(__name__)
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Processing chromosome {chrom_num}/{total_chroms}: {chrom}")
    logger.info(f"{'='*60}")
    logger.info(f"  Transcripts on {chrom}: {len(chrom_tx_ids)}")
    
    # Initialize metrics for this chromosome
    chrom_metrics = initialize_metrics()
    
    # Process the chromosome
    chrom_consensus = process_chromosome(
        chrom, chrom_tx_ids, tx_dict, support_df,
        mrna_metrics_df, cds_metrics_df, eval_df, tm_eval_df,
        ref_df, alignment_source_map, args, chrom_metrics, readthrough_gene_set,
        ref_gene_coords, genes_with_overlaps_in_ref, gp_attrs_transcript_biotypes, gp_attrs_gene_biotypes
    )
    
    logger.info(f"  ✓ Selected {len(chrom_consensus)} consensus transcripts for {chrom}")
    
    return chrom_consensus, chrom_metrics


if __name__ == "__main__":
    main()

