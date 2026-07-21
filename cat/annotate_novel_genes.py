#!/usr/bin/env python3
"""
Annotate novel gene predictions with homology-based descriptions.

"""

import argparse
import collections
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile

import pandas as pd
from pathlib import Path

logger = logging.getLogger(__name__)

NOVEL_CLASS = 'putative_novel'


# ── Reference protein extraction ─────────────────────────────────────────────

def build_tx_to_gene_name_map(gp_attrs_file):
    """
    Parse the gff3ToGenePred attrs file and return a dict of
    transcript_id -> gene_name (common name, e.g. 'BRCA2').

    File format (tab-separated, no header):
        transcript_id   attribute_name   attribute_value
    """
    tx_to_gene = {}
    with open(gp_attrs_file) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            tx_id, attr, value = parts[0], parts[1], parts[2]
            if attr == 'gene_name' and value and value not in ('N/A', '.', 'none', 'None'):
                tx_to_gene[tx_id] = value
    logger.info(f"Built transcript→gene-name map: {len(tx_to_gene)} entries")
    return tx_to_gene


def extract_ref_proteins(ref_gtf, ref_fasta, gp_attrs_file, output_fasta):
    """
    Use gffread to translate the reference annotation into proteins, then
    rename each header from '{transcript_id}' to '{gene_name}|{transcript_id}'
    so DIAMOND hits carry the gene common name.

    Returns the transcript_id -> gene_name mapping used for header rewriting.
    """
    tx_to_gene = build_tx_to_gene_name_map(gp_attrs_file)

    with tempfile.NamedTemporaryFile(suffix='.faa', delete=False) as tmp:
        tmp_fasta = tmp.name

    try:
        cmd = ['gffread', ref_gtf, '-g', ref_fasta, '-y', tmp_fasta]
        logger.info(f"Extracting reference proteins: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0 or not Path(tmp_fasta).exists():
            raise RuntimeError(f"gffread failed:\n{result.stderr}")

        count = 0
        with open(tmp_fasta) as inf, open(output_fasta, 'w') as outf:
            current_seq_lines = []
            current_header = None

            def flush(header, seq_lines):
                seq = ''.join(seq_lines).replace('.', '')
                if header and seq:
                    outf.write(header + '\n')
                    # Write in 60-char lines
                    for i in range(0, len(seq), 60):
                        outf.write(seq[i:i+60] + '\n')

            for line in inf:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    flush(current_header, current_seq_lines)
                    tx_id = line[1:].split()[0]
                    gene_name = tx_to_gene.get(tx_id, tx_id)
                    current_header = f'>{gene_name}|{tx_id}'
                    current_seq_lines = []
                    count += 1
                else:
                    current_seq_lines.append(line)
            flush(current_header, current_seq_lines)

        logger.info(f"Wrote {count} reference protein sequences to {output_fasta}")
    finally:
        if Path(tmp_fasta).exists():
            os.unlink(tmp_fasta)

    return tx_to_gene


# ── DIAMOND database & search ─────────────────────────────────────────────────

def build_diamond_db(protein_fasta, db_prefix, threads):
    """Build a DIAMOND protein database (writes {db_prefix}.dmnd)."""
    cmd = [
        'diamond', 'makedb',
        '--in', protein_fasta,
        '--db', db_prefix,
        '--threads', str(threads),
        '--quiet',
    ]
    logger.info(f"Building DIAMOND database from {protein_fasta}")
    subprocess.run(cmd, check=True)


def run_diamond_blastp(query_fasta, db_path, output_tsv,
                       threads, evalue, min_identity, min_query_cover):
    """
    Run DIAMOND blastp in sensitive mode with permissive thresholds.

    Output format (tabular):
        qseqid sseqid pident qcovhsp evalue bitscore
    """
    cmd = [
        'diamond', 'blastp',
        '--query', query_fasta,
        '--db', db_path,
        '--out', output_tsv,
        '--outfmt', '6',
        'qseqid', 'sseqid', 'pident', 'qcovhsp', 'evalue', 'bitscore',
        '--threads', str(threads),
        '--evalue', str(evalue),
        '--id', str(min_identity),
        '--query-cover', str(min_query_cover),
        '--max-target-seqs', '1',
        '--sensitive',
        '--quiet',
    ]
    logger.info(f"Running DIAMOND blastp (evalue≤{evalue}, id≥{min_identity}%, "
                f"qcov≥{min_query_cover}%)")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"DIAMOND blastp failed:\n{result.stderr}")


def parse_diamond_hits(diamond_tsv):
    """
    Parse DIAMOND tabular output and return the best hit (by bitscore) per
    query transcript ID.

    Returns dict: transcript_id -> {'sseqid', 'pident', 'qcovhsp', 'evalue', 'bitscore'}
    """
    hits = {}
    if not Path(diamond_tsv).exists() or Path(diamond_tsv).stat().st_size == 0:
        logger.warning("DIAMOND output is empty — no hits found")
        return hits

    try:
        df = pd.read_csv(
            diamond_tsv, sep='\t', header=None,
            names=['qseqid', 'sseqid', 'pident', 'qcovhsp', 'evalue', 'bitscore']
        )
        df = df.sort_values('bitscore', ascending=False)
        for _, row in df.iterrows():
            qid = row['qseqid']
            if qid not in hits:
                hits[qid] = row.to_dict()
    except pd.errors.EmptyDataError:
        pass

    logger.info(f"Parsed {len(hits)} best hits from DIAMOND output")
    return hits


# ── Novel protein extraction ──────────────────────────────────────────────────

def extract_novel_proteins(consensus_protein_fasta, novel_tx_ids, output_fasta):
    """
    Write a FASTA containing only proteins for transcripts in novel_tx_ids.
    Returns the number of sequences written.
    """
    count = 0
    write_seq = False
    with open(consensus_protein_fasta) as inf, open(output_fasta, 'w') as outf:
        for line in inf:
            if line.startswith('>'):
                tx_id = line[1:].split()[0]
                write_seq = tx_id in novel_tx_ids
                if write_seq:
                    outf.write(line)
                    count += 1
            elif write_seq:
                outf.write(line)

    logger.info(f"Extracted {count} novel protein sequences for DIAMOND query")
    return count


# ── Description assignment ────────────────────────────────────────────────────

def get_gene_name_from_subject(subject_id):
    """
    Extract gene common name from a DIAMOND subject ID.
    Expected format: '{gene_name}|{transcript_id}'
    Falls back to the full subject_id if the separator is absent.
    """
    if '|' in subject_id:
        return subject_id.split('|', 1)[0]
    return subject_id


def build_description(gene_name):
    """Return the description string for a novel gene with a given best hit."""
    return f"paralog of {gene_name}"


# ── gp_info update ────────────────────────────────────────────────────────────

def update_gp_info(gp_info_file, output_gp_info, tx_descriptions, tx_novel_class=None,
                   drop_tx_ids=None):
    """
    Add 'novel_gene_description' and 'novel_class' columns to the gp_info TSV.

    tx_descriptions: dict of transcript_id -> description string (or 'N/A')
    tx_novel_class:  dict of transcript_id -> 'paralog' | 'lineage_specific'
        for novel transcripts (others get 'N/A'). Distinguishes novel genes that
        are duplicate copies of a known (reference) gene from those with no
        reference homolog (candidate lineage-specific genes). Both are RETAINED;
        this only labels them.
    drop_tx_ids:     optional set of transcript_ids to remove entirely (excess
        novel paralog copies from the copy-number cap).
    """
    tx_novel_class = tx_novel_class or {}
    drop_tx_ids = drop_tx_ids or set()
    df = pd.read_csv(gp_info_file, sep='\t', index_col=[0, 1])
    if drop_tx_ids:
        keep = ~df.index.get_level_values('transcript_id').isin(drop_tx_ids)
        df = df[keep]
    tx_ids = df.index.get_level_values('transcript_id')
    df.insert(
        len(df.columns),
        'novel_gene_description',
        [tx_descriptions.get(tx_id, 'N/A') for tx_id in tx_ids]
    )
    df.insert(
        len(df.columns),
        'novel_class',
        [tx_novel_class.get(tx_id, 'N/A') for tx_id in tx_ids]
    )
    with open(output_gp_info, 'w') as fh:
        df.to_csv(fh, sep='\t', na_rep='N/A')

    n_assigned = sum(1 for v in tx_descriptions.values() if v != 'N/A')
    n_para = sum(1 for v in tx_novel_class.values() if v == 'paralog')
    n_ls = sum(1 for v in tx_novel_class.values() if v == 'lineage_specific')
    logger.info(f"Updated gp_info: {n_assigned} transcripts assigned descriptions, "
                f"novel_class = {n_para} paralog / {n_ls} lineage_specific "
                f"(written to {output_gp_info})")


# ── GFF3 update ───────────────────────────────────────────────────────────────

def _parse_gff3_attrs(attrs_str):
    """Return an ordered dict of key->value from a GFF3 attributes string."""
    d = collections.OrderedDict()
    for item in attrs_str.split(';'):
        item = item.strip()
        if not item:
            continue
        if '=' in item:
            k, v = item.split('=', 1)
            d[k] = v
        else:
            d[item] = ''
    return d


def _render_gff3_attrs(d):
    """Render an ordered dict back to a GFF3 attributes string."""
    return ';'.join(f"{k}={v}" if v != '' else k for k, v in d.items())


def update_gff3(input_gff3, output_gff3, gene_descriptions, tx_descriptions,
                drop_gene_ids=None, drop_tx_ids=None):
    """
    Add a 'description' attribute to gene and transcript records in the GFF3
    for any entry that has a novel gene description.

    gene_descriptions: dict gene_id -> description
    tx_descriptions:  dict transcript_id -> description
    drop_gene_ids/drop_tx_ids: optional sets of gene/transcript ids to remove
        entirely (excess novel paralog copies). Any feature whose ID or Parent is
        in these sets is skipped, cascading to child exon/CDS records.
    """
    drop_gene_ids = drop_gene_ids or set()
    drop_tx_ids = drop_tx_ids or set()
    n_gene = 0
    n_tx = 0
    with open(input_gff3) as inf, open(output_gff3, 'w') as outf:
        for line in inf:
            if line.startswith('#') or line.strip() == '':
                outf.write(line)
                continue

            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                outf.write(line)
                continue

            feature = parts[2]
            attrs = _parse_gff3_attrs(parts[8])

            if drop_gene_ids or drop_tx_ids:
                _id = attrs.get('ID', '')
                _parents = attrs.get('Parent', '').split(',') if attrs.get('Parent') else []
                if _id in drop_gene_ids or _id in drop_tx_ids or \
                   any(p in drop_gene_ids or p in drop_tx_ids for p in _parents):
                    continue

            description = None

            if feature == 'gene':
                gene_id = attrs.get('ID') or attrs.get('gene_id')
                if gene_id:
                    description = gene_descriptions.get(gene_id)
                    if description:
                        n_gene += 1

            elif feature in ('transcript', 'mRNA'):
                tx_id = attrs.get('ID')
                if tx_id:
                    description = tx_descriptions.get(tx_id)
                    if description:
                        n_tx += 1

            if description:
                # URL-encode characters that would break GFF3 attribute parsing
                desc_enc = description.replace('%', '%25').replace(
                    '=', '%3D').replace(';', '%3B').replace(',', '%2C')
                attrs['description'] = desc_enc
                parts[8] = _render_gff3_attrs(attrs)

            outf.write('\t'.join(parts) + '\n')

    logger.info(f"Updated GFF3: {n_gene} gene records, {n_tx} transcript records "
                f"annotated with descriptions")


# ── genePred filter ───────────────────────────────────────────────────────────

def filter_gp(input_gp, output_gp, drop_tx_ids):
    """Copy the genePred, dropping rows whose transcript id (column 0) is in
    ``drop_tx_ids`` (excess novel paralog copies). Keeps the structural .gp
    consistent with the capped gp_info / GFF3."""
    drop_tx_ids = drop_tx_ids or set()
    kept = 0
    dropped = 0
    with open(input_gp) as inf, open(output_gp, 'w') as outf:
        for line in inf:
            tx_id = line.split('\t', 1)[0]
            if tx_id in drop_tx_ids:
                dropped += 1
                continue
            outf.write(line)
            kept += 1
    logger.info(f"Filtered genePred: kept {kept}, dropped {dropped} transcript rows "
                f"(written to {output_gp})")


# ── Main ──────────────────────────────────────────────────────────────────────

def update_metrics_json(in_path, out_path, dropped_mode_counts, n_dropped_genes,
                        n_dropped_tx, fams_capped, cap):
    """Copy the consensus metrics JSON to out_path, correcting count-based metrics
    for the novel paralog copies removed by the cap so downstream plots reflect the
    final annotated set. Only counts we can adjust exactly are touched; quality
    distributions (coverage/identity/support/completeness, which are reference-
    transcript metrics unaffected by dropping novel paralogs) are left as-is."""
    import json
    with open(in_path) as fh:
        m = json.load(fh)
    tm = m.get('Transcript Modes')
    if isinstance(tm, dict):
        for mode, cnt in dropped_mode_counts.items():
            if isinstance(tm.get(mode), (int, float)):
                tm[mode] = max(0, int(tm[mode]) - int(cnt))
    m['Novel paralog cap'] = {
        'max_copies_per_family': int(cap),
        'families_capped': int(fams_capped),
        'dropped_genes': int(n_dropped_genes),
        'dropped_transcripts': int(n_dropped_tx),
    }
    with open(out_path, 'w') as fh:
        json.dump(m, fh, indent=1)


def main():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    parser = argparse.ArgumentParser(description=__doc__)

    # Input
    parser.add_argument('--consensus-protein-fasta', required=True,
                        help='Protein FASTA from generate_consensus '
                             '({genome}_consensus_protein.fasta).')
    parser.add_argument('--consensus-gp-info', required=True,
                        help='gp_info TSV from generate_consensus '
                             '({genome}_consensus.gp_info).')
    parser.add_argument('--consensus-gff3', required=True,
                        help='GFF3 from generate_consensus '
                             '({genome}_consensus.gff3).')
    parser.add_argument('--consensus-gp', default=None,
                        help='genePred from generate_consensus '
                             '({genome}_consensus.gp). Required with --output-gp so the '
                             'structural set stays consistent when paralog copies are capped.')
    parser.add_argument('--ref-gtf', required=True,
                        help='Reference annotation GTF '
                             '(reference/{ref_genome}.gtf from prepare_reference_files).')
    parser.add_argument('--ref-fasta', required=True,
                        help='Reference genome FASTA '
                             '(genome_files/{ref_genome}.fa).')
    parser.add_argument('--ref-gp-attrs', required=True,
                        help='Reference gp_attrs file for transcript→gene-name mapping '
                             '(reference/{ref_genome}.gp_attrs).')

    # Output
    parser.add_argument('--output-gff3', required=True,
                        help='Output GFF3 with novel gene descriptions added.')
    parser.add_argument('--output-gp-info', required=True,
                        help='Output gp_info TSV with novel_gene_description column.')
    parser.add_argument('--output-gp', default=None,
                        help='Output genePred with capped paralog copies removed '
                             '(requires --consensus-gp).')
    parser.add_argument('--consensus-metrics-json', default=None,
                        help='Metrics JSON from generate_consensus '
                             '({genome}_consensus.json). If given with --output-metrics-json, '
                             'count-based metrics (Transcript Modes) are corrected for the '
                             'paralog copies removed here so downstream plots reflect the '
                             'final annotated gene set.')
    parser.add_argument('--output-metrics-json', default=None,
                        help='Output metrics JSON reflecting the post-cap gene set '
                             '(requires --consensus-metrics-json).')

    # Novel paralog copy cap
    parser.add_argument('--novel-paralog-max-copies', type=int, default=0,
                        help='Cap the number of novel protein-coding copies per source gene '
                             'family (family = the "paralog of GENE" reference gene from DIAMOND), '
                             'keeping the longest-protein copies. Removes the protein multi-mapping '
                             'artifact where one protein (e.g. MIPOL1, L1TD1) seeds hundreds of '
                             'dispersed "novel paralog" predictions. Lineage-specific novels (no '
                             'reference homolog) are never capped. 0 disables [0].')

    # DIAMOND parameters (permissive defaults)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--evalue', type=float, default=1e-3,
                        help='E-value cutoff for DIAMOND blastp (default: 1e-3).')
    parser.add_argument('--min-identity', type=float, default=20.0,
                        help='Minimum %% identity for hits (default: 20).')
    parser.add_argument('--min-query-cover', type=float, default=20.0,
                        help='Minimum %% query coverage for hits (default: 20).')

    args = parser.parse_args()

    # ── Step 1: identify novel transcripts ───────────────────────────────────
    logger.info("Reading consensus gp_info...")
    gp_info_df = pd.read_csv(args.consensus_gp_info, sep='\t', index_col=[0, 1])

    tc_col = 'transcript_class'
    if tc_col not in gp_info_df.columns:
        logger.error(f"Column '{tc_col}' not found in gp_info — nothing to annotate.")
        shutil.copy(args.consensus_gff3, args.output_gff3)
        shutil.copy(args.consensus_gp_info, args.output_gp_info)
        if args.output_gp and args.consensus_gp:
            shutil.copy(args.consensus_gp, args.output_gp)
        if args.output_metrics_json and args.consensus_metrics_json:
            shutil.copy(args.consensus_metrics_json, args.output_metrics_json)
        sys.exit(0)

    novel_mask = gp_info_df[tc_col] == NOVEL_CLASS
    novel_rows = gp_info_df[novel_mask]
    novel_tx_ids = set(novel_rows.index.get_level_values('transcript_id'))

    logger.info(f"Found {len(novel_tx_ids)} novel transcripts (class='{NOVEL_CLASS}')")

    if not novel_tx_ids:
        logger.info("No novel transcripts — copying inputs to outputs unchanged.")
        shutil.copy(args.consensus_gff3, args.output_gff3)
        shutil.copy(args.consensus_gp_info, args.output_gp_info)
        if args.output_gp and args.consensus_gp:
            shutil.copy(args.consensus_gp, args.output_gp)
        if args.output_metrics_json and args.consensus_metrics_json:
            shutil.copy(args.consensus_metrics_json, args.output_metrics_json)
        sys.exit(0)

    # ── Remainder runs inside a temp directory ────────────────────────────────
    with tempfile.TemporaryDirectory(prefix='novel_annot_') as tmpdir:
        tmpdir = Path(tmpdir)

        # Step 2: extract reference proteins and build DIAMOND DB
        ref_prot_fasta = tmpdir / 'ref_proteins.faa'
        logger.info("Extracting reference proteins...")
        extract_ref_proteins(
            args.ref_gtf, args.ref_fasta, args.ref_gp_attrs, str(ref_prot_fasta)
        )

        if not ref_prot_fasta.exists() or ref_prot_fasta.stat().st_size == 0:
            logger.error("Reference protein FASTA is empty — cannot annotate novel genes.")
            shutil.copy(args.consensus_gff3, args.output_gff3)
            shutil.copy(args.consensus_gp_info, args.output_gp_info)
            if args.output_gp and args.consensus_gp:
                shutil.copy(args.consensus_gp, args.output_gp)
            if args.output_metrics_json and args.consensus_metrics_json:
                shutil.copy(args.consensus_metrics_json, args.output_metrics_json)
            sys.exit(1)

        diamond_db = str(tmpdir / 'ref_proteins')
        build_diamond_db(str(ref_prot_fasta), diamond_db, args.threads)

        # Step 3: extract novel gene protein sequences
        novel_prot_fasta = tmpdir / 'novel_proteins.faa'
        n_extracted = extract_novel_proteins(
            args.consensus_protein_fasta, novel_tx_ids, str(novel_prot_fasta)
        )

        if n_extracted == 0:
            logger.warning("None of the novel transcripts have CDS protein sequences "
                           "— descriptions will remain 'N/A'.")
            shutil.copy(args.consensus_gff3, args.output_gff3)
            shutil.copy(args.consensus_gp_info, args.output_gp_info)
            if args.output_gp and args.consensus_gp:
                shutil.copy(args.consensus_gp, args.output_gp)
            sys.exit(0)

        # Step 4: DIAMOND blastp
        diamond_out = str(tmpdir / 'diamond_hits.tsv')
        run_diamond_blastp(
            str(novel_prot_fasta), diamond_db, diamond_out,
            args.threads, args.evalue, args.min_identity, args.min_query_cover
        )

        # Step 5: parse hits → per-transcript descriptions
        hits = parse_diamond_hits(diamond_out)

        tx_descriptions = {}
        for tx_id, hit in hits.items():
            gene_name = get_gene_name_from_subject(str(hit['sseqid']))
            tx_descriptions[tx_id] = build_description(gene_name)

        # Step 6: derive per-gene descriptions
        # Build gene_id -> [transcript_ids] from novel rows
        gene_to_txs = collections.defaultdict(list)
        for gene_id, tx_id in novel_rows.index:
            gene_to_txs[gene_id].append(tx_id)

        # For each gene, pick the description from the first transcript that has one
        gene_descriptions = {}
        for gene_id, tx_list in gene_to_txs.items():
            for tx_id in tx_list:
                if tx_id in tx_descriptions:
                    gene_descriptions[gene_id] = tx_descriptions[tx_id]
                    break

        logger.info(f"Assigned descriptions: {len(tx_descriptions)} transcripts, "
                    f"{len(gene_descriptions)} genes")

        # Classify each novel gene: 'paralog' if its protein hits a reference gene
        # (duplicate/expansion of a known gene), else 'lineage_specific' (no
        # reference homolog). Applied at gene level and propagated to all its
        # isoforms. Neither is dropped -- this only labels them for honest stats.
        tx_novel_class = {}
        for gene_id, tx_list in gene_to_txs.items():
            desc = gene_descriptions.get(gene_id)
            cls = 'paralog' if (desc and desc.startswith('paralog')) else 'lineage_specific'
            for tx_id in tx_list:
                tx_novel_class[tx_id] = cls

        # Step 6b: cap novel paralog copies per source gene family. A few "sticky"
        # proteins multi-map to hundreds of dispersed repeat/pseudogene loci, each
        # becoming a "novel paralog" (e.g. MIPOL1 -> >1000 copies). These pass every
        # quality filter (coverage/identity/exon-count/expression); only copy-number
        # per family separates them. Keep the longest-protein copies per family.
        drop_gene_ids = set()
        drop_tx_ids = set()
        fams_capped = 0
        cap = int(args.novel_paralog_max_copies or 0)
        if cap > 0:
            prot_len = {}
            cur = None
            with open(args.consensus_protein_fasta) as fh:
                for line in fh:
                    if line.startswith('>'):
                        tid = line[1:].split()[0]
                        cur = tid if tid in novel_tx_ids else None
                        if cur:
                            prot_len[cur] = 0
                    elif cur:
                        prot_len[cur] += len(line.strip())
            fam_to_genes = collections.defaultdict(list)
            for gene_id, desc in gene_descriptions.items():
                if not (desc and desc.startswith('paralog of ')):
                    continue  # lineage-specific / undescribed -> never capped
                fam = desc[len('paralog of '):]
                glen = max((prot_len.get(tx, 0) for tx in gene_to_txs.get(gene_id, [])),
                           default=0)
                fam_to_genes[fam].append((glen, gene_id))
            fams_capped = 0
            for fam, glist in fam_to_genes.items():
                if len(glist) <= cap:
                    continue
                glist.sort(reverse=True)  # longest protein first
                for _, gid in glist[cap:]:
                    drop_gene_ids.add(gid)
                    drop_tx_ids.update(gene_to_txs.get(gid, []))
                fams_capped += 1
            logger.info(
                f"Novel paralog copy cap (<= {cap}/family): capped {fams_capped} families, "
                f"dropped {len(drop_gene_ids)} novel paralog genes "
                f"({len(drop_tx_ids)} transcripts)")

        # Step 7: update gp_info (dropping capped copies)
        update_gp_info(args.consensus_gp_info, args.output_gp_info,
                       tx_descriptions, tx_novel_class, drop_tx_ids=drop_tx_ids)

        # Step 8: update GFF3 (dropping capped copies)
        update_gff3(args.consensus_gff3, args.output_gff3, gene_descriptions, tx_descriptions,
                    drop_gene_ids=drop_gene_ids, drop_tx_ids=drop_tx_ids)

        # Step 9: emit a structurally consistent genePred if requested
        if args.output_gp:
            if args.consensus_gp:
                filter_gp(args.consensus_gp, args.output_gp, drop_tx_ids)
            else:
                logger.warning("--output-gp given without --consensus-gp; skipping genePred filter.")

        # Step 10: emit a post-cap metrics JSON so downstream plots reflect the
        # final annotated gene set (correct Transcript Modes counts for the
        # capped paralog copies we removed here).
        if args.output_metrics_json:
            if args.consensus_metrics_json:
                dropped_mode_counts = {}
                if drop_tx_ids and 'alignment_mode' in gp_info_df.columns:
                    lvl = gp_info_df.index.get_level_values('transcript_id')
                    mask = lvl.isin(drop_tx_ids)
                    dropped_mode_counts = (
                        gp_info_df.loc[mask, 'alignment_mode'].value_counts().to_dict())
                update_metrics_json(
                    args.consensus_metrics_json, args.output_metrics_json,
                    dropped_mode_counts, len(drop_gene_ids), len(drop_tx_ids),
                    fams_capped, cap)
                logger.info(
                    f"Wrote post-cap metrics JSON (adjusted modes: {dropped_mode_counts}) "
                    f"to {args.output_metrics_json}")
            else:
                logger.warning("--output-metrics-json given without --consensus-metrics-json; "
                               "skipping metrics update.")

    logger.info("Novel gene annotation complete.")


if __name__ == '__main__':
    main()
