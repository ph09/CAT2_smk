import tools.bio
import argparse
import collections
import json
import logging
import sys
import time
import warnings
import pandas as pd
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

DENOVO_PREFIXES = ('augPB-', 'strg-')

def _is_denovo(aln_id):
    """Check if an alignment ID belongs to a de novo mode (augPB or strg)."""
    return isinstance(aln_id, str) and aln_id.startswith(DENOVO_PREFIXES)

# Suppress warnings for better performance
warnings.filterwarnings('ignore')
warnings.filterwarnings('ignore', category=FutureWarning, module='pandas')
warnings.filterwarnings('ignore', category=DeprecationWarning, module='pandas')

logger = logging.getLogger(__name__)
ID_TEMPLATE = '{genome:.10}_{tag_type}{unique_id:07d}'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    
    # Run the main consensus logic
    metrics = generate_consensus(args)

    # Write the final metrics file
    with open(args.metrics_json, 'w') as outf:
        json.dump(metrics, outf, indent=4)
    logger.info(f"Successfully generated consensus gene set and metrics for {args.genome}.")


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
    parser.add_argument("--intron-rnaseq-support", type=float, default=5.0, help="Percent of introns that must be supported by RNA-seq.")
    parser.add_argument("--exon-rnaseq-support", type=float, default=5.0, help="Percent of exons supported by RNA-seq.")
    parser.add_argument("--intron-annot-support", type=float, default=20.0, help="Percent of introns supported by reference annotation.")
    parser.add_argument("--exon-annot-support", type=float, default=20.0, help="Percent of exons supported by reference annotation.")
    parser.add_argument("--original-intron-support", type=float, default=80.0, help="Percent of original introns that must be preserved.")
    parser.add_argument("--in-species-rna-support-only", action="store_true", help="Use in-species RNA-seq support only, ignoring cross-species evidence.")
    parser.add_argument("--filter-overlapping-genes", action="store_true", help="Filter out overlapping CDS intervals from different genes.")
    parser.add_argument("--overlapping-ignore-bases", type=int, default=0, help="Number of bases to ignore when clustering for overlap.")
    parser.add_argument("--txTM-min-coverage", type=float, default=0.0, help="Minimum alignment coverage (0-100) required for txTM transcripts. Transcripts below this threshold (with explicit 0 or very low coverage in the metrics DB) are filtered as likely gene fragments. Default: 0.0 (only filter explicit zeros).")
    
    # --- De Novo Parameters ---
    parser.add_argument("--denovo-tx-modes", nargs='*', default=[], help="List of de novo modes to consider (e.g., augPB).")
    parser.add_argument("--denovo-num-introns", type=int, default=1, help="A de novo isoform must have at least this many introns.")
    parser.add_argument("--denovo-splice-support", type=float, default=10.0, help="Percent of de novo splices that must be RNA-seq supported.")
    parser.add_argument("--denovo-exon-support", type=float, default=10.0, help="Percent of de novo exons that must be RNA-seq supported.")
    parser.add_argument("--denovo-ignore-novel-genes", action="store_true", help="If set, only incorporate de novo transcripts as novel isoforms, not novel genes.")
    parser.add_argument("--denovo-only-novel-genes", action="store_true", help="If set, only incorporate de novo transcripts if they are novel genes.")
    parser.add_argument("--denovo-allow-unsupported", action="store_true", help="Allow de novo transcripts with novel splices that lack RNA-seq support.")
    parser.add_argument("--denovo-allow-bad-annot-or-tm", action="store_true", help="Allow de novo models flagged as 'badAnnotOrTm'.")
    parser.add_argument("--denovo-allow-novel-ends", action="store_true", help="Allow de novo models with novel 5' or 3' ends.")
    parser.add_argument("--denovo-novel-end-distance", type=int, default=100)

    # --- PacBio Parameters ---
    parser.add_argument("--require-pacbio-support", action="store_true", help="If set, remove any consensus transcript not validated by Iso-Seq data.")
    parser.add_argument("--hints-db-has-rnaseq", action="store_true", help="Flag if the hints DB contains RNA-seq, for tagging purposes.")


def generate_consensus(args):
    """
    Main consensus finding logic.

    :param args: Argument namespace from luigi
    """
    # Create a mapping from alignment IDs to their source files
    alignment_source_map = {}
    
    # Load genePreds and track their sources
    for gp_file in args.gp_list:
        for t in tools.transcripts.gene_pred_iterator(gp_file):
            if 'transMap_pairwise' in gp_file:
                alignment_source_map[t.name] = 'transMap_pairwise'
            elif 'transMap' in gp_file:
                alignment_source_map[t.name] = 'transMap'
            elif 'txTM' in gp_file:
                alignment_source_map[t.name] = 'txTM'
            elif 'augTMR_pairwise' in gp_file:
                alignment_source_map[t.name] = 'augTMR_pairwise'
            elif 'augTM_pairwise' in gp_file:
                alignment_source_map[t.name] = 'augTM_pairwise'
            elif 'augTMR' in gp_file:
                alignment_source_map[t.name] = 'augTMR'
            elif 'augTM' in gp_file:
                alignment_source_map[t.name] = 'augTM'
            elif 'augMP' in gp_file:
                alignment_source_map[t.name] = 'augMP'
            elif 'augPB' in gp_file:
                alignment_source_map[t.name] = 'augPB'
            elif '_strg' in gp_file:
                alignment_source_map[t.name] = 'strg'
    
    # Set the global mapping for use by alignment_type function
    tools.nameConversions.set_alignment_source_map(alignment_source_map)
    
    # load all genePreds
    tx_dict = tools.transcripts.load_gps(args.gp_list)
    # load reference annotation information
    ref_df = tools.sqlInterface.load_annotation(args.ref_db_path)
    ref_biotype_counts = collections.Counter(ref_df.TranscriptBiotype)
    coding_count = ref_biotype_counts['protein_coding']
    non_coding_count = sum(y for x, y in ref_biotype_counts.items() if x != 'protein_coding')
    # gene transcript map to iterate over so that we capture missing gene information
    gene_biotype_map = tools.sqlInterface.get_gene_biotype_map(args.ref_db_path)
    transcript_biotype_map = tools.sqlInterface.get_transcript_biotype_map(args.ref_db_path)
    # load transMap evaluation data
    tm_eval_df = load_transmap_evals(args.db_path)
    # pd.set_option('max_columns', None)
    # Determine which modes are available for this genome
    txms = ['transMap', 'augTM', 'augTMR', 'augMP', 'txTM']
    
    # Check if augPB data exists for this genome by looking for augPB transcripts
    has_augpb = any(_is_denovo(tx_id) for tx_id in tx_dict.keys())
    if has_augpb:
        txms.append('augPB')
    
    # Check if augMP data exists
    has_augmp = any(tx_id.startswith('augMP-') for tx_id in tx_dict.keys())
    if not has_augmp and 'augMP' in txms:
        txms.remove('augMP')
    
    tx_modes = [x for x in txms if x in ['transMap', 'augTM', 'augTMR', 'augMP', 'txTM', 'augPB', 'strg']]
    # denovo and augMP don't have metrics tables, so exclude them from metrics loading
    tx_modes_with_metrics = [x for x in tx_modes if x not in ['augPB', 'strg', 'augMP']]
    # Get tm df but don't drop AlignmentId yet
    tm_eval = tools.sqlInterface.load_alignment_evaluation(args.db_path)
    tm_filter_eval = tools.sqlInterface.load_filter_evaluation(args.db_path)
    tm_eval_df = pd.merge(tm_eval, tm_filter_eval, on=['TranscriptId', 'AlignmentId'])
    
    # Create base support_df from evaluation tables
    support_df = tm_eval_df.filter(['GeneId', 'TranscriptId', 'AlignmentId'], axis=1)
    
    missing_transcripts = []
    for tx_id, tx in tx_dict.items():
        if tx_id not in support_df['AlignmentId'].values:
            base_tx_id = tools.nameConversions.strip_alignment_numbers(tx_id)
            gene_name_from_gp = tx.name2
            
            import re
            base_gene_name = re.sub(r'_\d+$', '', gene_name_from_gp) if gene_name_from_gp else None
            
            # Try to find gene ID from reference using base name
            gene_rows = ref_df[ref_df['TranscriptId'] == base_tx_id]
            if len(gene_rows) > 0:
                gene_id = gene_rows.iloc[0]['GeneId']
            else:
                if base_tx_id.endswith('_cp'):
                    gene_id = base_tx_id[:-3]  # Remove _cp suffix
                else:
                    gene_id = base_gene_name or f'UNKNOWN_GENE_{base_tx_id}'
            
            missing_transcripts.append({
                'GeneId': gene_id,
                'TranscriptId': base_tx_id,
                'AlignmentId': tx_id
            })
    
    if missing_transcripts:
        missing_df = pd.DataFrame(missing_transcripts)
        support_df = pd.concat([support_df, missing_df], ignore_index=True)

    def dummy_support_exons(aln_id):
        num_exon_frames = len(tx_dict[aln_id].exon_frames)
        # For denovo transcripts, give them full exon support (PacBio/StringTie-validated)
        if _is_denovo(aln_id):
            return [1] * num_exon_frames
        else:
            return [0] * num_exon_frames

    def dummy_support_introns(aln_id):
        num_intron_frames = len(tx_dict[aln_id].exon_frames) - 1
        # For denovo transcripts, give them full intron support (PacBio/StringTie-validated)
        if _is_denovo(aln_id):
            return [1] * num_intron_frames
        else:
            return [0] * num_intron_frames

    support_df['AllSpeciesIntronRnaSupport'] = support_df.apply(lambda x: dummy_support_introns(x['AlignmentId']), axis=1)
    support_df['AllSpeciesExonRnaSupport'] = support_df.apply(lambda x: dummy_support_exons(x['AlignmentId']), axis=1)
    support_df['IntronRnaSupport'] = support_df.apply(lambda x: dummy_support_introns(x['AlignmentId']), axis=1)
    support_df['ExonRnaSupport'] = support_df.apply(lambda x: dummy_support_exons(x['AlignmentId']), axis=1)
    # Annotation support should be 0 for augPB (novel splices/exons)
    support_df['IntronAnnotSupport'] = support_df.apply(lambda x: [0] * (len(tx_dict[x['AlignmentId']].exon_frames) - 1), axis=1)
    support_df['CdsAnnotSupport'] = support_df.apply(lambda x: [0] * len(tx_dict[x['AlignmentId']].exon_frames), axis=1)
    support_df['ExonAnnotSupport'] = support_df.apply(lambda x: [0] * len(tx_dict[x['AlignmentId']].exon_frames), axis=1)
    support_df['IntronAnnotSupportPercent'] = 0.0
    support_df['ExonAnnotSupportPercent'] = 0.0
    support_df['CdsAnnotSupportPercent'] = 0.0
    
    # Set support percentages differently for augPB transcripts
    def set_support_percentage(row, column_name, high_value, low_value):
        if _is_denovo(row['AlignmentId']):
            return high_value
        else:
            return low_value
    
    support_df['ExonRnaSupportPercent'] = support_df.apply(lambda x: set_support_percentage(x, 'ExonRnaSupportPercent', 80.0, 0.0), axis=1)
    support_df['IntronRnaSupportPercent'] = support_df.apply(lambda x: set_support_percentage(x, 'IntronRnaSupportPercent', 80.0, 0.0), axis=1)
    support_df['AllSpeciesExonRnaSupportPercent'] = support_df.apply(lambda x: set_support_percentage(x, 'AllSpeciesExonRnaSupportPercent', 80.0, 0.0), axis=1)
    support_df['AllSpeciesIntronRnaSupportPercent'] = support_df.apply(lambda x: set_support_percentage(x, 'AllSpeciesIntronRnaSupportPercent', 80.0, 0.0), axis=1)
    tm_eval_df = tm_eval_df.drop('AlignmentId', axis=1)

    # load the alignment metrics data (excluding augPB which doesn't have metrics tables)
    mrna_metrics_df = pd.concat([load_metrics_from_db(args.db_path, tx_mode, 'mRNA') for tx_mode in tx_modes_with_metrics])
    cds_metrics_df = pd.concat([load_metrics_from_db(args.db_path, tx_mode, 'CDS') for tx_mode in tx_modes_with_metrics])
    eval_df = pd.concat([load_evaluations_from_db(args.db_path, tx_mode) for tx_mode in tx_modes_with_metrics]).reset_index()
    coding_df, non_coding_df = combine_and_filter_dfs(tx_dict, support_df, mrna_metrics_df, cds_metrics_df, tm_eval_df,
                                                      ref_df, eval_df, args.intron_rnaseq_support,
                                                      args.exon_rnaseq_support, args.intron_annot_support,
                                                      args.exon_annot_support, args.original_intron_support,
                                                      args.in_species_rna_support_only, alignment_source_map=alignment_source_map,
                                                      txTM_min_coverage=args.txTM_min_coverage)
    if len(coding_df) + len(non_coding_df) == 0:
        raise RuntimeError('No transcripts pass filtering for species {}. '
                           'Consider lowering requirements. Please see the manual.'.format(args.genome))
    elif len(coding_df) == 0 and coding_count > 0:
        logger.warning('No protein coding transcripts pass filtering for species {}. '
                       'Consider lowering requirements. Please see the manual.'.format(args.genome))
    elif len(non_coding_df) == 0 and non_coding_count > 0:
        logger.warning('No non-coding transcripts pass filtering for species {}. '
                       'Consider lowering requirements. Please see the manual.'.format(args.genome))
    scored_coding_df, scored_non_coding_df = score_filtered_dfs(coding_df, non_coding_df,
                                                                args.in_species_rna_support_only)
    scored_df = merge_scored_dfs(scored_coding_df, scored_non_coding_df)
    
    txTM_aln_ids = [aln_id for aln_id, source in alignment_source_map.items() if source == 'txTM']
    txTM_cnv_aln_ids = [aln_id for aln_id in txTM_aln_ids if '_' in aln_id and aln_id.split('_')[-1].isdigit()]
    txTM_base_aln_ids = [aln_id for aln_id in txTM_aln_ids if aln_id not in txTM_cnv_aln_ids]
    
    is_txTM_cnv = scored_df['AlignmentId'].isin(txTM_cnv_aln_ids)
    is_txTM_base = scored_df['AlignmentId'].isin(txTM_base_aln_ids)
    is_other = ~(is_txTM_cnv | is_txTM_base)
    
    txTM_cnv_df = scored_df[is_txTM_cnv].copy()
    txTM_base_df = scored_df[is_txTM_base].copy()
    other_df = scored_df[is_other].copy()
    if len(txTM_cnv_df) > 0:
        best_cnv = txTM_cnv_df.groupby('AlignmentId')['TranscriptScore'].transform("max") == txTM_cnv_df['TranscriptScore']
        best_cnv_df = txTM_cnv_df[best_cnv].reset_index()
        logger.info(f"  TxTM CNV: kept {len(best_cnv_df)} transcripts (grouped by AlignmentId)")
    else:
        best_cnv_df = pd.DataFrame()
    
    if len(txTM_base_df) > 0:
        best_txTM_base = txTM_base_df.groupby('TranscriptId')['TranscriptScore'].transform("max") == txTM_base_df['TranscriptScore']
        best_txTM_base_df = txTM_base_df[best_txTM_base].reset_index()
        logger.info(f"  TxTM base: kept {len(best_txTM_base_df)} transcripts (grouped by TranscriptId, independent)")
    else:
        best_txTM_base_df = pd.DataFrame()
    
    # For other sources (transMap, augTM, etc): group by TranscriptId
    if len(other_df) > 0:
        best_other = other_df.groupby('TranscriptId')['TranscriptScore'].transform("max") == other_df['TranscriptScore']
        best_other_df = other_df[best_other].reset_index()
        logger.info(f"  Other sources: kept {len(best_other_df)} transcripts (grouped by TranscriptId)")
    else:
        best_other_df = pd.DataFrame()
    
    # Combine all results
    best_df = pd.concat([best_cnv_df, best_txTM_base_df, best_other_df], ignore_index=True)

    metrics = {'Transcript Missing': collections.Counter(),
               'Gene Missing': collections.Counter(),
               'Transcript Modes': collections.Counter(),  # coding only
               'Duplicate transcripts': collections.Counter(),
               'Discarded by strand resolution': 0,
               'Coverage': collections.defaultdict(list),
               'Identity': collections.defaultdict(list),
               'Splice Support': collections.defaultdict(list),
               'Exon Support': collections.defaultdict(list),
               'Original Introns': collections.defaultdict(list),
               'Splice Annotation Support': collections.defaultdict(list),
               'Exon Annotation Support': collections.defaultdict(list),
               'IsoSeq Transcript Validation': collections.Counter()}

    # we can keep track of missing stuff now
    for gene_biotype, tx_df in best_df.groupby('GeneBiotype'):
        biotype_genes = {gene_id for gene_id, b in gene_biotype_map.items() if b == gene_biotype}
        metrics['Gene Missing'][gene_biotype] = len(biotype_genes) - len(set(tx_df.GeneId))
    for tx_biotype, tx_df in best_df.groupby('TranscriptBiotype'):
        biotype_txs = {gene_id for gene_id, b in transcript_biotype_map.items() if b == tx_biotype}
        metrics['Transcript Missing'][tx_biotype] = len(biotype_txs) - len(set(tx_df.TranscriptId))

    consensus_dict = {}

    for aln_id, s in best_df.groupby('AlignmentId'):
        # Get the gene_id from the first row
        first_row = s.iloc[0]
        gene_id = first_row['GeneId']
        
        aln_id, m = incorporate_tx(s, gene_id, metrics, args.hints_db_has_rnaseq)
        consensus_dict[aln_id] = m

    # if we ran in either denovo mode, load those data and detect novel genes
    if len(args.denovo_tx_modes) > 0:
        metrics['denovo'] = {}
        for tx_mode in args.denovo_tx_modes:
            if args.denovo_allow_bad_annot_or_tm == True:
                metrics['denovo'][tx_mode] = {'Possible paralog': 0, 'Poor alignment': 0, 'Putative novel': 0,
                                          'Possible fusion': 0, 'Putative novel isoform': 0, 'Bad annot or tm': 0}
            else:
                metrics['denovo'][tx_mode] = {'Possible paralog': 0, 'Poor alignment': 0, 'Putative novel': 0,
                                          'Possible fusion': 0, 'Putative novel isoform': 0}
        
        denovo_dict = find_novel(args.db_path, tx_dict, consensus_dict, ref_df, metrics, gene_biotype_map,
                                 args.denovo_num_introns, args.in_species_rna_support_only,
                                 args.denovo_tx_modes, args.denovo_splice_support, args.denovo_exon_support,
                                 args.denovo_ignore_novel_genes, args.denovo_novel_end_distance,
                                 args.denovo_allow_unsupported, args.denovo_allow_bad_annot_or_tm,
                                 args.denovo_only_novel_genes, args.denovo_allow_novel_ends)
        consensus_dict.update(denovo_dict)

    # perform final filtering steps
    deduplicated_consensus = deduplicate_consensus(consensus_dict, tx_dict, metrics)
    deduplicated_strand_resolved_consensus = resolve_opposite_strand(deduplicated_consensus, tx_dict, metrics)

    if any(m in args.denovo_tx_modes for m in ('augPB', 'strg')):
        deduplicated_strand_resolved_consensus = validate_pacbio_splices(deduplicated_strand_resolved_consensus,
                                                                         args.db_path, tx_dict, metrics,
                                                                         args.require_pacbio_support)

    if args.filter_overlapping_genes == True:
        gene_resolved_consensus = resolve_overlapping_cds_intervals(args.overlapping_ignore_bases,
                                                                    deduplicated_strand_resolved_consensus, tx_dict)
    else:
        gene_resolved_consensus = deduplicated_strand_resolved_consensus

    # sort by genomic interval for prettily increasing numbers
    final_consensus = sorted(gene_resolved_consensus,
                             key=lambda tx_attrs: (tx_dict[tx_attrs[0]].chromosome, tx_dict[tx_attrs[0]].start))

    # calculate final gene set completeness
    calculate_completeness(final_consensus, metrics)
    if 'augTM' in tx_modes or 'augTMR' in tx_modes:
        calculate_improvement_metrics(final_consensus, scored_df, tm_eval_df, support_df, metrics)
    calculate_indel_metrics(final_consensus, eval_df, metrics)
    consensus_gene_dict = write_consensus_gps(args.consensus_gp, args.consensus_gp_info,
                                              final_consensus, tx_dict, args.genome)
    write_consensus_gff3(consensus_gene_dict, args.consensus_gff3)
    write_consensus_fastas(consensus_gene_dict, args.consensus_fasta, args.protein_fasta, args.fasta)

    return metrics


def load_transmap_evals(db_path):
    tm_eval = tools.sqlInterface.load_alignment_evaluation(db_path)
    tm_filter_eval = tools.sqlInterface.load_filter_evaluation(db_path)
    tm_pw_filter_eval = tools.sqlInterface.load_pairwise_filter_evaluation(db_path)
    
    # Merge regular transMap filter eval with evaluation metrics
    tm_eval_df = pd.merge(tm_eval, tm_filter_eval, on=['TranscriptId', 'AlignmentId'], how='outer')
    
    # Also merge pairwise transMap filter eval if it exists
    if len(tm_pw_filter_eval) > 0:
        # Merge pairwise filter eval with evaluation metrics (may not have metrics for all pairwise alignments)
        tm_pw_merged = pd.merge(tm_eval, tm_pw_filter_eval, on=['TranscriptId', 'AlignmentId'], how='right')
        # Combine both regular and pairwise, keeping all unique combinations
        tm_eval_df = pd.concat([tm_eval_df, tm_pw_merged], ignore_index=True).drop_duplicates(subset=['TranscriptId', 'AlignmentId'], keep='first')
    
    return tm_eval_df.drop('AlignmentId', axis=1)


def calculate_vector_support(s, resolve_nan=None, num_digits=4):
    return 100 * tools.mathOps.format_ratio(len([x for x in s if x > 0]), len(s), resolve_nan=resolve_nan,
                                            num_digits=num_digits)




def load_metrics_from_db(db_path, tx_mode, aln_mode):
    session = tools.sqlInterface.start_session(db_path)
    metrics_table = tools.sqlInterface.tables[aln_mode][tx_mode]['metrics']
    metrics_df = tools.sqlInterface.load_metrics(metrics_table, session)
    # unstack flattens the long-form data structure
    metrics_df = metrics_df.set_index(['AlignmentId', 'classifier']).unstack('classifier')
    metrics_df.columns = [col[1] for col in metrics_df.columns]
    metrics_df = metrics_df.reset_index()
    numeric_cols = ['AlnCoverage', 'AlnGoodness', 'AlnIdentity', 'PercentUnknownBases']
    optional_cols = numeric_cols + [
        'AdjStart', 'AdjStop', 'OriginalIntrons', 'ProperOrf', 'ValidStart', 'ValidStop',
    ]
    for col in optional_cols:
        if col not in metrics_df.columns:
            metrics_df[col] = pd.NA
    metrics_df[numeric_cols] = metrics_df[numeric_cols].apply(pd.to_numeric, errors='coerce')
    metrics_df['OriginalIntrons'] = metrics_df['OriginalIntrons'].fillna('').astype(str)
    metrics_df['OriginalIntrons'] = [list(map(int, x)) if len(x[0]) > 0 else [] for x in
                                     metrics_df['OriginalIntrons'].str.split(',').tolist()]
    metrics_df['OriginalIntronsPercent'] = metrics_df['OriginalIntrons'].apply(calculate_vector_support, resolve_nan=1)
    session.close()
    return metrics_df


def load_evaluations_from_db(db_path, tx_mode):
    def aggfunc(s):
        if s.value_CDS.any():
            c = set(s[s.value_CDS > 0].name)
        else:
            c = set(s[s.value_mRNA > 0].name)
        cols = ['Frameshift', 'CodingInsertion', 'CodingDeletion', 'CodingMult3Indel']
        return pd.Series(('CodingDeletion' in c or 'CodingInsertion' in c,
                          'CodingInsertion' in c, 'CodingDeletion' in c,
                          'CodingMult3Deletion' in c or 'CodingMult3Insertion' in c), index=cols)

    session = tools.sqlInterface.start_session(db_path)
    cds_table = tools.sqlInterface.tables['CDS'][tx_mode]['evaluation']
    mrna_table = tools.sqlInterface.tables['mRNA'][tx_mode]['evaluation']
    cds_df = tools.sqlInterface.load_evaluation(cds_table, session)
    mrna_df = tools.sqlInterface.load_evaluation(mrna_table, session)
    cds_df = cds_df.set_index('AlignmentId')
    mrna_df = mrna_df.set_index('AlignmentId')
    merged = mrna_df.reset_index().merge(cds_df.reset_index(), how='outer', on=['AlignmentId', 'name'],
                                         suffixes=['_mRNA', '_CDS'])
    eval_df = merged.groupby('AlignmentId').apply(aggfunc, include_groups=False)
    return eval_df


def load_alt_names(db_path, denovo_tx_modes):
    session = tools.sqlInterface.start_session(db_path)
    r = []
    for tx_mode in denovo_tx_modes:
        try:
            table = tools.sqlInterface.tables['alt_names'][tx_mode]
            
            # Check if the table exists in the database
            inspector = inspect(session.bind)
            if inspector.has_table(table.__tablename__):
                alt_data = tools.sqlInterface.load_alternatives(table, session)
                if not alt_data.empty:
                    r.append(alt_data)
        except Exception as e:
            logger.warning(f"Error loading alternative names for {tx_mode}: {e}")
    
    if not r:
        # If no valid data, create an empty dataframe with required columns
        empty_df = pd.DataFrame(columns=['TranscriptId', 'AssignedGeneId', 'AlternativeGeneIds', 'ResolutionMethod'])
        # rename TranscriptId to AlignmentId
        empty_df.columns = [x if x != 'TranscriptId' else 'AlignmentId' for x in empty_df.columns]
        return empty_df
    
    df = pd.concat(r)
    df.columns = [x if x != 'TranscriptId' else 'AlignmentId' for x in df.columns]
    return df


def combine_and_filter_dfs(tx_dict, support_df, mrna_metrics_df, cds_metrics_df, tm_eval_df, ref_df, eval_df,
                           intron_rnaseq_support, exon_rnaseq_support, intron_annot_support, exon_annot_support,
                           original_intron_support, in_species_rna_support_only, denovo_tx_modes=None, alignment_source_map=None,
                           txTM_min_coverage=0.0):
    import logging
    logger = logging.getLogger(__name__)
    import time
    
    start_time = time.time()
    logger.info(f"Starting combine_and_filter_dfs with {len(support_df)} transcripts...")
    
    if alignment_source_map:
        txTM_alns = [a for a, s in alignment_source_map.items() if s == 'txTM']
        txTM_in_support = support_df[support_df['AlignmentId'].isin(txTM_alns)]
    
    logger.info("  Merging with reference data...")
    
    # Merge on versionless TranscriptId to avoid mismatches like ENST000001.1 vs ENST000001
    if 'TranscriptId' in ref_df.columns:
        ref_df = ref_df.copy()
        ref_df['TranscriptId_base'] = ref_df['TranscriptId'].astype(str).str.replace(r"\.[0-9]+$", "", regex=True)
    if 'TranscriptId' in support_df.columns:
        support_df = support_df.copy()
        support_df['TranscriptId_base'] = support_df['TranscriptId'].astype(str).str.replace(r"\.[0-9]+$", "", regex=True)
    
    
    support_ref_df = pd.merge(support_df, ref_df, left_on='TranscriptId_base', right_on='TranscriptId_base', how='left', suffixes=('', '_ref'))
    # Keep the original GeneId from support_df, but update it if it's missing
    if 'GeneId_ref' in support_ref_df.columns:
        support_ref_df['GeneId'] = support_ref_df['GeneId'].fillna(support_ref_df['GeneId_ref'])
        support_ref_df = support_ref_df.drop('GeneId_ref', axis=1)
    
    # Prefer reference GeneName if available
    if 'GeneName_ref' in support_ref_df.columns:
        support_ref_df['GeneName'] = support_ref_df['GeneName'].fillna(support_ref_df['GeneName_ref'])
        support_ref_df = support_ref_df.drop('GeneName_ref', axis=1)
    
    # Drop helper key
    if 'TranscriptId_base' in support_ref_df.columns:
        support_ref_df = support_ref_df.drop('TranscriptId_base', axis=1)

    missing_names = support_ref_df['GeneName'].isna().sum()
    if missing_names > 0:
        logger.warning(f"  Warning: {missing_names} transcripts have no GeneName after reference merge")
    
    if alignment_source_map:
        txTM_alns = [a for a, s in alignment_source_map.items() if s == 'txTM']
        txTM_in_support_ref = support_ref_df[support_ref_df['AlignmentId'].isin(txTM_alns)]
        txTM_pc_genes = set()
        import re
        for idx, row in txTM_in_support_ref.iterrows():
            if row.get('TranscriptBiotype') == 'protein_coding' or row.get('GeneBiotype') == 'protein_coding':
                gene_id = row.get('GeneId', 'UNKNOWN')
                base_gene_id = re.sub(r'_\d+$', '', gene_id)
                txTM_pc_genes.add(base_gene_id)
    
    logger.info(f"  After ref merge: {len(support_ref_df)} rows ({time.time() - start_time:.1f}s)")

    is_augpb = support_ref_df['AlignmentId'].apply(_is_denovo)
    augpb_missing_ref = support_ref_df['TranscriptBiotype'].isna() & is_augpb
    denovo_missing_ref = augpb_missing_ref
    
    if denovo_missing_ref.any():
        support_ref_df.loc[denovo_missing_ref, 'TranscriptBiotype'] = 'protein_coding'  # Default for de novo
        support_ref_df.loc[denovo_missing_ref, 'GeneBiotype'] = 'protein_coding'
        # Fill other required columns with defaults
        for col in ['TranscriptName', 'GeneName']:
            if col in support_ref_df.columns:
                support_ref_df.loc[denovo_missing_ref, col] = support_ref_df.loc[denovo_missing_ref, 'TranscriptId']
    
    support_ref_tm_df = pd.merge(support_ref_df, tm_eval_df, on=['GeneId', 'TranscriptId'], how='left')
    logger.info(f"  After tm_eval merge: {len(support_ref_tm_df)} rows ({time.time() - start_time:.1f}s)")

    is_augpb_eval = support_ref_tm_df['AlignmentId'].isin(support_ref_df.loc[is_augpb, 'AlignmentId'])
    augpb_missing_eval = is_augpb_eval & support_ref_tm_df[tm_eval_df.columns[-1]].isna()
    denovo_missing_eval = augpb_missing_eval
    
    if denovo_missing_eval.any():
        for col in tm_eval_df.columns:
            if col not in ['GeneId', 'TranscriptId'] and col in support_ref_tm_df.columns:
                if support_ref_tm_df[col].dtype == 'bool':
                    support_ref_tm_df.loc[denovo_missing_eval, col] = False
                elif support_ref_tm_df[col].dtype in ['int64', 'float64']:
                    support_ref_tm_df.loc[denovo_missing_eval, col] = 0
                else:
                    support_ref_tm_df.loc[denovo_missing_eval, col] = None
    
    if alignment_source_map:
        txTM_before_filter = support_ref_tm_df[support_ref_tm_df['AlignmentId'].isin(txTM_alns)]
        txTM_pc_before = set()
        for idx, row in txTM_before_filter.iterrows():
            if row.get('TranscriptBiotype') == 'protein_coding' or row.get('GeneBiotype') == 'protein_coding':
                gene_id = row.get('GeneId', 'UNKNOWN')
                base_gene_id = re.sub(r'_\d+$', '', gene_id)
                txTM_pc_before.add(base_gene_id)
    
    support_ref_tm_df = support_ref_tm_df[support_ref_tm_df.AlignmentId.isin(tx_dict.keys())]
    
    if alignment_source_map:
        txTM_after_filter = support_ref_tm_df[support_ref_tm_df['AlignmentId'].isin(txTM_alns)]
        txTM_pc_after = set()
        for idx, row in txTM_after_filter.iterrows():
            if row.get('TranscriptBiotype') == 'protein_coding' or row.get('GeneBiotype') == 'protein_coding':
                gene_id = row.get('GeneId', 'UNKNOWN')
                base_gene_id = re.sub(r'_\d+$', '', gene_id)
                txTM_pc_after.add(base_gene_id)
    
    logger.info(f"  After filtering: {len(support_ref_tm_df)} transcripts ({time.time() - start_time:.1f}s)")
    
    if alignment_source_map:
        txTM_mask = support_ref_tm_df['AlignmentId'].isin(txTM_alns)
        fixed_count = 0
        
        for idx in support_ref_tm_df[txTM_mask].index:
            aln_id = support_ref_tm_df.at[idx, 'AlignmentId']
            if aln_id in tx_dict:
                tx_obj = tx_dict[aln_id]
                if hasattr(tx_obj, 'cds_size') and tx_obj.cds_size > 0:
                    current_biotype = support_ref_tm_df.at[idx, 'TranscriptBiotype']
                    if pd.isna(current_biotype) or current_biotype != 'protein_coding':
                        support_ref_tm_df.at[idx, 'TranscriptBiotype'] = 'protein_coding'
                        support_ref_tm_df.at[idx, 'GeneBiotype'] = 'protein_coding'
                        fixed_count += 1
        
        logger.info(f"  Fixed biotype for {fixed_count} txTM transcripts with CDS (now protein_coding)")

    # augMP transcripts may not be present in the reference DB (miniprot templates / missing ref match),
    # which can leave TranscriptBiotype/GeneBiotype null and cause them to be treated as non-coding.
    # If the projected model has a CDS, treat it as protein_coding for filtering/scoring.
    augmp_fixed = 0
    augmp_mask = support_ref_tm_df['AlignmentId'].apply(tools.nameConversions.aln_id_is_augustus_mp)
    for idx in support_ref_tm_df[augmp_mask].index:
        aln_id = support_ref_tm_df.at[idx, 'AlignmentId']
        if aln_id in tx_dict:
            tx_obj = tx_dict[aln_id]
            if getattr(tx_obj, 'cds_size', 0) > 0:
                current_biotype = support_ref_tm_df.at[idx, 'TranscriptBiotype']
                if pd.isna(current_biotype) or current_biotype != 'protein_coding':
                    support_ref_tm_df.at[idx, 'TranscriptBiotype'] = 'protein_coding'
                    support_ref_tm_df.at[idx, 'GeneBiotype'] = 'protein_coding'
                    augmp_fixed += 1
    if augmp_fixed > 0:
        logger.info(f"  Fixed biotype for {augmp_fixed} augMP transcripts with CDS (now protein_coding)")
    
    coding_df = support_ref_tm_df[support_ref_tm_df.TranscriptBiotype == 'protein_coding']
    non_coding_df = support_ref_tm_df[support_ref_tm_df.TranscriptBiotype != 'protein_coding']
    
    if alignment_source_map:
        txTM_coding = coding_df[coding_df['AlignmentId'].isin(txTM_alns)]
        txTM_pc_after_split = set()
        import re
        for idx, row in txTM_coding.iterrows():
            gene_id = row.get('GeneId', 'UNKNOWN')
            base_gene_id = re.sub(r'_\d+$', '', gene_id)
            txTM_pc_after_split.add(base_gene_id)
    
    logger.info(f"  Split into {len(coding_df)} coding, {len(non_coding_df)} non-coding ({time.time() - start_time:.1f}s)")
    
    metrics_df = pd.merge(mrna_metrics_df, cds_metrics_df, on='AlignmentId', how='left', suffixes=['_mRNA', '_CDS'])
    
    # Use left join to preserve augPB transcripts that may not have metrics
    coding_df = pd.merge(coding_df, metrics_df, on='AlignmentId', how='left')
    
    if alignment_source_map:
        txTM_after_metrics = coding_df[coding_df['AlignmentId'].isin(txTM_alns)]
        txTM_pc_after_metrics = set()
        for idx, row in txTM_after_metrics.iterrows():
            gene_id = row.get('GeneId', 'UNKNOWN')
            base_gene_id = re.sub(r'_\d+$', '', gene_id)
            txTM_pc_after_metrics.add(base_gene_id)
    
    logger.info(f"  After metrics merge: {len(coding_df)} coding transcripts ({time.time() - start_time:.1f}s)")
    
    if alignment_source_map:
        txTM_aln_ids = [aln_id for aln_id, source in alignment_source_map.items() if source == 'txTM']
        txTM_cnv_alns = [aln_id for aln_id in txTM_aln_ids if '_' in aln_id and aln_id.split('_')[-1].isdigit()]
        txTM_cnv_mask = coding_df['AlignmentId'].isin(txTM_cnv_alns)
        
        if txTM_cnv_mask.any():
            logger.info(f"  Found {txTM_cnv_mask.sum()} txTM CNV transcripts with missing metrics")
            
            cnv_to_base = {}
            for aln_id in txTM_cnv_alns:
                base_id = '_'.join(aln_id.split('_')[:-1]) if '_' in aln_id else aln_id
                cnv_to_base[aln_id] = base_id
            
            for idx in coding_df[txTM_cnv_mask].index:
                cnv_aln_id = coding_df.loc[idx, 'AlignmentId']
                base_aln_id = cnv_to_base.get(cnv_aln_id)
                
                if base_aln_id and base_aln_id in metrics_df['AlignmentId'].values:
                    # Get metrics from base transcript
                    base_metrics_row = metrics_df[metrics_df['AlignmentId'] == base_aln_id]
                    
                    if len(base_metrics_row) > 0:
                        base_metrics = base_metrics_row.iloc[0]
                        
                        # Copy metrics to CNV transcript (only if currently NaN)
                        for col in base_metrics.index:
                            if col != 'AlignmentId' and col in coding_df.columns:
                                try:
                                    cnv_value = coding_df.at[idx, col]
                                    # Skip if CNV already has a value
                                    if not pd.isna(cnv_value):
                                        continue
                                    
                                    # Get base value
                                    base_value = base_metrics[col]
                                    
                                    # Handle different value types
                                    import numpy as np
                                    if isinstance(base_value, (list, tuple, np.ndarray)):
                                        # For list/array columns (like support vectors), copy as-is
                                        coding_df.at[idx, col] = base_value
                                    elif hasattr(base_value, 'item'):
                                        # For numpy scalars, extract the Python value
                                        coding_df.at[idx, col] = base_value.item()
                                    else:
                                        # For regular scalars
                                        coding_df.at[idx, col] = base_value
                                except Exception as e:
                                    # Skip problematic columns
                                    continue
            
    
    augpb_mask = coding_df['AlignmentId'].isin(support_ref_df.loc[is_augpb, 'AlignmentId'])
    
    if alignment_source_map:
        all_txTM_alns = [aln_id for aln_id, source in alignment_source_map.items() if source == 'txTM']
        txTM_mask = coding_df['AlignmentId'].isin(all_txTM_alns)
        
        has_coverage = coding_df['AlnCoverage_mRNA'].notna() & (coding_df['AlnCoverage_mRNA'] > 0)
        txTM_missing_metrics = txTM_mask & ~has_coverage
    else:
        txTM_missing_metrics = pd.Series([False] * len(coding_df), index=coding_df.index)
    
    # augMP now has REAL PSL-derived metrics from miniprot (see Snakefile rules
    # ``miniprot_paf_to_genepred`` → ``generate_augMP_psl`` → ``store_psl_metrics_augMP``).
    # We deliberately do NOT include augMP in the "fill missing metrics with 100"
    # shortcut here: augMP records that lack a miniprot PSL row (extremely rare)
    # should stay NaN and be filtered honestly downstream, rather than being
    # falsely promoted to perfect coverage / identity.
    fill_mask = augpb_mask | txTM_missing_metrics

    if fill_mask.any():
        logger.info(
            f"  Filling missing metrics for {fill_mask.sum()} transcripts "
            f"({augpb_mask.sum()} augPB, {txTM_missing_metrics.sum()} txTM)"
        )
        for col in ['AlnCoverage_mRNA', 'AlnIdentity_mRNA', 'AlnCoverage_CDS', 'AlnIdentity_CDS']:
            if col in coding_df.columns:
                coding_df.loc[fill_mask, col] = coding_df.loc[fill_mask, col].fillna(100.0).astype(float)
        for col in ['AlnGoodness_mRNA', 'AlnGoodness_CDS']:
            if col in coding_df.columns:
                coding_df.loc[fill_mask, col] = coding_df.loc[fill_mask, col].fillna(1.0).astype(float)
        for col in ['OriginalIntronsPercent_mRNA', 'OriginalIntronsPercent_CDS']:
            if col in coding_df.columns:
                coding_df.loc[fill_mask, col] = coding_df.loc[fill_mask, col].fillna(100.0).astype(float)
    
    coding_df = pd.merge(coding_df, eval_df, on='AlignmentId', how='left')
    logger.info(f"  After eval merge: {len(coding_df)} coding transcripts ({time.time() - start_time:.1f}s)")
    
    coding_df = coding_df.copy()
    non_coding_df = non_coding_df.copy()
    
    coding_df['OriginalIntronsPercent_mRNA'] = coding_df['OriginalIntronsPercent_mRNA'].fillna(100).astype(float)
    coding_df['OriginalIntronsPercent_CDS'] = coding_df['OriginalIntronsPercent_CDS'].fillna(100).astype(float)
    non_coding_df['TransMapOriginalIntronsPercent'] = non_coding_df['TransMapOriginalIntronsPercent'].fillna(100).astype(float)
    coding_df['Frameshift'] = coding_df['Frameshift'].fillna(False).astype(bool)
    
    if alignment_source_map:
        txTM_aln_ids = [aln_id for aln_id, source in alignment_source_map.items() if source == 'txTM']
        is_txTM_coding = coding_df['AlignmentId'].isin(txTM_aln_ids)
        txTM_coding_df = coding_df[is_txTM_coding].copy()
        non_txTM_coding_df = coding_df[~is_txTM_coding].copy()
        
        logger.info(f"  Filtering: {len(txTM_coding_df)} txTM coding transcripts, {len(non_txTM_coding_df)} non-txTM coding transcripts (with filters)")
        if len(txTM_coding_df) > 0:
            # Filter txTM transcripts whose coverage was explicitly computed as <= txTM_min_coverage.
            # NaN coverage (no metrics in DB) is treated as trusted (fillna gives 100.0 above).
            # Explicit 0 coverage means BLAT found no alignment matches → likely a gene fragment.
            txTM_cov_filter = txTM_coding_df['AlnCoverage_mRNA'].fillna(100.0) > txTM_min_coverage
            n_before = len(txTM_coding_df)
            txTM_coding_df = txTM_coding_df[txTM_cov_filter]
            n_filtered = n_before - len(txTM_coding_df)
            if n_filtered > 0:
                logger.info(f"  Filtered {n_filtered} txTM coding transcripts with coverage <= {txTM_min_coverage}% (likely gene fragments)")
    else:
        txTM_coding_df = pd.DataFrame()
        non_txTM_coding_df = coding_df.copy()

    if in_species_rna_support_only == True:
        filt = ((non_txTM_coding_df.OriginalIntronsPercent_mRNA >= original_intron_support) &
                (non_txTM_coding_df.IntronAnnotSupportPercent >= intron_annot_support) &
                (non_txTM_coding_df.IntronRnaSupportPercent >= intron_rnaseq_support) &
                (non_txTM_coding_df.ExonAnnotSupportPercent >= exon_annot_support) &
                (non_txTM_coding_df.ExonRnaSupportPercent >= exon_rnaseq_support))
    else:
        filt = ((non_txTM_coding_df.OriginalIntronsPercent_mRNA >= original_intron_support) &
                (non_txTM_coding_df.IntronAnnotSupportPercent >= intron_annot_support) &
                (non_txTM_coding_df.AllSpeciesIntronRnaSupportPercent >= intron_rnaseq_support) &
                (non_txTM_coding_df.ExonAnnotSupportPercent >= exon_annot_support) &
                (non_txTM_coding_df.AllSpeciesExonRnaSupportPercent >= exon_rnaseq_support))

    # augMP does not have evaluation/metrics rows in the DB, so support columns can be 0/NaN even
    # for real loci. Allow augMP coding transcripts through this filter; downstream conflict
    # resolution + scoring still determine whether they survive in the final gene set.
    augmp_non_txTM = non_txTM_coding_df['AlignmentId'].apply(tools.nameConversions.aln_id_is_augustus_mp)
    filt = filt | augmp_non_txTM
    non_txTM_coding_df = non_txTM_coding_df[filt]
    
    coding_df = pd.concat([txTM_coding_df, non_txTM_coding_df], ignore_index=True)
    
    if alignment_source_map:
        is_txTM_noncoding = non_coding_df['AlignmentId'].isin(txTM_aln_ids)
        txTM_noncoding_df = non_coding_df[is_txTM_noncoding].copy()
        non_txTM_noncoding_df = non_coding_df[~is_txTM_noncoding].copy()
        if len(txTM_noncoding_df) > 0 and 'AlnCoverage_mRNA' in txTM_noncoding_df.columns:
            txTM_nc_cov_filter = txTM_noncoding_df['AlnCoverage_mRNA'].fillna(100.0) > txTM_min_coverage
            n_before = len(txTM_noncoding_df)
            txTM_noncoding_df = txTM_noncoding_df[txTM_nc_cov_filter]
            n_filtered = n_before - len(txTM_noncoding_df)
            if n_filtered > 0:
                logger.info(f"  Filtered {n_filtered} txTM non-coding transcripts with coverage <= {txTM_min_coverage}%")
    else:
        txTM_noncoding_df = pd.DataFrame()
        non_txTM_noncoding_df = non_coding_df.copy()

    if in_species_rna_support_only == True:
        filt = ((non_txTM_noncoding_df.TransMapOriginalIntronsPercent >= original_intron_support) &
                (non_txTM_noncoding_df.IntronAnnotSupportPercent >= intron_annot_support) &
                (non_txTM_noncoding_df.IntronRnaSupportPercent >= intron_rnaseq_support) &
                (non_txTM_noncoding_df.ExonAnnotSupportPercent >= exon_annot_support) &
                (non_txTM_noncoding_df.ExonRnaSupportPercent >= exon_rnaseq_support))
    else:
        filt = ((non_txTM_noncoding_df.TransMapOriginalIntronsPercent >= original_intron_support) &
                (non_txTM_noncoding_df.IntronAnnotSupportPercent >= intron_annot_support) &
                (non_txTM_noncoding_df.AllSpeciesIntronRnaSupportPercent >= intron_rnaseq_support) &
                (non_txTM_noncoding_df.ExonAnnotSupportPercent >= exon_annot_support) &
                (non_txTM_noncoding_df.AllSpeciesExonRnaSupportPercent >= exon_rnaseq_support))
    non_txTM_noncoding_df = non_txTM_noncoding_df[filt]
    
    non_coding_df = pd.concat([txTM_noncoding_df, non_txTM_noncoding_df], ignore_index=True)

    logger.info(f"  After filtering: {len(coding_df)} coding, {len(non_coding_df)} non-coding transcripts")
    logger.info(f"  combine_and_filter_dfs completed in {time.time() - start_time:.1f}s")
    return coding_df, non_coding_df


def score_filtered_dfs(coding_df, non_coding_df, in_species_rna_support_only):
    import logging
    logger = logging.getLogger(__name__)
    import time
    
    start_time = time.time()
    
    coding_df_copy = coding_df.copy() if len(coding_df) > 0 else coding_df
    non_coding_df_copy = non_coding_df.copy() if len(non_coding_df) > 0 else non_coding_df
    
    if len(coding_df_copy) > 0:
        # For coding transcripts
        aln_id = coding_df_copy['AlnIdentity_CDS'].fillna(0)
        aln_cov = coding_df_copy['AlnCoverage_CDS'].fillna(0)
        orig_intron = coding_df_copy['OriginalIntronsPercent_mRNA'].fillna(0)
        if in_species_rna_support_only:
            rna_support = coding_df_copy['ExonRnaSupportPercent'] + coding_df_copy['IntronRnaSupportPercent']
        else:
            rna_support = coding_df_copy['AllSpeciesExonRnaSupportPercent'] + coding_df_copy['AllSpeciesIntronRnaSupportPercent']
        coding_df_copy['TranscriptScore'] = (aln_id + aln_cov + 
                                              coding_df_copy['IntronAnnotSupportPercent'] + 
                                              coding_df_copy['ExonAnnotSupportPercent'] + 
                                              orig_intron + rna_support)
    
    if len(non_coding_df_copy) > 0:
        aln_id = non_coding_df_copy['TransMapIdentity'].fillna(0)
        aln_cov = non_coding_df_copy['TransMapCoverage'].fillna(0)
        orig_intron = non_coding_df_copy['TransMapOriginalIntronsPercent'].fillna(0)
        if in_species_rna_support_only:
            rna_support = non_coding_df_copy['ExonRnaSupportPercent'] + non_coding_df_copy['IntronRnaSupportPercent']
        else:
            rna_support = non_coding_df_copy['AllSpeciesExonRnaSupportPercent'] + non_coding_df_copy['AllSpeciesIntronRnaSupportPercent']
        non_coding_df_copy['TranscriptScore'] = (aln_id + aln_cov + 
                                                  non_coding_df_copy['IntronAnnotSupportPercent'] + 
                                                  non_coding_df_copy['ExonAnnotSupportPercent'] + 
                                                  orig_intron + rna_support)
    
    logger.info(f"  score_filtered_dfs completed in {time.time() - start_time:.1f}s")
    return coding_df_copy, non_coding_df_copy


def merge_scored_dfs(scored_coding_df, scored_non_coding_df):
    scored_non_coding_df = scored_non_coding_df.copy()
    for m in ['Coverage', 'Identity', 'Goodness']:
        scored_non_coding_df['Aln' + m + '_mRNA'] = scored_non_coding_df['TransMap' + m]
    merged_df = pd.concat([scored_non_coding_df, scored_coding_df])
    return merged_df


def validate_pacbio_splices(deduplicated_strand_resolved_consensus, db_path, tx_dict, metrics, require_pacbio_support):
    """
    Tag transcripts as having PacBio support.
    If users passed the --require-pacbio-support, remove any transcript which does not have support.
    """
    iso_txs = tools.sqlInterface.load_isoseq_txs(db_path)
    tx_ids, _ = list(zip(*deduplicated_strand_resolved_consensus))
    txs = [tx_dict[tx_id] for tx_id in tx_ids]
    
    if not iso_txs:
        pb_resolved_consensus = []
        for tx_id, d in deduplicated_strand_resolved_consensus:
            if _is_denovo(tx_id):
                d['pacbio_isoform_supported'] = True
                metrics['IsoSeq Transcript Validation'][True] += 1
                pb_resolved_consensus.append([tx_id, d])
            elif require_pacbio_support == False:
                d['pacbio_isoform_supported'] = False
                metrics['IsoSeq Transcript Validation'][False] += 1
                pb_resolved_consensus.append([tx_id, d])
        return pb_resolved_consensus
    
    clustered = tools.transcripts.cluster_txs(txs + iso_txs)
    divided_clusters = tools.transcripts.divide_clusters(clustered, tx_ids)
    subset_matches = tools.transcripts.calculate_subset_matches(divided_clusters)
    # invert the subset_matches to extract all validated tx_ids
    validated_ids = set()
    for tx_list in subset_matches.values():
        for tx in tx_list:
            validated_ids.add(tx.name)
    # begin resolving
    pb_resolved_consensus = []
    for tx_id, d in deduplicated_strand_resolved_consensus:
        if tx_id in validated_ids or _is_denovo(tx_id):
            d['pacbio_isoform_supported'] = True
            metrics['IsoSeq Transcript Validation'][True] += 1
            pb_resolved_consensus.append([tx_id, d])
        elif require_pacbio_support == False:
            d['pacbio_isoform_supported'] = False
            metrics['IsoSeq Transcript Validation'][False] += 1
            pb_resolved_consensus.append([tx_id, d])
    return pb_resolved_consensus


def incorporate_tx(best_rows, gene_id, metrics, hints_db_has_rnaseq):
    best_series = best_rows.iloc[0]
    is_augpb = _is_denovo(best_series.AlignmentId)
    is_denovo = is_augpb
    
    clean_gene_id = gene_id
    if gene_id and gene_id.endswith('_cp'):
        clean_gene_id = gene_id[:-3]
    
    d = {'source_transcript': best_series.TranscriptId if not is_denovo else 'N/A',
         'source_transcript_name': best_series.TranscriptName,
         'source_gene': clean_gene_id if not is_denovo else None,
         'score': int(10 * round(best_series.AlnGoodness_mRNA, 3)) if pd.notna(best_series.AlnGoodness_mRNA) else 100,
         'gene_biotype': best_series.GeneBiotype if not is_denovo else 'unknown_likely_coding',
         'transcript_biotype': best_series.TranscriptBiotype if not is_denovo else 'unknown_likely_coding',
         'alignment_id': str(best_series.AlignmentId),
         'frameshift': str(best_series.get('Frameshift', None)) if pd.notna(best_series.get('Frameshift', None)) else 'N/A',
         'exon_annotation_support': ','.join(map(str, best_series.ExonAnnotSupport)),
         'intron_annotation_support': ','.join(map(str, best_series.IntronAnnotSupport)),
         'transcript_class': 'putative_novel' if is_denovo else 'ortholog',
         'valid_start': bool(best_series.ValidStart) if pd.notna(best_series.ValidStart) else False,
         'valid_stop': bool(best_series.ValidStop) if pd.notna(best_series.ValidStop) else False,
         'adj_start': best_series.AdjStart_mRNA if pd.notna(best_series.AdjStart_mRNA) else 0,
         'adj_stop': best_series.AdjStop_mRNA if pd.notna(best_series.AdjStop_mRNA) else 0,
         'proper_orf': bool(best_series.ProperOrf) if pd.notna(best_series.ProperOrf) else False,
         'source_type': best_series.get('SourceType', 'unknown') if 'SourceType' in best_series else 'unknown',
         'AdditionalSources': best_series.get('AdditionalSources', '') if 'AdditionalSources' in best_series else ''}
    extra_tags = best_series.ExtraTags
    if pd.notna(extra_tags) and isinstance(extra_tags, str):
        for key, val in tools.misc.parse_gff_attr_line(extra_tags).items():
            d[key] = val
    
    if is_denovo:
        d['novel_5p_cap'] = True
        d['novel_poly_a'] = True  
        d['pacbio_isoform_supported'] = True  # augPB transcripts are by definition PacBio-supported
    
    if _is_denovo(best_series.AlignmentId):
        d['pacbio_isoform_supported'] = True
    
    if hasattr(best_series, 'ExonRnaSupport') and best_series.ExonRnaSupport is not None:
        d['exon_rna_support'] = ','.join(map(str, best_series.ExonRnaSupport))
    if hasattr(best_series, 'IntronRnaSupport') and best_series.IntronRnaSupport is not None:
        d['intron_rna_support'] = ','.join(map(str, best_series.IntronRnaSupport))
    if best_series.Paralogy is not None:
        d['paralogy'] = best_series.Paralogy
    if best_series.UnfilteredParalogy is not None:
        d['unfiltered_paralogy'] = best_series.UnfilteredParalogy
    if best_series.GeneAlternateLoci is not None:
        d['gene_alternate_contigs'] = best_series.GeneAlternateLoci
    if best_series.CollapsedGeneIds is not None:
        d['collapsed_gene_ids'] = best_series.CollapsedGeneIds
    if best_series.CollapsedGeneNames is not None:
        d['collapsed_gene_names'] = best_series.CollapsedGeneNames
    if best_series.PossibleSplitGeneLocations is not None:
        d['possible_split_gene_locations'] = best_series.PossibleSplitGeneLocations
    if best_series.GeneName is not None:
        if _is_denovo(best_series.AlignmentId):
            d['source_gene_common_name'] = None
        else:
            d['source_gene_common_name'] = best_series.GeneName
    # add information to the overall metrics
    coverage = best_series.AlnCoverage_mRNA if pd.notna(best_series.AlnCoverage_mRNA) else 0
    identity = best_series.AlnIdentity_mRNA if pd.notna(best_series.AlnIdentity_mRNA) else 0
    intron_rna_support = best_series.IntronRnaSupportPercent if pd.notna(best_series.IntronRnaSupportPercent) else 0
    exon_rna_support = best_series.ExonRnaSupportPercent if pd.notna(best_series.ExonRnaSupportPercent) else 0
    intron_annot_support = best_series.IntronAnnotSupportPercent if pd.notna(best_series.IntronAnnotSupportPercent) else 0
    exon_annot_support = best_series.ExonAnnotSupportPercent if pd.notna(best_series.ExonAnnotSupportPercent) else 0
    original_introns = best_series.OriginalIntronsPercent_mRNA if pd.notna(best_series.OriginalIntronsPercent_mRNA) else 0
    
    metrics['Coverage'][best_series.TranscriptBiotype].append(coverage)
    metrics['Identity'][best_series.TranscriptBiotype].append(identity)
    metrics['Splice Support'][best_series.TranscriptBiotype].append(intron_rna_support)
    metrics['Exon Support'][best_series.TranscriptBiotype].append(exon_rna_support)
    metrics['Splice Annotation Support'][best_series.TranscriptBiotype].append(intron_annot_support)
    metrics['Exon Annotation Support'][best_series.TranscriptBiotype].append(exon_annot_support)
    metrics['Original Introns'][best_series.TranscriptBiotype].append(original_introns)
    
    
    return best_series.AlignmentId, d


def evaluate_ties(best_rows):
    return ','.join(sorted(set([tools.nameConversions.alignment_type(x) for x in best_rows.AlignmentId])))


def find_novel(db_path, tx_dict, consensus_dict, ref_df, metrics, gene_biotype_map, denovo_num_introns,
               in_species_rna_support_only, denovo_tx_modes, denovo_splice_support, denovo_exon_support,
               denovo_ignore_novel_genes, denovo_novel_end_distance, denovo_allow_unsupported,
               denovo_allow_bad_annot_or_tm, denovo_only_novel_genes, denovo_allow_novel_ends):
    existing_splices = set()
    existing_5p = collections.defaultdict(set)
    existing_3p = collections.defaultdict(set)
    
    # Extract all splices, 5' and 3' ends we have already seen
    for consensus_tx in consensus_dict:
        tx_obj = tx_dict[consensus_tx]
        existing_splices.update(tx_obj.intron_intervals)
        existing_5p[tx_obj.chromosome].add(tx_obj.get_5p_interval())
        existing_3p[tx_obj.chromosome].add(tx_obj.get_3p_interval())
    
    existing_5p_sorted = {}
    existing_3p_sorted = {}
    existing_5p_starts = {}  # Pre-computed start positions for binary search
    existing_3p_starts = {}
    for chrom in existing_5p:
        existing_5p_sorted[chrom] = sorted(existing_5p[chrom], key=lambda i: i.start)
        existing_3p_sorted[chrom] = sorted(existing_3p[chrom], key=lambda i: i.start)
        existing_5p_starts[chrom] = [i.start for i in existing_5p_sorted[chrom]]
        existing_3p_starts[chrom] = [i.start for i in existing_3p_sorted[chrom]]
    
    def is_novel(s):
        """
        Determine if this transcript is possibly novel. If it is assigned a gene ID, pass this off to
         is_novel_supported()
        """
        try:
            if s.AssignedGeneId is not None:
                result = is_novel_supported(s)
                return result if result is not None else None
            if denovo_allow_bad_annot_or_tm == False and s.ResolutionMethod == 'badAnnotOrTm':
                return None
            elif s.ResolutionMethod == 'ambiguousOrFusion' and s.IntronRnaSupportPercent != 100:
                return None
            # validate the support level
            intron = s.IntronRnaSupportPercent if in_species_rna_support_only else s.AllSpeciesIntronRnaSupportPercent
            exon = s.ExonRnaSupportPercent if in_species_rna_support_only else s.AllSpeciesExonRnaSupportPercent
            if intron < denovo_splice_support or exon < denovo_exon_support:
                return None
            # if we previously flagged this as ambiguousOrFusion, propagate this tag
            if s.ResolutionMethod == 'ambiguousOrFusion':
                return 'possible_fusion'
            elif s.ResolutionMethod == 'badAnnotOrTm':
                return 'bad_annot_or_tm'
            # if we have alternatives, this is not novel but could be a gene family expansion
            elif s.AlternativeGeneIds is not None:
                return 'possible_paralog'
            # this may be a poor mapping
            elif bool(s.ExonAnnotSupportPercent > 0 or s.CdsAnnotSupportPercent > 0 or s.IntronAnnotSupportPercent > 0):
                return 'poor_alignment'
            # this is looking pretty novel, could still be a mapping problem in a complex region though
            else:
                return 'putative_novel'
        except Exception as e:
            # In case of any error, return None and continue
            logger.warning(f"Error in is_novel for alignment {s.AlignmentId}: {e}")
            return None

    def is_novel_supported(s):
        """Is this PB transcript with an assigned gene ID supported and have a novel splice?"""
        denovo_tx_obj = tx_dict[s.AlignmentId]
        
        if len(denovo_tx_obj.intron_intervals) < denovo_num_introns:
            return None
        elif in_species_rna_support_only and s.ExonRnaSupportPercent <= denovo_exon_support or \
                        s.IntronRnaSupportPercent <= denovo_splice_support:
            return None
        elif in_species_rna_support_only == False and s.AllSpeciesExonRnaSupportPercent <= denovo_exon_support or \
                        s.AllSpeciesIntronRnaSupportPercent <= denovo_splice_support:
            return None
        
        # look for splices that are not supported by the reference annotation
        # these splices may or may not be supported by RNA-seq based on the denovo_allow_unsupported flag
        new_supported_splices = set()
        intron_vector = s.IntronRnaSupport if in_species_rna_support_only else s.AllSpeciesIntronRnaSupport
        
        for intron, rna in zip(*[denovo_tx_obj.intron_intervals, intron_vector]):
            if intron not in existing_splices:
                if denovo_allow_unsupported == True or rna > 0:
                    new_supported_splices.add(intron)
        
        if len(new_supported_splices) == 0:
            # For augPB transcripts with assigned gene IDs, even if they don't have novel splices,
            # they represent novel isoform structures supported by PacBio data
            if _is_denovo(s.AlignmentId) and s.AssignedGeneId is not None:
                return 'putative_novel_isoform'
            return None
        
        # if any splices are both not supported by annotation and supported by RNA, call this as novel
        novel_and_supported = []
        for i, annot in zip(*[denovo_tx_obj.intron_intervals, s.IntronAnnotSupport]):
            if annot == 0 and i in new_supported_splices:
                novel_and_supported.append(i)
        
        if any(annot == 0 and i in new_supported_splices for i, annot in zip(*[denovo_tx_obj.intron_intervals,
                                                                               s.IntronAnnotSupport])):
            metrics['Transcript Modes'][tools.nameConversions.alignment_type(s.AlignmentId)] += 1
            tx_class = 'putative_novel_isoform'
        # if any splices are new, and supported by RNA-seq call this poor alignment
        else:
            tx_class = 'poor_alignment'
        return str(tx_class)  # Ensure we return a string, not a Series

    def has_novel_ends(s):
        denovo_tx_obj = tx_dict[s.AlignmentId]
        five_p = denovo_tx_obj.get_5p_interval()
        three_p = denovo_tx_obj.get_3p_interval()
        
        chrom = denovo_tx_obj.chromosome
        five_p_matches = interval_not_within_wiggle_room_sorted(
            existing_5p_sorted.get(chrom, []),
            existing_5p_starts.get(chrom, []),
            five_p, denovo_novel_end_distance)
        three_p_matches = interval_not_within_wiggle_room_sorted(
            existing_3p_sorted.get(chrom, []),
            existing_3p_starts.get(chrom, []),
            three_p, denovo_novel_end_distance)
        return pd.Series([five_p_matches, three_p_matches])
    
    def interval_not_within_wiggle_room_sorted(sorted_intervals, starts_list, query_interval, wiggle_room):
        if not sorted_intervals:
            return True
        
        query_start = query_interval.start
        query_stop = query_interval.stop
        search_margin = 2 * wiggle_room
        
        # Find the first interval that could overlap with our search range using pre-computed starts
        import bisect
        left_idx = bisect.bisect_left(starts_list, query_start - search_margin)
        
        # Only check intervals that could be within range (limited window)
        for i in range(max(0, left_idx - 5), min(len(sorted_intervals), left_idx + 50)):
            target_interval = sorted_intervals[i]
            
            # Early exit: if this interval starts way after our query, no more matches possible
            if target_interval.start > query_stop + search_margin:
                break
            
            # Check symmetric separation
            try:
                separation = sum(query_interval.symmetric_separation(target_interval))
                if separation <= search_margin:
                    return False  # Found a match, not novel
            except (TypeError, AttributeError):
                continue
        
        return True  # No match found, is novel

    # Create minimal support dataframe for de novo transcripts
    consensus_augpb = [aln_id for aln_id in consensus_dict.keys() if _is_denovo(aln_id)]
    
    if consensus_augpb and any(m in denovo_tx_modes for m in ('augPB', 'strg')):
        # Create a minimal dataframe with necessary columns for Augustus PB transcripts
        augpb_rows = []
        for aln_id in consensus_augpb:
            tx_obj = tx_dict[aln_id]
            # For augPB transcripts: high RNA support (PacBio-guided), low annotation support (novel)
            num_introns = len(tx_obj.intron_intervals)
            num_exons = len(tx_obj.exon_frames)
            
            augpb_rows.append({
                'AlignmentId': aln_id,
                # Annotation support - low for novel features
                'IntronAnnotSupport': [0] * num_introns,
                'ExonAnnotSupport': [0] * num_exons,
                'CdsAnnotSupport': [0] * num_exons,
                # RNA support - high for PacBio-supported transcripts
                'ExonRnaSupport': [1] * num_exons,
                'IntronRnaSupport': [1] * num_introns,
                'AllSpeciesExonRnaSupport': [1] * num_exons,
                'AllSpeciesIntronRnaSupport': [1] * num_introns
            })
        
        denovo_support_df = pd.DataFrame(augpb_rows)
        
        # Calculate support percentages (will be 0 for dummy data)
        cols = ['IntronAnnotSupport', 'ExonAnnotSupport', 'CdsAnnotSupport',
                'ExonRnaSupport', 'IntronRnaSupport',
                'AllSpeciesExonRnaSupport', 'AllSpeciesIntronRnaSupport']
        for col in cols:
            denovo_support_df[col + 'Percent'] = denovo_support_df[col].apply(calculate_vector_support, resolve_nan=1)
    else:
        return {}
    
    # remove the TranscriptId and GeneId columns so they can be populated by others
    if 'TranscriptId' in denovo_support_df.columns:
        denovo_support_df = denovo_support_df.drop(['GeneId', 'TranscriptId'], axis=1)
    
    # load the alignment metrics data
    denovo_alt_names = load_alt_names(db_path, denovo_tx_modes)
    
    # Merge with alternative names data (left join to preserve all support data)
    if len(denovo_alt_names) > 0:
        denovo_df = pd.merge(denovo_support_df, denovo_alt_names, on='AlignmentId', how='left')
    else:
        # If no alternative names data, add empty columns
        denovo_df = denovo_support_df.copy()
        denovo_df['AssignedGeneId'] = None
        denovo_df['AlternativeGeneIds'] = None
        denovo_df['ResolutionMethod'] = None
    
    # Assign gene IDs to augPB transcripts based on overlap with consensus transcripts
    def assign_gene_by_overlap(denovo_df, consensus_dict, tx_dict):
        start_time = time.time()
        
        assigned_count = 0
        consensus_by_chrom_strand = collections.defaultdict(list)
        
        for consensus_tx_id, consensus_info in consensus_dict.items():
            if consensus_tx_id not in tx_dict:
                continue
            
            consensus_tx = tx_dict[consensus_tx_id]
            gene_id = consensus_info.get('source_gene')
            
            if gene_id:
                # Pre-compute bounds for each consensus transcript
                start = min(e.start for e in consensus_tx.exon_intervals)
                end = max(e.stop for e in consensus_tx.exon_intervals)
                key = (consensus_tx.chromosome, consensus_tx.strand)
                consensus_by_chrom_strand[key].append((start, end, gene_id))
        
        # Sort each group by start position for efficient searching
        for key in consensus_by_chrom_strand:
            consensus_by_chrom_strand[key].sort()
        
        
        augpb_assignments = []
        
        for idx in range(len(denovo_df)):
            row = denovo_df.iloc[idx]
            
            # Skip if already has assigned gene ID
            if pd.notna(row.get('AssignedGeneId')):
                continue
                
            aln_id = row['AlignmentId']
            
            # Focus on augPB transcripts that need gene assignment
            if not _is_denovo(str(aln_id)):
                continue
                
            if aln_id not in tx_dict:
                continue
                
            denovo_tx = tx_dict[aln_id]
            denovo_start = min(e.start for e in denovo_tx.exon_intervals)
            denovo_end = max(e.stop for e in denovo_tx.exon_intervals)
            denovo_length = denovo_end - denovo_start
            
            # Get candidate consensus transcripts from spatial index
            key = (denovo_tx.chromosome, denovo_tx.strand)
            candidates = consensus_by_chrom_strand.get(key, [])
            
            # Find best overlapping gene
            best_overlap = 0
            best_gene_id = None
            
            for cons_start, cons_end, gene_id in candidates:
                # Early exit: if consensus starts after denovo ends, no more overlaps
                if cons_start >= denovo_end:
                    break
                
                # Skip if consensus ends before denovo starts
                if cons_end <= denovo_start:
                    continue
                
                # Calculate overlap
                overlap_start = max(denovo_start, cons_start)
                overlap_end = min(denovo_end, cons_end)
                
                if overlap_start < overlap_end:
                    overlap_length = overlap_end - overlap_start
                    overlap_fraction = overlap_length / denovo_length if denovo_length > 0 else 0
                    
                    if overlap_fraction > best_overlap:
                        best_overlap = overlap_fraction
                        best_gene_id = gene_id
            
            # Collect assignment if found with sufficient overlap
            if best_gene_id and best_overlap > 0.1:  # Require at least 10% overlap
                augpb_assignments.append((idx, best_gene_id))
                assigned_count += 1
        
        for idx, gene_id in augpb_assignments:
            denovo_df.at[idx, 'AssignedGeneId'] = gene_id
            denovo_df.at[idx, 'ResolutionMethod'] = 'overlap_assignment'
        
        elapsed = time.time() - start_time
        
        return denovo_df
    
    # Apply gene assignment by overlap
    denovo_df = assign_gene_by_overlap(denovo_df, consensus_dict, tx_dict)
    
    common_name_map = dict(list(zip(ref_df.GeneId, ref_df.GeneName)))

    denovo_df['CommonName'] = [common_name_map.get(x, None) for x in denovo_df.AssignedGeneId]
    denovo_df['GeneBiotype'] = [gene_biotype_map.get(x, None) for x in denovo_df.AssignedGeneId]

    # if we have an external reference, try to incorporate those names as well
    if 'exRef' in denovo_tx_modes:
        def add_exref_ids(s):
            if s.AlignmentId in exref_common_name_map:
                # if we have an assigned gene ID, defer the gene biotype to that but retain transcript biotype
                if s.AssignedGeneId is None:
                    return pd.Series([exref_common_name_map[s.AlignmentId], exref_gene_biotype_map[s.AlignmentId]])
                else:
                    return pd.Series([exref_common_name_map[s.AlignmentId], s.GeneBiotype])
            # pass along the original data
            return pd.Series([s.CommonName, s.GeneBiotype])
        exref_annot = tools.sqlInterface.load_annotation(db_path)
        exref_common_name_map = dict(list(zip(exref_annot.TranscriptId, exref_annot.GeneName)))
        exref_gene_biotype_map = dict(list(zip(exref_annot.TranscriptId, exref_annot.GeneBiotype)))
        denovo_df[['CommonName', 'GeneBiotype']] = denovo_df.apply(add_exref_ids, axis=1)
        exref_annot = exref_annot.set_index('TranscriptId')

    # apply the novel finding functions
    start_time = time.time()
    denovo_df['TranscriptClass'] = [is_novel(denovo_df.iloc[i]) for i in range(len(denovo_df))]
    
    start_time = time.time()
    novel_ends_results = [has_novel_ends(denovo_df.iloc[i]) for i in range(len(denovo_df))]
    denovo_df['Novel5pCap'] = [r[0] for r in novel_ends_results]
    denovo_df['NovelPolyA'] = [r[1] for r in novel_ends_results]
    
    # Update TranscriptClass for putative novel isoforms based on novel ends
    if denovo_allow_novel_ends:
        novel_isoform_mask = (
            (denovo_df['TranscriptClass'].isna()) & 
            (denovo_df['Novel5pCap'] | denovo_df['NovelPolyA'])
        )
        denovo_df.loc[novel_isoform_mask, 'TranscriptClass'] = 'putative_novel_isoform'
    # types of transcripts for later
    denovo_df['TranscriptMode'] = [tools.nameConversions.alignment_type(aln_id) for aln_id in denovo_df.AlignmentId]
    # filter out non-novel as well as fusions
    filtered_denovo_df = denovo_df[(~denovo_df.TranscriptClass.isnull())]
    filtered_denovo_df = filtered_denovo_df[filtered_denovo_df.TranscriptClass != 'possible_fusion']
    # fill in missing fields for novel loci
    filtered_denovo_df = filtered_denovo_df.copy()
    filtered_denovo_df['GeneBiotype'] = filtered_denovo_df['GeneBiotype'].fillna('unknown_likely_coding').astype(str)
    # filter out novel if requested by user
    if denovo_ignore_novel_genes == True:
        filtered_denovo_df = filtered_denovo_df[(filtered_denovo_df.TranscriptClass == 'possible_paralog') |
                                                (filtered_denovo_df.TranscriptClass == 'putative_novel_isoform')]
    elif denovo_only_novel_genes == True:
        filtered_denovo_df = filtered_denovo_df[~((filtered_denovo_df.TranscriptClass == 'possible_paralog') |
                                                (filtered_denovo_df.TranscriptClass == 'putative_novel_isoform'))]

    # construct aln_id -> features map to return
    denovo_tx_dict = {}
    for _, s in filtered_denovo_df.iterrows():
        aln_id = s.AlignmentId
        tx_mode = s.TranscriptMode
        denovo_tx_dict[aln_id] = {'source_gene': s.AssignedGeneId,
                                  'transcript_class': s.TranscriptClass,
                                  'novel_5p_cap': s.Novel5pCap,
                                  'novel_poly_a': s.NovelPolyA,
                                  'transcript_biotype': s.GeneBiotype,
                                  'gene_biotype': s.GeneBiotype,
                                  'intron_rna_support': ','.join(map(str, s.IntronRnaSupport)),
                                  'exon_rna_support': ','.join(map(str, s.ExonRnaSupport)),
                                  'exon_annotation_support': ','.join(map(str, s.ExonAnnotSupport)),
                                  'intron_annotation_support': ','.join(map(str, s.IntronAnnotSupport)),
                                  'alignment_id': aln_id,
                                  'source_gene_common_name': s.CommonName}

        # bring in extra tags for exRef
        if tools.nameConversions.aln_id_is_exref(aln_id):
            tags = exref_annot.loc[aln_id].ExtraTags
            if len(tags) > 0:
                for key, val in tools.misc.parse_gff_attr_line(tags).items():
                    # only add new tags if they don't already exist (and if they're informative)
                    if key not in denovo_tx_dict[aln_id] and val != "N/A":
                        denovo_tx_dict[aln_id][key] = val

        # record some metrics
        metrics['denovo'][tx_mode][s.TranscriptClass.replace('_', ' ').capitalize()] += 1
        metrics['Transcript Modes'][tx_mode] += 1
        metrics['Splice Support']['unknown_likely_coding'].append(s.IntronRnaSupportPercent)
        metrics['Exon Support']['unknown_likely_coding'].append(s.ExonRnaSupportPercent)

    # record how many of each type we threw out
    for tx_mode, df in denovo_df.groupby('TranscriptMode'):
        metrics['denovo'][tx_mode]['Discarded'] = len(df[(~df.TranscriptClass.isnull()) | (~df.Novel5pCap.isnull())
                                                         | (~df.NovelPolyA.isnull())])
    return denovo_tx_dict


def deduplicate_consensus(consensus_dict, tx_dict, metrics):
    """
    In the process of consensus building, we may find that we have ended up with more than one transcript for a gene
    that are actually identical. Remove these, picking the best based on their score, favoring the transcript
    whose biotype matches the parent.
    """
    def resolve_duplicate(tx_list, consensus_dict):
        biotype_txs = [tx for tx in tx_list if
                       consensus_dict[tx].get('gene_biotype', None) == consensus_dict[tx].get('transcript_biotype',
                                                                                              None)]
        if len(biotype_txs) > 0:
            tx_list = biotype_txs
        sorted_scores = sorted([[tx, consensus_dict[tx].get('score', 0)] for tx in tx_list],
                               key=lambda tx_s1: -tx_s1[1])
        return sorted_scores[0][0]

    def add_duplicate_field(best_tx, tx_list, consensus_dict, deduplicated_consensus):
        deduplicated_consensus[best_tx] = consensus_dict[best_tx]
        tx_list = [tools.nameConversions.strip_alignment_numbers(aln_id) for aln_id in tx_list]
        best_tx_base = tools.nameConversions.strip_alignment_numbers(best_tx)
        deduplicated_consensus[best_tx]['alternative_source_transcripts'] = ','.join(set(tx_list) - {best_tx_base})

    # build a dictionary mapping duplicates making use of hashing intervals
    # For gene families, also consider genomic location to avoid over-deduplication
    duplicates = collections.defaultdict(list)
    for aln_id in consensus_dict:
        tx = tx_dict[aln_id]
        
        is_txTM_cnv = False
        import re
        copy_match = re.search(r'_(\d+)$', aln_id)
        if copy_match:
            copy_num = int(copy_match.group(1))
            if copy_num > 0:
                is_txTM_cnv = True
        
        if is_txTM_cnv:
            location_key = (tx.chromosome, tx.start, tx.stop, aln_id)
            interval_key = (frozenset(tx.exon_intervals), location_key)
        else:
            location_key = (tx.chromosome, tx.start, tx.stop)
            interval_key = (frozenset(tx.exon_intervals), location_key)
        
        duplicates[interval_key].append(aln_id)

    # begin iterating
    deduplicated_consensus = {}
    for tx_list in duplicates.values():
        if len(tx_list) > 1:
            metrics['Duplicate transcripts'][len(tx_list)] += 1
            best_tx = resolve_duplicate(tx_list, consensus_dict)
            add_duplicate_field(best_tx, tx_list, consensus_dict, deduplicated_consensus)
        else:
            tx_id = tx_list[0]
            deduplicated_consensus[tx_id] = consensus_dict[tx_id]

    return deduplicated_consensus


def resolve_opposite_strand(deduplicated_consensus, tx_dict, metrics):
    """
    Resolves situations where multiple transcripts of the same gene are on opposite strands. Does so by looking for
    the largest sum of scores.
    """
    gene_dict = collections.defaultdict(list)
    deduplicated_strand_resolved_consensus = []
    for tx_id, attrs in deduplicated_consensus.items():
        tx_obj = tx_dict[tx_id]
        # don't try to resolve novel genes
        source_gene = attrs['source_gene']
        if source_gene is not None:
            gene_dict[source_gene].append([tx_obj, attrs])
        else:
            deduplicated_strand_resolved_consensus.append([tx_obj.name, attrs])

    for gene in gene_dict:
        tx_objs, attrs = list(zip(*gene_dict[gene]))
        if len(set(tx_obj.strand for tx_obj in tx_objs)) > 1:
            strand_scores = collections.Counter()
            for tx_obj, attrs in gene_dict[gene]:
                strand_scores[tx_obj.strand] += attrs.get('score', 0)
            best_strand = sorted(strand_scores.items())[1][0]
            for tx_obj, attrs in gene_dict[gene]:
                if tx_obj.strand == best_strand:
                    deduplicated_strand_resolved_consensus.append([tx_obj.name, attrs])
                else:
                    metrics['Discarded by strand resolution'] += 1
        else:
            deduplicated_strand_resolved_consensus.extend([[tx_obj.name, attrs] for tx_obj, attrs in gene_dict[gene]])
    return deduplicated_strand_resolved_consensus

def resolve_overlapping_cds_intervals(overlapping_ignore_bases, deduplicated_strand_resolved_consensus, tx_dict):
    attr_df = []
    with tools.fileOps.TemporaryFilePath() as tmp_gp, tools.fileOps.TemporaryFilePath() as tmp_clustered:
        with open(tmp_gp, 'w') as outf:
            for tx_id, attrs in deduplicated_strand_resolved_consensus:
                tx_obj = tx_dict[tx_id]
                tools.fileOps.print_row(outf, tx_obj.get_gene_pred())
                attr_df.append([tx_id, attrs['transcript_class'], attrs['gene_biotype'],
                                attrs.get('source_gene', tx_obj.name2), attrs.get('score', None)])
        # cluster
        cmd = ['clusterGenes', '-cds', f'-ignoreBases={overlapping_ignore_bases}',
               tmp_clustered, 'no', tmp_gp]
        tools.procOps.run_proc(cmd)
        cluster_df = pd.read_csv(tmp_clustered, sep='\t')
    attr_df = pd.DataFrame(attr_df, columns=['transcript_id', 'transcript_class', 'gene_biotype', 'gene_id', 'score'])
    m = attr_df.merge(cluster_df, left_on='transcript_id', right_on='gene')  # gene is transcript ID

    to_remove = set()  # list of transcript IDs to remove
    for cluster_id, group in m.groupby('#cluster'):
        if len(set(group['gene_id'])) > 1:
            if 'unknown_likely_coding' in set(group['gene_biotype']):  # pick longest ORF
                orfs = {tx_id: tx_dict[tx_id].cds_size for tx_id in group['transcript_id']}
                best_tx = sorted(iter(orfs.items()), key=lambda x: x[1])[-1][0]
                tx_df = group[group.transcript_id == best_tx].iloc[0]
                best_gene = tx_df.gene_id
            else:  # pick highest average score
                avg_scores = group[['gene_id', 'score']].groupby('gene_id', as_index=False).mean()
                best_gene = avg_scores.sort_values('score', ascending=False).iloc[0]['gene_id']
            to_remove.update(set(group[group.gene_id != best_gene].transcript_id))

    return [[tx_id, attrs] for tx_id, attrs in deduplicated_strand_resolved_consensus if tx_id not in to_remove]

def calculate_completeness(final_consensus, metrics):
    """calculates final completeness to make arithmetic easier"""
    genes = collections.defaultdict(set)
    txs = collections.Counter()
    for aln_id, c in final_consensus:
        # don't count novel transcripts towards completeness
        if tools.nameConversions.aln_id_is_pb(aln_id):
            continue
        genes[c['gene_biotype']].add(c['source_gene'])
        txs[c['transcript_biotype']] += 1
    genes = {biotype: len(gene_list) for biotype, gene_list in genes.items()}
    metrics['Completeness'] = {'Gene': genes, 'Transcript': txs}


def calculate_improvement_metrics(final_consensus, scored_df, tm_eval_df, support_df, metrics):
    """For coding transcripts, how much did we improve the metrics?"""
    tm_df = tm_eval_df.reset_index()[['TransMapOriginalIntronsPercent', 'TranscriptId']]
    support_df_subset = support_df[support_df['AlignmentId'].apply(tools.nameConversions.aln_id_is_transmap)]
    support_df_subset = support_df_subset[['TranscriptId', 'IntronAnnotSupportPercent', 'IntronRnaSupportPercent']]
    tm_df = pd.merge(tm_df, support_df_subset, on='TranscriptId')
    df = pd.merge(tm_df, scored_df.reset_index(), on='TranscriptId', suffixes=['TransMap', ''])
    df = df.drop_duplicates(subset='AlignmentId')  # why do I need to do this?
    df = df.set_index('AlignmentId')
    metrics['Evaluation Improvement'] = {'changes': [], 'unchanged': 0}
    for aln_id, c in final_consensus:
        if c['transcript_biotype'] != 'protein_coding':
            continue
        elif aln_id.startswith('exRef-') or _is_denovo(aln_id):
            continue
        elif aln_id.startswith('ENST'):  # transMap uses Ensembl IDs
            metrics['Evaluation Improvement']['unchanged'] += 1
            continue
        # Check if the alignment ID exists in the DataFrame before trying to access it
        if aln_id not in df.index:
            continue
        tx_s = df.loc[aln_id]
        metrics['Evaluation Improvement']['changes'].append([tx_s.TransMapOriginalIntronsPercent,
                                                             tx_s.IntronAnnotSupportPercentTransMap,
                                                             tx_s.IntronRnaSupportPercentTransMap,
                                                             tx_s.OriginalIntronsPercent_mRNA,
                                                             tx_s.IntronAnnotSupportPercent,
                                                             tx_s.IntronRnaSupportPercent,
                                                             tx_s.TransMapGoodness,
                                                             tx_s.AlnGoodness_mRNA])


def calculate_indel_metrics(final_consensus, eval_df, metrics):
    """How many transcripts in the final consensus have indels? How many did we have in transMap?"""
    if len(eval_df) == 0:  # edge case where no transcripts hit indel filters
        metrics['transMap Indels'] = {}
        metrics ['Consensus Indels'] = {}
        return
    eval_df_transmap = eval_df[eval_df['AlignmentId'].apply(tools.nameConversions.aln_id_is_transmap)]
    tm_vals = eval_df_transmap.set_index('AlignmentId').sum(axis=0)
    tm_vals = 100.0 * tm_vals / len(set(eval_df_transmap.index))
    metrics['transMap Indels'] = tm_vals.to_dict()
    consensus_ids = set(list(zip(*final_consensus))[0])
    consensus_vals = eval_df[eval_df['AlignmentId'].isin(consensus_ids)].set_index('AlignmentId').sum(axis=0)
    consensus_vals = 100.0 * consensus_vals / len(final_consensus)
    metrics['Consensus Indels'] = consensus_vals.to_dict()

def write_consensus_gps(consensus_gp, consensus_gp_info, final_consensus, tx_dict, genome, gene_offset=0):
    genes_seen = collections.defaultdict(dict)
    gene_count = gene_offset
    consensus_gene_dict = DefaultOrderedDict(lambda: DefaultOrderedDict(list))  # used to make gff3 next
    gp_infos = []

    # Write the GPD file directly
    with open(consensus_gp, 'w') as out_gp:
        for tx_count, (tx, attrs) in enumerate(final_consensus, 1):
            attrs = attrs.copy()
            tx_obj = tx_dict[tx]
            name = ID_TEMPLATE.format(genome=genome, tag_type='T', unique_id=tx_count)
            score = int(round(attrs.get('score', 0)))
            source_gene = attrs['source_gene']
            original_source_gene = source_gene  # Save for debugging
            alignment_id = attrs.get('alignment_id', '')
            transcript_class = attrs.get('transcript_class', 'ortholog')
            
            # Handle both None and string "None" from GFF3 conversion
            # Use alignment_id as fallback to ensure unique grouping
            if source_gene is None or source_gene == 'None' or source_gene == 'N/A':
                if alignment_id:
                    source_gene = f"unknown_{alignment_id}"
                else:
                    source_gene = f"unknown_{tx}"
            
            import re
            is_txTM_cnv = False
            if alignment_id:
                copy_match = re.search(r'_(\d+)$', alignment_id)
                if copy_match and int(copy_match.group(1)) > 0:
                    is_txTM_cnv = True
            
            if is_txTM_cnv:
                unique_gene_key = (source_gene, alignment_id)
            else:
                # Regular gene, use source_gene as key
                unique_gene_key = source_gene
            
            if unique_gene_key not in genes_seen[tx_obj.chromosome]:
                gene_count += 1
                genes_seen[tx_obj.chromosome][unique_gene_key] = gene_count
            gene_id = genes_seen[tx_obj.chromosome][unique_gene_key]
            name2 = ID_TEMPLATE.format(genome=genome, tag_type='G', unique_id=gene_id)
            out_gp.write('\t'.join(tx_obj.get_gene_pred(name=name, name2=name2, score=score)) + '\n')
            attrs['transcript_id'] = name
            attrs['gene_id'] = name2
            gp_infos.append(attrs)
            consensus_gene_dict[tx_obj.chromosome][name2].append([tx_obj, attrs])

    # Create DataFrame of GP infos
    gp_info_df = pd.DataFrame(gp_infos).set_index(['gene_id', 'transcript_id'])
    if 'alternative_source_transcripts' not in gp_info_df.columns:
        gp_info_df['alternative_source_transcripts'] = ['N/A'] * len(gp_info_df)

    # Write consensus GP info
    with open(consensus_gp_info, 'w') as outf:
        gp_info_df.to_csv(outf, sep='\t', na_rep='N/A')

    return consensus_gene_dict

def write_consensus_gff3(consensus_gene_dict, consensus_gff3):
    """
    Write the consensus set in gff3 format
    """
    def convert_attrs(attrs, id_field):
        import re
        attrs['ID'] = id_field
        if 'score' in attrs:
            score = 10 * attrs['score']
            del attrs['score']
        else:
            score = '.'
        if 'source_gene_common_name' in attrs and isinstance(attrs['source_gene_common_name'], str):
            attrs['Name'] = attrs['source_gene_common_name']
        else:
            attrs['Name'] = attrs['gene_id']
        attrs_str = []
        for key, val in attrs.items():
            # Handle NaN values properly
            if pd.isna(val) or val is None:
                val = 'None'
            else:
                val = str(val)
                # Replace 'nan' string with 'None' for better GFF3 compatibility
                if val.lower() == 'nan':
                    val = 'None'
                # Strip _cp suffixes from alignment_id and normalized_transcript_id in final output
                # (these are internal tracking IDs that shouldn't be exposed to users)
                if key in ['alignment_id', 'normalized_transcript_id']:
                    val = re.sub(r'_cp\d+$', '', val)
            val = val.replace('=', '%3D').replace(';', '%3B')
            key = key.replace('=', '%3D').replace(';', '%3B')
            attrs_str.append(f"{key}={val}")
        return score, ";".join(attrs_str)

    def find_feature_support(attrs, feature, i):
        try:
            vals = list(map(bool, attrs[feature].split(',')))
        except KeyError:
            return 'N/A'
        return vals[i] if i < len(vals) else False

    def get_gene_name(attrs):
        return attrs.get("source_gene_common_name") or attrs.get("gene_id")

    def generate_gene_record(chrom, tx_objs, gene_id, attrs_list):
        def find_all_sources(attrs_list):
            """Collect all sources from all transcripts under this gene by looking at alignment_id"""
            sources = set()
            for attrs in attrs_list:
                # Get source from alignment_id
                if 'alignment_id' in attrs and pd.notna(attrs['alignment_id']):
                    aln_id = attrs['alignment_id']
                    source = tools.nameConversions.alignment_type(aln_id)
                    sources.add(source)
                
                # Also add source_type if available
                if 'source_type' in attrs and pd.notna(attrs['source_type']):
                    sources.add(attrs['source_type'])
                
                # Add additional sources (from overlap resolution)
                if 'AdditionalSources' in attrs and pd.notna(attrs['AdditionalSources']) and attrs['AdditionalSources'] != '':
                    sources.update(attrs['AdditionalSources'].split(','))
            
            # Sort to ensure consistent ordering: txTM first, then others alphabetically
            sorted_sources = sorted(sources, key=lambda x: (x != 'txTM', x))
            return ','.join(sorted_sources)

        intervals = {iv for tx in tx_objs for iv in tx.exon_intervals}
        intervals = sorted(intervals)
        strand = tx_objs[0].strand
        attrs = {k: attrs_list[0][k] for k in [
            'source_gene_common_name', 'source_gene',
            'alternative_source_transcripts', 'gene_alternate_contigs',
            'gene_id', 'collapsed_gene_ids', 'collapsed_gene_names'
        ] if k in attrs_list[0]}
        
        # Select gene_biotype from highest-priority source
        # CRITICAL: Exclude augPB transcripts from biotype determination - they inherit biotype, not set it
        # TxTM/transMap/augTM/augMP transcripts contribute; augMP orthologs with CDS can
        # still force protein_coding below if no other transcript set it.
        # Priority: protein_coding > lncRNA/other known biotypes > unknown_likely_coding
        gene_biotype = attrs_list[0].get('gene_biotype', 'unknown')
        # Fix: if gene_biotype is None or N/A, set to 'unknown'
        if gene_biotype is None or (isinstance(gene_biotype, float) and pd.isna(gene_biotype)) or gene_biotype == 'N/A':
            gene_biotype = 'unknown'
        for attr in attrs_list:
            # Skip denovo transcripts - they should not influence gene biotype
            aln_id = attr.get('alignment_id', '')
            if _is_denovo(aln_id):
                continue
            
            biotype = attr.get('gene_biotype', 'unknown')
            # Fix: if biotype is None or N/A, skip it
            if biotype is None or (isinstance(biotype, float) and pd.isna(biotype)) or biotype == 'N/A':
                continue
            if biotype == 'protein_coding':
                gene_biotype = biotype
                break
            elif biotype not in ['unknown_likely_coding', 'unknown'] and gene_biotype in ['unknown_likely_coding', 'unknown']:
                gene_biotype = biotype
        # If augMP still has a full CDS as an ortholog model, call the gene protein_coding so
        # miniprot-rescued loci count like other projection tracks in summary stats.
        if gene_biotype != 'protein_coding':
            for tx_obj, attr in zip(tx_objs, attrs_list):
                if attr.get('alignment_mode') != 'augMP':
                    continue
                if attr.get('transcript_class') != 'ortholog':
                    continue
                if getattr(tx_obj, 'cds_size', 0) > 3:
                    gene_biotype = 'protein_coding'
                    break
        attrs['gene_biotype'] = gene_biotype
        
        attrs['gene_name'] = get_gene_name(attrs)
        attrs['gene_sources'] = find_all_sources(attrs_list)
        score, attrs_field = convert_attrs(attrs, gene_id)
        return [chrom, 'CAT', 'gene', intervals[0].start + 1, intervals[-1].stop, score, strand, '.', attrs_field]

    def generate_transcript_record(chrom, tx_obj, attrs):
        tx_id = attrs['transcript_id']
        attrs['Parent'] = attrs['gene_id']
        attrs['transcript_name'] = tx_id
        if 'gene_name' not in attrs:
            attrs['gene_name'] = get_gene_name(attrs)
        score, attrs_field = convert_attrs(attrs, tx_id)
        yield [chrom, 'CAT', 'transcript', tx_obj.start + 1, tx_obj.stop, score, tx_obj.strand, '.', attrs_field]
        attrs.pop('frameshift', None)
        for line in generate_intron_exon_records(chrom, tx_obj, tx_id, attrs):
            yield line
        if tx_obj.cds_size > 3:
            for line in generate_start_stop_codon_records(chrom, tx_obj, tx_id, attrs):
                yield line

    def generate_intron_exon_records(chrom, tx_obj, tx_id, attrs):
        attrs['Parent'] = tx_id
        cds_i = 0
        for i, (exon, exon_frame) in enumerate(zip(tx_obj.exon_intervals, tx_obj.exon_frames)):
            attrs['rna_support'] = find_feature_support(attrs, 'exon_rna_support', i)
            attrs['reference_support'] = find_feature_support(attrs, 'exon_annotation_support', i)
            score, attrs_field = convert_attrs(attrs, f'exon:{tx_id}:{i}')
            yield [chrom, 'CAT', 'exon', exon.start + 1, exon.stop, score, exon.strand, '.', attrs_field]
            cds_interval = exon.intersection(tx_obj.coding_interval)
            if cds_interval:
                score, attrs_field = convert_attrs(attrs, f'CDS:{tx_id}:{cds_i}')
                cds_i += 1
                yield [chrom, 'CAT', 'CDS', cds_interval.start + 1, cds_interval.stop, score, exon.strand, tools.transcripts.convert_frame(exon_frame), attrs_field]
        for i, intron in enumerate(tx_obj.intron_intervals):
            if len(intron) == 0:
                continue
            attrs['rna_support'] = find_feature_support(attrs, 'intron_rna_support', i)
            attrs['reference_support'] = find_feature_support(attrs, 'intron_annotation_support', i)
            score, attrs_field = convert_attrs(attrs, f'intron:{tx_id}:{i}')
            yield [chrom, 'CAT', 'intron', intron.start + 1, intron.stop, score, intron.strand, '.', attrs_field]

    def generate_start_stop_codon_records(chrom, tx_obj, tx_id, attrs):
        if attrs.get('valid_start'):
            score, attrs_field = convert_attrs(attrs, f'start_codon:{tx_id}')
            for iv in tx_obj.get_start_intervals():
                yield [chrom, 'CAT', 'start_codon', iv.start + 1, iv.stop, score, tx_obj.strand, iv.data, attrs_field]
        if attrs.get('valid_stop'):
            score, attrs_field = convert_attrs(attrs, f'stop_codon:{tx_id}')
            for iv in tx_obj.get_stop_intervals():
                yield [chrom, 'CAT', 'stop_codon', iv.start + 1, iv.stop, score, tx_obj.strand, iv.data, attrs_field]

    # Main writing logic: open file directly
    with open(consensus_gff3, 'w') as out_gff3:
        out_gff3.write('##gff-version 3\n')
        for chrom in sorted(consensus_gene_dict):
            for gene_id, tx_list in consensus_gene_dict[chrom].items():
                tx_objs, attrs_list = zip(*tx_list)
                # Write gene record first
                gene_row = generate_gene_record(chrom, tx_objs, gene_id, attrs_list)
                tools.fileOps.print_rows(out_gff3, [gene_row])
                
                # Write each transcript and all its child features together
                for tx_obj, attrs in tx_list:
                    tx_rows = list(generate_transcript_record(chrom, tx_obj, attrs))
                    # Sort child features by start position within each transcript
                    for row in sorted(tx_rows, key=lambda r: r[3]):
                        tools.fileOps.print_rows(out_gff3, [row])


def write_consensus_fastas(consensus_gene_dict, consensus_fasta, consensus_protein_fasta, fasta):
    """Write FASTA records for both transcripts and proteins"""
    seq_dict = tools.bio.get_sequence_dict(fasta)
    with open(consensus_fasta, 'w') as cfa, open(consensus_protein_fasta, 'w') as cpfa:
        for chrom in sorted(consensus_gene_dict):
            for gene_id, tx_list in consensus_gene_dict[chrom].items():
                for tx_obj, attrs in tx_list:
                    mrna_seq = tx_obj.get_mrna(seq_dict)
                    tools.bio.write_fasta(cfa, attrs['transcript_id'], mrna_seq)
                    if tx_obj.cds_size > 0:
                        protein_seq = tx_obj.get_protein_sequence(seq_dict)
                    elif str(attrs.get('transcript_class', '')).startswith('putative_novel'):
                        # For novel predictions without CDS annotations (e.g., StringTie),
                        # treat the full spliced transcript as coding sequence.
                        protein_seq = tools.bio.translate_sequence(mrna_seq)
                    else:
                        protein_seq = ''
                    if protein_seq:
                        tools.bio.write_fasta(cpfa, attrs['transcript_id'], protein_seq)

if __name__ == "__main__":
    main()
