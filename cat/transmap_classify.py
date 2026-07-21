import argparse
import bisect
import collections
import sqlite3

import pandas as pd

import tools.bio
import tools.mathOps
import tools.nameConversions
import tools.psl
import tools.tm2hints
import tools.transcripts


def generate_evaluation_df(args):
    """
    Generates a DataFrame containing classification metrics for transMap alignments.

    :param args: An argparse.Namespace object containing all file paths.
    :return: A pandas.DataFrame with the evaluation results.
    """
    # Load all necessary files
    psl_dict = tools.psl.get_alignment_dict(args.filtered_tm_psl)
    ref_psl_dict = tools.psl.get_alignment_dict(args.ref_psl)
    gp_dict = tools.transcripts.get_gene_pred_dict(args.filtered_tm_gp)
    ref_gp_dict = tools.transcripts.get_gene_pred_dict(args.annotation_gp)
    fasta = tools.bio.get_sequence_dict(args.fasta)

    # Pre-calculate synteny scores
    synteny_scores = synteny(ref_gp_dict, gp_dict)

    # Generate a record for each classifier for each transcript
    records = []
    for aln_id, tx in gp_dict.items():
        aln = psl_dict[aln_id]
        tx_id = tools.nameConversions.strip_alignment_numbers(aln_id)
        ref_aln = ref_psl_dict[tx_id]
        gene_id = ref_gp_dict[tx_id].name2
        
        # Add a dictionary of results for this transcript to the list of records
        records.extend([
            {'AlignmentId': aln_id, 'TranscriptId': tx_id, 'GeneId': gene_id, 'classifier': 'AlnExtendsOffContig', 'value': aln_extends_off_contig(aln)},
            {'AlignmentId': aln_id, 'TranscriptId': tx_id, 'GeneId': gene_id, 'classifier': 'AlnPartialMap', 'value': alignment_partial_map(aln)},
            {'AlignmentId': aln_id, 'TranscriptId': tx_id, 'GeneId': gene_id, 'classifier': 'AlnAbutsUnknownBases', 'value': aln_abuts_unknown_bases(tx, fasta)},
            {'AlignmentId': aln_id, 'TranscriptId': tx_id, 'GeneId': gene_id, 'classifier': 'PercentN', 'value': aln.percent_n},
            {'AlignmentId': aln_id, 'TranscriptId': tx_id, 'GeneId': gene_id, 'classifier': 'TransMapCoverage', 'value': 100 * aln.coverage},
            {'AlignmentId': aln_id, 'TranscriptId': tx_id, 'GeneId': gene_id, 'classifier': 'TransMapIdentity', 'value': 100 * aln.identity},
            {'AlignmentId': aln_id, 'TranscriptId': tx_id, 'GeneId': gene_id, 'classifier': 'TransMapGoodness', 'value': 100 * (1 - aln.badness)},
            {'AlignmentId': aln_id, 'TranscriptId': tx_id, 'GeneId': gene_id, 'classifier': 'TransMapOriginalIntronsPercent', 'value': percent_original_introns(aln, tx, ref_aln)},
            {'AlignmentId': aln_id, 'TranscriptId': tx_id, 'GeneId': gene_id, 'classifier': 'Synteny', 'value': synteny_scores[aln_id]},
            {'AlignmentId': aln_id, 'TranscriptId': tx_id, 'GeneId': gene_id, 'classifier': 'ValidStart', 'value': tools.transcripts.has_start_codon(fasta, tx)},
            {'AlignmentId': aln_id, 'TranscriptId': tx_id, 'GeneId': gene_id, 'classifier': 'ValidStop', 'value': tools.transcripts.has_stop_codon(fasta, tx)},
            {'AlignmentId': aln_id, 'TranscriptId': tx_id, 'GeneId': gene_id, 'classifier': 'ProperOrf', 'value': tx.cds_size % 3 == 0}
        ])
    
    expected_columns = ['AlignmentId', 'TranscriptId', 'GeneId', 'classifier', 'value']
    df = pd.DataFrame(columns=expected_columns) if not records else pd.DataFrame(records)
    if 'value' in df.columns:
        df['value'] = pd.to_numeric(df['value'], errors='coerce')
    return df


def main():
    """
    Main entry point for the transMap classification script. Parses arguments,
    runs the classification, and writes the results to a SQLite database.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    # Input files
    parser.add_argument("--filtered-tm-psl", required=True, help="Filtered transMap PSL file.")
    parser.add_argument("--filtered-tm-gp", required=True, help="Filtered transMap genePred file.")
    parser.add_argument("--ref-psl", required=True, help="Reference 'fake' PSL file.")
    parser.add_argument("--annotation-gp", required=True, help="Reference annotation genePred file.")
    parser.add_argument("--fasta", required=True, help="Target genome FASTA file.")
    # Output database
    parser.add_argument("--db-path", required=True, help="Path to the SQLite database to write results to.")
    parser.add_argument("--resolved-df", help="Optional path to save evaluation DataFrame as a Pickle file.")
    
    args = parser.parse_args()
    
    # Generate the evaluation dataframe
    df = generate_evaluation_df(args)
    
    # Write the dataframe to the specified database table
    con = sqlite3.connect(args.db_path)
    df.to_sql('TransMapEvaluation', con, if_exists='replace', index=False)
    con.close()
    df.to_pickle(args.resolved_df)


# Classifier functions remain unchanged
def aln_extends_off_contig(aln):
    return aln.t_start == 0 and aln.q_start != 0 or aln.t_end == aln.t_size and aln.q_end != aln.q_size

def alignment_partial_map(aln):
    return aln.q_size != aln.q_end - aln.q_start

def aln_abuts_unknown_bases(tx, fasta):
    chrom = tx.chromosome
    for exon in tx.exon_intervals:
        if exon.start > 0 and fasta[chrom][exon.start - 1] == 'N':
            return True
        if exon.stop < len(fasta[chrom]) and fasta[chrom][exon.stop] == 'N':
            return True
    return False

def synteny(ref_gp_dict, gp_dict):
    def get_sorted_gene_intervals(tx_dict):
        gene_intervals = collections.defaultdict(list)
        for tx in tx_dict.values():
            gene_intervals[tx.chromosome].append(tx.interval)
        
        sorted_intervals = {}
        for chrom, intervals in gene_intervals.items():
            merged = tools.intervals.gap_merge_intervals(intervals, float('inf'))
            for i in merged:
                i.data = {tx.name2 for tx in tx_dict.values() if tx.interval.overlap(i)}
            sorted_intervals[chrom] = sorted(merged)
        return sorted_intervals

    tm_chrom_intervals = get_sorted_gene_intervals(gp_dict)
    ref_chrom_intervals = get_sorted_gene_intervals(ref_gp_dict)
    ref_gene_map = {gene_id: interval for chrom in ref_chrom_intervals for interval in ref_chrom_intervals[chrom] for gene_id in interval.data}

    scores = {}
    for tx in gp_dict.values():
        target_intervals = tm_chrom_intervals.get(tx.chromosome, [])
        target_pos = bisect.bisect_left(target_intervals, tx.interval)
        target_genes = {g for i in target_intervals[max(0, target_pos - 5):target_pos + 6] for g in i.data}

        ref_interval = ref_gene_map.get(tx.name2)
        if not ref_interval:
            scores[tx.name] = 0
            continue
        
        ref_intervals = ref_chrom_intervals.get(ref_interval.chromosome, [])
        ref_pos = bisect.bisect_left(ref_intervals, ref_interval)
        reference_genes = {g for i in ref_intervals[max(0, ref_pos - 5):ref_pos + 6] for g in i.data}
        
        scores[tx.name] = len(reference_genes & target_genes)
    return scores

def percent_original_introns(aln, tx, ref_aln, fuzz_distance=7):
    ref_starts = tools.tm2hints.fix_ref_q_starts(ref_aln)
    supported_introns = sum(1 for intron in tx.intron_intervals if tools.tm2hints.is_fuzzy_intron(intron, aln, ref_starts, fuzz_distance))
    return 100 * tools.mathOps.format_ratio(supported_introns, len(tx.intron_intervals))

if __name__ == "__main__":
    main()
