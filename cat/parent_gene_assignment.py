import argparse
import collections
import itertools
import sqlite3

import pandas as pd

import tools.fileOps
import tools.intervals
import tools.mathOps
import tools.transcripts
from tools.defaultOrderedDict import DefaultOrderedDict


def generate_parent_df(args):
    """
    Main logic for assigning parental genes.

    :param args: An argparse.Namespace object containing all file paths and parameters.
    :return: A pandas.DataFrame with the parent gene assignments.
    """
    # Load transcript data from genePred files
    filtered_transmap_dict = tools.transcripts.get_gene_pred_dict(args.filtered_tm_gp)
    unfiltered_transmap_dict = tools.transcripts.get_gene_pred_dict(args.unfiltered_tm_gp)
    denovo_dict = tools.transcripts.get_gene_pred_dict(args.denovo_gp)

    # Check if de novo file is empty
    if not denovo_dict:
        print(f"Warning: No de novo gene predictions found in {args.denovo_gp}")
        # Return empty DataFrame with proper schema
        df = pd.DataFrame(columns=['TranscriptId', 'AssignedGeneId', 'AlternativeGeneIds', 'ResolutionMethod'])
        df = df.astype(str)
        return df

    # Dictionaries mapping chromosomes to transcripts for faster lookups
    tm_chrom_dict = create_chrom_dict(unfiltered_transmap_dict, args.chrom_sizes)
    denovo_chrom_dict = create_chrom_dict(denovo_dict)

    # Determine which transMap transcripts were filtered out
    filtered_ids = set(unfiltered_transmap_dict.keys()) - set(filtered_transmap_dict.keys())

    # Begin parent gene assignment process
    records = []
    for chrom, denovo_txs_on_chrom in denovo_chrom_dict.items():
        tm_txs_on_chrom = tm_chrom_dict.get(chrom, {})
        for denovo_tx in denovo_txs_on_chrom.values():
            # Find all unfiltered transMap transcripts that overlap this de novo transcript
            unfiltered_overlapping_txs = find_tm_overlaps(denovo_tx, tm_txs_on_chrom)
            
            # Find which of those are in the filtered set
            filtered_overlapping_txs = {tx for tx in unfiltered_overlapping_txs if tx.name not in filtered_ids}
            filtered_gene_ids = {tx.name2 for tx in filtered_overlapping_txs}

            resolved_name, resolution_method = None, None
            if len(filtered_gene_ids) > 1:
                # Attempt to resolve if multiple parent genes are found
                resolved_name, resolution_method = resolve_multiple_genes(denovo_tx, filtered_overlapping_txs,
                                                                          args.min_distance)
            elif len(filtered_gene_ids) == 1:
                # Exactly one parent gene found
                resolved_name = list(filtered_gene_ids)[0]

            # Identify alternative parent genes from the unfiltered set
            alternative_gene_ids = {tx.name2 for tx in unfiltered_overlapping_txs} - {resolved_name}
            
            # Convert alternative gene IDs to string, handling empty case
            alternative_genes_str = ','.join(sorted(alternative_gene_ids)) if alternative_gene_ids else ''
            
            records.append({'TranscriptId': denovo_tx.name,
                            'AssignedGeneId': resolved_name if resolved_name else '',
                            'AlternativeGeneIds': alternative_genes_str,
                            'ResolutionMethod': resolution_method if resolution_method else ''})

    # Create DataFrame with proper handling of empty case
    if records:
        df = pd.DataFrame(records)
        # Ensure all columns are string type to avoid SQL issues
        for col in df.columns:
            df[col] = df[col].astype(str)
    else:
        # Create empty DataFrame with proper column structure
        df = pd.DataFrame(columns=['TranscriptId', 'AssignedGeneId', 'AlternativeGeneIds', 'ResolutionMethod'])
        df = df.astype(str)
    
    return df


def main():
    """
    Main entry point for the parent assignment script. Parses arguments,
    runs the assignment logic, and writes results to a database.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    # Input files
    parser.add_argument("--filtered-tm-gp", required=True, help="Filtered transMap genePred file.")
    parser.add_argument("--unfiltered-tm-gp", required=True, help="Unfiltered transMap genePred file.")
    parser.add_argument("--chrom-sizes", required=True, help="Chromosome sizes file.")
    parser.add_argument("--denovo-gp", required=True, help="De novo genePred file (e.g., from AugustusPB).")
    # Output database info
    parser.add_argument("--db-path", required=True, help="Path to the SQLite database to write to.")
    parser.add_argument("--table-name", required=True, help="Name of the table to store results in.")
    # Parameters
    parser.add_argument("--min-distance", type=float, default=0.9, help="Minimum asymmetric distance difference to resolve multi-gene assignments.")
    
    args = parser.parse_args()
    
    # Generate the assignment dataframe
    df = generate_parent_df(args)
    
    # Debug: Print DataFrame info
    print(f"DataFrame shape: {df.shape}")
    print(f"DataFrame columns: {df.columns.tolist()}")
    print(f"DataFrame dtypes:\n{df.dtypes}")
    if not df.empty:
        print(f"First few rows:\n{df.head()}")
    else:
        print("DataFrame is empty!")
    
    # Write the results to the specified database table
    con = sqlite3.connect(args.db_path)
    try:
        # Clean the DataFrame before writing to SQL
        # Replace None values with empty strings for string columns
        df_clean = df.copy()
        for col in ['AssignedGeneId', 'AlternativeGeneIds', 'ResolutionMethod']:
            if col in df_clean.columns:
                df_clean[col] = df_clean[col].fillna('')
        
        # Ensure column names are valid SQL identifiers
        df_clean.columns = [col.replace('-', '_').replace(' ', '_') for col in df_clean.columns]
        
        df_clean.to_sql(args.table_name, con, if_exists='replace', index=False)
        print(f"Successfully wrote {len(df_clean)} records to database table '{args.table_name}'")
    except Exception as e:
        print(f"Error writing to database: {e}")
        print(f"DataFrame info: {df.info()}")
        raise
    finally:
        con.close()

def create_chrom_dict(tx_dict, chrom_sizes=None):
    """
    Split up a dictionary of Transcript objects by chromosome. Add in extra chromosomes based on a sizes file
    """
    chrom_dict = collections.defaultdict(dict)
    for tx_id, tx in tx_dict.items():
        chrom_dict[tx.chromosome][tx_id] = tx
    if chrom_sizes is not None:
        for chrom, size in tools.fileOps.iter_lines(chrom_sizes):
            if chrom not in chrom_dict:
                chrom_dict[chrom] = {}
    return chrom_dict


def find_tm_overlaps(denovo_tx, tm_tx_dict, cutoff=100):
    """Find overlap with transMap transcripts first on a genomic scale then an exonic scale"""
    r = DefaultOrderedDict(int)
    for tx in tm_tx_dict.values():
        for tx_exon in tx.exon_intervals:
            for denovo_exon in denovo_tx.exon_intervals:
                i = tx_exon.intersection(denovo_exon)
                if i is not None:
                    r[tx] += len(i)
    return [tx_id for tx_id, num_bases in r.items() if num_bases >= cutoff]


def resolve_multiple_genes(denovo_tx, overlapping_tm_txs, min_distance):
    """
    Resolve multiple assignments based on the following rules:
    """
    # use Jaccard metric to determine if the problem lies with transMap or annotation
    tm_txs_by_gene = tools.transcripts.group_transcripts_by_name2(overlapping_tm_txs)
    tm_jaccards = [find_highest_gene_jaccard(x, y) for x, y in itertools.combinations(list(tm_txs_by_gene.values()), 2)]
    if any(x > 0.001 for x in tm_jaccards):
        return None, 'badAnnotOrTm'
    # calculate asymmetric difference for this prediction
    scores = collections.defaultdict(list)
    for tx in overlapping_tm_txs:
        scores[tx.name2].append(calculate_asymmetric_closeness(denovo_tx, tx))
    best_scores = {gene_id: max(scores[gene_id]) for gene_id in scores}
    high_score = max(best_scores.values())
    if all(high_score - x >= min_distance for x in best_scores.values() if x != high_score):
        best = sorted(iter(best_scores.items()), key=lambda gene_id_score: gene_id_score[1])[-1][0]
        return best, 'rescued'
    else:
        return None, 'ambiguousOrFusion'


def find_highest_gene_jaccard(gene_list_a, gene_list_b):
    """
    Calculates the overall distance between two sets of transcripts by finding their distinct exonic intervals and then
    measuring the Jaccard distance.
    """
    def find_interval(gene_list):
        gene_intervals = set()
        for tx in gene_list:
            gene_intervals.update(tx.exon_intervals)
        gene_intervals = tools.intervals.gap_merge_intervals(gene_intervals, 0)
        return gene_intervals

    a_interval = find_interval(gene_list_a)
    b_interval = find_interval(gene_list_b)
    return tools.intervals.calculate_bed12_jaccard(a_interval, b_interval)


def calculate_asymmetric_closeness(denovo_tx, tm_tx):
    """
    Calculates the asymmetric closeness between two transcripts. This allows for denovo predictions that are subsets
    of existing annotations to get more weight.

    closeness = length(intersection) / length(denovo) in chromosome space
    """
    intersection = denovo_tx.interval.intersection(tm_tx.interval)
    if intersection is None:
        return 0
    return tools.mathOps.format_ratio(len(intersection), len(denovo_tx.interval))

if __name__ == "__main__":
    main()
