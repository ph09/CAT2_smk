import argparse
import collections
import hashlib
import json
import os
import tempfile
import shutil
import logging
from copy import deepcopy
from pathlib import Path
import pandas as pd

import tools.fileOps
import tools.intervals
import tools.mathOps
import tools.nameConversions
import tools.procOps
import tools.psl
import tools.sqlInterface
import tools.transcripts

pd.options.mode.chained_assignment = None


def run_filtering(tm_psl, ref_psl, tm_gp, db_path, psl_tgt, global_near_best, filter_overlapping_genes,
                    overlapping_ignore_bases, json_tgt, annotate_extra_paralogs, temp_dir=None,
                    min_cover=0.1, min_span=0.2, max_ref_span=5, paralog_rescue_min_coverage=0.5):
    """
    Main filtering logic. Wrapped in a function to be called by main().

    The recall-affecting cutoffs (``min_cover``, ``min_span``, ``max_ref_span``,
    ``paralog_rescue_min_coverage``) are parameterised so the Snakemake pipeline
    can loosen them (e.g. via its high_recall preset) without editing this file.
    """
    # Use provided temp directory or create one
    if temp_dir is None:
        temp_dir = tempfile.mkdtemp()
        cleanup_temp = True
    else:
        cleanup_temp = False
    
    unfiltered = tools.psl.get_alignment_dict(tm_psl)
    unfiltered_tx_dict = tools.transcripts.get_gene_pred_dict(tm_gp)
    ref_psl_dict = tools.psl.get_alignment_dict(ref_psl)
    size_filtered, num_too_long = ref_span(unfiltered, ref_psl_dict, max_span=max_ref_span)
    tmp_size_filtered = os.path.join(temp_dir, 'size_filtered.psl')
    with open(tmp_size_filtered, 'w') as outf:
        for aln in size_filtered.values():
            tools.fileOps.print_row(outf, aln.psl_string())
    transcript_gene_map = tools.sqlInterface.get_transcript_gene_map(db_path)
    transcript_biotype_map = tools.sqlInterface.get_transcript_biotype_map(db_path)
    gene_biotype_map = tools.sqlInterface.get_gene_biotype_map(db_path)
    annotation_df = tools.sqlInterface.load_annotation(db_path)
    gene_name_map = dict(list(zip(annotation_df.GeneId, annotation_df.GeneName)))
    
    def hash_aln(aln, aln_id):
        m = hashlib.sha256()
        for l in [aln.t_name, aln.t_start, aln.t_end, aln.matches, aln.mismatches, aln.block_count,
                  tuple(aln.t_starts), tuple(aln.q_starts), tuple(aln.block_sizes), aln_id]:
            m.update(str(l).encode('utf-8'))
        return m.hexdigest()

    unfiltered_hash_table = {}
    q_size_map = {}
    mismatched_q_sizes = collections.defaultdict(set)
    for aln_id, aln in size_filtered.items():
        stripped_id = tools.nameConversions.strip_alignment_numbers(aln_id)
        unfiltered_hash_table[hash_aln(aln, stripped_id)] = aln_id
        existing_q_size = q_size_map.get(stripped_id)
        if existing_q_size is None:
            q_size_map[stripped_id] = aln.q_size
        elif existing_q_size != aln.q_size:
            mismatched_q_sizes[stripped_id].update([existing_q_size, aln.q_size])
            if aln.q_size > existing_q_size:
                q_size_map[stripped_id] = aln.q_size
    if mismatched_q_sizes:
        logging.warning(
            "Normalizing inconsistent qSize values for %d transcript(s) (examples: %s)",
            len(mismatched_q_sizes),
            ', '.join(sorted(mismatched_q_sizes.keys())[:5])
        )
    
    # Use local temp files instead of TemporaryFilePath
    local_tmp = os.path.join(temp_dir, 'local_tmp.psl')
    strip_tmp = os.path.join(temp_dir, 'strip_tmp.psl')
    try:
        with open(strip_tmp, 'w') as outf:
            for rec in size_filtered.values():
                rec = deepcopy(rec)
                stripped_name = tools.nameConversions.strip_alignment_numbers(rec.q_name)
                rec.q_name = stripped_name
                if stripped_name in q_size_map:
                    rec.q_size = q_size_map[stripped_name]
                tools.fileOps.print_row(outf, rec.psl_string())
        cmd = ['pslCDnaFilter', '-globalNearBest={}'.format(global_near_best),
               '-minCover={}'.format(min_cover), '-verbose=0',
               '-minSpan={}'.format(min_span), strip_tmp, '/dev/stdout']
        tools.procOps.run_proc(cmd, stdout=local_tmp)
        filtered_alns = list(tools.psl.psl_iterator(local_tmp))
    finally:
        # Clean up these temp files
        for f in [local_tmp, strip_tmp]:
            if os.path.exists(f):
                os.remove(f)
    
    global_best = {unfiltered[unfiltered_hash_table[hash_aln(aln, aln.q_name)]] for aln in filtered_alns}

    def lookup_tx(aln):
        """Look up transcript in genePred, trying q_name and fallbacks (stripped, version variations)."""
        q_name = aln.q_name
        if q_name in unfiltered_tx_dict:
            return unfiltered_tx_dict[q_name]
        stripped = tools.nameConversions.strip_alignment_numbers(q_name)
        if stripped in unfiltered_tx_dict:
            return unfiltered_tx_dict[stripped]
        # Try without version suffix (NM_001146706.2 -> NM_001146706)
        import re
        base = re.sub(r'\.[0-9]+$', '', q_name)
        if base in unfiltered_tx_dict:
            return unfiltered_tx_dict[base]
        return None

    global_best_txs = []
    missing_tx_ids = []
    for aln in global_best:
        tx = lookup_tx(aln)
        if tx is not None:
            global_best_txs.append(tx)
        else:
            missing_tx_ids.append(aln.q_name)
    if missing_tx_ids:
        logging.warning(
            "Skipping %d global-best alignments not found in genePred (examples: %s)",
            len(missing_tx_ids),
            ', '.join(sorted(missing_tx_ids)[:5])
        )
        # Exclude these from global_best so downstream uses consistent set
        global_best = {aln for aln in global_best if lookup_tx(aln) is not None}

    global_best_ids = {x.name for x in global_best_txs}
    grouped = tools.psl.group_alignments_by_qname(global_best)
    unfiltered_grouped = tools.psl.group_alignments_by_qname(iter(unfiltered.values()))
    metrics = {'Paralogy': collections.defaultdict(lambda: collections.Counter()),
               'UnfilteredParalogy': collections.defaultdict(lambda: collections.Counter())}
    paralogy_df = []
    for tx_id, alns in grouped.items():
        biotype = transcript_biotype_map[tx_id]
        putative_paralogs = ','.join(sorted([x.q_name for x in alns if x.q_name not in global_best_ids]))
        all_alns = ','.join(sorted([x.q_name for x in unfiltered_grouped[tx_id] if x.q_name not in global_best_ids]))
        paralogy_df.append([tx_id, putative_paralogs, all_alns])
        metrics['Paralogy'][biotype][len(alns)] += 1
        metrics['UnfilteredParalogy'][biotype][len(unfiltered_grouped[tx_id])] += 1

    paralogy_df = pd.DataFrame(paralogy_df, columns=['TranscriptId', 'Paralogy', 'UnfilteredParalogy'])
    
    tmp_verbose = os.path.join(temp_dir, 'verbose.txt')
    try:
        cmd = ['pslCDnaFilter', '-verbose=5', tmp_size_filtered, '/dev/stdout']
        tools.procOps.run_proc(cmd, stderr=tmp_verbose, stdout='/dev/null')
        scores = parse_verbose(tmp_verbose)
    finally:
        if os.path.exists(tmp_verbose):
            os.remove(tmp_verbose)
    
    global_best_by_gene = tools.transcripts.group_transcripts_by_name2(global_best_txs)
    coding_genes = {gene_id for gene_id, tx_list in global_best_by_gene.items()
                    if any(x.cds_size > 0 for x in tx_list)}

    # Use local temp files for clustering
    coding_tmp = os.path.join(temp_dir, 'coding_tmp.txt')
    noncoding_tmp = os.path.join(temp_dir, 'noncoding_tmp.txt')
    coding_clusters = os.path.join(temp_dir, 'coding_clusters.gp')
    noncoding_clusters = os.path.join(temp_dir, 'noncoding_clusters.gp')
    try:
            with open(coding_clusters, 'w') as out_coding, open(noncoding_clusters, 'w') as out_noncoding:
                for tx in global_best_txs:
                    if tx.name2 in coding_genes:
                        tools.fileOps.print_row(out_coding, tx.get_gene_pred())
                    else:
                        tools.fileOps.print_row(out_noncoding, tx.get_gene_pred())
            cmd = ['clusterGenes', '-cds', f'-ignoreBases={overlapping_ignore_bases}',
                   coding_tmp, 'no', coding_clusters]
            tools.procOps.run_proc(cmd)
            cmd = ['clusterGenes', f'-ignoreBases={overlapping_ignore_bases}',
                   noncoding_tmp, 'no', noncoding_clusters]
            tools.procOps.run_proc(cmd)
            coding_clustered = pd.read_csv(coding_tmp, sep='\t')
            noncoding_clustered = pd.read_csv(noncoding_tmp, sep='\t')
    finally:
        # Clean up temp files
        for f in [coding_tmp, noncoding_tmp, coding_clusters, noncoding_clusters]:
            if os.path.exists(f):
                os.remove(f)

    metrics['Gene Family Collapse'] = collections.defaultdict(lambda: collections.Counter())
    coding_merged_df, coding_collapse_filtered = filter_clusters(coding_clustered, transcript_gene_map,
                                                                 gene_name_map, scores, metrics, gene_biotype_map,
                                                                 filter_overlapping_genes, annotate_extra_paralogs)
    noncoding_merged_df, noncoding_collapse_filtered = filter_clusters(noncoding_clustered, transcript_gene_map,
                                                                       gene_name_map, scores, metrics, gene_biotype_map,
                                                                       filter_overlapping_genes, annotate_extra_paralogs)

    merged_collapse_filtered = pd.concat([coding_collapse_filtered, noncoding_collapse_filtered])
    merged_df = pd.concat([coding_merged_df, noncoding_merged_df])
    high_cov_ids = {x.q_name for x in unfiltered.values() if x.coverage > paralog_rescue_min_coverage and x.q_name in scores}
    high_cov_ids -= set(merged_collapse_filtered.gene)  # gene is alignment ID
    putative_rescue_txs = {tx for aln_id, tx in unfiltered_tx_dict.items() if aln_id in high_cov_ids}
    unfiltered_by_gene = tools.transcripts.group_transcripts_by_name2(putative_rescue_txs)

    rescued_txs = []
    for gene_id, group in merged_collapse_filtered.groupby('gene_id'):
        assert len(set(group['#cluster'])) == 1
        tx_intervals = []
        for _, s in group.iterrows():
            tx_intervals.append(tools.intervals.ChromosomeInterval(s.chrom, s.txStart, s.txEnd, s.strand))
        tx_intervals = tools.intervals.hull_of_intervals(tx_intervals)
        assert tx_intervals is not None
        gene_interval = tx_intervals[0]
        for tx in unfiltered_by_gene[gene_id]:
            if tx.interval.overlap(gene_interval):
                rescued_txs.append(tx.name)
    combined_txs = rescued_txs + list(merged_collapse_filtered.gene)
    combined_tx_df = pd.DataFrame(combined_txs, columns=['AlignmentId'])
    combined_tx_df['score'] = [scores[x] for x in combined_tx_df.AlignmentId]
    combined_tx_df['TranscriptId'] = [tools.nameConversions.strip_alignment_numbers(x) for x in combined_tx_df.AlignmentId]
    combined_tx_df['GeneId'] = [transcript_gene_map[x] for x in combined_tx_df.TranscriptId]
    combined_tx_df = combined_tx_df.sort_values('score')
    combined_tx_df = combined_tx_df.groupby('TranscriptId', as_index=False).first()
    # Ensure GeneId columns are the same type (string) before merging
    # This is especially important when merged_df is empty (common in transMap_pairwise)
    combined_tx_df['GeneId'] = combined_tx_df['GeneId'].astype(str)
    if not merged_df.empty:
        merged_df['GeneId'] = merged_df['GeneId'].astype(str)
    resolved_df = combined_tx_df.merge(merged_df, on='GeneId', how='left')
    resolved_df = resolved_df.drop('score', axis=1)
    with psl_tgt.open('w') as outf:
        for aln_id in resolved_df.AlignmentId:
            aln = unfiltered[aln_id]
            tools.fileOps.print_row(outf, aln.psl_string())
    
    resolved_df, split_gene_metrics = resolve_split_genes(tmp_size_filtered, transcript_gene_map,
                                                          resolved_df, unfiltered_tx_dict, q_size_map)
    # Ensure TranscriptId columns are the same type (string) before merging
    resolved_df['TranscriptId'] = resolved_df['TranscriptId'].astype(str)
    if not paralogy_df.empty:
        paralogy_df['TranscriptId'] = paralogy_df['TranscriptId'].astype(str)
    resolved_df = resolved_df.merge(paralogy_df, on='TranscriptId')
    metrics['Split Genes'] = split_gene_metrics

    os.remove(tmp_size_filtered)
    tools.fileOps.ensure_file_dir(json_tgt)
    with json_tgt.open('w') as outf:
        json.dump(metrics, outf, indent=4)
    
    # Clean up temp directory if we created it
    if cleanup_temp and os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)

    return resolved_df.set_index(['GeneId', 'TranscriptId'])

def ref_span(aln_dict, ref_aln_dict, max_span=5):
    grouped = collections.defaultdict(list)
    for aln_id, aln in aln_dict.items():
        grouped[tools.nameConversions.strip_alignment_numbers(aln_id)].append(aln)
    r = collections.OrderedDict()
    missing_refs = 0
    missing_tx_ids = set()
    for tx_id, aln_list in grouped.items():
        ref_aln = ref_aln_dict.get(tx_id)
        if ref_aln is None:
            missing_refs += len(aln_list)
            missing_tx_ids.add(tx_id)
            continue
        ref_size = ref_aln.t_end - ref_aln.t_start
        ref_cutoff = ref_size * max_span
        alns = [aln for aln in aln_list if aln.t_end - aln.t_start <= ref_cutoff]
        for aln in alns:
            r[aln.q_name] = aln
    if missing_refs:
        logging.warning(
            "ref_span: skipped %d alignments lacking reference entries (examples: %s)",
            missing_refs,
            ', '.join(sorted(list(missing_tx_ids))[:5])
        )
    return r, len(aln_dict) - len(r)


def parse_stats(stats):
    stats = pd.read_csv(stats, sep='\t', names=['mode', 'seqs', 'alns'], index_col=0)
    stats.index = [x.replace(' ', '') for x in stats.index]
    stats = stats.T
    stats_dict = {}
    if 'dropminCover:' in stats:
        stats_dict['Coverage Filter'] = int(stats['dropminCover:'].alns)
    else:
        stats_dict['Coverage Filter'] = 0
    if 'dropminSpan:' in stats:
        stats_dict['Min Span Distance'] = int(stats['dropminSpan:'].alns)
    else:
        stats_dict['Min Span Distance'] = 0
    if 'dropglobalBest:' in stats:
        stats_dict['Paralog Filter'] = int(stats['dropglobalBest:'].alns)
    else:
        stats_dict['Paralog Filter'] = 0
    return stats_dict


def parse_verbose(verbose):
    scores = {}
    for l in open(verbose):
        if l.startswith('align'):
            l = l.split()
            aln_id = l[-3].rsplit(':', 1)[0].split(']')[1]
            score = l[5]
            score = float(score.split('=')[1])
            scores[aln_id] = score
    return scores


def find_best_group(group, key):
    avg_scores = group[[key, 'scores']].groupby(key, as_index=False).mean()
    if abs(avg_scores.sort_values('scores', ascending=False)['scores'][0] - avg_scores.sort_values('scores', ascending=False)['scores'][1]) < 0.01:
        return [avg_scores.sort_values('scores', ascending=False).iloc[0][key], avg_scores.sort_values('scores', ascending=False).iloc[1][key]]
    else:
        return [avg_scores.sort_values('scores', ascending=False).iloc[0][key]]


def construct_alt_loci(group, best_cluster):
    intervals = collections.defaultdict(list)
    for cluster_id, x in group.set_index('#cluster').iterrows():
        if cluster_id != best_cluster:
            intervals[x.chrom].append(tools.intervals.ChromosomeInterval(x.chrom, x.txStart, x.txEnd, '.'))
    merged_intervals = []
    for chrom, i in intervals.items():
        merged_intervals.extend(tools.intervals.gap_merge_intervals(i, 1000))
    return ','.join('{}:{}-{}'.format(x.chromosome, x.start, x.stop) for x in merged_intervals)


def filter_clusters(clustered, transcript_gene_map, gene_name_map, scores, metrics, gene_biotype_map,
                    filter_overlapping_genes, annotate_extra_paralogs):
    clustered['gene_id'] = [transcript_gene_map[tools.nameConversions.strip_alignment_numbers(x)] for x in clustered.gene]
    clustered['scores'] = [scores[x] for x in clustered.gene]
    cluster_done = []
    to_remove_alns = set()
    alt_loci = []
    for gene_id, group in clustered.groupby('gene_id'):
        if len(set(group['#cluster'])) > 1:
            best_clusters_list = find_best_group(group, '#cluster')
            best_clusters = [int(best_cluster) for best_cluster in best_clusters_list]
            best = best_clusters[0]
            if len(best_clusters) > 1:
                if best_clusters[0] in cluster_done:
                    best = best_clusters[1]
            cluster_done.append(best)
            alt_loci.append([gene_id, construct_alt_loci(group, best)])
            bad_clusters= group[group['#cluster'].isin(set(group['#cluster']) - {best})]
            to_remove_alns.update(set(bad_clusters['gene']))

    if len(alt_loci) > 0:
        paralog_df = pd.DataFrame(alt_loci, columns=['GeneId', 'GeneAlternateLoci'])
    else:
        paralog_df = pd.DataFrame(columns=['GeneId', 'GeneAlternateLoci'])
    paralog_filtered = clustered[~clustered['gene'].isin(to_remove_alns)]
    genes_to_remove = set()
    collapsed_genes = []
    for cluster_id, group in paralog_filtered.groupby('#cluster'):
        if len(set(group['gene_id'])) > 1:
            best_genes = find_best_group(group, 'gene_id')
            collapsed_gene_ids = set(group.gene_id)
            gene_biotypes = [gene_biotype_map[best_gene] for best_gene in best_genes]
            for gene_biotype in gene_biotypes:
                metrics['Gene Family Collapse'][gene_biotype][len(collapsed_gene_ids)] += 1
            collapsed_gene_names = {gene_name_map[x] for x in collapsed_gene_ids}
            genes_to_remove.update(collapsed_gene_ids)
            for best_gene in best_genes:
                collapsed_genes.append([best_gene, ','.join(collapsed_gene_ids), ','.join(collapsed_gene_names)])
    if filter_overlapping_genes == True:
        collapse_filtered = paralog_filtered[~paralog_filtered['gene_id'].isin(genes_to_remove)]
    else:
        collapse_filtered = paralog_filtered
    if len(collapsed_genes) > 0:
        collapsed_df = pd.DataFrame(collapsed_genes, columns=['GeneId', 'CollapsedGeneIds', 'CollapsedGeneNames'])
    else:
        collapsed_df = pd.DataFrame(columns=['GeneId', 'CollapsedGeneIds', 'CollapsedGeneNames'])
    # Ensure GeneId is string type in both DataFrames before merging
    # For empty DataFrames, we need to set dtype explicitly by creating a dummy row
    if paralog_df.empty:
        # Create with explicit string dtype by adding and removing a dummy row
        paralog_df = pd.DataFrame([['', '']], columns=['GeneId', 'GeneAlternateLoci']).astype(str).iloc[0:0]
    else:
        paralog_df['GeneId'] = paralog_df['GeneId'].astype(str)
    if collapsed_df.empty:
        # Create with explicit string dtype by adding and removing a dummy row
        collapsed_df = pd.DataFrame([['', '', '']], columns=['GeneId', 'CollapsedGeneIds', 'CollapsedGeneNames']).astype(str).iloc[0:0]
    else:
        collapsed_df['GeneId'] = collapsed_df['GeneId'].astype(str)
    merged_df = collapsed_df.merge(paralog_df, how='outer', on='GeneId')
    return merged_df, collapse_filtered

def find_split_genes(gene_id, g, resolved_interval, split_gene_data):
    intervals = collections.defaultdict(list)
    for aln in g:
        ref_i = tools.intervals.ChromosomeInterval(tools.nameConversions.strip_alignment_numbers(aln.q_name),
                                                   aln.q_start, aln.q_end, '.')
        tgt_i = tools.intervals.ChromosomeInterval(aln.t_name, aln.t_start, aln.t_end, aln.strand)
        intervals[ref_i].append(tgt_i)
    merged_intervals = tools.intervals.union_of_intervals(list(intervals.keys()))
    if len(merged_intervals) > 1:
        alt_intervals = collections.defaultdict(list)
        for interval_list in intervals.values():
            for i in interval_list:
                if not i.overlap(resolved_interval):
                    alt_intervals[i.chromosome].append(i)
        r = []
        for chrom, interval_list in alt_intervals.items():
            r.extend(tools.intervals.gap_merge_intervals(interval_list, 0))
        if len(alt_intervals) == 1 and list(alt_intervals.keys())[0] == resolved_interval.chromosome:
            split_gene_data['intra'].add(gene_id)
        else:
            split_gene_data['contig'].add(gene_id)
        if len(r) == 0:
            return None
        return ','.join(['{}:{}-{}'.format(i.chromosome, i.start, i.stop) for i in r])
    else:
        return None


def resolve_split_genes(tmp_size_filtered, transcript_gene_map, resolved_df, unfiltered_tx_dict, q_size_map=None):
    # Note: This function also needs temp_dir parameter if called
    local_tmp = tempfile.mktemp(suffix='.psl')
    stripped_tmp = tempfile.mktemp(suffix='.psl')
    try:
        if q_size_map is None:
            q_size_map = {}
        with open(stripped_tmp, 'w') as outf:
            for rec in tools.psl.psl_iterator(tmp_size_filtered):
                stripped_name = tools.nameConversions.strip_alignment_numbers(rec.q_name)
                rec.q_name = stripped_name
                max_q_size = q_size_map.get(stripped_name)
                if max_q_size is None or rec.q_size > max_q_size:
                    max_q_size = rec.q_size
                    q_size_map[stripped_name] = max_q_size
                rec.q_size = max_q_size
                tools.fileOps.print_row(outf, rec.psl_string())
        cmd = ['pslCDnaFilter', '-localNearBest=0.05',
               '-minCover=0.1', '-verbose=0',
               '-minSpan=0.2', stripped_tmp, '/dev/stdout']
        tools.procOps.run_proc(cmd, stdout=local_tmp)
        filtered_alns = list(tools.psl.psl_iterator(local_tmp))
    finally:
        for f in [local_tmp, stripped_tmp]:
            if os.path.exists(f):
                os.remove(f)
    
    resolved_ids = set(resolved_df.TranscriptId)
    filtered_alns = [x for x in filtered_alns if x.q_name in resolved_ids]
    grouped = tools.psl.group_alignments_by_qname(filtered_alns, strip=False)
    tx_intervals = {tx_id: unfiltered_tx_dict[aln_id].interval for
                    tx_id, aln_id in zip(resolved_df.TranscriptId, resolved_df.AlignmentId)}

    split_r = []
    split_gene_data = {'contig': set(), 'intra': set()}
    for tx_id, g in grouped.items():
        gene_id = transcript_gene_map[tx_id]
        split_r.append([tx_id, find_split_genes(gene_id, g, tx_intervals[tx_id], split_gene_data)])
    if len(split_r) > 0:
        split_df = pd.DataFrame(split_r, columns=['TranscriptId', 'PossibleSplitGeneLocations'])
    else:
        # Create empty DataFrame with correct string dtype
        split_df = pd.DataFrame([['', None]], columns=['TranscriptId', 'PossibleSplitGeneLocations']).astype({'TranscriptId': str}).iloc[0:0]
    # Ensure TranscriptId columns are the same type (string) before merging
    split_df['TranscriptId'] = split_df['TranscriptId'].astype(str)
    resolved_df['TranscriptId'] = resolved_df['TranscriptId'].astype(str)
    merged = split_df.merge(resolved_df, on='TranscriptId')
    split_gene_metrics = {'Number of contig split genes': len(split_gene_data['contig']),
                          'Number of intra-contig split genes': len(split_gene_data['intra'])}

    return merged, split_gene_metrics

def main():
    """
    Main entry point for the transMap filtering script.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tm-psl", required=True, help="Input transMap PSL file.")
    parser.add_argument("--ref-psl", required=True, help="Reference 'fake' PSL file.")
    parser.add_argument("--tm-gp", required=True, help="Input transMap genePred file.")
    parser.add_argument("--db-path", required=True, help="Path to reference genome database.")
    parser.add_argument("--filtered-psl", required=True, help="Output path for filtered PSL.")
    parser.add_argument("--filtered-gp", required=True, help="Output path for filtered genePred.")
    parser.add_argument("--metrics-json", required=True, help="Output path for filtering metrics JSON.")
    parser.add_argument("--resolved-df", required=True, help="Output path for resolved DataFrame (Pickle).")
    parser.add_argument("--global-near-best", type=float, default=0.1, help="globalNearBest value for pslCDnaFilter.")
    parser.add_argument("--filter-overlapping-genes", action='store_true', help="Filter out overlapping genes from collapsed clusters.")
    parser.add_argument("--overlapping-ignore-bases", type=int, default=0, help="Allowed overlap for clusterGenes.")
    parser.add_argument("--annotate-extra-paralogs", action='store_true', help="Annotate paralogs found outside the global near best set.")
    parser.add_argument("--min-cover", type=float, default=0.1, help="pslCDnaFilter -minCover (min fraction of query aligned). Lower to keep more alignments.")
    parser.add_argument("--min-span", type=float, default=0.2, help="pslCDnaFilter -minSpan (min fraction of query span). Lower to keep more alignments.")
    parser.add_argument("--max-ref-span", type=float, default=5, help="Drop alignments whose target span exceeds this multiple of the reference span. Raise to keep more.")
    parser.add_argument("--paralog-rescue-min-coverage", type=float, default=0.5, help="Min coverage (0-1) for rescuing a transcript collapsed by overlapping-gene filtering. Lower to rescue more paralogs.")
    parser.add_argument("--temp-dir", default=None, help="Directory for temporary files (default: auto-create in output directory)")

    args = parser.parse_args()
    
    # Create temp directory in the output directory if not specified
    if args.temp_dir is None:
        output_dir = os.path.dirname(os.path.abspath(args.filtered_psl))
        args.temp_dir = os.path.join(output_dir, f'.filter_transmap_temp_{os.getpid()}')
        os.makedirs(args.temp_dir, exist_ok=True)
        cleanup_temp_dir = True
    else:
        cleanup_temp_dir = False

    # The run_filtering function returns a dataframe that needs to be used.
    try:
        resolved_df = run_filtering(
            tm_psl=args.tm_psl,
            ref_psl=args.ref_psl,
            tm_gp=args.tm_gp,
            db_path=args.db_path,
            psl_tgt=Path(args.filtered_psl),
            global_near_best=args.global_near_best,
            filter_overlapping_genes=args.filter_overlapping_genes,
            overlapping_ignore_bases=args.overlapping_ignore_bases,
            json_tgt=Path(args.metrics_json),
            annotate_extra_paralogs=args.annotate_extra_paralogs,
            temp_dir=args.temp_dir,
            min_cover=args.min_cover,
            min_span=args.min_span,
            max_ref_span=args.max_ref_span,
            paralog_rescue_min_coverage=args.paralog_rescue_min_coverage,
        )
        resolved_df.to_pickle(args.resolved_df)
    finally:
        # Clean up temp directory if we created it
        if cleanup_temp_dir and os.path.exists(args.temp_dir):
            shutil.rmtree(args.temp_dir)

if __name__ == "__main__":
    main()
