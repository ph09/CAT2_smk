#!/usr/bin/env python3
"""Cluster-array-job-based parent-gene assignment 
"""

import argparse
import collections
import itertools
import os
import pickle
import sqlite3
import sys
from pathlib import Path

import pandas as pd

import tools.fileOps
import tools.intervals
import tools.mathOps
import tools.transcripts
from tools.defaultOrderedDict import DefaultOrderedDict

from cat.scheduler import get_scheduler


# ──────────────────────────────────────────────────────────────────────────────
# Worker invoked inside each array task
# ──────────────────────────────────────────────────────────────────────────────


def process_chromosome_worker(chrom, denovo_txs_pkl, tm_txs_pkl, filtered_ids_pkl, min_distance, output_file):
    """Process a single chromosome (called by one array task)."""
    with open(denovo_txs_pkl, "rb") as f:
        denovo_txs_on_chrom = pickle.load(f)
    with open(tm_txs_pkl, "rb") as f:
        tm_txs_on_chrom = pickle.load(f)
    with open(filtered_ids_pkl, "rb") as f:
        filtered_ids = pickle.load(f)

    records = []
    for denovo_tx in denovo_txs_on_chrom.values():
        unfiltered_overlapping_txs = find_tm_overlaps(denovo_tx, tm_txs_on_chrom)

        filtered_overlapping_txs = {tx for tx in unfiltered_overlapping_txs if tx.name not in filtered_ids}
        filtered_gene_ids = {tx.name2 for tx in filtered_overlapping_txs}

        resolved_name, resolution_method = None, None
        if len(filtered_gene_ids) > 1:
            resolved_name, resolution_method = resolve_multiple_genes(denovo_tx, filtered_overlapping_txs, min_distance)
        elif len(filtered_gene_ids) == 1:
            resolved_name = list(filtered_gene_ids)[0]

        alternative_gene_ids = {tx.name2 for tx in unfiltered_overlapping_txs} - {resolved_name}
        alternative_genes_str = ','.join(sorted(alternative_gene_ids)) if alternative_gene_ids else ''

        records.append({
            'TranscriptId': denovo_tx.name,
            'AssignedGeneId': resolved_name if resolved_name else '',
            'AlternativeGeneIds': alternative_genes_str,
            'ResolutionMethod': resolution_method if resolution_method else '',
        })

    with open(output_file, "wb") as f:
        pickle.dump(records, f)

    print(f"Processed chromosome {chrom}: {len(records)} transcripts")
    return len(records)


# ──────────────────────────────────────────────────────────────────────────────
# Job script + submission
# ──────────────────────────────────────────────────────────────────────────────


def _build_body(scheduler, work_dir, sentinel_dir, min_distance):
    """Bash body run by every array task.

    1-based at the scheduler level; IDX shift to 0-based for Python list
    lookup against chromosomes.pkl.
    """
    task_var = scheduler.task_id_env()
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    preamble = scheduler.script_preamble(conda_env=os.environ.get("CONDA_DEFAULT_ENV", "cat") or "cat")
    sentinel = scheduler.trap_sentinel(sentinel_dir)
    return f"""{preamble}
{sentinel}

TASK_ID="${{{task_var}:-1}}"
IDX=$((TASK_ID - 1))
echo "Processing chromosome index $IDX at $(date)"

python3 - <<PYEOF
import os, pickle, sys
sys.path.insert(0, {repo_root!r})
from cat.parent_gene_assignment_cluster import process_chromosome_worker

idx = int(os.environ.get('{task_var}', '1')) - 1
with open({str(work_dir / 'chromosomes.pkl')!r}, 'rb') as f:
    chromosomes = pickle.load(f)

chrom = chromosomes[idx]
print(f'Processing chromosome: {{chrom}}')

process_chromosome_worker(
    chrom,
    {str(work_dir / 'chromosomes')!r} + f'/{{chrom}}_denovo.pkl',
    {str(work_dir / 'chromosomes')!r} + f'/{{chrom}}_tm.pkl',
    {str(work_dir / 'filtered_ids.pkl')!r},
    {min_distance},
    {str(work_dir / 'results')!r} + f'/result_{{chrom}}.pkl',
)
PYEOF

echo "Finished chromosome index $IDX at $(date)"
"""


def run_cluster_parent_assignment(args):
    """Coordinate cluster array job submission across chromosomes."""
    scheduler = get_scheduler(args.execution_mode, _scheduler_config(args))

    print("Loading transcript data...")
    filtered_transmap_dict = tools.transcripts.get_gene_pred_dict(args.filtered_tm_gp)
    unfiltered_transmap_dict = tools.transcripts.get_gene_pred_dict(args.unfiltered_tm_gp)
    denovo_dict = tools.transcripts.get_gene_pred_dict(args.denovo_gp)

    if not denovo_dict:
        print(f"Warning: No de novo gene predictions found in {args.denovo_gp}")
        df = pd.DataFrame(columns=['TranscriptId', 'AssignedGeneId', 'AlternativeGeneIds', 'ResolutionMethod']).astype(str)
        write_to_database(df, args.db_path, args.table_name)
        return

    tm_chrom_dict = create_chrom_dict(unfiltered_transmap_dict, args.chrom_sizes)
    denovo_chrom_dict = create_chrom_dict(denovo_dict)
    filtered_ids = set(unfiltered_transmap_dict.keys()) - set(filtered_transmap_dict.keys())

    db_dir = os.path.dirname(os.path.abspath(args.db_path))
    work_dir = Path(db_dir) / f".parent_assignment_{os.path.basename(args.db_path).replace('.db', '')}_{args.table_name}"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "chromosomes").mkdir(exist_ok=True)
    (work_dir / "results").mkdir(exist_ok=True)
    (work_dir / "cluster_logs").mkdir(exist_ok=True)
    sentinel_dir = work_dir / "sentinels"
    sentinel_dir.mkdir(exist_ok=True)

    print(f"Work directory: {work_dir}")

    try:
        chromosomes = list(denovo_chrom_dict.keys())
        with open(work_dir / "chromosomes.pkl", "wb") as f:
            pickle.dump(chromosomes, f)

        with open(work_dir / "filtered_ids.pkl", "wb") as f:
            pickle.dump(filtered_ids, f)

        for chrom in chromosomes:
            with open(work_dir / f"chromosomes/{chrom}_denovo.pkl", "wb") as f:
                pickle.dump(denovo_chrom_dict[chrom], f)
            with open(work_dir / f"chromosomes/{chrom}_tm.pkl", "wb") as f:
                pickle.dump(tm_chrom_dict.get(chrom, {}), f)

        num_chroms = len(chromosomes)
        print(f"Processing {num_chroms} chromosomes via {scheduler.name} array jobs")

        log_out, log_err = scheduler.array_log_paths(work_dir / "cluster_logs", "parent")
        header = scheduler.header(
            job_name=f"parent_assign_{args.table_name}",
            cpus=args.cpus,
            mem=f"{args.memory}G",
            walltime=args.time,
            log_out=log_out,
            log_err=log_err,
            partition=args.partition,
            queue=args.partition,
            array=(1, num_chroms),
            max_concurrent=args.max_jobs,
        )
        body = _build_body(scheduler, work_dir, sentinel_dir, args.min_distance)
        script_path = scheduler.write_script(header + body, work_dir / "cluster_job.sh")

        job_id = scheduler.submit(script_path)
        print(f"Submitted {scheduler.name} array job {job_id}")

        result = scheduler.wait(
            job_id,
            num_tasks=num_chroms,
            timeout_s=args.timeout_hours * 3600,
            sentinel_dir=sentinel_dir,
        )
        if not result.ok:
            raise RuntimeError(f"{scheduler.name} job {job_id} failed: {result.detail}")
        print(f"Job {job_id} succeeded: {result.completed}/{result.total} tasks")

        all_records = []
        for chrom in chromosomes:
            result_file = work_dir / f"results/result_{chrom}.pkl"
            if result_file.exists():
                with open(result_file, "rb") as f:
                    records = pickle.load(f)
                    all_records.extend(records)
            else:
                print(f"Warning: Missing result file for chromosome {chrom}")

        if all_records:
            df = pd.DataFrame(all_records)
            for col in df.columns:
                df[col] = df[col].astype(str)
        else:
            df = pd.DataFrame(columns=['TranscriptId', 'AssignedGeneId', 'AlternativeGeneIds', 'ResolutionMethod']).astype(str)

        print(f"Completed processing: {len(df)} total transcript assignments")
        write_to_database(df, args.db_path, args.table_name)

    finally:
        if args.cleanup:
            import shutil
            try:
                shutil.rmtree(work_dir)
                print(f"Cleaned up work directory: {work_dir}")
            except Exception as e:  # noqa: BLE001
                print(f"Warning: Could not clean up work directory {work_dir}: {e}")


def _scheduler_config(args):
    """Build a per-backend config dict from CLI args (see classify_cluster._scheduler_config)."""
    return {
        "cluster": {
            "slurm": {
                "partition": args.partition,
                "exclude_nodes": args.exclude_nodes or "",
                "module_load": args.module_load or "",
            },
            "sge": {
                "queue": args.partition,
                "parallel_env": args.sge_parallel_env,
                "memory_flag": args.sge_memory_flag,
                "hostname_exclude": args.exclude_nodes or "",
                "module_load": args.module_load or "",
            },
        }
    }


# ──────────────────────────────────────────────────────────────────────────────
# Helpers unchanged from parent_gene_assignment_slurm.py
# ──────────────────────────────────────────────────────────────────────────────


def write_to_database(df, db_path, table_name):
    """Write *df* to *db_path* under *table_name* (replace semantics)."""
    con = sqlite3.connect(db_path)
    try:
        df_clean = df.copy()
        for col in ['AssignedGeneId', 'AlternativeGeneIds', 'ResolutionMethod']:
            if col in df_clean.columns:
                df_clean[col] = df_clean[col].fillna('')
        df_clean.columns = [col.replace('-', '_').replace(' ', '_') for col in df_clean.columns]
        df_clean.to_sql(table_name, con, if_exists='replace', index=False)
        print(f"Successfully wrote {len(df_clean)} records to database table '{table_name}'")
    finally:
        con.close()


def create_chrom_dict(tx_dict, chrom_sizes=None):
    """Group *tx_dict* by chromosome; ensure every chrom_sizes entry has a (possibly empty) bucket."""
    chrom_dict = collections.defaultdict(dict)
    for tx_id, tx in tx_dict.items():
        chrom_dict[tx.chromosome][tx_id] = tx
    if chrom_sizes is not None:
        for chrom, size in tools.fileOps.iter_lines(chrom_sizes):
            if chrom not in chrom_dict:
                chrom_dict[chrom] = {}
    return chrom_dict


def find_tm_overlaps(denovo_tx, tm_tx_dict, cutoff=100):
    """Find transMap transcripts that overlap *denovo_tx* by at least *cutoff* bases."""
    r = DefaultOrderedDict(int)
    denovo_start = denovo_tx.start
    denovo_stop = denovo_tx.stop

    for tx in tm_tx_dict.values():
        # Genomic-overlap pre-filter is much faster than per-exon checking.
        if tx.stop <= denovo_start or tx.start >= denovo_stop:
            continue
        for tx_exon in tx.exon_intervals:
            for denovo_exon in denovo_tx.exon_intervals:
                i = tx_exon.intersection(denovo_exon)
                if i is not None:
                    r[tx] += len(i)
    return [tx_id for tx_id, num_bases in r.items() if num_bases >= cutoff]


def resolve_multiple_genes(denovo_tx, overlapping_tm_txs, min_distance):
    """Disambiguate the parent gene when *denovo_tx* overlaps several transMap genes."""
    tm_txs_by_gene = tools.transcripts.group_transcripts_by_name2(overlapping_tm_txs)
    tm_jaccards = [find_highest_gene_jaccard(x, y) for x, y in itertools.combinations(list(tm_txs_by_gene.values()), 2)]
    if any(x > 0.001 for x in tm_jaccards):
        return None, 'badAnnotOrTm'
    scores = collections.defaultdict(list)
    for tx in overlapping_tm_txs:
        scores[tx.name2].append(calculate_asymmetric_closeness(denovo_tx, tx))
    best_scores = {gene_id: max(scores[gene_id]) for gene_id in scores}
    high_score = max(best_scores.values())
    if all(high_score - x >= min_distance for x in best_scores.values() if x != high_score):
        best = sorted(iter(best_scores.items()), key=lambda gene_id_score: gene_id_score[1])[-1][0]
        return best, 'rescued'
    return None, 'ambiguousOrFusion'


def find_highest_gene_jaccard(gene_list_a, gene_list_b):
    """Jaccard between exonic footprints of two groups of transcripts."""
    def find_interval(gene_list):
        gene_intervals = set()
        for tx in gene_list:
            gene_intervals.update(tx.exon_intervals)
        return tools.intervals.gap_merge_intervals(gene_intervals, 0)

    a_interval = find_interval(gene_list_a)
    b_interval = find_interval(gene_list_b)
    return tools.intervals.calculate_bed12_jaccard(a_interval, b_interval)


def calculate_asymmetric_closeness(denovo_tx, tm_tx):
    intersection = denovo_tx.interval.intersection(tm_tx.interval)
    if intersection is None:
        return 0
    return tools.mathOps.format_ratio(len(intersection), len(denovo_tx.interval))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--filtered-tm-gp", required=True)
    parser.add_argument("--unfiltered-tm-gp", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--denovo-gp", required=True)
    parser.add_argument("--db-path", required=True)
    parser.add_argument("--table-name", required=True)
    parser.add_argument("--min-distance", type=float, default=0.9)
    parser.add_argument("--execution-mode", choices=("auto", "slurm", "sge", "local"), default="auto")
    parser.add_argument("--partition", default="",
                        help="SLURM partition or SGE queue.")
    parser.add_argument("--exclude-nodes", default="")
    parser.add_argument("--module-load", default="")
    parser.add_argument("--sge-parallel-env", default="smp")
    parser.add_argument("--sge-memory-flag", default="h_vmem")
    parser.add_argument("--memory", type=int, default=8)
    parser.add_argument("--cpus", type=int, default=1)
    parser.add_argument("--time", default="12:00:00")
    parser.add_argument("--max-jobs", type=int, default=20)
    parser.add_argument("--timeout-hours", type=int, default=24)
    parser.add_argument("--cleanup", action="store_true")

    args = parser.parse_args()
    run_cluster_parent_assignment(args)


if __name__ == "__main__":
    main()
