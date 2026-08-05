#!/usr/bin/env python3
"""Cluster-array-job-based transcript evaluation for the CAT pipeline.
"""

import argparse
import collections
import os
import pickle
import sys
import time
from pathlib import Path

import pandas as pd

import tools.bio
import tools.mathOps
import tools.psl
import tools.sqlInterface
import tools.transcripts

from cat.scheduler import get_scheduler

# Distance allowed between intron locations to be considered equivalent.
FUZZ_DISTANCE = 7


# ──────────────────────────────────────────────────────────────────────────────
# Worker invoked inside each array task
# ──────────────────────────────────────────────────────────────────────────────


def process_psl_chunk_worker(chunk_file, ref_tx_dict_file, tx_dict_file,
                              tx_biotype_map_file, fasta_path, aln_mode, output_file):
    """Process a single chunk of PSL records (called by one array task)."""
    with open(chunk_file, "rb") as f:
        psl_chunk = pickle.load(f)
    with open(ref_tx_dict_file, "rb") as f:
        ref_tx_dict = pickle.load(f)
    with open(tx_dict_file, "rb") as f:
        tx_dict = pickle.load(f)
    with open(tx_biotype_map_file, "rb") as f:
        tx_biotype_map = pickle.load(f)
    seq_dict = tools.bio.get_sequence_dict(fasta_path)

    metrics_records = []
    eval_records = []

    for psl in psl_chunk:
        try:
            ref_tx = ref_tx_dict[psl.t_name]
            tx = tx_dict[psl.q_name]
            biotype = tx_biotype_map[psl.t_name]

            original_intron_vector = calculate_original_intron_vector(ref_tx, tx, psl, aln_mode)
            adj_start, adj_stop = find_adj_start_stop(tx, seq_dict)

            metrics_records.append([ref_tx.name2, ref_tx.name, tx.name, 'AlnCoverage', 100 * psl.target_coverage])
            metrics_records.append([ref_tx.name2, ref_tx.name, tx.name, 'AlnIdentity', 100 * psl.identity])
            metrics_records.append([ref_tx.name2, ref_tx.name, tx.name, 'AlnGoodness', 100 * (1 - psl.badness)])
            metrics_records.append([ref_tx.name2, ref_tx.name, tx.name, 'PercentUnknownBases', psl.percent_n])
            metrics_records.append([ref_tx.name2, ref_tx.name, tx.name, 'OriginalIntrons', original_intron_vector])
            metrics_records.append([ref_tx.name2, ref_tx.name, tx.name, 'ValidStart', tools.transcripts.has_start_codon(seq_dict, tx)])
            metrics_records.append([ref_tx.name2, ref_tx.name, tx.name, 'ValidStop', tools.transcripts.has_stop_codon(seq_dict, tx)])
            metrics_records.append([ref_tx.name2, ref_tx.name, tx.name, 'ProperOrf', tx.cds_size % 3 == 0])
            metrics_records.append([ref_tx.name2, ref_tx.name, tx.name, 'AdjStart', adj_start])
            metrics_records.append([ref_tx.name2, ref_tx.name, tx.name, 'AdjStop', adj_stop])

            eval_records.extend(find_indels(tx, psl, aln_mode))
            if biotype == 'protein_coding':
                line = in_frame_stop(tx, seq_dict)
                if line is not None:
                    eval_records.append(line)
        except Exception as e:  # noqa: BLE001
            print(f"Error processing transcript {psl.q_name}: {e}")
            continue

    with open(output_file, "wb") as f:
        pickle.dump({'metrics': metrics_records, 'evaluation': eval_records}, f)

    return len(metrics_records), len(eval_records)


def chunk_list(lst, chunk_size):
    """Split *lst* into chunks of size *chunk_size*."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]


# ──────────────────────────────────────────────────────────────────────────────
# Build the per-mode job script and submit it via Scheduler
# ──────────────────────────────────────────────────────────────────────────────


def _build_body(scheduler, work_dir, aln_mode, sentinel_dir):
    """Bash body that the array tasks all execute.

    Note: array indexing is 1-based at the scheduler level. The body shifts
    to 0-based ``IDX`` for indexing into chunk_${IDX}.pkl, matching the
    historical filename convention.
    """
    task_var = scheduler.task_id_env()
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    preamble = scheduler.script_preamble(conda_env=os.environ.get("CONDA_DEFAULT_ENV", "cat") or "cat")
    sentinel = scheduler.trap_sentinel(sentinel_dir)
    return f"""{preamble}
{sentinel}

TASK_ID="${{{task_var}:-1}}"
IDX=$((TASK_ID - 1))
echo "Processing chunk $IDX at $(date)"

python3 - <<PYEOF
import os, sys
sys.path.insert(0, {repo_root!r})
from cat.classify_cluster import process_psl_chunk_worker

idx = int(os.environ.get('{task_var}', '1')) - 1
fasta_path = open({str(work_dir / 'fasta_path.txt')!r}).read().strip()
process_psl_chunk_worker(
    {str(work_dir / 'chunks')!r} + f'/chunk_{{idx}}.pkl',
    {str(work_dir / 'ref_tx_dict.pkl')!r},
    {str(work_dir / 'tx_dict.pkl')!r},
    {str(work_dir / 'tx_biotype_map.pkl')!r},
    fasta_path,
    {aln_mode!r},
    {str(work_dir / 'results')!r} + f'/result_{{idx}}.pkl',
)
PYEOF

echo "Finished chunk $IDX at $(date)"
"""


def process_mode_with_cluster(tx_mode, path_dict, ref_tx_dict, tx_biotype_map,
                              fasta_path, args, results, scheduler):
    """Process one transcript mode using cluster array jobs."""
    print(f"Processing mode: {tx_mode}")

    tx_dict = tools.transcripts.get_gene_pred_dict(path_dict['gp'])

    for aln_mode in ['CDS', 'mRNA']:
        psl_path = path_dict.get(aln_mode)
        psl_list = list(tools.psl.psl_iterator(psl_path))

        print(f"  Processing {aln_mode} alignment mode with {len(psl_list)} records")

        if len(psl_list) == 0:
            print("    No PSL records found, skipping")
            metrics_tbl_name = tools.sqlInterface.tables[aln_mode][tx_mode]['metrics'].__tablename__
            eval_tbl_name = tools.sqlInterface.tables[aln_mode][tx_mode]['evaluation'].__tablename__

            mc_df = pd.DataFrame(columns=['GeneId', 'TranscriptId', 'AlignmentId', 'classifier', 'value']).set_index('AlignmentId')
            ec_df = pd.DataFrame(columns=['AlignmentId', 'chromosome', 'start', 'stop', 'name', 'score',
                                          'strand', 'thickStart', 'thickStop', 'rgb', 'blockCount',
                                          'blockSizes', 'blockStarts']).set_index('AlignmentId')
            results.append((metrics_tbl_name, mc_df))
            results.append((eval_tbl_name, ec_df))
            continue

        db_dir = os.path.dirname(os.path.abspath(args.db_path))
        work_dir = Path(db_dir) / f".classify_{tx_mode}_{aln_mode}_{os.path.basename(args.db_path).replace('.db', '')}"
        work_dir.mkdir(parents=True, exist_ok=True)
        (work_dir / "chunks").mkdir(exist_ok=True)
        (work_dir / "results").mkdir(exist_ok=True)
        (work_dir / "cluster_logs").mkdir(exist_ok=True)
        sentinel_dir = work_dir / "sentinels"
        sentinel_dir.mkdir(exist_ok=True)

        print(f"    Work directory: {work_dir}")

        try:
            with open(work_dir / "ref_tx_dict.pkl", "wb") as f:
                pickle.dump(ref_tx_dict, f)
            with open(work_dir / "tx_dict.pkl", "wb") as f:
                pickle.dump(tx_dict, f)
            with open(work_dir / "tx_biotype_map.pkl", "wb") as f:
                pickle.dump(tx_biotype_map, f)
            with open(work_dir / "fasta_path.txt", "w") as f:
                f.write(os.path.abspath(fasta_path))

            chunks = list(chunk_list(psl_list, args.chunk_size))
            print(f"    Split into {len(chunks)} chunks of size ~{args.chunk_size}")

            for chunk_id, chunk in enumerate(chunks):
                with open(work_dir / f"chunks/chunk_{chunk_id}.pkl", "wb") as f:
                    pickle.dump(chunk, f)

            # 1-based array; chunk filenames stay 0-based (IDX shift in body).
            num_chunks = len(chunks)
            log_out, log_err = scheduler.array_log_paths(work_dir / "cluster_logs", "eval")
            header = scheduler.header(
                job_name=f"evaluate_{tx_mode}_{aln_mode}",
                cpus=args.cpus,
                mem=f"{args.memory}G",
                walltime=args.time,
                log_out=log_out,
                log_err=log_err,
                partition=args.partition,
                queue=args.partition,
                array=(1, num_chunks),
                max_concurrent=args.max_jobs,
            )
            body = _build_body(scheduler, work_dir, aln_mode, sentinel_dir)
            script_path = scheduler.write_script(header + body, work_dir / f"cluster_job_{aln_mode}.sh")

            job_id = scheduler.submit(script_path)
            print(f"    Submitted {scheduler.name} array job {job_id}")

            result = scheduler.wait(
                job_id,
                num_tasks=num_chunks,
                timeout_s=args.timeout_hours * 3600,
                sentinel_dir=sentinel_dir,
            )
            if not result.ok:
                raise RuntimeError(
                    f"{scheduler.name} job {job_id} failed for {tx_mode}/{aln_mode}: {result.detail}"
                )
            print(f"    Job {job_id} succeeded: {result.completed}/{result.total} tasks")

            all_metrics = []
            all_eval = []
            for chunk_id in range(num_chunks):
                result_file = work_dir / f"results/result_{chunk_id}.pkl"
                if result_file.exists():
                    with open(result_file, "rb") as f:
                        chunk_results = pickle.load(f)
                        all_metrics.extend(chunk_results['metrics'])
                        all_eval.extend(chunk_results['evaluation'])
                else:
                    print(f"    Warning: Missing result file for chunk {chunk_id}")

            columns = ['GeneId', 'TranscriptId', 'AlignmentId', 'classifier', 'value']
            mc_df = pd.DataFrame(all_metrics, columns=columns).sort_values(columns).set_index('AlignmentId')

            columns = ['AlignmentId', 'chromosome', 'start', 'stop', 'name', 'score', 'strand',
                       'thickStart', 'thickStop', 'rgb', 'blockCount', 'blockSizes', 'blockStarts']
            ec_df = pd.DataFrame(all_eval, columns=columns).sort_values(columns).set_index('AlignmentId')

            metrics_tbl_name = tools.sqlInterface.tables[aln_mode][tx_mode]['metrics'].__tablename__
            eval_tbl_name = tools.sqlInterface.tables[aln_mode][tx_mode]['evaluation'].__tablename__
            results.append((metrics_tbl_name, mc_df))
            results.append((eval_tbl_name, ec_df))

            print(f"    Completed: {len(mc_df)} metrics records, {len(ec_df)} evaluation records")

        finally:
            if args.cleanup:
                import shutil
                try:
                    shutil.rmtree(work_dir)
                    print(f"    Cleaned up work directory: {work_dir}")
                except Exception as e:  # noqa: BLE001
                    print(f"    Warning: Could not clean up work directory {work_dir}: {e}")


def run_cluster_classification(args):
    """Coordinate cluster array job submission across modes."""
    scheduler = get_scheduler(args.execution_mode, _scheduler_config(args))

    ref_tx_dict = tools.transcripts.get_gene_pred_dict(args.annotation_gp)
    tx_biotype_map = tools.sqlInterface.get_transcript_biotype_map(args.ref_db_path)

    transcript_modes = collections.defaultdict(dict)
    for mode, gp_path, mrna_path, cds_path in args.mode_files:
        transcript_modes[mode]['gp'] = gp_path
        transcript_modes[mode]['mRNA'] = mrna_path
        transcript_modes[mode]['CDS'] = cds_path

    results = []
    for tx_mode, path_dict in transcript_modes.items():
        process_mode_with_cluster(tx_mode, path_dict, ref_tx_dict, tx_biotype_map,
                                  args.fasta, args, results, scheduler)

    with open(args.resolved_df, "wb") as f:
        pickle.dump(results, f)

    print(f"Successfully saved results to {args.resolved_df}")
    return results


def _scheduler_config(args):
    """Build a per-backend config dict from CLI args.

    Only the keys the scheduler reads at construction time matter
    (partition / queue, exclude_nodes, hostname_exclude, parallel_env,
    memory_flag, module_load). Higher-level args (cpus, memory, time, array
    geometry) are passed directly to header() at submission time.
    """
    cfg = {"cluster": {}}
    cfg["cluster"]["slurm"] = {
        "partition": args.partition,
        "exclude_nodes": args.exclude_nodes or "",
        "module_load": args.module_load or "",
    }
    cfg["cluster"]["sge"] = {
        "queue": args.partition,
        "parallel_env": args.sge_parallel_env,
        "memory_flag": args.sge_memory_flag,
        "hostname_exclude": args.exclude_nodes or "",
        "module_load": args.module_load or "",
    }
    return cfg


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation-gp", required=True)
    parser.add_argument("--ref-db-path", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--db-path", required=True)
    parser.add_argument("--resolved-df", required=True)
    parser.add_argument("--mode-files", nargs=4, action='append', required=True,
                        metavar=("MODE", "INPUT_GP", "MRNA_PSL", "CDS_PSL"))
    parser.add_argument("--execution-mode", choices=("auto", "slurm", "sge", "local"), default="auto")
    parser.add_argument("--partition", default="",
                        help="SLURM partition or SGE queue. Both backends read this.")
    parser.add_argument("--exclude-nodes", default="",
                        help="Comma list (SLURM) or SGE-native '!h1&!h2' expression.")
    parser.add_argument("--module-load", default="")
    parser.add_argument("--sge-parallel-env", default="smp",
                        help="SGE parallel environment name (site-specific).")
    parser.add_argument("--sge-memory-flag", default="h_vmem",
                        help="SGE memory resource flag (h_vmem vs mem_free vs s_vmem).")
    parser.add_argument("--memory", type=int, default=8)
    parser.add_argument("--cpus", type=int, default=2)
    parser.add_argument("--time", default="01:00:00")
    parser.add_argument("--max-jobs", type=int, default=25,
                        help="Maximum concurrent array tasks (SLURM %% suffix / SGE -tc).")
    parser.add_argument("--timeout-hours", type=int, default=24)
    parser.add_argument("--chunk-size", type=int, default=500)
    parser.add_argument("--cleanup", action="store_true")

    args = parser.parse_args()
    run_cluster_classification(args)


# ──────────────────────────────────────────────────────────────────────────────
# Per-transcript helpers (unchanged from classify_slurm.py)
# ──────────────────────────────────────────────────────────────────────────────


def calculate_original_intron_vector(ref_tx, tx, psl, aln_mode):
    if len(ref_tx.intron_intervals) == 0:
        return None
    ref_introns = get_intron_coordinates(ref_tx, aln_mode)
    tgt_introns = []
    for intron in get_intron_coordinates(tx, aln_mode):
        p = psl.query_coordinate_to_target(intron)
        if p is not None:
            tgt_introns.append(p)
    if len(tgt_introns) == 0:
        return ','.join(['0'] * len(ref_tx.intron_intervals))
    intron_vector = []
    for ref_intron in ref_introns:
        closest = tools.mathOps.find_closest(tgt_introns, ref_intron)
        if closest - FUZZ_DISTANCE < ref_intron < closest + FUZZ_DISTANCE:
            intron_vector.append(1)
        else:
            intron_vector.append(0)
    return ','.join(map(str, intron_vector))


def in_frame_stop(tx, fasta):
    for start_pos, stop_pos, codon in tx.codon_iterator(fasta):
        if tools.bio.translate_sequence(codon) == '*':
            bed = tx.get_bed(new_start=start_pos, new_stop=stop_pos, rgb='135,78,191', name='InFrameStop')
            return [tx.name] + bed


def find_adj_start_stop(tx, fasta):
    for start_pos, stop_pos, codon in tx.codon_iterator(fasta):
        if tools.bio.translate_sequence(codon) == '*':
            if tx.strand == '-':
                start = start_pos
                stop = tx.thick_stop
            else:
                stop = stop_pos
                start = tx.thick_start
            return start, stop
    return tx.thick_start, tx.thick_stop


def find_indels(tx, psl, aln_mode):
    def convert_coordinates_to_chromosome(left_pos, right_pos, coordinate_fn, strand):
        left_chrom_pos = coordinate_fn(left_pos)
        if left_chrom_pos is None:
            return None, None
        assert left_chrom_pos is not None
        right_chrom_pos = coordinate_fn(right_pos)
        if right_chrom_pos is None:
            assert aln_mode == "CDS"
            right_chrom_pos = coordinate_fn(tx.cds_size - 1)
        assert right_chrom_pos is not None
        if strand == '-':
            left_chrom_pos, right_chrom_pos = right_chrom_pos, left_chrom_pos
        assert right_chrom_pos >= left_chrom_pos
        return left_chrom_pos, right_chrom_pos

    def parse_indel(left_pos, right_pos, coordinate_fn, tx, offset, gap_type):
        left_chrom_pos, right_chrom_pos = convert_coordinates_to_chromosome(left_pos, right_pos, coordinate_fn,
                                                                            tx.strand)
        if left_chrom_pos is None or right_chrom_pos is None:
            assert aln_mode == 'CDS'
            return None

        if left_chrom_pos > tx.thick_start and right_chrom_pos < tx.thick_stop:
            indel_type = 'CodingMult3' if offset % 3 == 0 else 'Coding'
        else:
            indel_type = 'NonCoding'

        new_bed = tx.get_bed(new_start=left_chrom_pos, new_stop=right_chrom_pos, rgb=offset,
                             name=''.join([indel_type, gap_type]))
        return [tx.name] + new_bed

    if aln_mode == 'CDS':
        coordinate_fn = tx.cds_coordinate_to_chromosome
    else:
        coordinate_fn = tx.mrna_coordinate_to_chromosome
    r = []
    q_pos = 0
    t_pos = 0
    for block_size, q_start, t_start in zip(*[psl.block_sizes, psl.q_starts[1:], psl.t_starts[1:]]):
        q_offset = q_start - block_size - q_pos
        t_offset = t_start - block_size - t_pos
        assert (q_offset >= 0 and t_offset >= 0)
        if q_offset != 0:
            left_pos = q_start - q_offset
            right_pos = q_start
            row = parse_indel(left_pos, right_pos, coordinate_fn, tx, q_offset, 'Insertion')
            if row is not None:
                r.append(row)
        if t_offset != 0:
            left_pos = right_pos = q_start
            row = parse_indel(left_pos, right_pos, coordinate_fn, tx, t_offset, 'Deletion')
            if row is not None:
                r.append(row)
        q_pos = q_start
        t_pos = t_start
    return r


def convert_cds_frame(tx):
    offset = tx.offset
    mod3 = (tx.cds_size - offset) % 3
    if tx.strand == '+':
        b = tx.get_bed(new_start=tx.thick_start + offset, new_stop=tx.thick_stop - mod3)
    else:
        b = tx.get_bed(new_start=tx.thick_start + mod3, new_stop=tx.thick_stop - offset)
    return tools.transcripts.Transcript(b)


def get_intron_coordinates(tx, aln_mode):
    if aln_mode == 'CDS':
        tx = convert_cds_frame(tx)
        introns = [tx.chromosome_coordinate_to_cds(tx.start + x) for x in tx.block_starts[1:]]
    else:
        introns = [tx.chromosome_coordinate_to_mrna(tx.start + x) for x in tx.block_starts[1:]]
    return [x for x in introns if x is not None]


if __name__ == "__main__":
    main()
