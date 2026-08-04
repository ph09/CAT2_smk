#!/usr/bin/env python3
"""Cluster-array-job-based transcript alignment
"""

import argparse
import collections
import gc
import logging
import os
import pickle
import sys
import tempfile
import time
from pathlib import Path

import tools.bio
import tools.dataOps
import tools.fileOps
import tools.nameConversions
import tools.parasail_wrapper
import tools.procOps
import tools.sqlInterface
import tools.transcripts

from cat.scheduler import get_scheduler


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────────────────
# Chunk + worker
# ──────────────────────────────────────────────────────────────────────────────


# Pairs whose longer sequence exceeds this get their own chunk so parasail's
# peak RSS (~30G+ for ~110 kb mRNAs) does not share a cgroup with 499 others.
_HEAVY_SEQ_BP = 100_000


def _mem_gb_for_chunks(chunks, base_mem_gb):
    """Return configured memory; never exceed the YAML / CLI ceiling.

    Heavy sequences are already isolated into one-pair chunks by
    ``chunk_transcripts``. We used to raise requests to 128G/256G here, which
    ignored ``slurm.rules.align_transcripts.mem`` and could leave jobs queued
    forever on smaller clusters. Log a warning instead so operators can raise
    the YAML limit if long transcripts actually OOM.
    """
    base = max(1, int(base_mem_gb))
    max_len = max(max(len(x[1]), len(x[3])) for chunk in chunks for x in chunk)
    if max_len >= 200_000 and base < 256:
        logger.warning(
            f"Longest aligned sequence is {max_len} bp but memory is capped at "
            f"{base}G by config; raise slurm.rules.align_transcripts.mem if jobs OOM"
        )
    elif max_len >= _HEAVY_SEQ_BP and base < 128:
        logger.warning(
            f"Longest aligned sequence is {max_len} bp but memory is capped at "
            f"{base}G by config; raise slurm.rules.align_transcripts.mem if jobs OOM"
        )
    return base


def chunk_transcripts(seq_list, chunk_size=500):
    """Split *seq_list* into chunks for SLURM array tasks.

    Very long sequence pairs are placed one per chunk; remaining pairs are
    batched at *chunk_size*. This keeps peak memory predictable: parasail
    global alignment of ~100 kb mRNAs can exceed 32G, while normal chunks
    should stay well below that.
    """
    heavy = []
    normal = []
    for item in seq_list:
        tx_len = len(item[1])
        ref_len = len(item[3])
        if max(tx_len, ref_len) >= _HEAVY_SEQ_BP:
            heavy.append([item])
        else:
            normal.append(item)

    chunks = list(heavy)
    for i in range(0, len(normal), chunk_size):
        chunks.append(normal[i:i + chunk_size])
    return chunks


def save_chunk(chunk, chunk_dir, chunk_id):
    chunk_file = os.path.join(chunk_dir, f"chunk_{chunk_id}.pkl")
    with open(chunk_file, "wb") as f:
        pickle.dump(chunk, f)
    return chunk_file


def load_chunk(chunk_file):
    with open(chunk_file, "rb") as f:
        return pickle.load(f)


def process_chunk_worker(chunk_file, output_file):
    """Process one chunk of (tx_id, tx_seq, ref_tx_id, ref_tx_seq) tuples."""
    chunk = load_chunk(chunk_file)
    n = 0
    with open(output_file, "w") as f:
        for i, (tx_id, tx_seq, ref_tx_id, ref_tx_seq) in enumerate(chunk):
            p = tools.parasail_wrapper.aln_nucleotides(
                tx_seq, tx_id, ref_tx_seq, ref_tx_id
            )
            psl_str = '\t'.join(p.psl_string())
            if psl_str:
                f.write(psl_str + '\n')
                n += 1
            # Parasail can retain large trace buffers; release between alignments.
            del p
            if (i + 1) % 25 == 0:
                gc.collect()
    return n


def get_alignment_sequences(transcript_dict, ref_transcript_dict, genome_fasta,
                            ref_genome_fasta, mode, max_ref_span_ratio=5.0):
    """Build (tx_id, tx_seq, ref_tx_id, ref_tx_seq) tuples for the given mode.

    Drops targets whose genomic span exceeds ``max_ref_span_ratio`` × the
    reference transcript span (same idea as ``filter_transmap.ref_span`` /
    ``tm_max_ref_span``). Those are almost always chimeric / mis-spliced
    junk and they dominate parasail memory.
    """
    assert mode in ['mRNA', 'CDS']
    sequences = []
    dropped_span = 0
    for tx_id, tx in transcript_dict.items():
        ref_tx_id = tools.nameConversions.alignment_id_to_ref_transcript_id(tx_id)
        ref_tx = ref_transcript_dict.get(ref_tx_id)
        if ref_tx is None:
            continue
        ref_span = int(ref_tx.stop) - int(ref_tx.start)
        tx_span = int(tx.stop) - int(tx.start)
        if ref_span > 0 and max_ref_span_ratio > 0 and tx_span > ref_span * float(max_ref_span_ratio):
            dropped_span += 1
            continue
        try:
            tx_seq = tx.get_mrna(genome_fasta) if mode == 'mRNA' else tx.get_cds(genome_fasta)
        except KeyError:
            logger.warning(f"Skipping {tx_id}: chromosome {tx.chromosome} missing in target FASTA")
            continue
        try:
            ref_tx_seq = ref_tx.get_mrna(ref_genome_fasta) if mode == 'mRNA' else ref_tx.get_cds(ref_genome_fasta)
        except KeyError:
            logger.warning(f"Skipping {tx_id}: chromosome {ref_tx.chromosome} missing in reference FASTA")
            continue
        # Parasail behaves poorly on very short sequences.
        if len(ref_tx_seq) > 20 and len(tx_seq) > 20:
            sequences.append((tx_id, tx_seq, ref_tx_id, ref_tx_seq))
    if dropped_span:
        logger.info(
            f"Dropped {dropped_span} {mode} transcripts with genomic span > "
            f"{max_ref_span_ratio:g}× reference"
        )
    return sequences


# ──────────────────────────────────────────────────────────────────────────────
# Cluster submission
# ──────────────────────────────────────────────────────────────────────────────


def _build_body(scheduler, chunk_dir, output_dir, sentinel_dir):
    """Bash body executed by every array task.

    Uses 1-based TASK_ID with a shift to 0-based IDX so chunk filenames
    remain ``chunk_0.pkl`` … ``chunk_{N-1}.pkl``.
    """
    task_var = scheduler.task_id_env()
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    preamble = scheduler.script_preamble(conda_env=os.environ.get("CONDA_DEFAULT_ENV", "cat") or "cat")
    sentinel = scheduler.trap_sentinel(sentinel_dir)
    return f"""{preamble}
{sentinel}

TASK_ID="${{{task_var}:-1}}"
IDX=$((TASK_ID - 1))

CHUNK_FILE="{chunk_dir}/chunk_${{IDX}}.pkl"
OUTPUT_FILE="{output_dir}/results/result_${{IDX}}.psl"

echo "Processing chunk $IDX at $(date)"
echo "Chunk file: $CHUNK_FILE"
echo "Output file: $OUTPUT_FILE"

python3 - <<PYEOF
import sys
sys.path.insert(0, {repo_root!r})
from cat.align_transcripts_cluster import process_chunk_worker
process_chunk_worker("$CHUNK_FILE", "$OUTPUT_FILE")
PYEOF

echo "Completed chunk $IDX at $(date)"
"""


def merge_results(result_dir, output_psl, expected_chunks=None):
    """Concatenate result_*.psl in result_dir into output_psl."""
    logger.info(f"Merging results from {result_dir} to {output_psl}")
    result_files = sorted(Path(result_dir).glob("result_*.psl"))

    if not result_files:
        logger.error(f"No result files found in {result_dir}")
        if expected_chunks:
            logger.error(f"Expected {expected_chunks} result files but found 0 - all alignment tasks failed!")
        Path(output_psl).touch()
        return 0

    total_lines = 0
    with open(output_psl, "w") as outf:
        for result_file in result_files:
            with open(result_file, "r") as inf:
                for line in inf:
                    if line.strip():
                        outf.write(line)
                        total_lines += 1

    logger.info(f"Merged {len(result_files)} result files into {output_psl}")
    logger.info(f"Total alignments written: {total_lines}")

    if expected_chunks and len(result_files) < expected_chunks:
        missing = expected_chunks - len(result_files)
        logger.error(f"WARNING: Expected {expected_chunks} result files but only found {len(result_files)}")
        logger.error(f"Missing {missing} result files ({missing/expected_chunks*100:.1f}% of chunks failed)")

    return total_lines


def robust_rmtree(path: str, max_retries: int = 6, initial_delay_seconds: float = 0.5) -> None:
    """Remove a directory tree, retrying on transient NFS-style errors.

    Logs but never raises; cleanup failures must not fail the overall workflow.
    """
    import shutil

    delay = initial_delay_seconds
    for attempt in range(1, max_retries + 1):
        try:
            shutil.rmtree(path)
            return
        except OSError as e:
            if attempt == max_retries:
                logger.warning(
                    f"Cleanup: failed to remove '{path}' after {max_retries} attempts: {e}. "
                    "Proceeding without removing temp directory."
                )
                return
            logger.warning(
                f"Cleanup: attempt {attempt}/{max_retries} to remove '{path}' failed: {e}. "
                f"Retrying in {delay:.1f}s"
            )
            time.sleep(delay)
            delay = min(delay * 2, 8.0)


# ──────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ──────────────────────────────────────────────────────────────────────────────


def run_cluster_alignment_pipeline(args):
    """Cluster-based alignment pipeline driver."""
    scheduler = get_scheduler(args.execution_mode, _scheduler_config(args))

    logger.info(f"Starting {scheduler.name}-based transcript alignment pipeline")
    logger.info(f"Target genome: {args.genome}")

    transcript_modes = collections.defaultdict(dict)
    for mode, gp_path, mrna_path, cds_path in args.mode_files:
        transcript_modes[mode]['gp'] = gp_path
        transcript_modes[mode]['mRNA'] = mrna_path
        transcript_modes[mode]['CDS'] = cds_path

    logger.info(f"Processing {len(transcript_modes)} alignment modes: {list(transcript_modes.keys())}")

    logger.info("Loading reference data...")
    tx_biotype_map = tools.sqlInterface.get_transcript_biotype_map(args.ref_db_path)
    ref_transcript_dict = tools.transcripts.get_gene_pred_dict(args.annotation_gp)

    logger.info("Loading FASTA files...")
    genome_fasta = tools.bio.get_sequence_dict(args.genome_fasta, upper=False)
    ref_genome_fasta = tools.bio.get_sequence_dict(args.ref_genome_fasta, upper=False)

    # Biotypes treated as protein-coding for alignment purposes.
    protein_coding_biotypes = {
        'protein_coding',
        'nonsense_mediated_decay',
        'non_stop_decay',
        'protein_coding_LoF',
        'protein_coding_CDS_not_defined',
    }

    for tx_mode in ['transMap', 'transMap_pairwise', 'augTM', 'augTM_pairwise',
                    'augTMR', 'augTMR_pairwise', 'augMP', 'txTM']:
        if tx_mode not in transcript_modes:
            continue

        logger.info(f"\nProcessing alignment mode: {tx_mode}")

        gp_path = transcript_modes[tx_mode]['gp']
        mrna_path = transcript_modes[tx_mode]['mRNA']
        cds_path = transcript_modes[tx_mode]['CDS']

        logger.info(f"Loading transcripts from {gp_path}")
        transcript_dict = tools.transcripts.get_gene_pred_dict(gp_path)

        filtered_transcripts = {}
        for aln_id, tx in transcript_dict.items():
            ref_id = tools.nameConversions.alignment_id_to_ref_transcript_id(aln_id)
            biotype = tx_biotype_map.get(ref_id)
            if biotype in protein_coding_biotypes or biotype is None:
                filtered_transcripts[aln_id] = tx

        if not filtered_transcripts and transcript_dict:
            filtered_transcripts = transcript_dict

        transcript_dict = filtered_transcripts
        logger.info(f"Processing {len(transcript_dict)} transcripts")

        for aln_mode, out_path in [('mRNA', mrna_path), ('CDS', cds_path)]:
            logger.info(f"\nProcessing {aln_mode} alignments for {tx_mode}")

            sequences = get_alignment_sequences(
                transcript_dict, ref_transcript_dict,
                genome_fasta, ref_genome_fasta, aln_mode,
                max_ref_span_ratio=args.max_ref_span,
            )

            if not sequences:
                logger.warning(f"No sequences to align for {tx_mode} {aln_mode}")
                Path(out_path).touch()
                continue

            logger.info(f"Extracted {len(sequences)} sequence pairs for alignment")

            # Working dirs sit beside the final output (on shared storage so
            # compute nodes can read them, regardless of backend).
            tmp_root = Path(out_path).parent / "_cluster_work"
            tmp_root.mkdir(parents=True, exist_ok=True)
            work_dir = Path(tempfile.mkdtemp(
                prefix=f"align_{args.genome}_{tx_mode}_{aln_mode}_",
                dir=str(tmp_root),
            ))
            chunk_dir = work_dir / "chunks"
            output_dir = work_dir / "output"
            result_dir = output_dir / "results"
            log_dir = output_dir / "cluster_logs"
            sentinel_dir = work_dir / "sentinels"

            for d in (chunk_dir, result_dir, log_dir, sentinel_dir):
                d.mkdir(parents=True, exist_ok=True)

            logger.info(f"Working directory: {work_dir}")

            chunks = chunk_transcripts(sequences, chunk_size=args.chunk_size)
            logger.info(f"Split into {len(chunks)} chunks (chunk size: {args.chunk_size})")
            for i, chunk in enumerate(chunks):
                save_chunk(chunk, str(chunk_dir), i)

            num_chunks = len(chunks)
            mem_gb = _mem_gb_for_chunks(chunks, args.memory)
            header = scheduler.header(
                job_name=f"align_{tx_mode}_{aln_mode}",
                cpus=args.cpus,
                mem=f"{mem_gb}G",
                walltime=args.time,
                log_out=str(log_dir / "align_%A_%a.out"),
                log_err=str(log_dir / "align_%A_%a.err"),
                partition=args.partition,
                queue=args.partition,
                array=(1, num_chunks),
                max_concurrent=args.max_jobs,
            )
            body = _build_body(scheduler, chunk_dir, output_dir, sentinel_dir)
            script_path = scheduler.write_script(header + body, work_dir / "run_alignment.sh")

            job_id = scheduler.submit(script_path)
            logger.info(f"Submitted {scheduler.name} array job {job_id}")

            result = scheduler.wait(
                job_id,
                num_tasks=num_chunks,
                timeout_s=args.timeout_hours * 3600,
                sentinel_dir=sentinel_dir,
            )
            if not result.ok:
                logger.error(f"{scheduler.name} job {job_id} failed: {result.detail}")
                sys.exit(1)
            logger.info(f"Job {job_id} succeeded: {result.completed}/{result.total} tasks")

            num_lines = merge_results(str(result_dir), out_path, expected_chunks=num_chunks)

            if num_lines == 0:
                logger.error(f"ERROR: No alignments produced for {tx_mode} {aln_mode}")
                logger.error(f"Expected ~{len(sequences)} alignments but got 0")
            elif num_lines < len(sequences) * 0.5:
                logger.warning(f"WARNING: Only {num_lines}/{len(sequences)} alignments produced ({num_lines/len(sequences)*100:.1f}%)")

            if args.cleanup:
                logger.info(f"Cleaning up temporary directory: {work_dir}")
                robust_rmtree(str(work_dir))

    logger.info(f"\n{scheduler.name}-based transcript alignment pipeline completed successfully")


def _scheduler_config(args):
    """Build a per-backend config dict from CLI args."""
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


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--ref-genome-fasta", required=True)
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--annotation-gp", required=True)
    parser.add_argument("--ref-db-path", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--mode-files", nargs=4, action='append', required=True,
                        metavar=("MODE", "INPUT_GP", "MRNA_PSL", "CDS_PSL"))

    parser.add_argument("--execution-mode", choices=("auto", "slurm", "sge", "local"), default="auto")
    parser.add_argument("--partition", default="",
                        help="SLURM partition or SGE queue.")
    parser.add_argument("--exclude-nodes", default="",
                        help=("Comma list (SLURM) or SGE-native '!h1&!h2' expression. "
                              "Single source of truth for node exclusion; no env-var override."))
    parser.add_argument("--module-load", default="")
    parser.add_argument("--sge-parallel-env", default="smp")
    parser.add_argument("--sge-memory-flag", default="h_vmem")
    parser.add_argument("--memory", type=int, default=8)
    parser.add_argument("--cpus", type=int, default=1)
    parser.add_argument("--time", default="02:00:00")
    parser.add_argument("--max-jobs", type=int, default=50)
    parser.add_argument("--timeout-hours", type=int, default=12)
    parser.add_argument("--chunk-size", type=int, default=500)
    parser.add_argument(
        "--max-ref-span", type=float, default=5.0,
        help=("Drop target transcripts whose genomic span exceeds this multiple "
              "of the reference transcript span (default: 5; same as tm_max_ref_span). "
              "Set <=0 to disable."),
    )
    parser.add_argument("--cleanup", action="store_true")

    args = parser.parse_args()
    run_cluster_alignment_pipeline(args)


if __name__ == "__main__":
    main()
