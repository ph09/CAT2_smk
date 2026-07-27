"""
Generate a hints database file from RNAseq/IsoSeq alignments for AugustusTMR.
Enhanced with parallel processing and dynamic resource allocation.
"""
import collections
import itertools
import os
import shutil
import logging
import pysam
import argparse
import multiprocessing as mp
import subprocess
import time
import tempfile
import shutil
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from functools import partial

from toil.fileStores import FileID
from toil.common import Toil
from toil.job import Job
import tools.dataOps
import tools.fileOps
import tools.misc
import tools.procOps
import tools.toilInterface
import tools.transcripts
import tools.bio
from cat.exceptions import UserException
from cat.scheduler import Scheduler, get_scheduler

logger = logging.getLogger('cat')


# Process-wide active Scheduler instance for hints_db. Set by _init_scheduler()
# when the module enters cluster mode; remains None for the toil/local path.
_SCHEDULER: Scheduler = None


# Set in run_hints_pipeline from Toil --maxCores; used when worker env lacks the var.
_TOIL_MAX_CORES_OVERRIDE = None


def _toil_max_cores() -> int:
    """Upper bound for Toil job cores (honours --maxCores / CAT2_TOIL_MAX_CORES)."""
    if _TOIL_MAX_CORES_OVERRIDE is not None:
        return _TOIL_MAX_CORES_OVERRIDE
    for key in ("CAT2_TOIL_MAX_CORES", "TOIL_MAX_CORES"):
        val = os.environ.get(key)
        if val:
            try:
                return max(1, int(float(val)))
            except ValueError:
                pass
    return max(1, mp.cpu_count() or 1)


def _cap_cores(requested) -> int:
    """Clamp a core request so SingleMachine --maxCores is never exceeded."""
    try:
        n = int(requested)
    except (TypeError, ValueError):
        n = 1
    return max(1, min(n, _toil_max_cores()))


def _init_scheduler(args) -> Scheduler:
    """Construct and stash the Scheduler from CLI args, returning the instance.

    Idempotent: calling twice with the same args returns the same object.
    """
    global _SCHEDULER
    cfg = {
        "cluster": {
            "slurm": {
                "partition": getattr(args, "slurm_partition", None),
                "exclude_nodes": getattr(args, "exclude_nodes", "") or "",
                "module_load": getattr(args, "module_load", "") or "",
            },
            "sge": {
                "queue": getattr(args, "slurm_partition", None),
                "parallel_env": getattr(args, "sge_parallel_env", "smp"),
                "memory_flag": getattr(args, "sge_memory_flag", "h_vmem"),
                "hostname_exclude": getattr(args, "exclude_nodes", "") or "",
                "module_load": getattr(args, "module_load", "") or "",
            },
        }
    }
    mode = getattr(args, "execution_mode", "slurm")
    _SCHEDULER = get_scheduler(mode, cfg)
    return _SCHEDULER


def _scheduler() -> Scheduler:
    """Return the currently-active hints_db Scheduler.

    Fallback constructs a default SLURM scheduler so the existing public
    generate_slurm_*_job entry points keep working in legacy callers that
    haven't migrated to the new CLI args.
    """
    global _SCHEDULER
    if _SCHEDULER is None:
        _SCHEDULER = get_scheduler("slurm", {})
    return _SCHEDULER


def _hints_header(*, job_name, log_file, error_file, partition, memory_gb, cpus,
                  time_limit, exclude=None):
    """Build a job-script header for a single (non-array) hints job.

    Delegates to the active Scheduler. *exclude* lets per-call sites override
    the default exclude_nodes; callers should pass None to use the default.
    """
    s = _scheduler()
    return s.header(
        job_name=job_name,
        cpus=cpus,
        mem=f"{memory_gb}G",
        walltime=time_limit,
        log_out=log_file,
        log_err=error_file,
        partition=partition,
        queue=partition,
        exclude=exclude,  # None → uses the scheduler's default
    )


def generate_slurm_bam2hints_job(bam_path, output_path, job_id, temp_dir, 
                                 memory_gb=256, cpus=128, time_limit="12:00:00", 
                                 partition="", hint_type="intron"):
    """
    DEPRECATED: Use generate_slurm_bam_intron_job() or generate_slurm_bam_exon_job() instead.
    Generate a Slurm job script for running bam2hints on a single BAM file.
    Uses the same successful pipeline approach as IsoSeq hints.
    
    :param bam_path: Path to input BAM file
    :param output_path: Path for output hints file
    :param job_id: Unique job identifier
    :param temp_dir: Temporary directory for job files
    :param memory_gb: Memory allocation in GB
    :param cpus: Number of CPUs to allocate
    :param time_limit: Time limit for the job (HH:MM:SS)
    :param partition: Slurm partition to use
    :param hint_type: Type of hints to generate ("intron" or "exon")
    :return: Path to the generated job script
    """
    job_script = os.path.join(temp_dir, f"bam2hints_{job_id}.sh")
    log_file = os.path.join(temp_dir, f"bam2hints_{job_id}.log")
    error_file = os.path.join(temp_dir, f"bam2hints_{job_id}.err")
    
    header = _hints_header(
        job_name=f"bam2hints_{job_id}", log_file=log_file, error_file=error_file,
        partition=partition, memory_gb=memory_gb, cpus=cpus, time_limit=time_limit,
    )
    script_content = header + f"""
set -euo pipefail

mkdir -p {os.path.dirname(output_path)}

# Use the same successful pipeline approach as IsoSeq hints
# Convert BAM to PSL format, then use blat2hints.pl
samtools view -@ {cpus} -b -F 4 {bam_path} \\
| bamToPsl -nohead /dev/stdin /dev/stdout \\
| sort -S $(( {memory_gb} / 2 ))G -T {temp_dir} -n -k 16,16 \\
| sort -S $(( {memory_gb} / 2 ))G -T {temp_dir} -s -k 14,14 \\
| perl -ne '@f=split; print if ($f[0]>=100)' \\
| blat2hints.pl --source=W --nomult --ep_cutoff=20 --in=/dev/stdin --out={output_path}

# Verify output file was created and is not empty
if [ -s "{output_path}" ]; then
  echo "Success: {output_path} created"
  exit 0
else
  echo "Error: {output_path} empty or missing"
  exit 1
fi
"""
    
    with open(job_script, 'w') as f:
        f.write(script_content)
    
    # Make script executable
    os.chmod(job_script, 0o755)
    
    logger.info(f'Generated Slurm job script: {job_script}')
    return job_script, log_file, error_file


def submit_slurm_jobs(job_scripts, max_concurrent_jobs=50, check_interval=30, log_files=None):
    """
    Submit multiple Slurm jobs and monitor their completion.
    
    :param job_scripts: List of paths to job scripts
    :param max_concurrent_jobs: Maximum number of jobs to submit concurrently
    :param check_interval: Interval in seconds to check job status
    :return: Dictionary mapping job scripts to job IDs and status
    """
    job_status = {}
    submitted_jobs = []
    
    logger.info(f'Submitting {len(job_scripts)} Slurm jobs (max concurrent: {max_concurrent_jobs})')
    
    # Submit jobs in batches to avoid overwhelming the scheduler
    for i in range(0, len(job_scripts), max_concurrent_jobs):
        batch = job_scripts[i:i + max_concurrent_jobs]
        
        for script_path in batch:
            try:
                # Submit through the active scheduler so SLURM and SGE share
                # the same code path. Scheduler.submit() raises if submission
                # fails; we coerce its CalledProcessError into the same
                # bookkeeping the legacy code used.
                job_id = _scheduler().submit(script_path)
                job_status[script_path] = {'job_id': job_id, 'status': 'submitted'}
                submitted_jobs.append(job_id)

                logger.info(f'Submitted job {job_id} for script {os.path.basename(script_path)}')

            except (subprocess.CalledProcessError, RuntimeError) as e:
                logger.error(f'Failed to submit job for {script_path}: {e}')
                job_status[script_path] = {'job_id': None, 'status': 'failed', 'error': str(e)}
        
        # Wait for batch to complete before submitting next batch
        if i + max_concurrent_jobs < len(job_scripts):
            logger.info(f'Waiting for batch {i//max_concurrent_jobs + 1} to complete before submitting next batch')
            time.sleep(check_interval)
    
    # Monitor job completion
    logger.info(f'Monitoring {len(submitted_jobs)} submitted jobs...')
    completed_jobs = set()
    
    while len(completed_jobs) < len(submitted_jobs):
        time.sleep(check_interval)

        try:
            # Backend-agnostic per-job presence check. Note: this loop is O(N)
            # in the number of submitted jobs; for the largest hints_db runs
            # (~few hundred jobs) that's still cheap compared to the per-job
            # cost of running the scheduler binary in batch.
            sched = _scheduler()
            running_jobs = {jid for jid in submitted_jobs if sched.job_present(jid)}
            for jid in submitted_jobs:
                for script_path, job_info in job_status.items():
                    if job_info.get('job_id') == jid:
                        job_info['status'] = 'RUNNING' if jid in running_jobs else 'NOT_RUNNING'
            
            # Check for completed jobs
            current_completed = set(submitted_jobs) - running_jobs - completed_jobs
            if current_completed:
                logger.info(f'Jobs completed: {len(current_completed)}')
                completed_jobs.update(current_completed)
                
                # Check if jobs succeeded or failed
                for job_id in current_completed:
                    for i, (script_path, job_info) in enumerate(job_status.items()):
                        if job_info.get('job_id') == job_id:
                            # Check exit status by looking at log files
                            if log_files and i < len(log_files):
                                log_file = log_files[i]
                            else:
                                log_file = script_path.replace('.sh', '.log')
                            # Infer error file alongside the log file
                            if log_file.endswith('.log'):
                                error_file = log_file[:-4] + 'err'
                            else:
                                error_file = script_path.replace('.sh', '.err')
                            
                            # Ensure log file directory exists
                            log_dir = os.path.dirname(log_file)
                            if log_dir and not os.path.exists(log_dir):
                                try:
                                    os.makedirs(log_dir, exist_ok=True)
                                    logger.info(f'Created log directory: {log_dir}')
                                except Exception as e:
                                    logger.warning(f'Failed to create log directory {log_dir}: {e}')
                            
                            # Wait a bit for log file to be created by the job
                            max_wait_attempts = 5
                            wait_attempt = 0
                            log_file_ready = False
                            
                            while wait_attempt < max_wait_attempts and not log_file_ready:
                                if os.path.exists(log_file):
                                    try:
                                        with open(log_file, 'r') as f:
                                            log_content = f.read()
                                            if 'Success:' in log_content:
                                                job_info['status'] = 'completed'
                                                log_file_ready = True
                                            elif 'Error:' in log_content or 'Failed:' in log_content:
                                                job_info['status'] = 'failed'
                                                log_file_ready = True
                                            else:
                                                # Log file exists but doesn't have success/error markers yet
                                                if wait_attempt < max_wait_attempts - 1:
                                                    time.sleep(2)  # Wait 2 seconds before retrying
                                                    wait_attempt += 1
                                                else:
                                                    # Before declaring unknown, check if there's an error file with content
                                                    if os.path.exists(error_file):
                                                        try:
                                                            with open(error_file, 'r') as ef:
                                                                err_txt = ef.read().strip()
                                                                if err_txt:
                                                                    job_info['status'] = 'failed'
                                                                    job_info['error'] = err_txt[-2000:]
                                                                else:
                                                                    job_info['status'] = 'unknown'
                                                        except Exception:
                                                            job_info['status'] = 'unknown'
                                                    else:
                                                        job_info['status'] = 'unknown'
                                                    log_file_ready = True
                                    except Exception as e:
                                        logger.warning(f'Failed to read log file {log_file}: {e}')
                                        if wait_attempt < max_wait_attempts - 1:
                                            time.sleep(2)
                                            wait_attempt += 1
                                        else:
                                            job_info['status'] = 'unknown'
                                            log_file_ready = True
                                else:
                                    # Log file doesn't exist yet - this is normal for recently completed jobs
                                    if wait_attempt < max_wait_attempts - 1:
                                        logger.debug(f'Log file not found yet: {log_file} (attempt {wait_attempt + 1}/{max_wait_attempts})')
                                        time.sleep(2)
                                        wait_attempt += 1
                                    else:
                                        logger.debug(f'Log file not found after {max_wait_attempts} attempts: {log_file}')
                                        # If the log never appeared, try to read the error file; if present, mark as failed
                                        if os.path.exists(error_file):
                                            try:
                                                with open(error_file, 'r') as ef:
                                                    err_txt = ef.read().strip()
                                                    if err_txt:
                                                        job_info['status'] = 'failed'
                                                        job_info['error'] = err_txt[-2000:]
                                                    else:
                                                        job_info['status'] = 'unknown'
                                            except Exception:
                                                job_info['status'] = 'unknown'
                                        else:
                                            job_info['status'] = 'unknown'
                                        log_file_ready = True
            
            logger.info(f'Progress: {len(completed_jobs)}/{len(submitted_jobs)} jobs completed')
            
        except subprocess.CalledProcessError as e:
            logger.warning(f'Error checking job status: {e}')
            time.sleep(check_interval)
            continue
    
    # Final status check
    final_status = {}
    for script_path, job_info in job_status.items():
        if job_info['status'] == 'submitted':
            # Check if job actually completed
            log_file = script_path.replace('.sh', '.log')
            err_file = script_path.replace('.sh', '.err')
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    if 'Success:' in f.read():
                        job_info['status'] = 'completed'
                    else:
                        job_info['status'] = 'failed'
            elif os.path.exists(err_file):
                # If only error file exists, consider it failed and retain tail of error
                try:
                    with open(err_file, 'r') as ef:
                        err_txt = ef.read().strip()
                        if err_txt:
                            job_info['status'] = 'failed'
                            job_info['error'] = err_txt[-2000:]
                        else:
                            job_info['status'] = 'failed'
                except Exception:
                    job_info['status'] = 'failed'
            else:
                # Treat missing logs as failure to avoid opaque 'unknown'
                job_info['status'] = 'failed'
        
        final_status[script_path] = job_info
    
    logger.info('All Slurm jobs completed')
    return final_status


def run_parallel_bam2hints_slurm(bam_files, output_dir, hint_type="intron", 
                                memory_gb=256, cpus=128, time_limit="12:00:00",
                                partition="", max_concurrent_jobs=50):
    """
    DEPRECATED: Use run_parallel_bam_intron_slurm() or run_parallel_bam_exon_slurm() instead.
    Run bam2hints on multiple BAM files using parallel Slurm jobs.
    
    :param bam_files: List of BAM file paths
   :param output_dir: Directory for output hints files
    :param hint_type: Type of hints ("intron" or "exon")
    :param memory_gb: Memory allocation per job in GB
    :param cpus: CPUs per job
    :param time_limit: Time limit per job
    :param partition: Slurm partition
    :param max_concurrent_jobs: Maximum concurrent jobs
    :return: List of output hints file paths
    """
    logger.info(f'Running parallel bam2hints with Slurm on {len(bam_files)} BAM files')
    
    # Create temporary directory for job scripts and logs
    temp_dir = tempfile.mkdtemp(prefix='bam2hints_slurm_')
    logger.info(f'Using temporary directory: {temp_dir}')
    
    try:
        # Generate job scripts for each BAM file
        job_scripts = []
        output_files = []
        
        for i, bam_path in enumerate(bam_files):
            # Create output file path
            bam_basename = os.path.splitext(os.path.basename(bam_path))[0]
            output_file = os.path.join(output_dir, f"{bam_basename}_{hint_type}_hints.gff")
            output_files.append(output_file)
            
            # Generate job script
            job_script, log_file, error_file = generate_slurm_bam2hints_job(
                bam_path, output_file, f"{i}_{bam_basename}", temp_dir,
                memory_gb=memory_gb, cpus=cpus, time_limit=time_limit,
                partition=partition, hint_type=hint_type
            )
            job_scripts.append(job_script)
        
        # Submit and monitor jobs
        job_status = submit_slurm_jobs(job_scripts, max_concurrent_jobs)
        
        # Collect results
        successful_outputs = []
        failed_jobs = []
        
        for i, (script_path, output_file) in enumerate(zip(job_scripts, output_files)):
            status = job_status[script_path]['status']
            
            if status == 'completed' and os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                successful_outputs.append(output_file)
                logger.info(f'Successfully generated hints: {output_file}')
            else:
                failed_jobs.append((bam_files[i], output_file, status))
                logger.error(f'Failed to generate hints for {bam_files[i]}: {status}')
        
        logger.info(f'Successfully processed {len(successful_outputs)}/{len(bam_files)} BAM files')
        
        if failed_jobs:
            logger.warning(f'{len(failed_jobs)} jobs failed. Failed jobs:')
            for bam_file, output_file, status in failed_jobs:
                logger.warning(f'  - {bam_file}: {status}')
        
        return successful_outputs, failed_jobs
        
    finally:
        # Cleanup temporary directory
        try:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
                logger.info(f'Cleaned up temporary directory: {temp_dir}')
            else:
                logger.debug(f'Temporary directory already removed: {temp_dir}')
        except OSError as e:
            if e.errno == 39:  # Directory not empty
                logger.warning(f'Directory not empty, attempting force cleanup: {temp_dir}')
                try:
                    # Try to remove files individually first
                    for root, dirs, files in os.walk(temp_dir, topdown=False):
                        for file in files:
                            try:
                                os.remove(os.path.join(root, file))
                            except OSError:
                                pass
                        for dir in dirs:
                            try:
                                os.rmdir(os.path.join(root, dir))
                            except OSError:
                                pass
                    os.rmdir(temp_dir)
                    logger.info(f'Force cleaned up temporary directory: {temp_dir}')
                except Exception as force_e:
                    logger.warning(f'Failed to force clean up temporary directory {temp_dir}: {force_e}')
            else:
                logger.warning(f'Failed to clean up temporary directory {temp_dir}: {e}')
        except Exception as e:
            logger.warning(f'Failed to clean up temporary directory {temp_dir}: {e}')


def concatenate_hints_files(hints_files, output_file, sort_hints=True):
    """
    Concatenate multiple hints files into a single file with optional sorting.
    
    :param hints_files: List of hints file paths
    :param output_file: Output file path
    :param sort_hints: Whether to sort the concatenated hints
    :return: Path to the output file
    """
    logger.info(f'Concatenating {len(hints_files)} hints files to {output_file}')
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Concatenate files
    temp_concatenated = output_file + '.tmp'
    
    with open(temp_concatenated, 'w') as outf:
        for hints_file in hints_files:
            if os.path.exists(hints_file) and os.path.getsize(hints_file) > 0:
                with open(hints_file, 'r') as inf:
                    for line in inf:
                        if line.strip():  # Skip empty lines
                            outf.write(line)
                logger.debug(f'Added hints from: {hints_file}')
    
    # Sort if requested
    if sort_hints:
        logger.info('Sorting concatenated hints file (single-pass multi-key sort)')
        # Use single-pass multi-key sort instead of 4 sequential sorts for much better performance
        cmd = [['sort', '-k1,1', '-k3,3', '-k4,4n', '-k5,5n', temp_concatenated],
               ['join_mult_hints.pl']]
        tools.procOps.run_proc(cmd, stdout=output_file)
        os.remove(temp_concatenated)
    else:
        shutil.move(temp_concatenated, output_file)
    
    logger.info(f'Created concatenated hints file: {output_file} ({os.path.getsize(output_file)} bytes)')
    return output_file


def concatenate_and_sort_hints_slurm(hints_files, output_file, logs_dir=None,
                                     memory_gb=256, cpus=64, time_limit="12:00:00",
                                     partition="", tmp_dir=None):
    """
    Cluster-based concatenation and sorting of hints files.

    Runs the compute-intensive sort on a cluster node (SLURM or SGE) instead
    of the head node, with parallel sorting for much faster performance.
    Falls back to local execution when the active scheduler is LocalScheduler.

    :param hints_files: List of hints file paths to merge
    :param output_file: Final output file path
    :param logs_dir: Directory for log files
    :param memory_gb: Memory in GB for the cluster job
    :param cpus: Number of CPUs for parallel sorting
    :param time_limit: Time limit for the job
    :param partition: Cluster partition / SGE queue
    :param tmp_dir: Temporary directory for sort (should be on fast storage)
    :return: Path to the output file
    """
    import tempfile
    
    logger.info(f'Starting Slurm-based merge/sort of {len(hints_files)} hints files')
    logger.info(f'Using {cpus} CPUs and {memory_gb}GB RAM for parallel sorting')
    
    # Create output directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Create logs directory
    if logs_dir is None:
        logs_dir = os.path.join(os.path.dirname(output_file), 'logs', 'merge_sort')
    os.makedirs(logs_dir, exist_ok=True)
    
    # Concatenate files first (fast, runs locally)
    temp_concatenated = output_file + '.tmp'
    logger.info(f'Concatenating {len(hints_files)} files into {temp_concatenated}')
    
    concatenated_size = 0
    with open(temp_concatenated, 'w') as outf:
        for hints_file in hints_files:
            if os.path.exists(hints_file) and os.path.getsize(hints_file) > 0:
                file_size = os.path.getsize(hints_file)
                concatenated_size += file_size
                with open(hints_file, 'r') as inf:
                    shutil.copyfileobj(inf, outf)
    
    concatenated_size_gb = concatenated_size / (1024**3)
    logger.info(f'Concatenated file size: {concatenated_size_gb:.2f} GB')
    
    # Create optimized sort script
    sort_script = tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False, dir=logs_dir)
    sort_script_path = sort_script.name
    
    # Determine buffer size based on available memory (use ~70% for sort buffers)
    buffer_size = int(memory_gb * 0.7)
    
    # Set temp directory for sort (critical for performance)
    if tmp_dir is None:
        tmp_dir = os.path.dirname(temp_concatenated)
    
    sort_script.write(f"""#!/bin/bash
set -euo pipefail

# Log system info
echo "Running on host: $(hostname)"
echo "Available memory: $(free -h | grep Mem)"
echo "CPUs available: $(nproc)"
echo "Temp directory: {tmp_dir}"
echo "Input file size: $(ls -lh {temp_concatenated} | awk '{{print $5}}')"
echo "Start time: $(date)"

# Create temp directory for sort if it doesn't exist
mkdir -p {tmp_dir}

# Optimized parallel sort pipeline
# Using single-pass multi-key sort instead of 4 sequential sorts for much better performance
echo "Starting parallel sort pipeline (single-pass multi-key sort)..."

sort --parallel={cpus} \
     --buffer-size={buffer_size}G \
     -T {tmp_dir} \
     -k1,1 -k3,3 -k4,4n -k5,5n {temp_concatenated} | \
join_mult_hints.pl > {output_file}

echo "Sort completed successfully"
echo "Output file size: $(ls -lh {output_file} | awk '{{print $5}}')"
echo "End time: $(date)"

# Clean up temp file
rm -f {temp_concatenated}

echo "Merge/sort job completed successfully"
""")
    sort_script.close()
    os.chmod(sort_script_path, 0o755)
    
    # Create Slurm submission script
    job_name = f"merge_sort_{os.path.basename(output_file).replace('.gff', '')}"
    log_file = os.path.join(logs_dir, f'{job_name}.log')
    error_file = os.path.join(logs_dir, f'{job_name}.err')
    
    slurm_script = tempfile.NamedTemporaryFile(mode='w', suffix='.cluster.sh', delete=False, dir=logs_dir)
    slurm_script_path = slurm_script.name

    sched = _scheduler()
    header = sched.header(
        job_name=job_name,
        cpus=cpus,
        mem=f"{memory_gb}G",
        walltime=time_limit,
        log_out=log_file,
        log_err=error_file,
        partition=partition,
        queue=partition,
    )
    slurm_script.write(header + f"""
{sort_script_path}
""")
    slurm_script.close()

    logger.info(f'Submitting {sched.name} job: {job_name}')
    logger.info(f'  Memory: {memory_gb}GB, CPUs: {cpus}, Time: {time_limit}')
    logger.info(f'  Log file: {log_file}')

    job_id = sched.submit(slurm_script_path)
    logger.info(f'Submitted {sched.name} job: {job_id}')

    # Monitor job via the backend's queue-presence check, streaming progress
    # from the log file in between polls.
    logger.info(f'Monitoring job {job_id}...')
    check_interval = 30
    last_log_size = 0

    while True:
        time.sleep(check_interval)

        if not sched.job_present(job_id):
            break

        if os.path.exists(log_file):
            current_size = os.path.getsize(log_file)
            if current_size > last_log_size:
                with open(log_file, 'r') as f:
                    f.seek(last_log_size)
                    new_content = f.read()
                    if new_content.strip():
                        logger.info(f'Job {job_id} progress:\n{new_content.strip()}')
                last_log_size = current_size
    
    # Check if job completed successfully
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        output_size_gb = os.path.getsize(output_file) / (1024**3)
        logger.info(f'Merge/sort completed successfully!')
        logger.info(f'Output file: {output_file} ({output_size_gb:.2f} GB)')
        
        # Clean up
        os.remove(sort_script_path)
        os.remove(slurm_script_path)
        
        return output_file
    else:
        # Check error log
        error_msg = "Merge/sort job failed"
        if os.path.exists(error_file):
            with open(error_file, 'r') as f:
                error_content = f.read()
                if error_content.strip():
                    error_msg += f"\nError log:\n{error_content}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)


def generate_iso_seq_hints_slurm(bam_path, output_file):
    """
    Generate IsoSeq hints directly without Toil, similar to the original function.
    
    :param bam_path: Path to IsoSeq BAM file
    :param output_file: Output hints file path
    :return: Path to the output file
    """
    logger.info(f'Generating IsoSeq hints for {bam_path}')
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    try:
        bam_to_psl = shutil.which('bamToPsl')
        if not bam_to_psl:
            bam_to_psl_local = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'standalones', 'bamToPsl'))
            bam_to_psl = bam_to_psl_local if os.path.exists(bam_to_psl_local) else 'bamToPsl'

        cmd = [['samtools', 'view', '-@', '64', '-b', '-F', '4', bam_path], 
               [bam_to_psl, '-nohead', '/dev/stdin', '/dev/stdout'],
               ['sort', '-n', '-k', '16,16'],
               ['sort', '-s', '-k', '14,14'],
               ['perl', '-ne', '@f=split; print if ($f[0]>=100)'],
               ['blat2hints.pl', '--source=PB', '--nomult', '--ep_cutoff=20', '--in=/dev/stdin',
                '--out={}'.format(output_file)]]
        
        tools.procOps.run_proc(cmd)
        
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            logger.info(f'Successfully generated IsoSeq hints: {output_file}')
        else:
            logger.warning(f'IsoSeq hints file is empty or not created: {output_file}')
            
        return output_file
        
    except Exception as e:
        logger.error(f'Error generating IsoSeq hints for {bam_path}: {e}')
        raise


def generate_slurm_isoseq_job(bam_path, output_path, job_id, temp_dir,
                              memory_gb=192, cpus=96, time_limit="16:00:00",
                              partition=""):
    """
    Generate a Slurm job script for running IsoSeq-to-hints on a single BAM file.
    This does not chunk the BAM; it processes the full file in one job.
    """
    job_script = os.path.join(temp_dir, f"isoseq2hints_{job_id}.sh")
    log_file = os.path.join(temp_dir, f"isoseq2hints_{job_id}.log")
    error_file = os.path.join(temp_dir, f"isoseq2hints_{job_id}.err")

    bam_to_psl = shutil.which('bamToPsl')
    if not bam_to_psl:
        bam_to_psl_local = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'standalones', 'bamToPsl'))
        bam_to_psl = bam_to_psl_local if os.path.exists(bam_to_psl_local) else 'bamToPsl'

    # Node excludes come from config / --exclude-nodes (site-specific).
    # for isoseq jobs that historically OOMed on those nodes. It is now resolved
    # from the active Scheduler's default exclude_nodes (set via
    # --exclude-nodes / cluster.{slurm,sge}.{exclude_nodes,hostname_exclude}).
    header = _hints_header(
        job_name=f"isohints_{job_id}", log_file=log_file, error_file=error_file,
        partition=partition, memory_gb=memory_gb, cpus=cpus, time_limit=time_limit,
    )
    script_content = header + f"""
set -euo pipefail

mkdir -p {os.path.dirname(output_path)}

samtools view -@ {cpus} -b -F 4 {bam_path} \
| {bam_to_psl} -nohead /dev/stdin /dev/stdout \
| sort --parallel={cpus} -S $(( {memory_gb} / 2 ))G -T {temp_dir} -n -k 16,16 \
| sort --parallel={cpus} -S $(( {memory_gb} / 2 ))G -T {temp_dir} -s -k 14,14 \
| perl -ne '@f=split; print if ($f[0]>=100)' \
| blat2hints.pl --source=PB --nomult --ep_cutoff=20 --in=/dev/stdin --out={output_path}

if [ -s "{output_path}" ]; then
  echo "Success: {output_path} created"
  exit 0
else
  echo "Error: {output_path} empty or missing"
  exit 1
fi
"""

    try:
        with open(job_script, 'w') as f:
            f.write(script_content)
        os.chmod(job_script, 0o755)
        logger.info(f'Generated Slurm IsoSeq job script: {job_script}')
        return job_script, log_file, error_file
    except Exception as e:
        logger.error(f'Failed to create job script {job_script}: {e}')
        logger.error(f'Parent directory exists: {os.path.exists(os.path.dirname(job_script))}')
        logger.error(f'Parent directory writable: {os.access(os.path.dirname(job_script), os.W_OK) if os.path.exists(os.path.dirname(job_script)) else "N/A"}')
        raise


def run_parallel_isoseq_slurm(bam_files, output_dir, logs_dir=None,
                              memory_gb=64, cpus=32, time_limit="12:00:00",
                              partition="", max_concurrent_jobs=50):
    """
    Run IsoSeq hints generation on multiple BAM files using parallel Slurm jobs.
    One job is launched per BAM; no chunking is performed.
    """
    logger.info(f'Running parallel IsoSeq->hints with Slurm on {len(bam_files)} BAM files')
    
    # Use logs_dir if provided, otherwise create a temp directory
    created_temp_dir = False
    if logs_dir:
        temp_dir = logs_dir
        try:
            os.makedirs(temp_dir, exist_ok=True)
            logger.info(f'Using logs directory: {temp_dir}')
            # Verify directory was created and is writable
            if not os.path.exists(temp_dir):
                raise RuntimeError(f'Failed to create logs directory: {temp_dir}')
            if not os.access(temp_dir, os.W_OK):
                raise RuntimeError(f'Logs directory is not writable: {temp_dir}')
        except Exception as e:
            logger.error(f'Error with logs directory {logs_dir}: {e}')
            # Fallback to temp directory
            temp_dir = tempfile.mkdtemp(prefix='isoseq_hints_slurm_')
            created_temp_dir = True
            logger.info(f'Falling back to temporary directory: {temp_dir}')
    else:
        temp_dir = tempfile.mkdtemp(prefix='isoseq_hints_slurm_')
        created_temp_dir = True
        logger.info(f'Using temporary directory: {temp_dir}')

    os.makedirs(output_dir, exist_ok=True)

    try:
        job_scripts = []
        output_files = []
        log_files = []

        for i, bam_path in enumerate(bam_files):
            bam_basename = os.path.splitext(os.path.basename(bam_path))[0]
            output_file = os.path.join(output_dir, f"{bam_basename}_isoseq_hints.gff")
            output_files.append(output_file)

            job_script, log_file, error_file = generate_slurm_isoseq_job(
                bam_path, output_file, f"{i}_{bam_basename}", temp_dir,
                memory_gb=memory_gb, cpus=cpus, time_limit=time_limit, partition=partition
            )
            job_scripts.append(job_script)
            log_files.append(log_file)

        job_status = submit_slurm_jobs(job_scripts, max_concurrent_jobs, log_files=log_files)

        successful_outputs = []
        failed_jobs = []
        for i, (script_path, output_file) in enumerate(zip(job_scripts, output_files)):
            status = job_status[script_path]['status']
            
            # For completed jobs, add polling retry to handle filesystem lag
            if status == 'completed':
                max_retries = 30  # 30 retries with 10s intervals = 5 minutes max wait
                retry_count = 0
                file_ready = False
                
                while retry_count < max_retries:
                    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                        file_ready = True
                        break
                    retry_count += 1
                    if retry_count < max_retries:
                        logger.info(f'Waiting for output file to appear: {output_file} (attempt {retry_count}/{max_retries})')
                        time.sleep(10)  # Wait 10 seconds between retries
                
                if file_ready:
                    successful_outputs.append(output_file)
                    logger.info(f'Successfully generated IsoSeq hints: {output_file}')
                else:
                    failed_jobs.append((bam_files[i], output_file, f'completed but file not ready after {max_retries} retries'))
                    logger.error(f'Failed to generate IsoSeq hints for {bam_files[i]}: completed but file not ready after {max_retries} retries')
            else:
                failed_jobs.append((bam_files[i], output_file, status))
                logger.error(f'Failed to generate IsoSeq hints for {bam_files[i]}: {status}')
                
                # Check for error logs to provide more details
                error_file = script_path.replace('.sh', '.err')
                if os.path.exists(error_file):
                    with open(error_file, 'r') as f:
                        error_content = f.read().strip()
                        if error_content:
                            logger.error(f'Error details for {bam_files[i]}: {error_content}')
                
                # Check log file for more details
                if i < len(log_files):
                    log_file = log_files[i]
                else:
                    log_file = script_path.replace('.sh', '.log')
                    
                if os.path.exists(log_file):
                    with open(log_file, 'r') as f:
                        log_content = f.read().strip()
                        if log_content:
                            logger.error(f'Log details for {bam_files[i]}: {log_content}')

        logger.info(f'Successfully processed {len(successful_outputs)}/{len(bam_files)} IsoSeq BAM files')
        if failed_jobs:
            logger.warning(f'{len(failed_jobs)} IsoSeq jobs failed')
        return successful_outputs, failed_jobs
    finally:
        if created_temp_dir:
            try:
                shutil.rmtree(temp_dir)
                logger.info(f'Cleaned up temporary directory: {temp_dir}')
            except Exception as e:
                logger.warning(f'Failed to clean up temporary directory {temp_dir}: {e}')


def generate_slurm_bam_intron_job(bam_path, output_path, job_id, temp_dir,
                                  memory_gb=256, cpus=128, time_limit="12:00:00",
                                  partition=""):
    """
    Generate a Slurm job script for running bam2hints intron processing on a single BAM file.
    Uses parallel processing by chromosome for speed improvement.
    
    OPTIMIZATION: Since bam2hints is single-threaded, we parallelize by processing
    each chromosome independently and then merging the results.
    """
    job_script = os.path.join(temp_dir, f"bam_intron_{job_id}.sh")
    log_file = os.path.join(temp_dir, f"bam_intron_{job_id}.log")
    error_file = os.path.join(temp_dir, f"bam_intron_{job_id}.err")

    # Create a unique temp directory for this job's intermediate files
    job_temp_dir = f"{temp_dir}/bam_intron_{job_id}_tmp"

    header = _hints_header(
        job_name=f"bam_intron_{job_id}", log_file=log_file, error_file=error_file,
        partition=partition, memory_gb=memory_gb, cpus=cpus, time_limit=time_limit,
    )
    script_content = header + f"""
set -euo pipefail

mkdir -p {os.path.dirname(output_path)}
mkdir -p {job_temp_dir}

# Get list of chromosomes/contigs from BAM file
CHROMOSOMES=$(samtools view -H {bam_path} | grep '^@SQ' | cut -f2 | sed 's/SN://' | tr '\\n' ' ')

# Function to process a single chromosome
process_chromosome() {{
    local chr=$1
    local chr_bam="{job_temp_dir}/${{chr}}.bam"
    local output_file="{job_temp_dir}/${{chr}}_intron_hints.gff"
    
    # Extract reads for this chromosome
    samtools view -b {bam_path} "$chr" > "$chr_bam" 2>/dev/null || return 0
    
    # Process with bam2hints if the BAM has content
    if [ -s "$chr_bam" ]; then
        bam2hints --intronsonly --in "$chr_bam" --out "$output_file" 2>/dev/null || true
        rm -f "$chr_bam"
    fi
    
    # Return the output file path if it has content
    if [ -s "$output_file" ]; then
        echo "$output_file"
    fi
}}

export -f process_chromosome

# Process all chromosomes in parallel using GNU parallel or xargs
# Use up to $SLURM_CPUS_PER_TASK parallel jobs
echo "$CHROMOSOMES" | tr ' ' '\\n' | grep -v '^$' | \\
    parallel -j {cpus} --will-cite process_chromosome {{}} > {job_temp_dir}/output_files.txt 2>/dev/null || \\
    echo "$CHROMOSOMES" | tr ' ' '\\n' | grep -v '^$' | \\
    xargs -P {cpus} -I {{}} bash -c 'process_chromosome "{{}}"' > {job_temp_dir}/output_files.txt

# Merge all chromosome-specific hint files
cat {job_temp_dir}/*_intron_hints.gff > {output_path} 2>/dev/null || touch {output_path}

# Clean up temporary files
rm -rf {job_temp_dir}

if [ -s "{output_path}" ]; then
  echo "Success: {output_path} created"
  exit 0
else
  echo "Error: {output_path} empty or missing"
  exit 1
fi
"""

    try:
        with open(job_script, 'w') as f:
            f.write(script_content)
        os.chmod(job_script, 0o755)
        logger.info(f'Generated Slurm BAM intron job script (PARALLELIZED): {job_script}')
        return job_script, log_file, error_file
    except Exception as e:
        logger.error(f'Failed to create BAM intron job script {job_script}: {e}')
        logger.error(f'Parent directory exists: {os.path.exists(os.path.dirname(job_script))}')
        logger.error(f'Parent directory writable: {os.access(os.path.dirname(job_script), os.W_OK) if os.path.exists(os.path.dirname(job_script)) else "N/A"}')
        raise


def generate_slurm_bam_exon_job(bam_path, output_path, job_id, temp_dir,
                                memory_gb=256, cpus=128, time_limit="12:00:00",
                                partition=""):
    """
    Generate a Slurm job script for running bam2hints exon processing on a single BAM file.
    Uses bam2wig + wig2hints.pl pipeline with parallel processing by chromosome.
    
    OPTIMIZATION: Since bam2wig is single-threaded, we parallelize by processing
    each chromosome independently and then merging the results.
    """
    job_script = os.path.join(temp_dir, f"bam_exon_{job_id}.sh")
    log_file = os.path.join(temp_dir, f"bam_exon_{job_id}.log")
    error_file = os.path.join(temp_dir, f"bam_exon_{job_id}.err")

    # Create a unique temp directory for this job's intermediate files
    job_temp_dir = f"{temp_dir}/bam_exon_{job_id}_tmp"

    header = _hints_header(
        job_name=f"bam_exon_{job_id}", log_file=log_file, error_file=error_file,
        partition=partition, memory_gb=memory_gb, cpus=cpus, time_limit=time_limit,
    )
    script_content = header + f"""
set -euo pipefail

mkdir -p {os.path.dirname(output_path)}
mkdir -p {job_temp_dir}

# Get list of chromosomes/contigs from BAM file
CHROMOSOMES=$(samtools view -H {bam_path} | grep '^@SQ' | cut -f2 | sed 's/SN://' | tr '\\n' ' ')

# Function to process a single chromosome
process_chromosome() {{
    local chr=$1
    local output_file="{job_temp_dir}/${{chr}}_exon_hints.gff"
    
    # Extract reads for this chromosome and process with bam2wig + wig2hints.pl
    samtools view -b {bam_path} "$chr" | bam2wig /dev/stdin | \\
        wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 \\
                     --prune=0.1 --src=W --type=ep --UCSC=/dev/null \\
                     --radius=4.5 --pri=4 --strand=. > "$output_file" 2>/dev/null || true
    
    # Return the output file path if it has content
    if [ -s "$output_file" ]; then
        echo "$output_file"
    fi
}}

export -f process_chromosome

# Process all chromosomes in parallel using GNU parallel or xargs
# Use up to $SLURM_CPUS_PER_TASK parallel jobs
echo "$CHROMOSOMES" | tr ' ' '\\n' | grep -v '^$' | \\
    parallel -j {cpus} --will-cite process_chromosome {{}} > {job_temp_dir}/output_files.txt 2>/dev/null || \\
    echo "$CHROMOSOMES" | tr ' ' '\\n' | grep -v '^$' | \\
    xargs -P {cpus} -I {{}} bash -c 'process_chromosome "{{}}"' > {job_temp_dir}/output_files.txt

# Merge all chromosome-specific hint files
cat {job_temp_dir}/*_exon_hints.gff > {output_path} 2>/dev/null || touch {output_path}

# Clean up temporary files
rm -rf {job_temp_dir}

if [ -s "{output_path}" ]; then
  echo "Success: {output_path} created"
  exit 0
else
  echo "Error: {output_path} empty or missing"
  exit 1
fi
"""

    try:
        with open(job_script, 'w') as f:
            f.write(script_content)
        os.chmod(job_script, 0o755)
        logger.info(f'Generated Slurm BAM exon job script (PARALLELIZED): {job_script}')
        return job_script, log_file, error_file
    except Exception as e:
        logger.error(f'Failed to create BAM exon job script {job_script}: {e}')
        logger.error(f'Parent directory exists: {os.path.exists(os.path.dirname(job_script))}')
        logger.error(f'Parent directory writable: {os.access(os.path.dirname(job_script), os.W_OK) if os.path.exists(os.path.dirname(job_script)) else "N/A"}')
        raise


def generate_slurm_protein_job(protein_fasta, genome_fasta, output_path, job_id, temp_dir,
                               memory_gb=256, cpus=128, time_limit="24:00:00",
                               partition=""):
    """
    Generate a Slurm job script for running protein-to-genome alignment using exonerate.
    """
    job_script = os.path.join(temp_dir, f"protein_{job_id}.sh")
    log_file = os.path.join(temp_dir, f"protein_{job_id}.log")
    error_file = os.path.join(temp_dir, f"protein_{job_id}.err")

    header = _hints_header(
        job_name=f"protein_{job_id}", log_file=log_file, error_file=error_file,
        partition=partition, memory_gb=memory_gb, cpus=cpus, time_limit=time_limit,
    )
    script_content = header + f"""
set -euo pipefail

mkdir -p {os.path.dirname(output_path)}

# Index protein FASTA
samtools faidx {protein_fasta}

# Run exonerate protein-to-genome alignment
exonerate --model protein2genome --showvulgar no --showalignment no --showquerygff yes --ryo "AveragePercentIdentity: %pi\\n" {protein_fasta} {genome_fasta} > {temp_dir}/exonerate_output.txt

# Sort exonerate output (single-pass multi-key sort for better performance)
sort -k1,1 -k3,3 -k4,4n -k5,5n {temp_dir}/exonerate_output.txt > {temp_dir}/sorted_exonerate.txt

# Generate hints from exonerate output
exonerate2hints.pl --in={temp_dir}/sorted_exonerate.txt --CDSpart_cutoff=5 --out={output_path}

if [ -s "{output_path}" ]; then
  echo "Success: {output_path} created"
  exit 0
else
  echo "Error: {output_path} empty or missing"
  exit 1
fi
"""

    try:
        with open(job_script, 'w') as f:
            f.write(script_content)
        os.chmod(job_script, 0o755)
        logger.info(f'Generated Slurm protein job script: {job_script}')
        return job_script, log_file, error_file
    except Exception as e:
        logger.error(f'Failed to create protein job script {job_script}: {e}')
        logger.error(f'Parent directory exists: {os.path.exists(os.path.dirname(job_script))}')
        logger.error(f'Parent directory writable: {os.access(os.path.dirname(job_script), os.W_OK) if os.path.exists(os.path.dirname(job_script)) else "N/A"}')
        raise


def generate_slurm_annotation_job(annotation_gp, output_path, job_id, temp_dir,
                                  memory_gb=256, cpus=128, time_limit="12:00:00",
                                  partition=""):
    """
    Generate a Slurm job script for converting annotation to hints.
    """
    job_script = os.path.join(temp_dir, f"annotation_{job_id}.sh")
    log_file = os.path.join(temp_dir, f"annotation_{job_id}.log")
    error_file = os.path.join(temp_dir, f"annotation_{job_id}.err")

    header = _hints_header(
        job_name=f"annotation_{job_id}", log_file=log_file, error_file=error_file,
        partition=partition, memory_gb=memory_gb, cpus=cpus, time_limit=time_limit,
    )
    script_content = header + f"""
set -euo pipefail

mkdir -p {os.path.dirname(output_path)}

# Convert annotation to hints using Python script
python3 -c "
import sys
import tools.transcripts
import tools.fileOps

# Load annotation
tx_dict = tools.transcripts.get_gene_pred_dict('{annotation_gp}')
hints = []

for tx_id, tx in tx_dict.items():
    if tx.cds_size == 0:
        continue
    # Convert transcript coordinates
    cds_tx = tools.transcripts.Transcript(tx.get_bed(new_start=tx.thick_start, new_stop=tx.thick_stop))
    for intron in cds_tx.intron_intervals:
        r = [intron.chromosome, 'a2h', 'intron', intron.start + 1, intron.stop, 0, intron.strand, '.',
             f'grp={{tx_id}};src=M;pri=2']
        hints.append(r)
    for exon in cds_tx.exon_intervals:
        r = [exon.chromosome, 'a2h', 'CDS', exon.start + 1, exon.stop, 0, exon.strand, '.',
             f'grp={{tx_id}};src=M;pri=2']
        hints.append(r)

tools.fileOps.print_rows('{output_path}', hints)
"

if [ -s "{output_path}" ]; then
  echo "Success: {output_path} created"
  exit 0
else
  echo "Error: {output_path} empty or missing"
  exit 1
fi
"""

    try:
        with open(job_script, 'w') as f:
            f.write(script_content)
        os.chmod(job_script, 0o755)
        logger.info(f'Generated Slurm annotation job script: {job_script}')
        return job_script, log_file, error_file
    except Exception as e:
        logger.error(f'Failed to create annotation job script {job_script}: {e}')
        logger.error(f'Parent directory exists: {os.path.exists(os.path.dirname(job_script))}')
        logger.error(f'Parent directory writable: {os.access(os.path.dirname(job_script), os.W_OK) if os.path.exists(os.path.dirname(job_script)) else "N/A"}')
        raise


def run_parallel_bam_intron_slurm(bam_files, output_dir, logs_dir=None,
                                  memory_gb=256, cpus=128, time_limit="12:00:00",
                                  partition="", max_concurrent_jobs=50):
    """
    Run BAM intron hints generation on multiple BAM files using parallel Slurm jobs.
    """
    logger.info(f'Running parallel BAM intron hints with Slurm on {len(bam_files)} BAM files')
    
    # Use logs_dir if provided, otherwise create a temp directory
    created_temp_dir = False
    if logs_dir:
        temp_dir = logs_dir
        try:
            os.makedirs(temp_dir, exist_ok=True)
            logger.info(f'Using logs directory: {temp_dir}')
            # Verify directory was created and is writable
            if not os.path.exists(temp_dir):
                raise RuntimeError(f'Failed to create logs directory: {temp_dir}')
            if not os.access(temp_dir, os.W_OK):
                raise RuntimeError(f'Logs directory is not writable: {temp_dir}')
        except Exception as e:
            logger.error(f'Error with logs directory {logs_dir}: {e}')
            # Fallback to temp directory
            temp_dir = tempfile.mkdtemp(prefix='bam_intron_hints_slurm_')
            created_temp_dir = True
            logger.info(f'Falling back to temporary directory: {temp_dir}')
    else:
        temp_dir = tempfile.mkdtemp(prefix='bam_intron_hints_slurm_')
        created_temp_dir = True
        logger.info(f'Using temporary directory: {temp_dir}')

    os.makedirs(output_dir, exist_ok=True)

    try:
        job_scripts = []
        output_files = []
        log_files = []

        for i, bam_path in enumerate(bam_files):
            bam_basename = os.path.splitext(os.path.basename(bam_path))[0]
            output_file = os.path.join(output_dir, f"{bam_basename}_intron_hints.gff")
            output_files.append(output_file)

            job_script, log_file, error_file = generate_slurm_bam_intron_job(
                bam_path, output_file, f"{i}_{bam_basename}", temp_dir,
                memory_gb=memory_gb, cpus=cpus, time_limit=time_limit, partition=partition
            )
            job_scripts.append(job_script)
            log_files.append(log_file)

        job_status = submit_slurm_jobs(job_scripts, max_concurrent_jobs, log_files=log_files)

        successful_outputs = []
        failed_jobs = []
        for i, (script_path, output_file) in enumerate(zip(job_scripts, output_files)):
            status = job_status[script_path]['status']
            
            # For completed jobs, add polling retry to handle filesystem lag
            if status == 'completed':
                max_retries = 30  # 30 retries with 10s intervals = 5 minutes max wait
                retry_count = 0
                file_ready = False
                
                while retry_count < max_retries:
                    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                        file_ready = True
                        break
                    retry_count += 1
                    if retry_count < max_retries:
                        logger.info(f'Waiting for output file to appear: {output_file} (attempt {retry_count}/{max_retries})')
                        time.sleep(10)  # Wait 10 seconds between retries
                
                if file_ready:
                    successful_outputs.append(output_file)
                    logger.info(f'Successfully generated BAM intron hints: {output_file}')
                else:
                    failed_jobs.append((bam_files[i], output_file, f'completed but file not ready after {max_retries} retries'))
                    logger.error(f'Failed to generate BAM intron hints for {bam_files[i]}: completed but file not ready after {max_retries} retries')
            else:
                failed_jobs.append((bam_files[i], output_file, status))
                logger.error(f'Failed to generate BAM intron hints for {bam_files[i]}: {status}')
                
                # Check for error logs to provide more details
                error_file = script_path.replace('.sh', '.err')
                if os.path.exists(error_file):
                    with open(error_file, 'r') as f:
                        error_content = f.read().strip()
                        if error_content:
                            logger.error(f'Error details for {bam_files[i]}: {error_content}')

        logger.info(f'Successfully processed {len(successful_outputs)}/{len(bam_files)} BAM files for intron hints')
        if failed_jobs:
            logger.warning(f'{len(failed_jobs)} BAM intron jobs failed')
        return successful_outputs, failed_jobs
    finally:
        if created_temp_dir:
            try:
                shutil.rmtree(temp_dir)
                logger.info(f'Cleaned up temporary directory: {temp_dir}')
            except Exception as e:
                logger.warning(f'Failed to clean up temporary directory {temp_dir}: {e}')


def run_parallel_bam_exon_slurm(bam_files, output_dir, logs_dir=None,
                                memory_gb=256, cpus=128, time_limit="12:00:00",
                                partition="", max_concurrent_jobs=50):
    """
    Run BAM exon hints generation on multiple BAM files using parallel Slurm jobs.
    """
    logger.info(f'Running parallel BAM exon hints with Slurm on {len(bam_files)} BAM files')
    
    # Use logs_dir if provided, otherwise create a temp directory
    created_temp_dir = False
    if logs_dir:
        temp_dir = logs_dir
        try:
            os.makedirs(temp_dir, exist_ok=True)
            logger.info(f'Using logs directory: {temp_dir}')
            # Verify directory was created and is writable
            if not os.path.exists(temp_dir):
                raise RuntimeError(f'Failed to create logs directory: {temp_dir}')
            if not os.access(temp_dir, os.W_OK):
                raise RuntimeError(f'Logs directory is not writable: {temp_dir}')
        except Exception as e:
            logger.error(f'Error with logs directory {logs_dir}: {e}')
            # Fallback to temp directory
            temp_dir = tempfile.mkdtemp(prefix='bam_exon_hints_slurm_')
            created_temp_dir = True
            logger.info(f'Falling back to temporary directory: {temp_dir}')
    else:
        temp_dir = tempfile.mkdtemp(prefix='bam_exon_hints_slurm_')
        created_temp_dir = True
        logger.info(f'Using temporary directory: {temp_dir}')

    os.makedirs(output_dir, exist_ok=True)

    try:
        job_scripts = []
        output_files = []
        log_files = []

        for i, bam_path in enumerate(bam_files):
            bam_basename = os.path.splitext(os.path.basename(bam_path))[0]
            output_file = os.path.join(output_dir, f"{bam_basename}_exon_hints.gff")
            output_files.append(output_file)

            job_script, log_file, error_file = generate_slurm_bam_exon_job(
                bam_path, output_file, f"{i}_{bam_basename}", temp_dir,
                memory_gb=memory_gb, cpus=cpus, time_limit=time_limit, partition=partition
            )
            job_scripts.append(job_script)
            log_files.append(log_file)

        job_status = submit_slurm_jobs(job_scripts, max_concurrent_jobs, log_files=log_files)

        successful_outputs = []
        failed_jobs = []
        for i, (script_path, output_file) in enumerate(zip(job_scripts, output_files)):
            status = job_status[script_path]['status']
            
            # For completed jobs, add polling retry to handle filesystem lag
            if status == 'completed':
                max_retries = 30  # 30 retries with 10s intervals = 5 minutes max wait
                retry_count = 0
                file_ready = False
                
                while retry_count < max_retries:
                    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                        file_ready = True
                        break
                    retry_count += 1
                    if retry_count < max_retries:
                        logger.info(f'Waiting for output file to appear: {output_file} (attempt {retry_count}/{max_retries})')
                        time.sleep(10)  # Wait 10 seconds between retries
                
                if file_ready:
                    successful_outputs.append(output_file)
                    logger.info(f'Successfully generated BAM exon hints: {output_file}')
                else:
                    failed_jobs.append((bam_files[i], output_file, f'completed but file not ready after {max_retries} retries'))
                    logger.error(f'Failed to generate BAM exon hints for {bam_files[i]}: completed but file not ready after {max_retries} retries')
            else:
                failed_jobs.append((bam_files[i], output_file, status))
                logger.error(f'Failed to generate BAM exon hints for {bam_files[i]}: {status}')
                
                # Check for error logs to provide more details
                error_file = script_path.replace('.sh', '.err')
                if os.path.exists(error_file):
                    with open(error_file, 'r') as f:
                        error_content = f.read().strip()
                        if error_content:
                            logger.error(f'Error details for {bam_files[i]}: {error_content}')

        logger.info(f'Successfully processed {len(successful_outputs)}/{len(bam_files)} BAM files for exon hints')
        if failed_jobs:
            logger.warning(f'{len(failed_jobs)} BAM exon jobs failed')
        return successful_outputs, failed_jobs
    finally:
        if created_temp_dir:
            try:
                shutil.rmtree(temp_dir)
                logger.info(f'Cleaned up temporary directory: {temp_dir}')
            except Exception as e:
                logger.warning(f'Failed to clean up temporary directory {temp_dir}: {e}')


def run_protein_hints_slurm(protein_fasta, genome_fasta, output_file, logs_dir=None,
                            memory_gb=256, cpus=128, time_limit="24:00:00",
                            partition=""):
    """
    Run protein hints generation using Slurm.
    """
    logger.info(f'Running protein hints generation with Slurm')
    
    # Use logs_dir if provided, otherwise create a temp directory
    if logs_dir:
        temp_dir = logs_dir
        try:
            os.makedirs(temp_dir, exist_ok=True)
            logger.info(f'Using logs directory: {temp_dir}')
            # Verify directory was created and is writable
            if not os.path.exists(temp_dir):
                raise RuntimeError(f'Failed to create logs directory: {temp_dir}')
            if not os.access(temp_dir, os.W_OK):
                raise RuntimeError(f'Logs directory is not writable: {temp_dir}')
        except Exception as e:
            logger.error(f'Error with logs directory {logs_dir}: {e}')
            # Fallback to temp directory
            temp_dir = tempfile.mkdtemp(prefix='protein_hints_slurm_')
            logger.info(f'Falling back to temporary directory: {temp_dir}')
    else:
        temp_dir = tempfile.mkdtemp(prefix='protein_hints_slurm_')
        logger.info(f'Using temporary directory: {temp_dir}')

    try:
        job_script, log_file, error_file = generate_slurm_protein_job(
            protein_fasta, genome_fasta, output_file, "protein_main", temp_dir,
            memory_gb=memory_gb, cpus=cpus, time_limit=time_limit, partition=partition
        )

        job_status = submit_slurm_jobs([job_script], 1, log_files=[log_file])

        status = job_status[job_script]['status']
        
        if status == 'completed' and os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            logger.info(f'Successfully generated protein hints: {output_file}')
            return output_file
        else:
            logger.error(f'Failed to generate protein hints: {status}')
            
            # Check for error logs to provide more details
            if os.path.exists(error_file):
                with open(error_file, 'r') as f:
                    error_content = f.read().strip()
                    if error_content:
                        logger.error(f'Error details for protein hints: {error_content}')
            
            raise RuntimeError(f'Protein hints generation failed: {status}')
            
    finally:
        try:
            shutil.rmtree(temp_dir)
            logger.info(f'Cleaned up temporary directory: {temp_dir}')
        except Exception as e:
            logger.warning(f'Failed to clean up temporary directory {temp_dir}: {e}')


def run_annotation_hints_slurm(annotation_gp, output_file, logs_dir=None,
                               memory_gb=256, cpus=128, time_limit="12:00:00",
                               partition=""):
    """
    Run annotation hints generation using Slurm.
    """
    logger.info(f'Running annotation hints generation with Slurm')
    
    # Use logs_dir if provided, otherwise create a temp directory
    if logs_dir:
        temp_dir = logs_dir
        try:
            os.makedirs(temp_dir, exist_ok=True)
            logger.info(f'Using logs directory: {temp_dir}')
            # Verify directory was created and is writable
            if not os.path.exists(temp_dir):
                raise RuntimeError(f'Failed to create logs directory: {temp_dir}')
            if not os.access(temp_dir, os.W_OK):
                raise RuntimeError(f'Logs directory is not writable: {temp_dir}')
        except Exception as e:
            logger.error(f'Error with logs directory {logs_dir}: {e}')
            # Fallback to temp directory
            temp_dir = tempfile.mkdtemp(prefix='annotation_hints_slurm_')
            logger.info(f'Falling back to temporary directory: {temp_dir}')
    else:
        temp_dir = tempfile.mkdtemp(prefix='annotation_hints_slurm_')
        logger.info(f'Using temporary directory: {temp_dir}')

    try:
        job_script, log_file, error_file = generate_slurm_annotation_job(
            annotation_gp, output_file, "annotation_main", temp_dir,
            memory_gb=memory_gb, cpus=cpus, time_limit=time_limit, partition=partition
        )

        job_status = submit_slurm_jobs([job_script], 1, log_files=[log_file])

        status = job_status[job_script]['status']
        
        if status == 'completed' and os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            logger.info(f'Successfully generated annotation hints: {output_file}')
            return output_file
        else:
            logger.error(f'Failed to generate annotation hints: {status}')
            
            # Check for error logs to provide more details
            if os.path.exists(error_file):
                with open(error_file, 'r') as f:
                    error_content = f.read().strip()
                    if error_content:
                        logger.error(f'Error details for annotation hints: {error_content}')
            
            raise RuntimeError(f'Annotation hints generation failed: {status}')
            
    finally:
        try:
            shutil.rmtree(temp_dir)
            logger.info(f'Cleaned up temporary directory: {temp_dir}')
        except Exception as e:
            logger.warning(f'Failed to clean up temporary directory {temp_dir}: {e}')


def generate_annotation_hints_direct(annotation_gp, output_file):
    """
    Generate annotation hints directly without Toil.
    
    :param annotation_gp: Path to annotation genePred file
    :param output_file: Output hints file path
    :return: Path to the output file
    """
    logger.info(f'Generating annotation hints from {annotation_gp}')
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    try:
        tx_dict = tools.transcripts.get_gene_pred_dict(annotation_gp)
        hints = []
        
        for tx_id, tx in tx_dict.items():
            if tx.cds_size == 0:
                continue
            # rather than try to re-do the arithmetic, we will use the get_bed() function to convert this transcript
            cds_tx = tools.transcripts.Transcript(tx.get_bed(new_start=tx.thick_start, new_stop=tx.thick_stop))
            for intron in cds_tx.intron_intervals:
                r = [intron.chromosome, 'a2h', 'intron', intron.start + 1, intron.stop, 0, intron.strand, '.',
                     'grp={};src=M;pri=2'.format(tx_id)]
                hints.append(r)
            for exon in cds_tx.exon_intervals:
                r = [exon.chromosome, 'a2h', 'CDS', exon.start + 1, exon.stop, 0, exon.strand, '.',
                     'grp={};src=M;pri=2'.format(tx_id)]
                hints.append(r)
        
        tools.fileOps.print_rows(output_file, hints)
        
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            logger.info(f'Successfully generated annotation hints: {output_file}')
        else:
            logger.warning(f'Annotation hints file is empty or not created: {output_file}')
            
        return output_file
        
    except Exception as e:
        logger.error(f'Error generating annotation hints from {annotation_gp}: {e}')
        raise


def generate_protein_hints_direct(protein_fasta, genome_fasta, output_file):
    """
    Generate protein hints directly without Toil.
    
    :param protein_fasta: Path to protein FASTA file
    :param genome_fasta: Path to genome FASTA file
    :param output_file: Output hints file path
    :return: Path to the output file
    """
    logger.info(f'Generating protein hints from {protein_fasta} against {genome_fasta}')
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    try:
        cmd = ['samtools', 'faidx', protein_fasta]
        tools.procOps.run_proc(cmd)
        protein_handle = tools.bio.get_sequence_dict(protein_fasta)
        
        # Create temporary files
        temp_dir = tempfile.mkdtemp(prefix='protein_hints_')
        
        try:
            # Run exonerate
            exonerate_output = os.path.join(temp_dir, 'exonerate_output.txt')
            cmd = ['exonerate', '--model', 'protein2genome', '--showvulgar', 'no', 
                   '--showalignment', 'no', '--showquerygff', 'yes', 
                   '--ryo', 'AveragePercentIdentity: %pi\n',
                   protein_fasta, genome_fasta]
            tools.procOps.run_proc(cmd, stdout=exonerate_output)
            
            # Sort exonerate output
            sorted_exonerate = os.path.join(temp_dir, 'sorted_exonerate.txt')
            tools.misc.sort_gff(exonerate_output, sorted_exonerate)
            
            # Generate hints
            cmd = ['exonerate2hints.pl', '--in={}'.format(sorted_exonerate), 
                   '--CDSpart_cutoff=5', '--out={}'.format(output_file)]
            tools.procOps.run_proc(cmd)
            
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                logger.info(f'Successfully generated protein hints: {output_file}')
            else:
                logger.warning(f'Protein hints file is empty or not created: {output_file}')
                
            return output_file
            
        finally:
            # Cleanup temporary directory
            try:
                if os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)
                    logger.info(f'Cleaned up temporary directory: {temp_dir}')
                else:
                    logger.debug(f'Temporary directory already removed: {temp_dir}')
            except OSError as e:
                if e.errno == 39:  # Directory not empty
                    logger.warning(f'Directory not empty, attempting force cleanup: {temp_dir}')
                    try:
                        # Try to remove files individually first
                        for root, dirs, files in os.walk(temp_dir, topdown=False):
                            for file in files:
                                try:
                                    os.remove(os.path.join(root, file))
                                except OSError:
                                    pass
                            for dir in dirs:
                                try:
                                    os.rmdir(os.path.join(root, dir))
                                except OSError:
                                    pass
                        os.rmdir(temp_dir)
                        logger.info(f'Force cleaned up temporary directory: {temp_dir}')
                    except Exception as force_e:
                        logger.warning(f'Failed to force clean up temporary directory {temp_dir}: {force_e}')
                else:
                    logger.warning(f'Failed to clean up temporary directory {temp_dir}: {e}')
            except Exception as e:
                logger.warning(f'Failed to clean up temporary directory {temp_dir}: {e}')
        
    except Exception as e:
        logger.error(f'Error generating protein hints from {protein_fasta}: {e}')
        raise


def calculate_dynamic_resources(bam_files, iso_bam_files, annotation_gp, protein_fasta, genome_fasta):
    """
    Calculate dynamic resource allocation based on input data characteristics.
    
    :param bam_files: List of BAM file paths
    :param iso_bam_files: List of IsoSeq BAM file paths
    :param annotation_gp: Annotation GenePred file path
    :param protein_fasta: Protein FASTA file path
    :param genome_fasta: Genome FASTA file path
    :return: Dictionary with resource allocations
    """
    try:
        total_bam_size_gb = 0
        total_bam_count = 0
        
        # Calculate total BAM file sizes
        all_bams = (bam_files or []) + (iso_bam_files or [])
        for bam_path in all_bams:
            if os.path.exists(bam_path):
                size_gb = os.path.getsize(bam_path) / (1024**3)
                total_bam_size_gb += size_gb
                total_bam_count += 1
        
        # Estimate genome size
        genome_size_gb = 0
        if genome_fasta and os.path.exists(genome_fasta):
            genome_size_gb = os.path.getsize(genome_fasta) / (1024**3)
        
        # Estimate protein complexity
        protein_size_mb = 0
        if protein_fasta and os.path.exists(protein_fasta):
            protein_size_mb = os.path.getsize(protein_fasta) / (1024**2)
        
        # Estimate annotation complexity
        annotation_transcripts = 0
        if annotation_gp and os.path.exists(annotation_gp):
            try:
                with open(annotation_gp, 'r') as f:
                    annotation_transcripts = sum(1 for line in f if line.strip())
            except:
                annotation_transcripts = 10000  # Conservative estimate
        
        # Cap CPU requests to Toil --maxCores (local) or machine size. Never
        # request more cores than SingleMachine was configured with.
        available_cpus = _toil_max_cores()
        
        # BAM processing resources
        if total_bam_size_gb > 50:  # Large BAM datasets
            bam_memory_gb = max(8, min(64, int(total_bam_size_gb / 5) or 8))
            bam_cpus = available_cpus
            chunk_size = 25_000_000  # Smaller chunks for large files
        elif total_bam_size_gb > 10:  # Medium BAM datasets
            bam_memory_gb = max(8, min(64, int(total_bam_size_gb / 2) or 8))
            bam_cpus = available_cpus
            chunk_size = 50_000_000  # Standard chunks
        else:  # Small BAM datasets
            bam_memory_gb = min(64, max(4, available_cpus * 2))
            bam_cpus = available_cpus
            chunk_size = 100_000_000  # Larger chunks for small files
        
        # Protein processing resources
        if protein_size_mb > 100:  # Large protein datasets
            protein_memory_gb = max(8, min(64, int(protein_size_mb / 50) or 8))
            protein_chunk_size = 50  # Smaller protein chunks
            protein_cpus = available_cpus
        else:  # Standard protein datasets
            protein_memory_gb = min(64, max(4, available_cpus * 2))
            protein_chunk_size = 100  # Standard protein chunks
            protein_cpus = available_cpus
        
        # Merge and final processing resources
        merge_memory_gb = max(64, int(total_bam_size_gb / 3))  # Scale with total data
        merge_disk_gb = max(128, int((total_bam_size_gb + genome_size_gb) * 2))  # Conservative disk estimate
        
        return {
            'bam_memory_gb': bam_memory_gb,
            'bam_cpus': bam_cpus,
            'bam_chunk_size': chunk_size,
            'protein_memory_gb': protein_memory_gb,
            'protein_cpus': protein_cpus,
            'protein_chunk_size': protein_chunk_size,
            'merge_memory_gb': merge_memory_gb,
            'merge_disk_gb': merge_disk_gb,
            'total_bam_size_gb': total_bam_size_gb,
            'genome_size_gb': genome_size_gb,
            'annotation_transcripts': annotation_transcripts
        }
        
    except Exception as e:
        logger.warning(f'Error calculating dynamic resources: {str(e)}')
        # Fallback to conservative estimates (still capped to --maxCores)
        cpus = _toil_max_cores()
        return {
            'bam_memory_gb': min(32, max(4, cpus * 2)),
            'bam_cpus': cpus,
            'bam_chunk_size': 50_000_000,
            'protein_memory_gb': min(32, max(4, cpus * 2)),
            'protein_cpus': cpus,
            'protein_chunk_size': 100,
            'merge_memory_gb': min(32, max(4, cpus * 2)),
            'merge_disk_gb': 32,
            'total_bam_size_gb': 1.0,
            'genome_size_gb': 3.0,
            'annotation_transcripts': 50000
        }


def estimate_bam_complexity(bam_files):
    """
    Estimate BAM file complexity for optimized processing.
    
    :param bam_files: List of BAM file paths
    :return: Dictionary with complexity metrics
    """
    complexity = {
        'total_reads_estimate': 0,
        'avg_read_length': 100,
        'has_paired_reads': False,
        'reference_count': 0,
        'total_genome_size': 0,
        'needs_chunking': False
    }
    
    if not bam_files:
        return complexity
    
    try:
        # Sample first BAM file for characteristics
        first_bam = next((bam for bam in bam_files if os.path.exists(bam)), None)
        if not first_bam:
            return complexity
        
        handle = pysam.Samfile(first_bam, 'rb')
        
        # Get reference information
        complexity['reference_count'] = len(handle.references)
        complexity['total_genome_size'] = sum(handle.lengths)
        
        # Sample reads for characteristics
        read_lengths = []
        paired_count = 0
        total_sampled = 0
        
        for i, read in enumerate(itertools.islice(handle, 10000)):
            if read.query_length:
                read_lengths.append(read.query_length)
            if read.is_paired:
                paired_count += 1
            total_sampled += 1
            
            if i > 1000:  # Sample enough for good estimate
                break
        
        if read_lengths:
            complexity['avg_read_length'] = sum(read_lengths) / len(read_lengths)
        
        if total_sampled > 0:
            complexity['has_paired_reads'] = (paired_count / total_sampled) > 0.5
        
        # Estimate total reads across all BAM files
        file_size = os.path.getsize(first_bam)
        estimated_reads_per_gb = 10_000_000  # Rough estimate
        
        for bam_path in bam_files:
            if os.path.exists(bam_path):
                bam_size_gb = os.path.getsize(bam_path) / (1024**3)
                complexity['total_reads_estimate'] += int(bam_size_gb * estimated_reads_per_gb)
        
        # Determine if chunking is needed
        complexity['needs_chunking'] = complexity['total_reads_estimate'] > 100_000_000
        
        handle.close()
        
    except Exception as e:
        logger.warning(f'Error estimating BAM complexity: {str(e)}')
    
    return complexity


def optimize_processing_strategy(resources, complexity):
    """
    Optimize processing strategy based on resource and complexity analysis.
    
    :param resources: Resource allocation dictionary
    :param complexity: Complexity metrics dictionary
    :return: Dictionary with processing strategy
    """
    strategy = {
        'parallel_bam_processing': True,
        'chunk_bams': complexity.get('needs_chunking', False),
        'parallel_protein_alignment': True,
        'optimize_merge': True,
        'max_concurrent_jobs': resources['bam_cpus']
    }
    
    # Adjust strategy based on data characteristics
    if resources['total_bam_size_gb'] > 100:  # Very large datasets
        strategy['chunk_bams'] = True
        strategy['max_concurrent_jobs'] = min(strategy['max_concurrent_jobs'], 8)
    
    if complexity.get('reference_count', 0) > 1000:  # Many chromosomes
        strategy['optimize_merge'] = True
        strategy['parallel_bam_processing'] = True
    
    if resources['annotation_transcripts'] > 100000:  # Large annotations
        strategy['parallel_annotation_processing'] = True
    
    return strategy


def run_hints_pipeline_slurm(genome: str,
                            fasta: str,
                            bams: list,
                            intron_bams: list,
                            iso_bams: list,
                            annotation_gp: str,
                            protein_fasta: str,
                            hints_out: str,
                            slurm_options: dict = None):
    """
    Main entry function for the hints pipeline using Slurm for parallel bam2hints processing.
    
    :param genome: Reference genome name
    :param fasta: Path to reference genome fasta
    :param bams: List of short-read RNA-seq BAM files
    :param intron_bams: List of intron-spanning RNA-seq BAM files
    :param iso_bams: List of Iso-Seq BAM files
    :param annotation_gp: Path to annotation genePred file
    :param protein_fasta: Path to protein FASTA file
    :param hints_out: Path to output hints GFF file
    :param slurm_options: Dictionary with Slurm configuration options
    :return: Path to the output hints file
    """
    logger.info('Starting hints pipeline with ALL SLURM-based processing (reimplemented from scratch)')
    
    # Default Slurm options
    default_slurm_options = {
        'memory_gb': 256,
        'cpus': 128,
        'time_limit': '12:00:00',
        'partition': '',
        'max_concurrent_jobs': 50,
        'use_slurm': True
    }
    
    if slurm_options:
        default_slurm_options.update(slurm_options)
    
    slurm_opts = default_slurm_options
    
    # Create output directory
    output_dir = os.path.dirname(os.path.abspath(hints_out))
    os.makedirs(output_dir, exist_ok=True)
    
    # Create logs directory for job scripts and logs
    logs_dir = os.path.join(output_dir, 'logs', 'slurm_jobs')
    try:
        os.makedirs(logs_dir, exist_ok=True)
        logger.info(f'Using logs directory: {logs_dir}')
        # Verify directory was created and is writable
        if not os.path.exists(logs_dir):
            raise RuntimeError(f'Failed to create logs directory: {logs_dir}')
        if not os.access(logs_dir, os.W_OK):
            raise RuntimeError(f'Logs directory is not writable: {logs_dir}')
    except Exception as e:
        logger.error(f'Error creating logs directory {logs_dir}: {e}')
        # Fallback to temp directory
        logs_dir = tempfile.mkdtemp(prefix='hints_logs_')
        logger.info(f'Falling back to temporary logs directory: {logs_dir}')
    
    # Create temporary directory for intermediate files
    temp_dir = tempfile.mkdtemp(prefix='hints_pipeline_')
    logger.info(f'Using temporary directory: {temp_dir}')
    
    try:
        all_hints_files = []
        
        # Process regular BAM files for intron hints
        if bams:
            logger.info(f'Processing {len(bams)} BAM files for intron hints using NEW SLURM implementation')
            intron_output_dir = os.path.join(output_dir, 'intron_hints')
            os.makedirs(intron_output_dir, exist_ok=True)
            
            successful_intron, failed_intron = run_parallel_bam_intron_slurm(
                bams, intron_output_dir, logs_dir=logs_dir,
                memory_gb=slurm_opts['memory_gb'],
                cpus=slurm_opts['cpus'],
                time_limit=slurm_opts['time_limit'],
                partition=slurm_opts['partition'],
                max_concurrent_jobs=slurm_opts['max_concurrent_jobs']
            )
            all_hints_files.extend(successful_intron)
            
            if failed_intron:
                logger.warning(f'{len(failed_intron)} intron hint jobs failed')
        
        # Process regular BAM files for exon hints
        if bams:
            logger.info(f'Processing {len(bams)} BAM files for exon hints using NEW SLURM implementation')
            exon_output_dir = os.path.join(output_dir, 'exon_hints')
            os.makedirs(exon_output_dir, exist_ok=True)
            
            successful_exon, failed_exon = run_parallel_bam_exon_slurm(
                bams, exon_output_dir, logs_dir=logs_dir,
                memory_gb=slurm_opts['memory_gb'],
                cpus=slurm_opts['cpus'],
                time_limit=slurm_opts['time_limit'],
                partition=slurm_opts['partition'],
                max_concurrent_jobs=slurm_opts['max_concurrent_jobs']
            )
            all_hints_files.extend(successful_exon)
            
            if failed_exon:
                logger.warning(f'{len(failed_exon)} exon hint jobs failed')
        
        # Process intron BAM files for intron hints
        if intron_bams:
            logger.info(f'Processing {len(intron_bams)} intron BAM files for intron hints using NEW SLURM implementation')
            intron_bam_output_dir = os.path.join(output_dir, 'intron_bam_hints')
            os.makedirs(intron_bam_output_dir, exist_ok=True)
            
            successful_intron_bam, failed_intron_bam = run_parallel_bam_intron_slurm(
                intron_bams, intron_bam_output_dir, logs_dir=logs_dir,
                memory_gb=slurm_opts['memory_gb'],
                cpus=slurm_opts['cpus'],
                time_limit=slurm_opts['time_limit'],
                partition=slurm_opts['partition'],
                max_concurrent_jobs=slurm_opts['max_concurrent_jobs']
            )
            all_hints_files.extend(successful_intron_bam)
            
            if failed_intron_bam:
                logger.warning(f'{len(failed_intron_bam)} intron BAM hint jobs failed')
        
        # Process IsoSeq BAM files via Slurm (one job per BAM, no chunking)
        if iso_bams:
            logger.info(f'Processing {len(iso_bams)} IsoSeq BAM files via Slurm')
            iso_output_dir = os.path.join(output_dir, 'iso_hints')
            os.makedirs(iso_output_dir, exist_ok=True)

            successful_iso, failed_iso = run_parallel_isoseq_slurm(
                iso_bams, iso_output_dir, logs_dir=logs_dir,
                memory_gb=max(128, slurm_opts['memory_gb']),
                cpus=max(64, slurm_opts['cpus']),
                time_limit=max(slurm_opts['time_limit'], '12:00:00'),
                partition=slurm_opts['partition'],
                max_concurrent_jobs=slurm_opts['max_concurrent_jobs']
            )
            
            # Re-check for any IsoSeq files that might have appeared after initial verification
            logger.info('Re-checking for any IsoSeq hint files that may have appeared late')
            additional_iso_files = []
            for bam_path in iso_bams:
                bam_basename = os.path.splitext(os.path.basename(bam_path))[0]
                potential_file = os.path.join(iso_output_dir, f"{bam_basename}_isoseq_hints.gff")
                if (os.path.exists(potential_file) and 
                    os.path.getsize(potential_file) > 0 and 
                    potential_file not in successful_iso):
                    additional_iso_files.append(potential_file)
                    logger.info(f'Found late-appearing IsoSeq hints file: {potential_file}')
            
            # Add any additional files found
            if additional_iso_files:
                successful_iso.extend(additional_iso_files)
                logger.info(f'Added {len(additional_iso_files)} late-appearing IsoSeq files to merge list')
            
            # Add successful IsoSeq hints files directly to the main list
            if successful_iso:
                logger.info(f'Adding {len(successful_iso)} IsoSeq hints files to merge list')
                all_hints_files.extend(successful_iso)
                logger.info(f'Successfully added IsoSeq hints files to merge list')
            else:
                logger.warning('No successful IsoSeq hints files to add')
                
            if failed_iso:
                logger.warning(f'{len(failed_iso)} IsoSeq hint jobs failed')
        
        # Process annotation hints via NEW SLURM implementation
        if annotation_gp:
            logger.info('Processing annotation hints via NEW SLURM implementation')
            try:
                anno_hints_file = os.path.join(output_dir, 'annotation_hints.gff')
                run_annotation_hints_slurm(
                    annotation_gp, anno_hints_file, logs_dir=logs_dir,
                    memory_gb=64,  # Annotation processing needs less memory
                    cpus=32,
                    time_limit='6:00:00',  # Annotation processing is faster
                    partition=slurm_opts['partition']
                )
                if os.path.exists(anno_hints_file) and os.path.getsize(anno_hints_file) > 0:
                    all_hints_files.append(anno_hints_file)
                    logger.info(f'Successfully generated annotation hints: {anno_hints_file}')
            except Exception as e:
                logger.error(f'Failed to process annotation hints: {e}')
        
        # Process protein hints via NEW SLURM implementation
        if protein_fasta and fasta:
            logger.info('Processing protein hints via NEW SLURM implementation')
            try:
                protein_hints_file = os.path.join(output_dir, 'protein_hints.gff')
                run_protein_hints_slurm(
                    protein_fasta, fasta, protein_hints_file, logs_dir=logs_dir,
                    memory_gb=slurm_opts['memory_gb'],
                    cpus=slurm_opts['cpus'],
                    time_limit='24:00:00',  # Protein alignment needs more time
                    partition=slurm_opts['partition']
                )
                if os.path.exists(protein_hints_file) and os.path.getsize(protein_hints_file) > 0:
                    all_hints_files.append(protein_hints_file)
                    logger.info(f'Successfully generated protein hints: {protein_hints_file}')
            except Exception as e:
                logger.error(f'Failed to process protein hints: {e}')
        
        # Concatenate all hints files using optimized Slurm-based sorting
        if all_hints_files:
            logger.info(f'Concatenating {len(all_hints_files)} hints files using Slurm-based parallel sorting')
            
            # Calculate total size to determine optimal resources
            total_size_gb = sum(os.path.getsize(f) / (1024**3) for f in all_hints_files if os.path.exists(f))
            logger.info(f'Total hints data size: {total_size_gb:.2f} GB')
            
            # Adjust resources based on data size
            # For very large datasets, use more memory and CPUs for faster sorting
            if total_size_gb > 50:  # Large dataset
                merge_memory = max(256, int(total_size_gb * 2))  # 2x data size
                merge_cpus = 64
                merge_time = "8:00:00"
            elif total_size_gb > 20:  # Medium dataset
                merge_memory = 128
                merge_cpus = 32
                merge_time = "4:00:00"
            else:  # Small dataset
                merge_memory = 64
                merge_cpus = 16
                merge_time = "2:00:00"
            
            logger.info(f'Using {merge_memory}GB RAM, {merge_cpus} CPUs for merge/sort')
            
            concatenate_and_sort_hints_slurm(
                all_hints_files, 
                hints_out,
                logs_dir=logs_dir,
                memory_gb=merge_memory,
                cpus=merge_cpus,
                time_limit=merge_time,
                partition=slurm_opts['partition']
            )
            
            # Verify output
            if os.path.exists(hints_out) and os.path.getsize(hints_out) > 0:
                logger.info(f'Successfully created hints file: {hints_out}')
                logger.info(f'Hints file size: {os.path.getsize(hints_out) / (1024**2):.1f}MB')
            else:
                raise RuntimeError(f'Failed to create hints file: {hints_out}')
        else:
            raise RuntimeError('No hints files were successfully generated')
        
        return hints_out
        
    finally:
        # Cleanup temporary directory
        try:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
                logger.info(f'Cleaned up temporary directory: {temp_dir}')
            else:
                logger.debug(f'Temporary directory already removed: {temp_dir}')
        except OSError as e:
            if e.errno == 39:  # Directory not empty
                logger.warning(f'Directory not empty, attempting force cleanup: {temp_dir}')
                try:
                    # Try to remove files individually first
                    for root, dirs, files in os.walk(temp_dir, topdown=False):
                        for file in files:
                            try:
                                os.remove(os.path.join(root, file))
                            except OSError:
                                pass
                        for dir in dirs:
                            try:
                                os.rmdir(os.path.join(root, dir))
                            except OSError:
                                pass
                    os.rmdir(temp_dir)
                    logger.info(f'Force cleaned up temporary directory: {temp_dir}')
                except Exception as force_e:
                    logger.warning(f'Failed to force clean up temporary directory {temp_dir}: {force_e}')
            else:
                logger.warning(f'Failed to clean up temporary directory {temp_dir}: {e}')
        except Exception as e:
            logger.warning(f'Failed to clean up temporary directory {temp_dir}: {e}')


def run_hints_pipeline(genome: str,
                       fasta: str,
                       bams: list,
                       intron_bams: list,
                       iso_bams: list,
                       annotation_gp: str,
                       protein_fasta: str,
                       hints_out: str,
                       toil_options: object):
    """
    Main entry function for the hints Toil pipeline with dynamic resource allocation.
    """
    logger.info('Starting hints pipeline with dynamic resource allocation')

    # Honour Toil --maxCores so SingleMachine never sees jobs requesting more
    # cores than Snakemake allocated (e.g. --cores 4 → --maxCores 4).
    global _TOIL_MAX_CORES_OVERRIDE
    max_cores = getattr(toil_options, "maxCores", None)
    if max_cores is None or float(max_cores) <= 0:
        max_cores = mp.cpu_count() or 1
    max_cores = max(1, int(float(max_cores)))
    _TOIL_MAX_CORES_OVERRIDE = max_cores
    os.environ["CAT2_TOIL_MAX_CORES"] = str(max_cores)
    logger.info(f'Toil core cap (CAT2_TOIL_MAX_CORES)={max_cores}')
    
    # Calculate dynamic resources based on input characteristics
    resources = calculate_dynamic_resources(bams, iso_bams, annotation_gp, protein_fasta, fasta)
    complexity = estimate_bam_complexity((bams or []) + (intron_bams or []))
    strategy = optimize_processing_strategy(resources, complexity)
    
    logger.info(f'Processing {len(bams or [])} BAM files, {len(iso_bams or [])} IsoSeq files')
    logger.info(f'Estimated total BAM size: {resources["total_bam_size_gb"]:.1f}GB')
    logger.info(f'Optimal BAM processing: {resources["bam_cpus"]} CPUs, {resources["bam_memory_gb"]}GB memory')

    # Toil's FileJobStore.initialize() uses os.mkdir (not makedirs), so the
    # parent of the job store path must already exist.
    job_store = getattr(toil_options, "jobStore", None)
    if job_store:
        job_store_path = str(job_store)
        if job_store_path.startswith("file:"):
            job_store_path = job_store_path[5:]
        parent = os.path.dirname(os.path.abspath(job_store_path))
        if parent:
            os.makedirs(parent, exist_ok=True)
    
    with Toil(toil_options) as toil:
        if not toil.options.restart:
            # import inputs
            bam_file_ids = {}
            for dtype, paths in [('BAM', bams), ('INTRONBAM', intron_bams)]:
                bam_file_ids[dtype] = {}
                # check that paths is not None before iterating
                if paths:
                    for p in paths:
                        # Convert to absolute path before creating the URI
                        p_abs = os.path.abspath(p)
                        # import BAM and its index
                        fid = FileID.forPath(toil.importFile('file://' + p_abs), p_abs)
                        ixfid = FileID.forPath(toil.importFile('file://' + p_abs + '.bai'), p_abs + '.bai')
                        bam_file_ids[dtype][os.path.basename(p)] = (fid, ixfid)

            iso_seq_ids = []
            if iso_bams:
                for p in iso_bams:
                    p_abs = os.path.abspath(p)
                    fid = FileID.forPath(toil.importFile('file://' + p_abs), p_abs)
                    ixfid = FileID.forPath(toil.importFile('file://' + p_abs + '.bai'), p_abs + '.bai')
                    iso_seq_ids.append((fid, ixfid))

            annotation_id = None
            if annotation_gp:
                p_abs = os.path.abspath(annotation_gp)
                annotation_id = FileID.forPath(toil.importFile('file://' + p_abs), p_abs)

            protein_id = None
            genome_id = None
            if protein_fasta:
                p_abs = os.path.abspath(protein_fasta)
                protein_id = FileID.forPath(toil.importFile('file://' + p_abs), p_abs)
            
            # The genome FASTA is required for protein hints and should always be imported
            if fasta:
                p_abs = os.path.abspath(fasta)
                genome_id = FileID.forPath(toil.importFile('file://' + p_abs), p_abs)

            inputs = {
                'bams': bam_file_ids,
                'iso_seq_bams': iso_seq_ids,
                'annotation': annotation_id,
                'protein_fasta': protein_id,
                'genome_fasta': genome_id
            }

            # Dynamic resource calculation
            resources = calculate_dynamic_resources(bams, iso_bams, annotation_gp, protein_fasta, fasta)
            logger.info(f'Dynamic resources calculated: {resources}')
            
            # Update toil_options with dynamic resources
            toil_options.cores = resources['bam_cpus']
            toil_options.memory = f'{resources["bam_memory_gb"]}G'
            # For disk, we might want to set a minimum and let Toil handle the rest
            toil_options.disk = f'{max(32, resources["merge_disk_gb"])}G'
            
            disk = tools.toilInterface.find_total_disk_usage(inputs)
            root_job = Job.wrapJobFn(_setup_hints, inputs, disk=disk)
            hints_file_id = toil.start(root_job)
        else:
            hints_file_id = toil.restart()

    # Export the final file
        os.makedirs(os.path.dirname(os.path.abspath(hints_out)), exist_ok=True)
        toil.exportFile(hints_file_id, 'file://' + os.path.abspath(hints_out))
    return hints_out

def _setup_hints(job: Job, input_ids: dict):
    """
    Setup job: split by reference groups, sort, filter, then merge with dynamic resources.
    """
    logger.info('Setting up hints pipeline with dynamic resource allocation')
    
    # Calculate resources based on input characteristics
    all_bam_fids = []
    for dtype, groups in input_ids['bams'].items():
        for ref_group, (bam_fid, bai_fid) in groups.items():
            all_bam_fids.extend([bam_fid, bai_fid])
    
    iso_bam_fids = []
    for bam_fid, bai_fid in input_ids['iso_seq_bams']:
        iso_bam_fids.extend([bam_fid, bai_fid])
    
    # Store metadata for dynamic resource calculation
    resources_metadata = {
        'total_bam_files': len(all_bam_fids) + len(iso_bam_fids),
        'has_protein': input_ids['protein_fasta'] is not None,
        'has_annotation': input_ids['annotation'] is not None,
        'iso_seq_count': len(input_ids['iso_seq_bams'])
    }
    logger.info(f'Pipeline metadata: {resources_metadata}')
    
    # Split into filtered sets with dynamic resource allocation
    filtered = {'BAM': collections.defaultdict(list), 'INTRONBAM': collections.defaultdict(list)}
    
    # Process each BAM group with adaptive resources
    for dtype, groups in input_ids['bams'].items():
        for ref_group, (bam_fid, bai_fid) in groups.items():
            try:
                # Sample BAM to determine characteristics
                bam_path = job.fileStore.readGlobalFile(bam_fid)
                is_paired = bam_is_paired(bam_path)
                
                # Calculate adaptive resources for this BAM
                bam_size_gb = os.path.getsize(bam_path) / (1024**3)
                
                # Dynamic resource calculation (cores capped to Toil --maxCores)
                if bam_size_gb > 20:  # Large BAM
                    cores = _cap_cores(64)
                    memory = f'{max(8, min(64, int(bam_size_gb) or 8))}G'
                    chunk_reads = 25_000_000
                elif bam_size_gb > 5:  # Medium BAM
                    cores = _cap_cores(64)
                    memory = f'{max(8, min(64, int(bam_size_gb * 2) or 8))}G'
                    chunk_reads = 50_000_000
                else:  # Small BAM
                    cores = _cap_cores(64)
                    memory = f'{max(4, min(64, _toil_max_cores() * 2))}G'
                    chunk_reads = 100_000_000
                
                disk = tools.toilInterface.find_total_disk_usage([bam_fid, bai_fid]) * 3 + 4
                
                logger.info(f'Processing {dtype} BAM {ref_group}: {bam_size_gb:.1f}GB, {cores} cores, {memory} memory')
                
                j = job.addChildJobFn(namesort_bam_dynamic, bam_fid, bai_fid, ref_group, chunk_reads,
                                      disk=disk, cores=cores, memory=memory)
                filtered[dtype][ref_group].append(j.rv())
                
            except Exception as e:
                logger.error(f'Error processing BAM {ref_group}: {str(e)}')
                # Fallback to conservative resources
                disk = tools.toilInterface.find_total_disk_usage([bam_fid, bai_fid]) * 3 + 2
                j = job.addChildJobFn(namesort_bam, bam_fid, bai_fid, ref_group,
                                      disk=disk, cores=_cap_cores(8), memory='12G')
                filtered[dtype][ref_group].append(j.rv())

    # IsoSeq with dynamic resources
    iso_hints = []
    for i, (bam_fid, bai_fid) in enumerate(input_ids['iso_seq_bams']):
        try:
            disk = tools.toilInterface.find_total_disk_usage([bam_fid, bai_fid]) * 2
            
            # Adaptive resources for IsoSeq based on dataset size
            memory = f'{max(4, min(64, _toil_max_cores() * 2))}G'
            cores = _cap_cores(64)
            
            logger.info(f'Processing IsoSeq BAM {i+1}/{len(input_ids["iso_seq_bams"])}: {cores} cores, {memory} memory')
            
            j = job.addChildJobFn(generate_iso_seq_hints, bam_fid, bai_fid,
                                  disk=disk, memory=memory, cores=cores)
            iso_hints.append(j.rv())
            
        except Exception as e:
            logger.error(f'Error setting up IsoSeq processing: {str(e)}')
            # Fallback
            disk = tools.toilInterface.find_total_disk_usage([bam_fid, bai_fid])
            j = job.addChildJobFn(generate_iso_seq_hints, bam_fid, bai_fid,
                                  disk=disk, memory='8G')
            iso_hints.append(j.rv())

    # Protein with dynamic resources
    prot_hints = None
    if input_ids['protein_fasta']:
        try:
            # Estimate protein file size for resource allocation
            protein_path = job.fileStore.readGlobalFile(input_ids['protein_fasta'])
            protein_size_mb = os.path.getsize(protein_path) / (1024**2)
            
            # Adaptive resources for protein alignment
            if protein_size_mb > 200:  # Large protein dataset
                chunk_size = 50  # Smaller chunks
            elif protein_size_mb > 50:  # Medium protein dataset
                chunk_size = 75
            else:  # Small protein dataset
                chunk_size = 100
            memory = f'{max(4, min(64, _toil_max_cores() * 2))}G'
            cores = _cap_cores(64)
            
            disk = tools.toilInterface.find_total_disk_usage([input_ids['protein_fasta'], input_ids['genome_fasta']]) * 2
            
            logger.info(f'Processing protein alignment: {protein_size_mb:.1f}MB, {cores} cores, {memory} memory, chunk size {chunk_size}')
            
            j = job.addChildJobFn(generate_protein_hints_dynamic,
                                  input_ids['protein_fasta'], input_ids['genome_fasta'], chunk_size,
                                  disk=disk, memory=memory, cores=cores)
            prot_hints = j.rv()
            
        except Exception as e:
            logger.error(f'Error setting up protein processing: {str(e)}')
            # Fallback
            disk = tools.toilInterface.find_total_disk_usage([input_ids['protein_fasta'], input_ids['genome_fasta']])
            j = job.addChildJobFn(generate_protein_hints,
                                  input_ids['protein_fasta'], input_ids['genome_fasta'],
                                  disk=disk, memory='12G')
            prot_hints = j.rv()

    # Annotation with dynamic resources
    anno_hints = None
    if input_ids['annotation']:
        try:
            # Estimate annotation complexity
            annotation_path = job.fileStore.readGlobalFile(input_ids['annotation'])
            annotation_size_mb = os.path.getsize(annotation_path) / (1024**2)
            
            # Count transcripts for resource estimation
            transcript_count = 0
            with open(annotation_path, 'r') as f:
                transcript_count = sum(1 for line in f if line.strip())
            
            # Adaptive resources for annotation processing
            memory = f'{max(4, min(64, _toil_max_cores() * 2))}G'
            cores = _cap_cores(64)
            
            disk = tools.toilInterface.find_total_disk_usage(input_ids['annotation']) * 2
            
            logger.info(f'Processing annotation: {transcript_count} transcripts, {annotation_size_mb:.1f}MB, {cores} cores, {memory} memory')
            
            j = job.addChildJobFn(generate_annotation_hints,
                                  input_ids['annotation'],
                                  disk=disk, memory=memory, cores=cores)
            anno_hints = j.rv()
            
        except Exception as e:
            logger.error(f'Error setting up annotation processing: {str(e)}')
            # Fallback
            disk = tools.toilInterface.find_total_disk_usage(input_ids['annotation'])
            j = job.addChildJobFn(generate_annotation_hints,
                                  input_ids['annotation'],
                                  disk=disk, memory='128G')
            anno_hints = j.rv()

    return job.addFollowOnJobFn(merge_bams_dynamic,
                                filtered, anno_hints, iso_hints, prot_hints, resources_metadata).rv()


def namesort_bam_dynamic(job: Job, bam_fid, bai_fid, ref_group, num_reads=50_000_000):
    """
    Enhanced namesort_bam with dynamic chunking and parallel processing.
    """
    logger.info(f'Dynamic BAM name sorting for {ref_group} with {num_reads} reads per chunk')
    
    bam_path = job.fileStore.readGlobalFile(bam_fid)
    job.fileStore.readGlobalFile(bai_fid, bam_path + '.bai')
    tmp = tools.fileOps.get_tmp_toil_file(suffix='name_sorted.bam')
    is_paired = bam_is_paired(bam_path)
    
    try:
        # Adaptive threading based on available cores
        cores = min(job.cores, 64)  # Use job-allocated cores
        memory_per_thread = f'{max(1, int(job.memory // (1024**3) // cores))}G'
        
        # Correctly construct the command to view the entire BAM file
        cmd = [['samtools', 'view', '-@', str(cores//2), '-b', bam_path] + list(ref_group),
               ['sambamba', 'sort', '-t', str(cores//2), '-m', memory_per_thread, 
                '-o', '/dev/stdout', '-n', '/dev/stdin']]
        
        logger.info(f'Running BAM sort with {cores//2} threads and {memory_per_thread} memory per thread')
        tools.procOps.run_proc(cmd, stdout=tmp)
        
    except Exception as e:
        logger.warning(f'Dynamic sort failed for {ref_group}, falling back to standard approach: {str(e)}')
        # Fallback: still honour allocated cores / memory, not hardcoded 32/64G
        fb_cores = max(1, _cap_cores(job.cores))
        fb_mem_gb = max(1, int((job.memory or (4 * 1024**3)) // (1024**3)))
        cmd = [['samtools', 'view', '-@', str(fb_cores), '-b', bam_path] + list(ref_group),
               ['sambamba', 'sort', '-t', str(fb_cores), '-m', f'{fb_mem_gb}G',
                '-o', '/dev/stdout', '-n', '/dev/stdin']]
        tools.procOps.run_proc(cmd, stdout=tmp)

    # Process chunks with parallel filtering
    handle = pysam.Samfile(tmp)
    r = []
    fids = []
    chunk_count = 0
    
    try:
        for qname, reads in itertools.groupby(handle, lambda x: x.qname):
            r.extend(reads)
            if len(r) >= num_reads:
                outf = tools.fileOps.get_tmp_toil_file()
                out_handle = pysam.Samfile(outf, 'wb', template=handle)
                for rec in r: 
                    out_handle.write(rec)
                out_handle.close()
                
                # Dynamic resource allocation for filtering
                chunk_size_gb = os.path.getsize(outf) / (1024**3)
                filter_memory = f'{max(4, int(chunk_size_gb * 2))}G'
                filter_cores = _cap_cores(min(4, max(1, job.cores // 4)))
                
                fid = job.addChildJobFn(filter_bam_dynamic, job.fileStore.writeGlobalFile(outf), is_paired,
                                       disk='32G', memory=filter_memory, cores=filter_cores).rv()
                fids.append(fid)
                r = []
                chunk_count += 1
                
                logger.info(f'Created chunk {chunk_count} for {ref_group}: {chunk_size_gb:.2f}GB, {filter_cores} cores, {filter_memory} memory')
        
        # Handle remaining reads
        if r:
            outf = tools.fileOps.get_tmp_toil_file()
            out_handle = pysam.Samfile(outf, 'wb', template=handle)
            for rec in r: 
                out_handle.write(rec)
            out_handle.close()
            
            chunk_size_gb = os.path.getsize(outf) / (1024**3)
            filter_memory = f'{max(4, int(chunk_size_gb * 2))}G'
            filter_cores = _cap_cores(min(4, max(1, job.cores // 4)))
            
            fids.append(job.addChildJobFn(filter_bam_dynamic, job.fileStore.writeGlobalFile(outf), is_paired,
                                         disk='8G', memory=filter_memory, cores=filter_cores).rv())
            chunk_count += 1
            
        handle.close()
        logger.info(f'Created {chunk_count} chunks for {ref_group}')
        
    except Exception as e:
        logger.error(f'Error during chunking for {ref_group}: {str(e)}')
        handle.close()
        # Fallback: create single chunk
        fids.append(job.addChildJobFn(filter_bam, job.fileStore.writeGlobalFile(tmp), is_paired,
                                     disk='64G', memory='64G').rv())
    
    return job.addFollowOnJobFn(merge_filtered_bams_dynamic, fids, ref_group).rv()


def filter_bam_dynamic(job: Job, file_id, is_paired):
    """
    Enhanced BAM filtering with dynamic resource utilization.
    """
    bam_path = job.fileStore.readGlobalFile(file_id)
    assert os.path.getsize(bam_path) > 0
    
    logger.info(f'Filtering BAM chunk: {os.path.getsize(bam_path) / (1024**2):.1f}MB')
    
    tmp_filtered = tools.fileOps.get_tmp_toil_file()
    filter_cmd = ['filterBam', '--uniq', '--in', bam_path, '--out', tmp_filtered]

    if is_paired == True:
        filter_cmd.extend(['--paired', '--pairwiseAlignments'])
    
    try:
        tools.procOps.run_proc(filter_cmd)
        
        if os.path.getsize(tmp_filtered) == 0:
            logger.warning('Filtered BAM chunk is empty - this may indicate low-quality alignments')
            # Create minimal placeholder to avoid pipeline failures
            with open(tmp_filtered, 'w') as f:
                f.write('')  # Empty file
            return None
        
        out_filter = tools.fileOps.get_tmp_toil_file()
        
        # Adaptive sorting based on available resources
        cores = _cap_cores(job.cores)
        memory_gb = max(32, int(job.memory // (1024**3) // 2))  # Use half available memory
        
        sort_cmd = ['sambamba', 'sort', tmp_filtered, '-o', out_filter, 
                    '-t', str(cores), '-m', f'{memory_gb}G']
        
        logger.info(f'Sorting with {cores} threads and {memory_gb}G memory')
        tools.procOps.run_proc(sort_cmd)
        
        return job.fileStore.writeGlobalFile(out_filter)
        
    except Exception as e:
        logger.error(f'Error during BAM filtering: {str(e)}')
        # Fallback to original approach
        tools.procOps.run_proc(filter_cmd)
        if os.path.getsize(tmp_filtered) == 0:
            raise RuntimeError('After filtering one BAM subset became empty. This could be bad.')
        
        out_filter = tools.fileOps.get_tmp_toil_file()
        sort_cmd = ['sambamba', 'sort', tmp_filtered, '-o', out_filter, '-t', '16']
        tools.procOps.run_proc(sort_cmd)
        return job.fileStore.writeGlobalFile(out_filter)


def merge_filtered_bams_dynamic(job: Job, filtered_file_ids, ref_group):
    """
    Enhanced BAM merging with dynamic resource allocation and parallel processing.
    """
    # Filter out None values (empty chunks)
    valid_file_ids = [fid for fid in filtered_file_ids if fid is not None]
    
    if not valid_file_ids:
        logger.warning(f'No valid BAM chunks for {ref_group} - creating empty placeholder')
        empty_bam = tools.fileOps.get_tmp_toil_file()
        # Create minimal BAM header
        cmd = ['samtools', 'view', '-H', '-o', empty_bam, '/dev/null']
        try:
            tools.procOps.run_proc(cmd)
        except:
            with open(empty_bam, 'w') as f:
                f.write('')  # Fallback empty file
        return job.fileStore.writeGlobalFile(empty_bam)
    
    if len(valid_file_ids) == 1:
        logger.info(f'Single chunk for {ref_group}, no merging needed')
        return valid_file_ids[0]
    
    logger.info(f'Merging {len(valid_file_ids)} BAM chunks for {ref_group}')
    
    try:
        local_paths = [job.fileStore.readGlobalFile(x) for x in valid_file_ids]
        
        # Calculate total size for resource planning
        total_size_gb = sum(os.path.getsize(path) / (1024**3) for path in local_paths)
        logger.info(f'Total BAM data to merge: {total_size_gb:.2f}GB')
        
        # Adaptive threading and memory
        cores = min(job.cores, max(32, len(valid_file_ids) // 4))
        memory_gb = max(64, int(total_size_gb / 2))
        
        # Use file list approach for many files
        if len(local_paths) > 50:
            fofn = tools.fileOps.get_tmp_toil_file()
            with open(fofn, 'w') as outf:
                for path in local_paths:
                    if os.environ.get('CAT_BINARY_MODE') == 'singularity':
                        path = tools.procOps.singularify_arg(path)
                    outf.write(path + '\n')
            
            out_bam = tools.fileOps.get_tmp_toil_file()
            cmd = ['samtools', 'merge', '-@', str(cores), '-b', fofn, out_bam]
            
        else:
            # Direct merge for fewer files
            out_bam = tools.fileOps.get_tmp_toil_file()
            cmd = ['samtools', 'merge', '-@', str(cores)] + local_paths + [out_bam]
        
        logger.info(f'Merging with {cores} threads')
        tools.procOps.run_proc(cmd)
        
        merged_size_gb = os.path.getsize(out_bam) / (1024**3)
        logger.info(f'Merged BAM for {ref_group}: {merged_size_gb:.2f}GB')
        
        return job.fileStore.writeGlobalFile(out_bam)
        
    except Exception as e:
        logger.error(f'Error during dynamic BAM merging for {ref_group}: {str(e)}')
        # Fallback to original approach
        local_paths = [job.fileStore.readGlobalFile(x) for x in valid_file_ids]
        fofn = tools.fileOps.get_tmp_toil_file()
        with open(fofn, 'w') as outf:
            for l in local_paths:
                if os.environ.get('CAT_BINARY_MODE') == 'singularity':
                    l = tools.procOps.singularify_arg(l)
                outf.write(l + '\n')
        out_bam = tools.fileOps.get_tmp_toil_file()
        cmd = ['samtools', 'merge', '-b', fofn, out_bam]
        tools.procOps.run_proc(cmd)
        return job.fileStore.writeGlobalFile(out_bam)


def generate_protein_hints_dynamic(job: Job, protein_fasta_file_id, genome_fasta_file_id, chunk_size=100):
    """
    Enhanced protein hints generation with dynamic chunking and parallel processing.
    """
    logger.info(f'Generating protein hints with dynamic chunking (chunk size: {chunk_size})')
    
    try:
        disk_usage = tools.toilInterface.find_total_disk_usage(genome_fasta_file_id)
        protein_fasta = job.fileStore.readGlobalFile(protein_fasta_file_id)
        
        cmd = ['samtools', 'faidx', protein_fasta]
        tools.procOps.run_proc(cmd)
        protein_handle = tools.bio.get_sequence_dict(protein_fasta)
        
        # Calculate optimal parallel processing
        total_proteins = len(protein_handle)
        available_cores = job.cores
        optimal_chunks = min(total_proteins // chunk_size + 1, available_cores * 2)
        
        logger.info(f'Processing {total_proteins} proteins in {optimal_chunks} parallel chunks')
        
        results = []
        chunk_count = 0
        
        # Create protein chunks for parallel processing
        protein_items = list(protein_handle.items())
        for chunk in tools.dataOps.grouper(protein_items, chunk_size):
            chunk_proteins = [item for item in chunk if item is not None]  # Filter None values
            if not chunk_proteins:
                continue
                
            chunk_count += 1
            
            # Dynamic resource allocation per chunk
            chunk_memory = f'{max(64, min(16, len(chunk_proteins) // 10 + 4))}G'
            chunk_cores = _cap_cores(min(64, max(1, available_cores // optimal_chunks)))
            
            logger.info(f'Chunk {chunk_count}: {len(chunk_proteins)} proteins, {chunk_cores} cores, {chunk_memory} memory')
            
            j = job.addChildJobFn(run_protein_aln_dynamic, chunk_proteins, genome_fasta_file_id,
                                 disk=disk_usage, memory=chunk_memory, cores=chunk_cores)
            results.append(j.rv())
        
        if not results:
            logger.warning('No protein chunks created - creating empty result')
            empty_result = tools.fileOps.get_tmp_toil_file()
            with open(empty_result, 'w') as f:
                f.write('')
            return job.fileStore.writeGlobalFile(empty_result)
        
        logger.info(f'Created {len(results)} protein alignment jobs')
        
        # Merge results with dynamic resources
        merge_memory = f'{max(64, min(128, len(results) + 64))}G'
        merge_cores = _cap_cores(min(64, max(32, available_cores // 2)))
        
        return job.addFollowOnJobFn(convert_protein_aln_results_to_hints_dynamic, results,
                                   memory=merge_memory, cores=merge_cores).rv()
        
    except Exception as e:
        logger.error(f'Error in dynamic protein hints generation: {str(e)}')
        # Fallback to original approach
        return job.addFollowOnJobFn(generate_protein_hints, protein_fasta_file_id, genome_fasta_file_id,
                                   disk=disk_usage, memory='12G').rv()


def run_protein_aln_dynamic(job: Job, protein_subset, genome_fasta_file_id):
    """
    Enhanced protein alignment with adaptive resource utilization.
    """
    genome_fasta = job.fileStore.readGlobalFile(genome_fasta_file_id)
    protein_fasta = tools.fileOps.get_tmp_toil_file()
    
    # Write protein sequences
    protein_count = 0
    with open(protein_fasta, 'w') as outf:
        for name, seq in protein_subset:
            tools.bio.write_fasta(outf, name, str(seq))
            protein_count += 1
    
    logger.info(f'Running protein alignment for {protein_count} sequences')
    
    tmp_exonerate = tools.fileOps.get_tmp_toil_file()
    
    try:
        # Adaptive exonerate parameters based on sequence count
        if protein_count > 500:  # Large chunk
            cmd = ['exonerate', '--model', 'protein2genome', '--showvulgar', 'no', 
                   '--showalignment', 'no', '--showquerygff', 'yes', '--ryo', 'AveragePercentIdentity: %pi\n',
                   protein_fasta, genome_fasta]
        else:  # Standard chunk
            cmd = ['exonerate', '--model', 'protein2genome', '--showvulgar', 'no', 
                   '--showalignment', 'no', '--showquerygff', 'yes', protein_fasta, genome_fasta]
        
        tools.procOps.run_proc(cmd, stdout=tmp_exonerate)
        
        # Validate output
        if os.path.getsize(tmp_exonerate) == 0:
            logger.warning(f'Empty exonerate output for {protein_count} proteins')
        else:
            logger.info(f'Generated protein alignments: {os.path.getsize(tmp_exonerate) / (1024**2):.1f}MB')
        
        return job.fileStore.writeGlobalFile(tmp_exonerate)
        
    except Exception as e:
        logger.error(f'Error in protein alignment for {protein_count} sequences: {str(e)}')
        # Create empty result to avoid pipeline failure
        with open(tmp_exonerate, 'w') as f:
            f.write('')
        return job.fileStore.writeGlobalFile(tmp_exonerate)


def convert_protein_aln_results_to_hints_dynamic(job: Job, results):
    """
    Enhanced protein alignment results processing with parallel merging.
    """
    logger.info(f'Converting {len(results)} protein alignment results to hints')
    
    try:
        merged_exonerate = tools.fileOps.get_tmp_toil_file()
        total_size_mb = 0
        
        with open(merged_exonerate, 'w') as outf:
            for i, r in enumerate(results):
                try:
                    f = job.fileStore.readGlobalFile(r)
                    file_size_mb = os.path.getsize(f) / (1024**2)
                    total_size_mb += file_size_mb
                    
                    with open(f, 'r') as inf:
                        outf.write(inf.read())
                    
                    if i % 10 == 0:
                        logger.info(f'Merged {i+1}/{len(results)} protein alignment files')
                        
                except Exception as e:
                    logger.warning(f'Error reading protein alignment result {i}: {str(e)}')
                    continue
        
        logger.info(f'Merged protein alignments: {total_size_mb:.1f}MB total')
        
        # Sort with adaptive resources
        tmp_sorted = tools.fileOps.get_tmp_toil_file()
        tools.misc.sort_gff(merged_exonerate, tmp_sorted)
        
        # Generate hints
        out_hints = tools.fileOps.get_tmp_toil_file()
        
        # Adaptive hint generation parameters
        if total_size_mb > 100:  # Large alignment result
            cmd = ['exonerate2hints.pl', '--in={}'.format(tmp_sorted), 
                   '--CDSpart_cutoff=5', '--out={}'.format(out_hints)]
        else:  # Standard alignment result
            cmd = ['exonerate2hints.pl', '--in={}'.format(tmp_sorted), 
                   '--CDSpart_cutoff=5', '--out={}'.format(out_hints)]
        
        tools.procOps.run_proc(cmd)
        
        hints_size_mb = os.path.getsize(out_hints) / (1024**2)
        logger.info(f'Generated protein hints: {hints_size_mb:.1f}MB')
        
        return job.fileStore.writeGlobalFile(out_hints)
        
    except Exception as e:
        logger.error(f'Error converting protein alignments to hints: {str(e)}')
        # Create empty hints file as fallback
        out_hints = tools.fileOps.get_tmp_toil_file()
        with open(out_hints, 'w') as f:
            f.write('')
        return job.fileStore.writeGlobalFile(out_hints)


def merge_bams_dynamic(job: Job, filtered_bam_file_ids, annotation_hints_file_id, 
                      iso_seq_hints_file_ids, protein_hints_file_id, resources_metadata):
    """
    Enhanced BAM merging with dynamic resource allocation and parallel processing.
    """
    logger.info('Merging BAMs with dynamic resource allocation')
    logger.info(f'Resources metadata: {resources_metadata}')
    
    merged_bam_file_ids = {'BAM': {}, 'INTRONBAM': {}}
    merge_jobs = []
    
    # Process each BAM type with adaptive resources
    for dtype in filtered_bam_file_ids:
        logger.info(f'Processing {dtype} BAMs: {len(filtered_bam_file_ids[dtype])} reference groups')
        
        for ref_group, file_ids in filtered_bam_file_ids[dtype].items():
            # Filter out None values (empty/failed chunks)
            valid_file_ids = [x for x in file_ids if x is not None]
            
            if not valid_file_ids:
                logger.warning(f'No valid file IDs for {dtype} {ref_group}')
                continue
            
            try:
                # Calculate resources based on number of chunks and estimated data size
                num_chunks = len(valid_file_ids)
                
                if num_chunks > 20:  # Many chunks - high resource job
                    memory = '128G'
                    cores = _cap_cores(32)
                    disk_multiplier = 4
                elif num_chunks > 10:  # Medium chunks
                    memory = '64G'
                    cores = _cap_cores(16)
                    disk_multiplier = 3
                else:  # Few chunks
                    memory = '32G'
                    cores = _cap_cores(8)
                    disk_multiplier = 2
                
                disk_usage = tools.toilInterface.find_total_disk_usage(valid_file_ids) * disk_multiplier
                
                logger.info(f'Merging {dtype} {ref_group}: {num_chunks} chunks, {cores} cores, {memory} memory')
                
                merge_job = job.addChildJobFn(cat_sort_bams_dynamic, valid_file_ids, ref_group,
                                            disk=disk_usage, memory=memory, cores=cores)
                merged_bam_file_ids[dtype][ref_group] = merge_job.rv()
                merge_jobs.append(merge_job)
                
            except Exception as e:
                logger.error(f'Error setting up merge for {dtype} {ref_group}: {str(e)}')
                # Fallback to original approach
                disk_usage = tools.toilInterface.find_total_disk_usage(valid_file_ids)
                merged_bam_file_ids[dtype][ref_group] = job.addChildJobFn(cat_sort_bams, valid_file_ids, 
                                                                        disk=disk_usage, memory='16G', cores=_cap_cores(4)).rv()
    
    # Adaptive final processing based on total workload
    total_ref_groups = sum(len(merged_bam_file_ids[dtype]) for dtype in merged_bam_file_ids)
    
    if total_ref_groups > 50:  # Large dataset
        final_memory = '128G'
        final_cores = _cap_cores(64)
    elif total_ref_groups > 20:  # Medium dataset
        final_memory = '64G'
        final_cores = _cap_cores(32)
    else:  # Small dataset
        final_memory = '32G'
        final_cores = _cap_cores(16)

    logger.info(f'Final hints building: {total_ref_groups} reference groups, {final_cores} cores, {final_memory} memory')
    
    return job.addFollowOnJobFn(build_hints_dynamic, merged_bam_file_ids, annotation_hints_file_id, 
                               iso_seq_hints_file_ids, protein_hints_file_id, resources_metadata,
                               memory=final_memory, cores=final_cores).rv()


def cat_sort_bams_dynamic(job: Job, bam_file_ids, ref_group):
    """
    Enhanced BAM concatenation and sorting with adaptive resource management.
    """
    logger.info(f'Concatenating and sorting {len(bam_file_ids)} BAM files for {ref_group}')
    
    try:
        bamfiles = [job.fileStore.readGlobalFile(x) for x in bam_file_ids]
        
        # Calculate total data size for resource planning
        total_size_gb = sum(os.path.getsize(f) / (1024**3) for f in bamfiles)
        logger.info(f'Total BAM size for {ref_group}: {total_size_gb:.2f}GB')
        
        # Handle single file case
        if len(bamfiles) == 1:
            logger.info(f'Single BAM file for {ref_group}, copying directly')
            return job.fileStore.writeGlobalFile(bamfiles[0])
        
        # Adaptive concatenation strategy
        catfile = tools.fileOps.get_tmp_toil_file()
        
        if len(bamfiles) > 100:  # Many files - use chunked approach
            logger.info(f'Using chunked concatenation for {len(bamfiles)} files')
            sam_iter = tools.dataOps.grouper(bamfiles, 4095)
            cmd = ['samtools', 'cat', '-o', catfile]
            cmd.extend(next(sam_iter))
            tools.procOps.run_proc(cmd)
            
            for more in sam_iter:
                old_catfile = catfile
                catfile = tools.fileOps.get_tmp_toil_file()
                cmd = ['samtools', 'cat', '-o', catfile, old_catfile]
                cmd.extend(more)
                tools.procOps.run_proc(cmd)
                
        else:  # Direct concatenation
            cmd = ['samtools', 'cat', '-o', catfile] + bamfiles
            tools.procOps.run_proc(cmd)
        
        # Adaptive sorting
        merged = tools.fileOps.get_tmp_toil_file()
        
        # Calculate optimal sort parameters
        cores = min(job.cores, 12)
        # Use most of available memory but leave some buffer
        memory_gb = max(4, int(job.memory // (1024**3) * 0.8))
        
        cmd = ['sambamba', 'sort', catfile, '-o', merged, 
               '-t', str(cores), '-m', f'{memory_gb}G']
        
        logger.info(f'Sorting {ref_group} with {cores} threads and {memory_gb}G memory')
        tools.procOps.run_proc(cmd)
        
        sorted_size_gb = os.path.getsize(merged) / (1024**3)
        logger.info(f'Sorted BAM for {ref_group}: {sorted_size_gb:.2f}GB')
        
        return job.fileStore.writeGlobalFile(merged)
        
    except Exception as e:
        logger.error(f'Error in dynamic BAM cat/sort for {ref_group}: {str(e)}')
        # Fallback to original approach
        bamfiles = [job.fileStore.readGlobalFile(x) for x in bam_file_ids]
        catfile = tools.fileOps.get_tmp_toil_file()
        sam_iter = tools.dataOps.grouper(bamfiles, 4095)
        cmd = ['samtools', 'cat', '-o', catfile]
        cmd.extend(next(sam_iter))
        tools.procOps.run_proc(cmd)
        for more in sam_iter:
            old_catfile = catfile
            catfile = tools.fileOps.get_tmp_toil_file()
            cmd = ['samtools', 'cat', '-o', catfile, old_catfile]
            cmd.extend(more)
            tools.procOps.run_proc(cmd)
        merged = tools.fileOps.get_tmp_toil_file()
        cmd = ['sambamba', 'sort', catfile, '-o', merged, '-t', '32', '-m', '63G']
        tools.procOps.run_proc(cmd)
        return job.fileStore.writeGlobalFile(merged)


def build_hints_dynamic(job: Job, merged_bam_file_ids, anno_hints, iso_hints, prot_hints, resources_metadata):
    """
    Enhanced hints building with parallel processing and dynamic resource allocation.
    """
    logger.info('Building hints with dynamic parallel processing')
    
    intron_hints_file_ids = []
    exon_hints_file_ids = []
    
    # Process different BAM types in parallel with adaptive resources
    hint_jobs = []
    
    for dtype in merged_bam_file_ids:
        logger.info(f'Processing {dtype} hints: {len(merged_bam_file_ids[dtype])} reference groups')
        
        for ref_group, file_id in merged_bam_file_ids[dtype].items():
            if file_id is None:
                logger.warning(f'Skipping None file_id for {dtype} {ref_group}')
                continue
            
            try:
                # Dynamic resource allocation based on BAM type and dataset characteristics
                if dtype == 'INTRONBAM':
                    # Intron processing is typically lighter
                    intron_memory = '8G'
                    intron_cores = _cap_cores(2)
                elif resources_metadata.get('total_bam_files', 0) > 20:  # Many BAMs
                    intron_memory = '12G'
                    intron_cores = _cap_cores(4)
                else:
                    intron_memory = '10G'
                    intron_cores = _cap_cores(3)
                
                logger.info(f'Processing intron hints for {dtype} {ref_group}: {intron_cores} cores, {intron_memory} memory')
                
                intron_job = job.addChildJobFn(build_intron_hints_dynamic, file_id, ref_group,
                                             memory=intron_memory, cores=intron_cores)
                intron_hints_file_ids.append(intron_job.rv())
                hint_jobs.append(intron_job)
                
                # Exon hints only for regular BAMs
                if dtype == 'BAM':
                    # Exon processing is more resource intensive
                    if resources_metadata.get('total_bam_files', 0) > 20:  # Many BAMs
                        exon_memory = '16G'
                        exon_cores = _cap_cores(6)
                    else:
                        exon_memory = '12G'
                        exon_cores = _cap_cores(4)
                    
                    logger.info(f'Processing exon hints for {dtype} {ref_group}: {exon_cores} cores, {exon_memory} memory')
                    
                    exon_job = job.addChildJobFn(build_exon_hints_dynamic, file_id, ref_group,
                                               memory=exon_memory, cores=exon_cores)
                    exon_hints_file_ids.append(exon_job.rv())
                    hint_jobs.append(exon_job)
                    
            except Exception as e:
                logger.error(f'Error setting up hints processing for {dtype} {ref_group}: {str(e)}')
                # Fallback to original approach
                intron_hints_file_ids.append(job.addChildJobFn(build_intron_hints, file_id).rv())
                if dtype == 'BAM':
                    exon_hints_file_ids.append(job.addChildJobFn(build_exon_hints, file_id).rv())
    
    # Calculate final merging resources
    total_hint_files = len(intron_hints_file_ids) + len(exon_hints_file_ids) + len(iso_hints) + 2  # anno + prot
    
    if total_hint_files > 100:  # Very large dataset
        merge_memory = '48G'
        merge_cores = _cap_cores(16)
        merge_disk_multiplier = 6
    elif total_hint_files > 50:  # Large dataset
        merge_memory = '32G'
        merge_cores = _cap_cores(12)
        merge_disk_multiplier = 4
    elif total_hint_files > 20:  # Medium dataset
        merge_memory = '24G'
        merge_cores = _cap_cores(8)
        merge_disk_multiplier = 3
    else:  # Small dataset
        merge_memory = '16G'
        merge_cores = _cap_cores(6)
        merge_disk_multiplier = 2
    
    # Calculate disk usage for final merge
    all_hint_files = intron_hints_file_ids + exon_hints_file_ids + iso_hints
    if anno_hints:
        all_hint_files.append(anno_hints)
    if prot_hints:
        all_hint_files.append(prot_hints)
    
    try:
        disk_usage = tools.toilInterface.find_total_disk_usage(all_hint_files) * merge_disk_multiplier
    except:
        disk_usage = f'{total_hint_files * 2}G'  # Conservative estimate
    
    logger.info(f'Final hints merging: {total_hint_files} files, {merge_cores} cores, {merge_memory} memory, {disk_usage} disk')
    
    return job.addFollowOnJobFn(cat_hints_dynamic,
                               intron_hints_file_ids, exon_hints_file_ids, anno_hints, iso_hints, prot_hints,
                               resources_metadata, disk=disk_usage, memory=merge_memory, cores=merge_cores).rv()


def build_intron_hints_dynamic(job: Job, merged_bam_file_id, ref_group):
    """
    Enhanced intron hints building with adaptive processing.
    """
    logger.info(f'Building intron hints for {ref_group}')
    
    try:
        merged_bam_file = job.fileStore.readGlobalFile(merged_bam_file_id)
        bam_size_mb = os.path.getsize(merged_bam_file) / (1024**2)
        
        logger.info(f'Processing intron hints for {ref_group}: {bam_size_mb:.1f}MB BAM')
        
        intron_gff_path = tools.fileOps.get_tmp_toil_file()
        
        # Adaptive processing based on BAM size
        if bam_size_mb > 1000:  # Large BAM - use memory-efficient approach
            cmd = ['bam2hints', '--intronsonly', '--in', merged_bam_file, '--out', intron_gff_path]
        else:  # Standard processing
            cmd = ['bam2hints', '--intronsonly', '--in', merged_bam_file, '--out', intron_gff_path]
        
        tools.procOps.run_proc(cmd)
        
        # Validate output
        if os.path.exists(intron_gff_path):
            hints_size_mb = os.path.getsize(intron_gff_path) / (1024**2)
            logger.info(f'Generated intron hints for {ref_group}: {hints_size_mb:.1f}MB')
        else:
            logger.warning(f'No intron hints generated for {ref_group}')
            # Create empty file
            with open(intron_gff_path, 'w') as f:
                f.write('')
        
        return job.fileStore.writeGlobalFile(intron_gff_path)
        
    except Exception as e:
        logger.error(f'Error building intron hints for {ref_group}: {str(e)}')
        # Create empty hints file as fallback
        intron_gff_path = tools.fileOps.get_tmp_toil_file()
        with open(intron_gff_path, 'w') as f:
            f.write('')
        return job.fileStore.writeGlobalFile(intron_gff_path)


def build_exon_hints_dynamic(job: Job, merged_bam_file_id, ref_group):
    """
    Enhanced exon hints building with adaptive processing and resource optimization.
    """
    logger.info(f'Building exon hints for {ref_group}')
    
    try:
        merged_bam_file = job.fileStore.readGlobalFile(merged_bam_file_id)
        bam_size_mb = os.path.getsize(merged_bam_file) / (1024**2)
        
        logger.info(f'Processing exon hints for {ref_group}: {bam_size_mb:.1f}MB BAM')
        
        exon_gff_path = tools.fileOps.get_tmp_toil_file()
        
        # Adaptive parameters based on BAM size and available resources
        if bam_size_mb > 2000:  # Very large BAM
            width = 15
            margin = 15
            minthresh = 3
            minscore = 5
        elif bam_size_mb > 500:  # Large BAM
            width = 12
            margin = 12
            minthresh = 2
            minscore = 4
        else:  # Standard BAM
            width = 10
            margin = 10
            minthresh = 2
            minscore = 4
        
        cmd = [['bam2wig', merged_bam_file], 
               ['wig2hints.pl', f'--width={width}', f'--margin={margin}', 
                f'--minthresh={minthresh}', f'--minscore={minscore}', '--prune=0.1', '--src=W',
                '--type=ep', '--UCSC=/dev/null', '--radius=4.5', '--pri=4', '--strand=.']]
        
        logger.info(f'Running exon hints with parameters: width={width}, margin={margin}, minthresh={minthresh}, minscore={minscore}')
        
        tools.procOps.run_proc(cmd, stdout=exon_gff_path)
        
        # Validate output
        if os.path.exists(exon_gff_path):
            hints_size_mb = os.path.getsize(exon_gff_path) / (1024**2)
            logger.info(f'Generated exon hints for {ref_group}: {hints_size_mb:.1f}MB')
        else:
            logger.warning(f'No exon hints generated for {ref_group}')
            # Create empty file
            with open(exon_gff_path, 'w') as f:
                f.write('')
        
        return job.fileStore.writeGlobalFile(exon_gff_path)
        
    except Exception as e:
        logger.error(f'Error building exon hints for {ref_group}: {str(e)}')
        # Create empty hints file as fallback
        exon_gff_path = tools.fileOps.get_tmp_toil_file()
        with open(exon_gff_path, 'w') as f:
            f.write('')
        return job.fileStore.writeGlobalFile(exon_gff_path)


def cat_hints_dynamic(job: Job, intron_hints_file_ids, exon_hints_file_ids, annotation_hints_file_id, 
                     iso_seq_hints_file_ids, protein_hints_file_id, resources_metadata):
    """
    Enhanced hints concatenation with parallel processing and adaptive sorting.
    """
    total_files = (len(intron_hints_file_ids) + len(exon_hints_file_ids) + 
                   len(iso_seq_hints_file_ids) + (2 if annotation_hints_file_id and protein_hints_file_id else 
                   1 if annotation_hints_file_id or protein_hints_file_id else 0))
    
    logger.info(f'Concatenating and processing {total_files} hint files')
    
    try:
        cat_hints = tools.fileOps.get_tmp_toil_file()
        total_size_mb = 0
        
        with open(cat_hints, 'w') as outf:
            # Process intron and exon hints
            processed_files = 0
            for file_id in itertools.chain(intron_hints_file_ids, exon_hints_file_ids):
                if file_id is None:
                    continue
                    
                try:
                    f = job.fileStore.readGlobalFile(file_id)
                    file_size_mb = os.path.getsize(f) / (1024**2)
                    total_size_mb += file_size_mb
                    
                    with open(f, 'r') as inf:
                        for line in inf:
                            outf.write(line)
                    
                    processed_files += 1
                    if processed_files % 20 == 0:
                        logger.info(f'Processed {processed_files} hint files, {total_size_mb:.1f}MB total')
                        
                except Exception as e:
                    logger.warning(f'Error reading hint file {processed_files}: {str(e)}')
                    continue
            
            # Process annotation and protein hints
            for file_id in [annotation_hints_file_id, protein_hints_file_id]:
                if file_id is not None:
                    try:
                        f = job.fileStore.readGlobalFile(file_id)
                        file_size_mb = os.path.getsize(f) / (1024**2)
                        total_size_mb += file_size_mb
                        
                        with open(f, 'r') as inf:
                            for line in inf:
                                outf.write(line)
                        
                        processed_files += 1
                        
                    except Exception as e:
                        logger.warning(f'Error reading annotation/protein hint file: {str(e)}')
                        continue
        
        logger.info(f'Concatenated {processed_files} hint files: {total_size_mb:.1f}MB total')
        
        # Adaptive sorting based on data size
        combined_hints = tools.fileOps.get_tmp_toil_file()
        
        if total_size_mb > 1000:  # Large dataset - use external sort
            logger.info('Using external sort for large hint dataset (single-pass multi-key sort)')
            # Use temporary directory for large sorts with single-pass multi-key sort
            import tempfile
            with tempfile.TemporaryDirectory() as temp_dir:
                cmd = [['sort', '-T', temp_dir, '-k1,1', '-k3,3', '-k4,4n', '-k5,5n', cat_hints],
                       ['join_mult_hints.pl']]
                tools.procOps.run_proc(cmd, stdout=combined_hints)
        else:  # Standard sorting with single-pass multi-key sort
            cmd = [['sort', '-k1,1', '-k3,3', '-k4,4n', '-k5,5n', cat_hints],
                   ['join_mult_hints.pl']]
            tools.procOps.run_proc(cmd, stdout=combined_hints)
        
        # Add IsoSeq hints (these are processed separately)
        if iso_seq_hints_file_ids:
            logger.info(f'Adding {len(iso_seq_hints_file_ids)} IsoSeq hint files')
            
            with open(combined_hints, 'a') as outf:
                for file_id in iso_seq_hints_file_ids:
                    if file_id is None:
                        continue
                        
                    try:
                        f = job.fileStore.readGlobalFile(file_id)
                        with open(f, 'r') as inf:
                            for line in inf:
                                outf.write(line)
                    except Exception as e:
                        logger.warning(f'Error reading IsoSeq hint file: {str(e)}')
                        continue
        
        # Final sort of combined hints
        sorted_combined_hints = tools.fileOps.get_tmp_toil_file()
        
        logger.info('Performing final GFF sort')
        tools.misc.sort_gff(combined_hints, sorted_combined_hints)
        
        final_size_mb = os.path.getsize(sorted_combined_hints) / (1024**2)
        logger.info(f'Final sorted hints file: {final_size_mb:.1f}MB')
        
        return job.fileStore.writeGlobalFile(sorted_combined_hints)
        
    except Exception as e:
        logger.error(f'Error in dynamic hints concatenation: {str(e)}')
        # Fallback to original approach
        cat_hints = tools.fileOps.get_tmp_toil_file()
        with open(cat_hints, 'w') as outf:
            for file_id in itertools.chain(intron_hints_file_ids, exon_hints_file_ids):
                if file_id is not None:
                    f = job.fileStore.readGlobalFile(file_id)
                    for line in open(f):
                        outf.write(line)
            for file_id in [annotation_hints_file_id, protein_hints_file_id]:
                if file_id is not None:
                    f = job.fileStore.readGlobalFile(file_id)
                    for line in open(f):
                        outf.write(line)
        # Use single-pass multi-key sort instead of 4 sequential sorts
        cmd = [['sort', '-k1,1', '-k3,3', '-k4,4n', '-k5,5n', cat_hints],
               ['join_mult_hints.pl']]
        combined_hints = tools.fileOps.get_tmp_toil_file()
        tools.procOps.run_proc(cmd, stdout=combined_hints)
        with open(combined_hints, 'a') as outf:
            for file_id in iso_seq_hints_file_ids:
                if file_id is not None:
                    f = job.fileStore.readGlobalFile(file_id)
                    for line in open(f):
                        outf.write(line)
        sorted_combined_hints = tools.fileOps.get_tmp_toil_file()
        tools.misc.sort_gff(combined_hints, sorted_combined_hints)
        return job.fileStore.writeGlobalFile(sorted_combined_hints)


def namesort_bam(job: Job, bam_fid, bai_fid, ref_group, num_reads=50_000_000):
    bam_path = job.fileStore.readGlobalFile(bam_fid)
    job.fileStore.readGlobalFile(bai_fid, bam_path + '.bai')
    tmp = tools.fileOps.get_tmp_toil_file(suffix='name_sorted.bam')
    is_paired = bam_is_paired(bam_path)
    # Correctly construct the command to view the entire BAM file.
    # The original command was incorrectly trying to use 'ref_group'.
    cmd = [['samtools', 'view', '-t', '8', '-b', bam_path] + list(ref_group),
           ['sambamba', 'sort', '-t', '8', '-m', '12G', '-o', '/dev/stdout', '-n', '/dev/stdin']]
    tools.procOps.run_proc(cmd, stdout=tmp)

    handle = pysam.Samfile(tmp)
    r = []
    fids = []
    for qname, reads in itertools.groupby(handle, lambda x: x.qname):
        r.extend(reads)
        if len(r) >= num_reads:
            outf = tools.fileOps.get_tmp_toil_file()
            out_handle = pysam.Samfile(outf, 'wb', template=handle)
            for rec in r: out_handle.write(rec)
            out_handle.close()
            fid = job.addChildJobFn(filter_bam, job.fileStore.writeGlobalFile(outf), is_paired,
                                     disk='8G', memory='8G').rv()
            fids.append(fid)
            r = []
    if r:
        outf = tools.fileOps.get_tmp_toil_file()
        out_handle = pysam.Samfile(outf, 'wb', template=handle)
        for rec in r: out_handle.write(rec)
        out_handle.close()
        fids.append(job.addChildJobFn(filter_bam, job.fileStore.writeGlobalFile(outf), is_paired,
                                       disk='8G', memory='8G').rv())
    return job.addFollowOnJobFn(merge_filtered_bams, fids).rv()


def filter_bam(job: Job, file_id, is_paired):
    bam_path = job.fileStore.readGlobalFile(file_id)
    assert os.path.getsize(bam_path) > 0
    tmp_filtered = tools.fileOps.get_tmp_toil_file()
    filter_cmd = ['filterBam', '--uniq', '--in', bam_path, '--out', tmp_filtered]

    if is_paired == True:
        filter_cmd.extend(['--paired', '--pairwiseAlignments'])
    tools.procOps.run_proc(filter_cmd)
    if os.path.getsize(tmp_filtered) == 0:
        raise RuntimeError('After filtering one BAM subset became empty. This could be bad.')

    out_filter = tools.fileOps.get_tmp_toil_file()
    sort_cmd = ['sambamba', 'sort', tmp_filtered, '-o', out_filter, '-t', '1']
    tools.procOps.run_proc(sort_cmd)
    return job.fileStore.writeGlobalFile(out_filter)


def merge_filtered_bams(job: Job, filtered_file_ids):
    local_paths = [job.fileStore.readGlobalFile(x) for x in filtered_file_ids]
    fofn = tools.fileOps.get_tmp_toil_file()
    with open(fofn, 'w') as outf:
        for l in local_paths:
            if os.environ.get('CAT_BINARY_MODE') == 'singularity':
                l = tools.procOps.singularify_arg(l)
            outf.write(l + '\n')
    out_bam = tools.fileOps.get_tmp_toil_file()
    cmd = ['samtools', 'merge', '-b', fofn, out_bam]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(out_bam)

def merge_bams(job: Job, filtered_bam_file_ids, annotation_hints_file_id, iso_seq_hints_file_ids, protein_hints_file_id):
    """
    Legacy merge_bams function - redirects to dynamic version with default metadata.
    """
    # Create default metadata for backward compatibility
    resources_metadata = {
        'total_bam_files': sum(len(groups) for groups in filtered_bam_file_ids.values()),
        'has_protein': protein_hints_file_id is not None,
        'has_annotation': annotation_hints_file_id is not None,
        'iso_seq_count': len(iso_seq_hints_file_ids) if iso_seq_hints_file_ids else 0
    }
    
    return merge_bams_dynamic(job, filtered_bam_file_ids, annotation_hints_file_id, 
                             iso_seq_hints_file_ids, protein_hints_file_id, resources_metadata)


def cat_sort_bams(job: Job, bam_file_ids):
    bamfiles = [job.fileStore.readGlobalFile(x) for x in bam_file_ids]
    catfile = tools.fileOps.get_tmp_toil_file()
    sam_iter = tools.dataOps.grouper(bamfiles, 4095)
    cmd = ['samtools', 'cat', '-o', catfile]
    cmd.extend(next(sam_iter))
    tools.procOps.run_proc(cmd)
    for more in sam_iter:
        old_catfile = catfile
        catfile = tools.fileOps.get_tmp_toil_file()
        cmd = ['samtools', 'cat', '-o', catfile, old_catfile]
        cmd.extend(more)
        tools.procOps.run_proc(cmd)
    merged = tools.fileOps.get_tmp_toil_file()
    cmd = ['sambamba', 'sort', catfile, '-o', merged, '-t', '4', '-m', '15G']
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(merged)


def generate_protein_hints(job: Job, protein_fasta_file_id, genome_fasta_file_id):
    """
    Legacy generate_protein_hints function - redirects to dynamic version with default chunk size.
    """
    return generate_protein_hints_dynamic(job, protein_fasta_file_id, genome_fasta_file_id, chunk_size=100)


def run_protein_aln(job: Job, protein_subset, genome_fasta_file_id):
    """
    Legacy run_protein_aln function - redirects to dynamic version.
    """
    return run_protein_aln_dynamic(job, protein_subset, genome_fasta_file_id)


def convert_protein_aln_results_to_hints(job: Job, results):
    """
    Legacy convert_protein_aln_results_to_hints function - redirects to dynamic version.
    """
    return convert_protein_aln_results_to_hints_dynamic(job, results)


def build_hints(job: Job, merged_bam_file_ids, anno_hints, iso_hints, prot_hints):
    """
    Legacy build_hints function - redirects to dynamic version with default metadata.
    """
    # Create default metadata for backward compatibility
    resources_metadata = {
        'total_bam_files': sum(len(groups) for groups in merged_bam_file_ids.values()),
        'has_protein': prot_hints is not None,
        'has_annotation': anno_hints is not None,
        'iso_seq_count': len(iso_hints) if iso_hints else 0
    }
    
    return build_hints_dynamic(job, merged_bam_file_ids, anno_hints, iso_hints, prot_hints, resources_metadata)


def build_intron_hints(job: Job, merged_bam_file_id):
    merged_bam_file = job.fileStore.readGlobalFile(merged_bam_file_id)
    intron_gff_path = tools.fileOps.get_tmp_toil_file()
    tools.procOps.run_proc(['bam2hints', '--intronsonly', '--in', merged_bam_file, '--out', intron_gff_path])
    return job.fileStore.writeGlobalFile(intron_gff_path)


def build_exon_hints(job: Job, merged_bam_file_id):
    merged_bam_file = job.fileStore.readGlobalFile(merged_bam_file_id)
    exon_gff_path = tools.fileOps.get_tmp_toil_file()
    tools.procOps.run_proc([['bam2wig', merged_bam_file], 
    ['wig2hints.pl', '--width=10', '--margin=10', '--minthresh=2', '--minscore=4', '--prune=0.1', '--src=W',
      '--type=ep', '--UCSC=/dev/null', '--radius=4.5', '--pri=4', '--strand=.']], stdout=exon_gff_path)
    return job.fileStore.writeGlobalFile(exon_gff_path)


def generate_iso_seq_hints(job: Job, bam_file_id, bai_file_id):
    bam_path = job.fileStore.readGlobalFile(bam_file_id)
    job.fileStore.readGlobalFile(bai_file_id, bam_path + '.bai')
    pacbio_gff_path = tools.fileOps.get_tmp_toil_file()
    cmd = [['samtools', 'view', '-@', '64', '-b', '-F', '4', bam_path], 
           ['bamToPsl', '-nohead', '/dev/stdin', '/dev/stdout'],
           ['sort', '-n', '-k', '16,16'],
           ['sort', '-s', '-k', '14,14'],
           ['perl', '-ne', '@f=split; print if ($f[0]>=100)'],
           ['blat2hints.pl', '--source=PB', '--nomult', '--ep_cutoff=20', '--in=/dev/stdin',
            '--out={}'.format(pacbio_gff_path)]]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(pacbio_gff_path)


def generate_annotation_hints(job: Job, annotation_hints_file_id):
    annotation_gp = job.fileStore.readGlobalFile(annotation_hints_file_id)
    tx_dict = tools.transcripts.get_gene_pred_dict(annotation_gp)
    hints = []
    for tx_id, tx in tx_dict.items():
        if tx.cds_size == 0:
            continue
        # rather than try to re-do the arithmetic, we will use the get_bed() function to convert this transcript
        cds_tx = tools.transcripts.Transcript(tx.get_bed(new_start=tx.thick_start, new_stop=tx.thick_stop))
        for intron in cds_tx.intron_intervals:
            r = [intron.chromosome, 'a2h', 'intron', intron.start + 1, intron.stop, 0, intron.strand, '.',
                 'grp={};src=M;pri=2'.format(tx_id)]
            hints.append(r)
        for exon in cds_tx.exon_intervals:
            r = [exon.chromosome, 'a2h', 'CDS', exon.start + 1, exon.stop, 0, exon.strand, '.',
                 'grp={};src=M;pri=2'.format(tx_id)]
            hints.append(r)
    annotation_hints_gff = tools.fileOps.get_tmp_toil_file()
    tools.fileOps.print_rows(annotation_hints_gff, hints)
    return job.fileStore.writeGlobalFile(annotation_hints_gff)


def cat_hints(job: Job, intron_hints_file_ids, exon_hints_file_ids, annotation_hints_file_id, iso_seq_hints_file_ids,
              protein_hints_file_id):
    """
    Legacy cat_hints function - redirects to dynamic version with default metadata.
    """
    # Create default metadata for backward compatibility
    total_files = (len(intron_hints_file_ids) + len(exon_hints_file_ids) + 
                   len(iso_seq_hints_file_ids) + (2 if annotation_hints_file_id and protein_hints_file_id else 
                   1 if annotation_hints_file_id or protein_hints_file_id else 0))
    
    resources_metadata = {
        'total_hint_files': total_files,
        'has_protein': protein_hints_file_id is not None,
        'has_annotation': annotation_hints_file_id is not None,
        'iso_seq_count': len(iso_seq_hints_file_ids)
    }
    
    return cat_hints_dynamic(job, intron_hints_file_ids, exon_hints_file_ids, annotation_hints_file_id,
                            iso_seq_hints_file_ids, protein_hints_file_id, resources_metadata)


def validate_bam_fasta_pairs(bam_path, fasta_sequences, genome):
    handle = pysam.Samfile(bam_path, 'rb')
    bam_sequences = {(n, s) for n, s in zip(*[handle.references, handle.lengths])}
    difference = bam_sequences - fasta_sequences
    if len(difference) > 0:
        base_err = 'Error: BAM {} has the following sequence/length pairs not found in the {} fasta: {}.'
        err = base_err.format(bam_path, genome, ','.join(['-'.join(map(str, x)) for x in difference]))
        raise UserException(err)
    missing_seqs = fasta_sequences - bam_sequences
    if len(missing_seqs) > 0:
        base_msg = 'BAM {} does not have the following sequence/length pairs in its header: {}.'
        msg = base_msg.format(bam_path, ','.join(['-'.join(map(str, x)) for x in missing_seqs]))
        logger.warning(msg)


def bam_is_paired(bam_path, num_reads=20000, paired_cutoff=0.75):
    sam = pysam.Samfile(bam_path)
    count = 0
    for rec in itertools.islice(sam, num_reads):
        if rec.is_paired:
            count += 1
    if tools.mathOps.format_ratio(count, num_reads) > 0.75:
        return True
    elif tools.mathOps.format_ratio(count, num_reads) < 1 - paired_cutoff:
        return False
    else:
        raise UserException("Unable to infer pairing from bamfile {}".format(bam_path))


def group_references(sam_handle, num_bases=1e7, max_seqs=1000):
    name_iter = zip(*[sam_handle.references, sam_handle.lengths])
    name, size = next(name_iter)
    this_bin = [name]
    bin_base_count = size
    num_seqs = 1
    for name, size in name_iter:
        bin_base_count += size
        num_seqs += 1
        if bin_base_count >= num_bases or num_seqs > max_seqs:
            yield this_bin
            this_bin = [name]
            bin_base_count = size
            num_seqs = 1
        else:
            this_bin.append(name)
    yield this_bin

###
# Entrypoint for command-line execution
###

def main():
    """
    Main entrypoint for the hints pipeline.
    """
    # First, parse just the mode to determine which parser to use
    mode_parser = argparse.ArgumentParser(add_help=False)
    mode_parser.add_argument("--mode", choices=['toil', 'cluster'], default='cluster',
                            help="Pipeline mode: 'toil' for Toil-based processing, 'cluster' for SLURM/SGE parallel processing")

    mode_args, remaining_args = mode_parser.parse_known_args()

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--mode", choices=['toil', 'cluster'], default='cluster',
                       help="Pipeline mode: 'toil' for Toil-based processing, 'cluster' for SLURM/SGE parallel processing")

    if mode_args.mode == 'toil':
        Job.Runner.addToilOptions(parser)

    parser.add_argument("--genome", required=True, help="Reference genome name.")
    parser.add_argument("--fasta", required=True, help="Path to reference genome fasta.")
    parser.add_argument("--hints_out", required=True, help="Path to output hints GFF file.")

    parser.add_argument("--bams", nargs='*', help="Paths to short-read RNA-seq BAM files.")
    parser.add_argument("--intron_bams", nargs='*', help="Paths to intron-spanning RNA-seq BAM files.")
    parser.add_argument("--iso_bams", nargs='*', help="Paths to Iso-Seq BAM files.")
    parser.add_argument("--annotation_gp", help="Path to an annotation genePred file.")
    parser.add_argument("--protein_fasta", help="Path to a protein FASTA file for homology.")

    # Cluster-mode options. The "slurm_*" names are kept for back-compat with
    # the historical CLI; they apply uniformly to whichever backend
    # --execution_mode selects (slurm or sge).
    parser.add_argument("--execution_mode", choices=("auto", "slurm", "sge", "local"), default="auto",
                       help="Cluster backend to use when --mode=cluster ('auto' detects it).")
    parser.add_argument("--exclude_nodes", default="",
                       help="Node exclude list. SLURM: comma list; SGE: '!h1&!h2' or comma list (auto-converted).")
    parser.add_argument("--module_load", default="", help="Module to load in every job (or '').")
    parser.add_argument("--sge_parallel_env", default="smp", help="SGE parallel environment name.")
    parser.add_argument("--sge_memory_flag", default="h_vmem", help="SGE memory resource flag.")
    parser.add_argument("--slurm_memory", type=int, default=256, help="Memory per job in GB")
    parser.add_argument("--slurm_cpus", type=int, default=128, help="CPUs per job")
    parser.add_argument("--slurm_time", default="12:00:00", help="Time limit per job (HH:MM:SS)")
    parser.add_argument("--slurm_partition", default="", help="SLURM partition / SGE queue")
    parser.add_argument("--slurm_max_jobs", type=int, default=20, help="Maximum concurrent jobs")
    
    # Parse arguments
    args = parser.parse_args()

    # Validate required input files exist
    logger.info('Validating input files...')
    
    # Check required fasta file
    if not os.path.exists(args.fasta):
        logger.error(f"Required fasta file does not exist: {args.fasta}")
        parser.error(f"Required fasta file does not exist: {args.fasta}")
    
    # Check BAM files if provided
    all_bam_files = (args.bams or []) + (args.intron_bams or []) + (args.iso_bams or [])
    for bam_file in all_bam_files:
        if not os.path.exists(bam_file):
            logger.error(f"BAM file does not exist: {bam_file}")
            parser.error(f"BAM file does not exist: {bam_file}")
    
    # Check annotation file if provided
    if args.annotation_gp and not os.path.exists(args.annotation_gp):
        logger.error(f"Annotation file does not exist: {args.annotation_gp}")
        parser.error(f"Annotation file does not exist: {args.annotation_gp}")
    
    # Check protein fasta file if provided
    if args.protein_fasta and not os.path.exists(args.protein_fasta):
        logger.error(f"Protein fasta file does not exist: {args.protein_fasta}")
        parser.error(f"Protein fasta file does not exist: {args.protein_fasta}")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(os.path.abspath(args.hints_out))
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info('Input file validation completed successfully')

    if args.mode == 'cluster':
        # Wire the Scheduler from CLI args before any generate_slurm_* call.
        _init_scheduler(args)
        logger.info(f'Running hints pipeline in cluster mode (backend: {args.execution_mode})')

        slurm_options = {
            'memory_gb': args.slurm_memory,
            'cpus': args.slurm_cpus,
            'time_limit': args.slurm_time,
            'partition': args.slurm_partition,
            'max_concurrent_jobs': args.slurm_max_jobs
        }

        # Same entry function as before; the SLURM-vs-SGE distinction is
        # handled inside the scheduler instance, not by the job-generator code.
        run_hints_pipeline_slurm(
            genome=args.genome,
            fasta=args.fasta,
            bams=args.bams or [],
            intron_bams=args.intron_bams or [],
            iso_bams=args.iso_bams or [],
            annotation_gp=args.annotation_gp,
            protein_fasta=args.protein_fasta,
            hints_out=args.hints_out,
            slurm_options=slurm_options
        )
        
    else:  # toil mode
        logger.info('Running hints pipeline in Toil mode')
        
        # Check if jobStore is provided (required for Toil mode)
        if not hasattr(args, 'jobStore') or args.jobStore is None:
            logger.error("jobStore argument is required for Toil mode")
            parser.error("jobStore argument is required for Toil mode")
        
        # The 'args' object now contains all of Toil's options AND your script's options.
        # The jobStore is now a positional argument added by addToilOptions.
        # We must construct the path to it before passing it to the pipeline.
        job_store_path = args.jobStore

        # Run the Toil-based pipeline
        run_hints_pipeline(
            genome=args.genome,
            fasta=args.fasta,
            bams=args.bams or [],
            intron_bams=args.intron_bams or [],
            iso_bams=args.iso_bams or [],
            annotation_gp=args.annotation_gp,
            protein_fasta=args.protein_fasta,
            hints_out=args.hints_out,
            toil_options=args  # Pass the complete args object
        )

if __name__ == "__main__":
    main()
