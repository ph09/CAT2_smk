#!/usr/bin/env python3
"""
Parallel Augustus PB Pipeline for CAT
"""

import argparse
import os
import sys
import subprocess
import tempfile
import logging
from pathlib import Path
import shutil
from typing import List, Tuple, Dict, Optional
import time
import collections

from cat.scheduler import get_scheduler


def _scheduler_from_args(args):
    """Construct a Scheduler from CLI args (mirrors augustus_parallel.py)."""
    cfg = {
        "cluster": {
            "slurm": {
                "partition": getattr(args, "slurm_partition", None),
                "exclude_nodes": getattr(args, "slurm_exclude_nodes", "") or "",
                "module_load": getattr(args, "module_load", "") or "",
            },
            "sge": {
                "queue": getattr(args, "slurm_partition", None),
                "parallel_env": getattr(args, "sge_parallel_env", "smp"),
                "memory_flag": getattr(args, "sge_memory_flag", "h_vmem"),
                "hostname_exclude": getattr(args, "slurm_exclude_nodes", "") or "",
                "module_load": getattr(args, "module_load", "") or "",
            },
        }
    }
    return get_scheduler(getattr(args, "execution_mode", "slurm"), cfg)


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class ParallelAugustusPB:
    """Main class for running parallel Augustus PB pipeline."""
    
    def __init__(self, args):
        self.args = args
        self._scheduler = _scheduler_from_args(args)
        self.args.genome_fasta = os.path.abspath(self.args.genome_fasta)
        self.args.chrom_sizes = os.path.abspath(self.args.chrom_sizes)
        self.args.pb_cfg = os.path.abspath(self.args.pb_cfg)
        self.args.hints_gff = os.path.abspath(self.args.hints_gff)
        self.args.raw_gtf = os.path.abspath(self.args.raw_gtf)
        self.args.gtf = os.path.abspath(self.args.gtf)
        self.args.gp = os.path.abspath(self.args.gp)
        self.args.work_dir = os.path.abspath(self.args.work_dir)

        self.work_dir = Path(os.path.abspath(args.work_dir))
        self.temp_dir = self.work_dir / "augustus_pb_parallel_temp"
        self.split_dir = self.temp_dir / "split"
        self.jobs_dir = self.temp_dir / "jobs"
        self.results_dir = self.temp_dir / "results"

        # Unique job prefix to disambiguate parallel runs.
        self.job_prefix = os.path.basename(self.args.work_dir)

        for dir_path in [self.temp_dir, self.split_dir, self.jobs_dir, self.results_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

    def validate_inputs(self) -> bool:
        """Validate all required input files exist."""
        required_files = [
            self.args.genome_fasta,
            self.args.chrom_sizes,
            self.args.pb_cfg,
            self.args.hints_gff
        ]
        
        for file_path in required_files:
            if not os.path.exists(file_path):
                logger.error(f"Required file not found: {file_path}")
                return False
        
        return True
    
    def split_genome_by_chromosome(self) -> bool:
        """Step 1: Split genome into chromosomes and get chromosome lengths."""
        logger.info("Step 1: Splitting genome by chromosome...")
        
        try:
            # Run splitMfasta.pl to split genome
            split_cmd = f"splitMfasta.pl {self.args.genome_fasta} --outputpath={self.split_dir}"
            
            logger.info(f"Running: {split_cmd}")
            result = subprocess.run(split_cmd, shell=True, capture_output=True, text=True, check=True)
            
            # Rename files to use chromosome names
            for split_file in self.split_dir.glob("*.split.*"):
                # Extract chromosome name from FASTA header
                with open(split_file, 'r') as f:
                    first_line = f.readline().strip()
                    if first_line.startswith('>'):
                        chrom_name = first_line[1:].split()[0]  # Get chromosome name
                        new_name = f"{chrom_name}.fa"
                        shutil.move(str(split_file), str(self.split_dir / new_name))
                        logger.info(f"Renamed {split_file.name} to {new_name}")
            
            # Generate ACGT content summary for each split chromosome
            summary_out = self.temp_dir / "summary.out"
            with open(summary_out, 'w') as f:
                for split_file in self.split_dir.glob("*.fa"):
                    summary_cmd = f"summarizeACGTcontent.pl {split_file}"
                    logger.info(f"Running: {summary_cmd}")
                    subprocess.run(summary_cmd, shell=True, stdout=f, stderr=subprocess.PIPE, text=True, check=True)
            
            logger.info("Genome splitting completed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error splitting genome: {e}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            return False
    
    def create_hints_location_file(self) -> bool:
        """Step 2: Create hints location file for each chromosome."""
        logger.info("Step 2: Creating hints location file...")
        
        try:
            summary_file = self.temp_dir / "summary.out"
            chr_lst_file = self.temp_dir / "chr.lst"
            
            # Parse summary.out to create chr.lst
            with open(summary_file, 'r') as f_in, open(chr_lst_file, 'w') as f_out:
                for line in f_in:
                    if "bases" in line:
                        # Extract chromosome name and length with robust format handling
                        parts = line.strip().split()
                        chrom_length = None
                        chrom_name = None
                        # Format A: "123456789 bases in chrX"
                        if len(parts) >= 4 and parts[1] == "bases" and parts[2] == "in":
                            chrom_length = parts[0]
                            chrom_name = parts[3]
                        # Format B: "210155 bases.\tchr16 ..."
                        elif len(parts) >= 3 and parts[1].startswith("bases"):
                            chrom_length = parts[0]
                            chrom_name = parts[2]
                        if chrom_length and chrom_name:
                            # Create hints file path (will be created later)
                            hints_file = self.temp_dir / f"{chrom_name}_hints.gff"
                            # Write to chr.lst format: chromosome_file\thints_file\t1\tlength
                            f_out.write(f"{self.split_dir}/{chrom_name}.fa\t{hints_file}\t1\t{chrom_length}\n")
            
            logger.info(f"Created chromosome list: {chr_lst_file}")
            return True
            
        except Exception as e:
            logger.error(f"Error creating hints location file: {e}")
            return False
    
    def create_chromosome_hints_files(self) -> bool:
        """Create PB hints files for each chromosome."""
        logger.info("Creating chromosome-specific PB hints files...")
        
        try:
            # Load PB hints
            hints_file = self.args.hints_gff
            hints = [x.split('\t') for x in open(hints_file) if 'src=PB' in x]
            
            if len(hints) == 0:
                logger.warning("No PB hints found.")
                return True
            
            logger.info(f'Found {len(hints)} PB hints')
            
            # Convert the start/stops to ints and break up by chromosome
            hints_by_chrom = collections.defaultdict(list)
            for h in hints:
                h[3] = int(h[3])
                h[4] = int(h[4])
                hints_by_chrom[h[0]].append(h)
            
            # Create hints for each chromosome
            for chrom, chrom_hints in hints_by_chrom.items():
                hints_file = self.temp_dir / f"{chrom}_hints.gff"
                
                with open(hints_file, 'w') as outf:
                    for hint in chrom_hints:
                        outf.write('\t'.join(map(str, hint)) + '\n')
                
                logger.info(f"Created hints file for chromosome {chrom}: {hints_file} ({len(chrom_hints)} hints)")
            
            return True
            
        except Exception as e:
            logger.error(f"Error creating chromosome hints files: {e}")
            return False
    
    def create_augustus_jobs(self) -> bool:
        """Step 3: Create Augustus job list using createAugustusJoblist.pl."""
        logger.info("Step 3: Creating Augustus PB jobs...")
        
        try:
            chr_lst_file = self.temp_dir / "chr.lst"
            jobs_lst_file = self.temp_dir / "jobs_PB.lst"
            
            # Build Augustus command for PB mode with parameters that createAugustusJoblist.pl will augment.
            # IMPORTANT: Do NOT include {seqfile}, --predictionStart/End, or --hintsfile placeholders here.
            aug_cmd = [
                "augustus",
                "--softmasking=1",
                "--allow_hinted_splicesites=atac",
                "--alternatives-from-evidence=1",
                "--UTR={}".format(int(self.args.utr)),
                "--extrinsicCfgFile={}".format(self.args.pb_cfg),
                "--species={}".format(self.args.species),
                "--/augustus/verbosity=0"
            ]
            
            aug_call = " ".join(aug_cmd)
            
            # Create job list using argument list (avoids fragile line continuations)
            create_jobs_cmd = [
                "createAugustusJoblist.pl",
                f"--sequences={chr_lst_file}",
                "--wrap=#",
                f"--overlap={self.args.overlap}",
                f"--chunksize={self.args.chunksize}",
                f"--outputdir={self.jobs_dir}",
                f"--joblist={jobs_lst_file}",
                f"--jobprefix={self.job_prefix}_aug_PB_",
                f"--command={aug_call}"
            ]
            
            logger.info(f"Running: {' '.join(map(str, create_jobs_cmd))}")
            # Run in temp_dir so job files are created there (createAugustusJoblist.pl creates jobs in CWD)
            result = subprocess.run(create_jobs_cmd, capture_output=True, text=True, check=True, cwd=str(self.temp_dir))
            
            logger.info(f"Created job list: {jobs_lst_file}")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating Augustus jobs: {e}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            return False

    def create_initial_setup_slurm_script(self) -> str:
        """Create cluster job script for initial setup (genome splitting and chr.lst creation)."""
        logger.info("Creating cluster script for initial setup (genome splitting and chr.lst)...")

        header = self._scheduler.header(
            job_name="augPB_setup",
            cpus=1,
            mem="16G",
            walltime=getattr(self.args, 'slurm_setup_time', '04:00:00'),
            log_out=f"{self.temp_dir}/augPB_setup_%j.out",
            log_err=f"{self.temp_dir}/augPB_setup_%j.err",
            partition=getattr(self.args, 'slurm_partition', 'medium'),
            queue=getattr(self.args, 'slurm_partition', 'medium'),
        )
        slurm_script = header + rf"""
set -euo pipefail

echo "Starting PB initial setup..."
echo "Split directory: {self.split_dir}"

# Step 1: Split genome by chromosome
splitMfasta.pl {self.args.genome_fasta} --outputpath={self.split_dir}

# Rename files to use chromosome names
for split_file in {self.split_dir}/*.split.*; do
    if [ -f "$split_file" ]; then
        chrom_name=$(head -n1 "$split_file" | sed 's/^>//' | awk '{{print $1}}')
        new_name="{self.split_dir}/${{chrom_name}}.fa"
        mv "$split_file" "$new_name"
        echo "Renamed $(basename $split_file) to $(basename $new_name)"
    fi
done

# Generate ACGT content summary
summary_out="{self.temp_dir}/summary.out"
> "$summary_out"
for split_file in {self.split_dir}/*.fa; do
    if [ -f "$split_file" ]; then
        summarizeACGTcontent.pl "$split_file" >> "$summary_out"
    fi
done

# Step 2: Create chromosome list (chr.lst)
chr_lst_file="{self.temp_dir}/chr.lst"
chrom_list="{self.temp_dir}/chromosomes.txt"
> "$chr_lst_file"
> "$chrom_list"

while IFS= read -r line; do
    if [[ "$line" == *"bases"* ]]; then
        # Extract chromosome name and length from summarizeACGTcontent.pl output
        parts=($line)
        if [[ ${{#parts[@]}} -ge 3 && "${{parts[1]}}" == "bases."* ]]; then
            chrom_length="${{parts[0]}}"
            chrom_name="${{parts[2]}}"  # The chromosome name is at position 2
            
            # Create hints file path (will be created by array job)
            hints_file="{self.temp_dir}/${{chrom_name}}_hints.gff"
            
            # Write to chr.lst format: chromosome_file	hints_file	1	length
            echo -e "{self.split_dir}/${{chrom_name}}.fa\\t${{hints_file}}\\t1\\t${{chrom_length}}" >> "$chr_lst_file"
            
            # Also save chromosome name for array job
            echo "$chrom_name" >> "$chrom_list"
        fi
    fi
done < "$summary_out"

echo "Created chromosome list: $chr_lst_file"
echo "Created chromosome names file: $chrom_list"
echo "PB initial setup completed successfully!"
"""

        slurm_script_path = self.temp_dir / "augPB_setup.slurm"
        with open(slurm_script_path, 'w') as f:
            f.write(slurm_script)

        logger.info(f"Created PB initial setup SLURM script: {slurm_script_path}")
        return str(slurm_script_path)

    def create_hints_generation_slurm_script(self, num_chromosomes: int, dependency_job_id: str = None) -> str:
        """Create cluster array job script for PB hints splitting (per chromosome)."""
        logger.info(f"Creating cluster array script for PB hints generation ({num_chromosomes} chromosomes)...")

        dependency = self._scheduler.depends_on_job_id(dependency_job_id) if dependency_job_id else None
        task_var = self._scheduler.task_id_env()
        header = self._scheduler.header(
            job_name="augPB_hints",
            cpus=1,
            mem=getattr(self.args, 'slurm_hints_mem', '8G'),
            walltime=getattr(self.args, 'slurm_hints_time', '04:00:00'),
            log_out=f"{self.temp_dir}/augPB_hints_%A_%a.out",
            log_err=f"{self.temp_dir}/augPB_hints_%A_%a.err",
            partition=getattr(self.args, 'slurm_partition', 'medium'),
            queue=getattr(self.args, 'slurm_partition', 'medium'),
            array=(1, num_chromosomes),
            max_concurrent=getattr(self.args, 'slurm_hints_concurrency', 10),
            dependency=dependency,
        )
        slurm_script = header + rf"""
set -euo pipefail

TASK_ID="${{{task_var}}}"
CHROM=$(sed -n "${{TASK_ID}}p" {self.temp_dir}/chromosomes.txt)
echo "Creating PB hints for chromosome: $CHROM"

awk -v c="$CHROM" -F"\t" 'BEGIN{{OFS="\t"}} $1==c && $0 !~ /^#/ && $9 ~ /src=PB/ {{print $0}}' {self.args.hints_gff} > {self.temp_dir}/${{CHROM}}_hints.gff

echo "Wrote hints: {self.temp_dir}/${{CHROM}}_hints.gff"
"""

        slurm_script_path = self.temp_dir / "augPB_hints_array.slurm"
        with open(slurm_script_path, 'w') as f:
            f.write(slurm_script)
        logger.info(f"Created PB hints generation cluster array script: {slurm_script_path}")
        return str(slurm_script_path)

    def create_joblist_generation_slurm_script(self, dependency_job_id: str = None) -> str:
        """Create cluster job script to generate PB Augustus job list."""
        logger.info("Creating cluster script for PB job list generation...")

        dependency = self._scheduler.depends_on_job_id(dependency_job_id) if dependency_job_id else None
        header = self._scheduler.header(
            job_name="augPB_joblist",
            cpus=1,
            mem="16G",
            walltime=getattr(self.args, 'slurm_setup_time', '04:00:00'),
            log_out=f"{self.temp_dir}/augPB_joblist_%j.out",
            log_err=f"{self.temp_dir}/augPB_joblist_%j.err",
            partition=getattr(self.args, 'slurm_partition', 'medium'),
            queue=getattr(self.args, 'slurm_partition', 'medium'),
            dependency=dependency,
        )
        slurm_script = header + rf"""
set -euo pipefail

chr_lst_file="{self.temp_dir}/chr.lst"
aug_cmd="augustus --softmasking=1 --allow_hinted_splicesites=atac --alternatives-from-evidence=1 --UTR={int(self.args.utr)} --extrinsicCfgFile={self.args.pb_cfg} --species={self.args.species} --/augustus/verbosity=0"

# Change to the temp directory to ensure createAugustusJoblist.pl creates job files here
cd "{self.temp_dir}"

createAugustusJoblist.pl --sequences="$chr_lst_file" --wrap="#" --overlap={self.args.overlap} --chunksize={self.args.chunksize} --outputdir="{self.jobs_dir}" --joblist="{self.temp_dir}/jobs_PB.lst" --jobprefix={self.job_prefix}_aug_PB_ --command="$aug_cmd"

echo "Created PB job list: {self.temp_dir}/jobs_PB.lst"
"""

        slurm_script_path = self.temp_dir / "augPB_joblist.slurm"
        with open(slurm_script_path, 'w') as f:
            f.write(slurm_script)
        logger.info(f"Created PB job list SLURM script: {slurm_script_path}")
        return str(slurm_script_path)
    
    def create_slurm_script(self, num_jobs: int = 30, dependency_job_id: str = None) -> str:
        """Step 4: Create cluster array job script for parallel execution."""
        logger.info("Step 4: Creating cluster script for PB mode...")

        dependency = self._scheduler.depends_on_job_id(dependency_job_id) if dependency_job_id else None
        task_var = self._scheduler.task_id_env()
        header = self._scheduler.header(
            job_name="augustus_PB",
            cpus=1,
            mem=getattr(self.args, 'slurm_jobs_mem', '32G'),
            walltime=getattr(self.args, 'slurm_jobs_time', '24:00:00'),
            log_out=f"{self.temp_dir}/augustus_PB_%A_%a.out",
            log_err=f"{self.temp_dir}/augustus_PB_%A_%a.err",
            partition=getattr(self.args, 'slurm_jobs_partition', 'long'),
            queue=getattr(self.args, 'slurm_jobs_partition', 'long'),
            array=(1, num_jobs),
            max_concurrent=getattr(self.args, 'slurm_jobs_concurrency', 10),
            dependency=dependency,
        )
        slurm_script = header + f"""
# Change to temp directory (job paths in the list are relative to this directory)
cd {self.temp_dir}

# Get the job file for this array task (1-based — works on both SLURM and SGE).
TASK_ID="${{{task_var}}}"
JOB_FILE=$(sed -n "${{TASK_ID}}p" jobs_PB.lst)

# Run the Augustus job
if [ -f "$JOB_FILE" ]; then
    echo "Running Augustus PB job: $JOB_FILE"
    bash "$JOB_FILE"
    echo "Completed job: $JOB_FILE"
else
    echo "Job file not found: $JOB_FILE"
    echo "Current directory: $(pwd)"
    echo "Looking for: $JOB_FILE"
    ls -la "$JOB_FILE" 2>&1 || true
    exit 1
fi
"""
        
        slurm_script_path = self.temp_dir / "augustus_PB.slurm"
        with open(slurm_script_path, 'w') as f:
            f.write(slurm_script)
        
        logger.info(f"Created SLURM script: {slurm_script_path}")
        return str(slurm_script_path)
    
    def submit_slurm_job(self, slurm_script_path: str, job_name: str = "job") -> str:
        """Submit a cluster job and return its backend job ID."""
        logger.info(f"Submitting cluster {job_name}...")
        try:
            job_id = self._scheduler.submit(slurm_script_path)
            logger.info(f"Submitted cluster {job_name}: {job_id}")
            return job_id
        except (subprocess.CalledProcessError, RuntimeError) as e:
            logger.error(f"Error submitting cluster {job_name}: {e}")
            return None

    def wait_for_slurm_job(self, job_id: str, job_name: str = "job", check_interval: int = 30) -> bool:
        """Wait for a cluster job to complete and verify success.

        Poll ``job_present()`` until the job drains, then ask the scheduler
        to validate the final state via ``verify_completed()`` (SLURM uses
        sacct; SGE is a no-op since qacct is non-portable). Returns False on
        failure so the caller can abort the dependency chain instead of
        submitting another job that will instantly land in
        ``DependencyNeverSatisfied``.
        """
        logger.info(f"Waiting for cluster {job_name} {job_id} to complete...")
        try:
            while True:
                if not self._scheduler.job_present(job_id):
                    logger.info(f"Cluster {job_name} {job_id} no longer in queue")
                    break
                time.sleep(check_interval)
        except KeyboardInterrupt:
            logger.warning(f"Interrupted while waiting for {job_name} {job_id}")
            return False

        result = self._scheduler.verify_completed(job_id)
        if not result.ok:
            logger.error(f"Cluster {job_name} {job_id} failed: {result.detail}")
            return False
        if result.detail:
            logger.info(f"Cluster {job_name} {job_id} verified: {result.detail}")
        return True

    def run_slurm_jobs(self, slurm_script_path: str) -> bool:
        """Submit and wait for a cluster job array.

        Returns True once the array drains from the queue. Per-task success
        is validated downstream via output-file inspection (portable across
        SLURM/SGE flavors).
        """
        logger.info(f"Submitting {self._scheduler.name} jobs...")
        try:
            job_id = self._scheduler.submit(slurm_script_path)
            logger.info(f"Submitted {self._scheduler.name} job array: {job_id}")

            logger.info(f"Waiting for {self._scheduler.name} jobs to complete...")
            while True:
                if not self._scheduler.job_present(job_id):
                    logger.info(f"{self._scheduler.name} jobs no longer in queue")
                    break
                time.sleep(60)

            result = self._scheduler.verify_completed(job_id)
            if not result.ok:
                logger.error(f"{self._scheduler.name} job array {job_id} failed: {result.detail}")
                return False
            logger.info(f"{self._scheduler.name} job array {job_id} completed")
            return True

        except (subprocess.CalledProcessError, RuntimeError) as e:
            logger.error(f"Error running cluster jobs: {e}")
            return False
    
    def run_local_jobs(self, jobs_file: str, num_cpus: int = None) -> bool:
        """Run Augustus jobs locally using multiprocessing."""
        import multiprocessing
        
        logger.info("Running Augustus jobs locally...")
        
        try:
            # Read job commands from file
            with open(jobs_file, 'r') as f:
                jobs = [line.strip() for line in f if line.strip()]
            
            logger.info(f"Total jobs to run: {len(jobs)}")
            
            # Determine number of CPUs to use
            if num_cpus is None:
                num_cpus = multiprocessing.cpu_count()
            
            logger.info(f"Using {num_cpus} CPUs for parallel execution")
            
            # Function to run a single job
            def run_job(job_cmd: str) -> tuple:
                try:
                    result = subprocess.run(job_cmd, shell=True, capture_output=True, text=True, check=True)
                    return (True, job_cmd, "")
                except subprocess.CalledProcessError as e:
                    error_msg = f"STDOUT: {e.stdout}\nSTDERR: {e.stderr}"
                    return (False, job_cmd, error_msg)
            
            # Run jobs in parallel using multiprocessing
            with multiprocessing.Pool(processes=num_cpus) as pool:
                results = pool.map(run_job, jobs)
            
            # Check results
            failed_jobs = [(job, error) for success, job, error in results if not success]
            
            if failed_jobs:
                logger.error(f"{len(failed_jobs)} jobs failed:")
                for job, error in failed_jobs[:5]:  # Show first 5 errors
                    logger.error(f"Job: {job}\nError: {error}")
                return False
            
            logger.info(f"All {len(jobs)} jobs completed successfully")
            return True
            
        except Exception as e:
            logger.error(f"Error running local jobs: {e}")
            return False
    
    def merge_augustus_output(self) -> str:
        """Step 5: Merge Augustus output by globbing produced GFFs from jobs directory."""
        logger.info("Step 5: Merging Augustus PB output...")

        try:
            raw_gtf = self.temp_dir / "augustus_PB_raw.gtf"
            raw_gtf_fofn = self.temp_dir / "augustus_PB_raw.fofn"

            # Discover GFF outputs directly in jobs directory
            gff_files = sorted([str(p) for p in (self.jobs_dir).glob('*.gff')])
            valid_files = []
            for gf in gff_files:
                if os.path.exists(gf) and os.path.getsize(gf) > 0:
                    valid_files.append(os.path.abspath(gf))
                else:
                    logger.warning(f"Skipping empty/missing GFF: {gf}")

            if not valid_files:
                logger.error("No Augustus PB GFF outputs found in jobs directory.")
                return None

            # Write FOFN with priority column for joingenes and raw concatenated GTF
            with open(raw_gtf_fofn, 'w') as fofn, open(raw_gtf, 'w') as raw:
                for gf in valid_files:
                    fofn.write(f"{gf}\t1\n")
                    with open(gf, 'r') as inf:
                        raw.write(inf.read())

            logger.info(f"Created raw Augustus PB output: {raw_gtf}")
            return str(raw_gtf)

        except Exception as e:
            logger.error(f"Error merging Augustus output: {e}")
            return None

    def save_intermediate_outputs(self) -> bool:
        """Optionally save intermediate PB outputs (GFF, ERR, joblist, job files)."""
        try:
            intermediate_dir = self.temp_dir / "intermediate_PB_outputs"
            intermediate_dir.mkdir(exist_ok=True)
            # Save GFFs
            for p in self.jobs_dir.glob('*.gff'):
                if p.exists() and p.is_file():
                    shutil.copy2(p, intermediate_dir / p.name)
            # Save job scripts and errors
            for p in self.jobs_dir.glob('*'):
                if p.suffix in ('.err', '.out') and p.is_file():
                    shutil.copy2(p, intermediate_dir / p.name)
            # Save job list
            jl = self.temp_dir / "jobs_PB.lst"
            if jl.exists():
                shutil.copy2(jl, intermediate_dir / jl.name)
            logger.info(f"Saved PB intermediate outputs to: {intermediate_dir}")
            return True
        except Exception as e:
            logger.warning(f"Could not save PB intermediate outputs: {e}")
            return False
    
    def run_joingenes(self, raw_gtf: str) -> str:
        """Run joingenes to merge overlapping predictions."""
        logger.info("Running joingenes to merge overlapping predictions...")
        
        try:
            join_genes_file = self.temp_dir / "augustus_PB_joined.gtf"
            join_genes_gp = self.temp_dir / "augustus_PB_joined.gp"
            raw_gtf_fofn = self.temp_dir / "augustus_PB_raw.fofn"
            
            # First pass: run joingenes
            tmp_join_genes_file = self.temp_dir / "augustus_PB_tmp_joined.gtf"
            # Ensure FOFN exists and is non-empty
            if (not os.path.exists(raw_gtf_fofn)) or os.path.getsize(raw_gtf_fofn) == 0:
                raise RuntimeError("FOFN for joingenes is missing or empty; cannot perform gene joining.")

            cmd = f'joingenes -f {raw_gtf_fofn} -o {tmp_join_genes_file}'
            
            logger.info(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, check=True)
            
            # Filter and format the joingenes output
            def filter_joingenes(injoingenes_file, out_joingenes_file):
                """Filter and format joingenes output."""
                import re
                matcher = re.compile("\tAUGUSTUS\t(exon|CDS|start_codon|stop_codon|tts|tss)\t")
                lines_written = 0
                with open(out_joingenes_file, "w") as ofh:
                    for l in open(injoingenes_file):
                        if matcher.search(l):
                            l = l.replace("jg", "augPB-")
                            ofh.write(l)
                            lines_written += 1
                return lines_written
            
            filtered_lines = filter_joingenes(tmp_join_genes_file, join_genes_file)
            if filtered_lines == 0:
                raise RuntimeError("joingenes produced zero final gene predictions after filtering.")
            logger.info(f'Joingenes produced {filtered_lines} final gene predictions')
            
            # Convert to GenePred format and back to GTF for proper formatting
            cmd = f'gtfToGenePred -genePredExt {join_genes_file} {join_genes_gp}'
            logger.info(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, check=True)
            
            # Convert back to GTF with proper formatting
            cmd = f'genePredToGtf file {join_genes_gp} -utr -honorCdsStat -source=augustusPB {join_genes_file}'
            logger.info(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, check=True)
            
            logger.info('Successfully completed Augustus PB gene joining and formatting')
            return str(join_genes_file)
            
        except Exception as e:
            logger.error(f"Error during gene joining: {e}")
            # Fail fast; do not proceed to downstream conversion on invalid inputs
            return None
    
    def process_final_output(self, final_gtf: str) -> bool:
        """Step 6: Process final output and create all required files."""
        logger.info("Step 6: Processing final Augustus PB output...")
        
        try:
            # Create output directories
            os.makedirs(os.path.dirname(self.args.gtf), exist_ok=True)
            os.makedirs(os.path.dirname(self.args.raw_gtf), exist_ok=True)
            os.makedirs(os.path.dirname(self.args.gp), exist_ok=True)
            
            # Copy final GTF to output location
            shutil.copy2(final_gtf, self.args.gtf)
            
            # Create raw GTF (copy from temp directory)
            raw_gtf_temp = self.temp_dir / "augustus_PB_raw.gtf"
            if raw_gtf_temp.exists():
                shutil.copy2(raw_gtf_temp, self.args.raw_gtf)
            
            # Convert GTF to GenePred
            gtf_to_gp_cmd = f'gtfToGenePred -genePredExt {self.args.gtf} {self.args.gp}'
            logger.info(f"Running: {gtf_to_gp_cmd}")
            subprocess.run(gtf_to_gp_cmd, shell=True, check=True)
            
            logger.info(f"Created final files:")
            logger.info(f"  GTF: {self.args.gtf}")
            logger.info(f"  Raw GTF: {self.args.raw_gtf}")
            logger.info(f"  GenePred: {self.args.gp}")
            
            return True
            
        except Exception as e:
            logger.error(f"Error processing final output: {e}")
            return False
    
    def cleanup_temp_files(self):
        """Clean up temporary files."""
        logger.info("Cleaning up temporary files...")
        
        try:
            if self.temp_dir.exists():
                shutil.rmtree(self.temp_dir)
            logger.info("Cleanup completed")
        except Exception as e:
            logger.warning(f"Error during cleanup: {e}")
    
    def run_pipeline(self) -> bool:
        """Run the complete parallel Augustus PB pipeline (mirrors TM pipeline stages)."""
        logger.info("Starting parallel Augustus PB pipeline...")

        try:
            if not self.validate_inputs():
                return False

            # Step 1: Run preprocessing (either on SLURM or locally)
            joblist_job_id = None  # Will be set if using SLURM preprocessing
            
            if not getattr(self.args, 'use_slurm_preprocessing', True):
                logger.info("Step 1: Running preprocessing steps locally...")
                if not self.split_genome_by_chromosome():
                    return False
                
                if not self.create_hints_location_file():
                    return False
                
                if not self.create_chromosome_hints_files():
                    return False
                
                if not self.create_augustus_jobs():
                    return False
            else:
                logger.info("Step 1: Running preprocessing steps on SLURM (3-stage parallel pipeline)...")
                
                # Stage 1: Initial setup (genome splitting + chr.lst creation)
                logger.info("Stage 1a: Genome splitting and chromosome list creation...")
                setup_script = self.create_initial_setup_slurm_script()
                setup_job_id = self.submit_slurm_job(setup_script, "PB initial setup")
                
                if not setup_job_id:
                    logger.error("Failed to submit PB initial setup job")
                    return False
                
                # Wait for setup to complete
                if not self.wait_for_slurm_job(setup_job_id, "PB initial setup", check_interval=30):
                    logger.error("PB initial setup job failed")
                    return False
                
                # Verify setup outputs and count chromosomes
                logger.info("Stage 1b: Verifying setup outputs and counting chromosomes...")
                chromosomes_file = self.temp_dir / "chromosomes.txt"
                expected_files = [
                    self.temp_dir / "chr.lst",
                    self.temp_dir / "summary.out",
                    chromosomes_file
                ]
                
                missing_files = [f for f in expected_files if not f.exists()]
                if missing_files:
                    logger.error("Setup completed but expected files are missing:")
                    for f in missing_files:
                        logger.error(f"  - {f}")
                    
                    # Check for SLURM output files
                    slurm_out_files = list(self.temp_dir.glob("augPB_setup_*.out"))
                    slurm_err_files = list(self.temp_dir.glob("augPB_setup_*.err"))
                    
                    if slurm_out_files:
                        logger.error(f"\nCheck SLURM output file: {slurm_out_files[0]}")
                    if slurm_err_files:
                        logger.error(f"Check SLURM error file: {slurm_err_files[0]}")
                    
                    return False
                
                # Count chromosomes for array job
                with open(chromosomes_file, 'r') as f:
                    chromosomes = [line.strip() for line in f if line.strip()]
                num_chromosomes = len(chromosomes)
                logger.info(f"Found {num_chromosomes} chromosomes to process: {chromosomes}")
                
                # Stage 2: Parallel hints generation (one job per chromosome)
                logger.info(f"Stage 2: Submitting parallel hints generation for {num_chromosomes} chromosomes...")
                hints_script = self.create_hints_generation_slurm_script(num_chromosomes, setup_job_id)
                hints_job_id = self.submit_slurm_job(hints_script, "PB hints generation array")
                
                if not hints_job_id:
                    logger.error("Failed to submit PB hints generation array job")
                    return False
                
                # Wait for all hints jobs to complete
                if not self.wait_for_slurm_job(hints_job_id, "PB hints generation array", check_interval=30):
                    logger.error("PB hints generation array job failed")
                    return False
                
                # Verify hints files were created
                logger.info("Verifying hints files were created...")
                missing_hints = []
                for chrom in chromosomes:
                    hints_file = self.temp_dir / f"{chrom}_hints.gff"
                    if not hints_file.exists():
                        missing_hints.append(chrom)
                        logger.error(f"Missing hints file for chromosome: {chrom}")
                
                if missing_hints:
                    logger.error(f"Hints generation failed for {len(missing_hints)} chromosomes: {missing_hints}")
                    # Show some error logs
                    slurm_err_files = list(self.temp_dir.glob("augPB_hints_*.err"))
                    if slurm_err_files:
                        logger.error(f"\nCheck SLURM error files in: {self.temp_dir}")
                        for err_file in slurm_err_files[:3]:  # Show first 3
                            logger.error(f"  - {err_file.name}")
                    return False
                
                logger.info(f"Successfully created hints files for all {num_chromosomes} chromosomes")
                
                # Stage 3: Job list generation
                logger.info("Stage 3: Creating Augustus job lists...")
                joblist_script = self.create_joblist_generation_slurm_script(hints_job_id)
                joblist_job_id = self.submit_slurm_job(joblist_script, "PB joblist generation")
                
                if not joblist_job_id:
                    logger.error("Failed to submit PB joblist generation job")
                    return False
                
                # Wait for joblist generation to complete
                if not self.wait_for_slurm_job(joblist_job_id, "PB joblist generation", check_interval=30):
                    logger.error("PB joblist generation job failed")
                    return False
                
                # Verify joblist outputs
                logger.info("Verifying joblist outputs...")
                expected_files = [self.temp_dir / "jobs_PB.lst"]
                
                missing_files = [f for f in expected_files if not f.exists()]
                if missing_files:
                    logger.error("Joblist generation completed but expected files are missing:")
                    for f in missing_files:
                        logger.error(f"  - {f}")
                    
                    # Check for SLURM output files
                    slurm_out_files = list(self.temp_dir.glob("augPB_joblist_*.out"))
                    slurm_err_files = list(self.temp_dir.glob("augPB_joblist_*.err"))
                    
                    if slurm_out_files:
                        logger.error(f"\nCheck SLURM output file: {slurm_out_files[0]}")
                    if slurm_err_files:
                        logger.error(f"Check SLURM error file: {slurm_err_files[0]}")
                    
                    return False
                
                logger.info("All preprocessing stages completed successfully!")

            # Step 2: Run Augustus jobs
            logger.info("Step 2: Running Augustus PB jobs...")
            
            # Count number of jobs
            jobs_file = self.temp_dir / "jobs_PB.lst"
            if not jobs_file.exists():
                logger.error(f"jobs_PB.lst not found: {jobs_file}")
                return False
            
            with open(jobs_file, 'r') as f:
                num_jobs = len([line for line in f if line.strip()])
            
            # Execute Augustus jobs
            if getattr(self.args, 'no_slurm_jobs', False):
                logger.info(f"Running {num_jobs} Augustus PB jobs locally...")
                if not self.run_local_jobs(str(jobs_file), getattr(self.args, 'num_cpus', None)):
                    return False
            else:
                logger.info(f"Running {num_jobs} Augustus PB jobs on SLURM...")
                # Use the last preprocessing job ID as dependency (joblist_job_id if SLURM, None if local)
                last_preprocessing_job_id = joblist_job_id if getattr(self.args, 'use_slurm_preprocessing', True) else None
                slurm_script = self.create_slurm_script(num_jobs, last_preprocessing_job_id)
                
                if not self.run_slurm_jobs(slurm_script):
                    return False

            # Save intermediate outputs for debugging
            if getattr(self.args, 'save_intermediate', True):
                logger.info("Saving intermediate PB outputs for debugging...")
                self.save_intermediate_outputs()
            else:
                logger.info("Skipping intermediate PB output saving (disabled)")

            # Step 3: Merge and post-process
            logger.info("Step 3: Merging and post-processing Augustus PB outputs...")
            raw_gtf = self.merge_augustus_output()
            if not raw_gtf:
                return False
            
            final_gtf = self.run_joingenes(raw_gtf)
            if not final_gtf:
                return False
            
            if not self.process_final_output(final_gtf):
                return False

            logger.info("Parallel Augustus PB pipeline completed successfully!")
            return True

        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            return False
        finally:
            if not getattr(self.args, 'keep_temp', True):
                self.cleanup_temp_files()


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Parallel Augustus PB Pipeline for CAT")
    
    # Input Files
    parser.add_argument("--genome_fasta", required=True, 
                       help="Genome FASTA file for Augustus gene prediction.")
    parser.add_argument("--chrom_sizes", required=True, 
                       help="Chromosome sizes file (tab-separated: chrom_name<tab>size).")
    parser.add_argument("--hints_gff", required=True, 
                       help="GFF file of extrinsic hints (will be filtered for src=PB only).")
    parser.add_argument("--pb_cfg", required=True, 
                       help="Augustus configuration file optimized for PacBio mode.")
    
    # Output Files
    parser.add_argument("--raw_gtf", required=True, 
                       help="Output path for raw (pre-joingenes) GTF predictions.")
    parser.add_argument("--gtf", required=True, 
                       help="Output path for final (post-joingenes) GTF with merged overlapping genes.")
    parser.add_argument("--gp", required=True, 
                       help="Output path for final GenePred format file.")
    
    # Augustus Parameters
    parser.add_argument("--species", required=True, 
                       help="Species parameter for Augustus (e.g., 'human', 'mouse', 'fly').")
    parser.add_argument("--utr", type=int, required=True, choices=[0, 1], 
                       help="UTR prediction parameter for Augustus (0=no UTRs, 1=predict UTRs).")
    
    # Chunking Parameters
    parser.add_argument("--chunksize", type=int, default=11000000, 
                       help="Size of genomic chunks for parallel processing (default: 11000000 bp).")
    parser.add_argument("--overlap", type=int, default=1000000, 
                       help="Overlap between genomic chunks (default: 1000000 bp).")
    
    # Pipeline options
    parser.add_argument("--work_dir", required=True,
                       help="Working directory for temporary files.")
    parser.add_argument("--keep_temp", action="store_true", default=True,
                       help="Keep temporary files for debugging (default: True).")
    parser.add_argument("--no_keep_temp", action="store_true",
                       help="Clean up temporary files after completion.")
    parser.add_argument("--use_slurm_preprocessing", action="store_true", default=True,
                       help="Use SLURM for preprocessing steps (genome splitting, hints creation, job list generation). Default: True")
    parser.add_argument("--no_slurm_preprocessing", action="store_true",
                       help="Disable SLURM preprocessing and run preprocessing steps locally.")
    parser.add_argument("--save_intermediate", action="store_true", default=True,
                       help="Save intermediate Augustus outputs for debugging. Default: True")
    parser.add_argument("--no_save_intermediate", action="store_true",
                       help="Disable saving intermediate Augustus outputs.")
    
    # Execution mode
    parser.add_argument("--no_slurm_jobs", action="store_true",
                       help="Run Augustus jobs locally instead of using SLURM (default: use SLURM).")
    parser.add_argument("--num_cpus", type=int, default=None,
                       help="Number of CPUs to use for local job execution (default: all available CPUs).")
    # SLURM resource configuration (applied when using SLURM mode)
    parser.add_argument("--slurm_partition", default="medium",
                       help="SLURM partition for preprocessing steps (setup/hints/joblist). Default: medium.")
    parser.add_argument("--slurm_jobs_partition", default="long",
                       help="SLURM partition for Augustus PB execution array jobs. Default: long.")
    parser.add_argument("--slurm_hints_mem", default="8G",
                       help="Memory per hints generation SLURM task. Default: 8G.")
    parser.add_argument("--slurm_jobs_mem", default="32G",
                       help="Memory per Augustus PB execution SLURM task. Default: 32G.")
    parser.add_argument("--slurm_setup_time", default="04:00:00",
                       help="Time limit for setup and joblist SLURM jobs. Default: 04:00:00.")
    parser.add_argument("--slurm_hints_time", default="04:00:00",
                       help="Time limit for hints generation SLURM array jobs. Default: 04:00:00.")
    parser.add_argument("--slurm_jobs_time", default="24:00:00",
                       help="Time limit for Augustus PB execution SLURM array jobs. Default: 24:00:00.")
    parser.add_argument("--slurm_hints_concurrency", type=int, default=10,
                       help="Max concurrent hints generation array tasks. Default: 10.")
    parser.add_argument("--slurm_jobs_concurrency", type=int, default=10,
                       help="Max concurrent Augustus PB execution array tasks. Default: 10.")
    parser.add_argument("--slurm_exclude_nodes", default="",
                       help="Node exclude list. SLURM: comma list; SGE: '!h1&!h2' or comma list (auto-converted).")
    # Cluster-backend options.
    parser.add_argument("--execution_mode", choices=("auto", "slurm", "sge", "local"), default="auto",
                       help="Cluster backend used for all submit/wait operations ('auto' detects it).")
    parser.add_argument("--module_load", default="",
                       help="Module to load at the top of every job script (or '').")
    parser.add_argument("--sge_parallel_env", default="smp",
                       help="SGE parallel environment name (site-specific; ignored on SLURM).")
    parser.add_argument("--sge_memory_flag", default="h_vmem",
                       help="SGE memory resource flag (h_vmem / mem_free / s_vmem).")

    args = parser.parse_args()
    
    # Handle SLURM preprocessing options
    if args.no_slurm_preprocessing:
        args.use_slurm_preprocessing = False
    
    # Handle intermediate output saving options
    if args.no_save_intermediate:
        args.save_intermediate = False
    
    # Handle no_keep_temp flag
    if args.no_keep_temp:
        args.keep_temp = False
    
    # Run pipeline
    pipeline = ParallelAugustusPB(args)
    success = pipeline.run_pipeline()
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
