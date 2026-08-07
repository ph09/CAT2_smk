#!/usr/bin/env python3
"""
Parallel Augustus Pipeline for CAT

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

# Import CAT tools (assume they're in the same directory structure)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tools.tm2hints
import tools.transcripts
import tools.psl
import tools.nameConversions
import tools.intervals
from tools.hintsDatabaseInterface import get_rnaseq_hints, reflect_hints_db

from cat.scheduler import get_scheduler


def _scheduler_from_args(args):
    """Build a Scheduler from the parsed CLI args.

    Honors --execution_mode (slurm/sge/local) plus the cluster-flavour knobs
    that each backend needs (parallel_env, memory_flag, hostname_exclude).
    """
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
                "memory_per_slot": getattr(args, "sge_memory_per_slot", True),
            },
        }
    }
    return get_scheduler(getattr(args, "execution_mode", "slurm"), cfg)

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def _run_augustus_job_script(job_path: str) -> tuple:
    """Run one Augustus job script (module-level for multiprocessing pickling)."""
    try:
        subprocess.run(['bash', job_path], capture_output=True, text=True, check=True)
        return (True, job_path, "")
    except subprocess.CalledProcessError as e:
        # Augustus usually logs to --errfile, not process stderr — surface that too.
        errfile_notes = []
        try:
            with open(job_path, 'r') as jf:
                for line in jf:
                    if '--errfile=' in line:
                        for part in line.split():
                            if part.startswith('--errfile='):
                                err_path = part.split('=', 1)[1]
                                if os.path.isfile(err_path) and os.path.getsize(err_path) > 0:
                                    with open(err_path, 'r') as ef:
                                        tail = ef.read()[-4000:]
                                    errfile_notes.append(f"ERRFILE ({err_path}):\n{tail}")
                                elif os.path.isfile(err_path):
                                    errfile_notes.append(f"ERRFILE ({err_path}): <empty>")
                                else:
                                    errfile_notes.append(f"ERRFILE missing: {err_path}")
        except OSError:
            pass
        error_msg = (
            f"exit={e.returncode}\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}\n"
            + ("\n".join(errfile_notes) if errfile_notes else "(no --errfile found in job script)")
        )
        return (False, job_path, error_msg)


class ParallelAugustus:
    """Main class for running parallel Augustus pipeline."""
    
    def __init__(self, args):
        self.args = args
        # Single scheduler instance shared across all submission sites in this pipeline.
        self._scheduler = _scheduler_from_args(args)

        # Normalize ALL paths to absolute paths to avoid SLURM relative-path issues
        self.args.genome_fasta = os.path.abspath(args.genome_fasta)
        self.args.coding_gp = os.path.abspath(args.coding_gp)
        self.args.filtered_tm_psl = os.path.abspath(args.filtered_tm_psl)
        self.args.ref_psl = os.path.abspath(args.ref_psl)
        self.args.annotation_gp = os.path.abspath(args.annotation_gp)
        self.args.tm_cfg = os.path.abspath(args.tm_cfg)
        self.args.miniprot_hints_gff = os.path.abspath(args.miniprot_hints_gff) if args.miniprot_hints_gff else None
        self.args.augustus_tm_gtf = os.path.abspath(args.augustus_tm_gtf)
        
        if args.augustus_tmr_gtf:
            self.args.augustus_tmr_gtf = os.path.abspath(args.augustus_tmr_gtf)
            self.args.tmr_cfg = os.path.abspath(args.tmr_cfg)
            self.args.augustus_hints_db = os.path.abspath(args.augustus_hints_db)
        
        # Normalize work_dir to an absolute path
        self.work_dir = Path(os.path.abspath(args.work_dir))
        self.temp_dir = self.work_dir / "augustus_parallel_temp"
        self.split_dir = self.temp_dir / "split"
        self.jobs_dir = self.temp_dir / "jobs"
        self.results_dir = self.temp_dir / "results"
        
        # Create directories
        for dir_path in [self.temp_dir, self.split_dir, self.jobs_dir, self.results_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

    def _wrap_cluster_header(self, header: str) -> str:
        """Append conda activation so cluster jobs use the active CAT env."""
        return (
            header
            + self._scheduler.script_preamble(
                conda_env=self._scheduler.resolve_conda_env(),
                set_strict=False,
            )
            + "\n"
        )

    def validate_inputs(self) -> bool:
        """Validate all required input files exist."""
        required_files = [
            self.args.genome_fasta,
            self.args.coding_gp,
            self.args.filtered_tm_psl,
            self.args.ref_psl,
            self.args.annotation_gp,
            self.args.tm_cfg
        ]
        
        # Add miniprot hints only if provided (augMP mode)
        if self.args.miniprot_hints_gff:
            required_files.append(self.args.miniprot_hints_gff)
        
        if self.args.augustus_tmr_gtf:
            required_files.extend([self.args.augustus_hints_db, self.args.tmr_cfg])
        
        for file_path in required_files:
            if not os.path.exists(file_path):
                logger.error(f"Required file not found: {file_path}")
                return False
        
        return True
    
    def split_genome_by_chromosome(self) -> bool:
        """Step 1: Split genome into chromosomes and get chromosome lengths."""
        logger.info("Step 1: Splitting genome by chromosome...")
        
        try:
            # Debug: Check input file exists and get size
            if not os.path.exists(self.args.genome_fasta):
                logger.error(f"Input genome file does not exist: {self.args.genome_fasta}")
                return False
            
            genome_size = os.path.getsize(self.args.genome_fasta)
            logger.info(f"Input genome file size: {genome_size} bytes")
            
            # Debug: Check if split directory exists, create if not
            if not self.split_dir.exists():
                logger.info(f"Creating split directory: {self.split_dir}")
                self.split_dir.mkdir(parents=True, exist_ok=True)
            else:
                logger.info(f"Split directory already exists: {self.split_dir}")
            
            # Run splitMfasta.pl to split genome
            split_cmd = [
                "splitMfasta.pl",
                self.args.genome_fasta,
                f"--outputpath={self.split_dir}"
            ]
            
            logger.info(f"Running: {' '.join(split_cmd)}")
            result = subprocess.run(split_cmd, capture_output=True, text=True, check=True)
            
            # Debug: Log splitMfasta.pl output
            logger.info(f"splitMfasta.pl stdout: {result.stdout}")
            if result.stderr:
                logger.info(f"splitMfasta.pl stderr: {result.stderr}")
            
            # Debug: List all files created in split directory
            split_files = list(self.split_dir.glob("*"))
            logger.info(f"Files created in split directory: {[f.name for f in split_files]}")
            
            # Rename files to use chromosome names
            renamed_count = 0
            for split_file in self.split_dir.glob("*.split.*"):
                # Extract chromosome name from FASTA header
                with open(split_file, 'r') as f:
                    first_line = f.readline().strip()
                    if first_line.startswith('>'):
                        chrom_name = first_line[1:].split()[0]  # Get chromosome name
                        new_name = f"{chrom_name}.fa"
                        new_path = self.split_dir / new_name
                        shutil.move(str(split_file), str(new_path))
                        logger.info(f"Renamed {split_file.name} to {new_name}")
                        renamed_count += 1
                        
                        # Debug: Check file size after rename
                        file_size = os.path.getsize(new_path)
                        logger.info(f"Chromosome {chrom_name} file size: {file_size} bytes")
            
            logger.info(f"Renamed {renamed_count} chromosome files")
            
            # Debug: List final chromosome files
            final_files = list(self.split_dir.glob("*.fa"))
            logger.info(f"Final chromosome files: {[f.name for f in final_files]}")
            
            # Generate ACGT content summary for each split chromosome
            summary_out = self.temp_dir / "summary.out"
            logger.info(f"Creating ACGT content summary: {summary_out}")
            
            with open(summary_out, 'w') as f:
                for split_file in self.split_dir.glob("*.fa"):
                    summary_cmd = ["summarizeACGTcontent.pl", str(split_file)]
                    logger.info(f"Running: {' '.join(summary_cmd)}")
                    summary_result = subprocess.run(summary_cmd, stdout=f, stderr=subprocess.PIPE, text=True, check=True)
                    if summary_result.stderr:
                        logger.info(f"summarizeACGTcontent.pl stderr for {split_file.name}: {summary_result.stderr}")
            
            # Debug: Check summary file was created and has content
            if summary_out.exists():
                summary_size = os.path.getsize(summary_out)
                logger.info(f"Summary file created with {summary_size} bytes")
                # Log first few lines of summary
                with open(summary_out, 'r') as f:
                    lines = f.readlines()[:5]
                    logger.info(f"Summary file preview: {lines}")
            else:
                logger.error("Summary file was not created!")
                return False
            
            logger.info("Genome splitting completed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error splitting genome: {e}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"Unexpected error in split_genome_by_chromosome: {e}")
            return False
    
    def create_hints_location_file(self) -> bool:
        """Step 2: Create hints location file for each chromosome."""
        logger.info("Step 2: Creating hints location file...")
        
        try:
            summary_file = self.temp_dir / "summary.out"
            chr_lst_file = self.temp_dir / "chr.lst"
            
            # Debug: Check if summary file exists
            if not summary_file.exists():
                logger.error(f"Summary file not found: {summary_file}")
                return False
            
            logger.info(f"Reading chromosome information from: {summary_file}")
            
            # Parse summary.out to create chr.lst
            chromosomes_found = 0
            with open(summary_file, 'r') as f_in, open(chr_lst_file, 'w') as f_out:
                for line_num, line in enumerate(f_in, 1):
                    logger.info(f"Processing line {line_num}: {line.strip()}")
                    
                    if "bases" in line:
                        # Handle multiple possible formats from summarizeACGTcontent.pl
                        # Format A: "123456789 bases in PAN011.chr7.haplotype1"
                        # Format B: "210155 bases.\tchr16 BASE COUNT ..."
                        parts = line.strip().split()
                        chrom_length = None
                        chrom_name = None

                        if len(parts) >= 4 and parts[1] == "bases" and parts[2] == "in":
                            # Format A
                            chrom_length = parts[0]
                            chrom_name = parts[3]
                            logger.info(f"Format A detected: {chrom_name} = {chrom_length} bases")
                        elif len(parts) >= 3 and parts[1].startswith("bases"):
                            # Format B (e.g., "bases.")
                            chrom_length = parts[0]
                            chrom_name = parts[2]
                            logger.info(f"Format B detected: {chrom_name} = {chrom_length} bases")
                        else:
                            logger.warning(f"Unexpected format in line {line_num}: {line.strip()}")

                        if chrom_length and chrom_name:
                            # Create hints file path (will be created later)
                            hints_file = self.temp_dir / f"{chrom_name}_hints.gff"
                            seq_file = self.split_dir / f"{chrom_name}.fa"
                            
                            # Debug: Check if sequence file exists
                            if not seq_file.exists():
                                logger.error(f"Sequence file not found: {seq_file}")
                                return False
                            
                            # Write to chr.lst format: chromosome_file\thints_file\t1\tlength
                            chr_line = f"{seq_file}\t{hints_file}\t1\t{chrom_length}\n"
                            f_out.write(chr_line)
                            logger.info(f"Added to chr.lst: {chr_line.strip()}")
                            chromosomes_found += 1
                        else:
                            logger.warning(f"Could not parse chromosome info from line {line_num}")
            
            logger.info(f"Found {chromosomes_found} chromosomes")
            
            # Debug: Verify chr.lst file was created
            if chr_lst_file.exists():
                chr_lst_size = os.path.getsize(chr_lst_file)
                logger.info(f"Created chromosome list with {chr_lst_size} bytes")
                # Log content of chr.lst
                with open(chr_lst_file, 'r') as f:
                    content = f.read()
                    logger.info(f"chr.lst content:\n{content}")
            else:
                logger.error("chr.lst file was not created!")
                return False
            
            logger.info(f"Created chromosome list: {chr_lst_file}")
            return True
            
        except Exception as e:
            logger.error(f"Error creating hints location file: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            return False
    
    def create_chromosome_hints_files(self) -> bool:
        """Create hints files for each chromosome."""
        logger.info("Creating chromosome-specific hints files...")
        
        try:
            # Use imported modules for hints processing
            
            # Load required data
            ref_psl_dict = tools.psl.get_alignment_dict(self.args.ref_psl)
            tm_psl_dict = tools.psl.get_alignment_dict(self.args.filtered_tm_psl)
            ref_tx_dict = tools.transcripts.get_gene_pred_dict(self.args.annotation_gp)
            tx_dict = tools.transcripts.get_gene_pred_dict(self.args.coding_gp)
            
            # augMP: Snakefile passes --miniprot_hints_gff only for augustus_run_mp (not augTM/TMR)
            is_augmp_mode = self.args.miniprot_hints_gff is not None
            
            # Load miniprot hints ONLY for augMP mode
            mp_hints = None
            if is_augmp_mode:
                with open(self.args.miniprot_hints_gff) as mpf:
                    mp_hints = mpf.read()
            
            # Load TMR hints if available (NOT for augMP mode - only for augTMR)
            tmr_hints_data = None
            if self.args.augustus_tmr_gtf and self.args.augustus_hints_db and not is_augmp_mode:
                try:
                    hints_db_file = self.args.augustus_hints_db
                    speciesnames, seqnames, hints, featuretypes, session = reflect_hints_db(hints_db_file)
                    tmr_hints_data = (speciesnames, seqnames, hints, featuretypes, session)
                except Exception as e:
                    logger.warning(f"Could not load TMR hints database: {e}")
            
            # Group transcripts by chromosome
            chrom_transcripts = {}
            for tx_id, tx in tx_dict.items():
                chrom = tx.chromosome
                if chrom not in chrom_transcripts:
                    chrom_transcripts[chrom] = []
                chrom_transcripts[chrom].append(tx_id)
            
            # Create hints for each chromosome
            for chrom, tx_ids in chrom_transcripts.items():
                hints_file = self.temp_dir / f"{chrom}_hints.gff"
                
                with open(hints_file, 'w') as outf:
                    for tx_id in tx_ids:
                        tx = tx_dict[tx_id]
                        ref_tx = ref_tx_dict.get(tools.nameConversions.remove_alignment_number(tx_id))
                        tm_psl = tm_psl_dict.get(tx_id)
                        ref_psl = ref_psl_dict.get(tools.nameConversions.remove_alignment_number(tx_id))
                        
                        if all([ref_tx, tm_psl, ref_psl]):
                            # Generate TM hints (only for augTM/TMR mode)
                            tm_hints = tools.tm2hints.tm_to_hints(tx, tm_psl, ref_psl)
                            outf.write(tm_hints + '\n')
                    
                    # Add miniprot hints (only for augMP mode)
                    if mp_hints:
                        outf.write(mp_hints)
                    
                    # Add TMR hints if available
                    if tmr_hints_data:
                        speciesnames, seqnames, hints, featuretypes, session = tmr_hints_data
                        # Get RNA-seq hints for this chromosome
                        # Note: This is simplified - in practice you might want to filter by position
                        rnaseq_hints = get_rnaseq_hints(
                            self.args.genome, chrom, 0, 999999999,  # Large range for now
                            speciesnames, seqnames, hints, featuretypes, session
                        )
                        outf.write(rnaseq_hints + '\n')
                
                logger.info(f"Created hints file for chromosome {chrom}: {hints_file}")

            # chr.lst lists a hints path for every genome chromosome; coding_gp may
            # only cover a subset (or be empty in sparse pairwise modes). Touch any
            # missing hints files so Augustus --hintsfile= does not fail open().
            chr_lst_file = self.temp_dir / "chr.lst"
            if chr_lst_file.exists():
                with open(chr_lst_file) as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            hints_path = Path(parts[1])
                            if not hints_path.exists():
                                hints_path.parent.mkdir(parents=True, exist_ok=True)
                                hints_path.touch()
                                logger.info(f"Created empty hints file (no transcripts): {hints_path}")
            
            if tmr_hints_data:
                session.close()
            
            return True
            
        except Exception as e:
            logger.error(f"Error creating chromosome hints files: {e}")
            return False
    
    def create_augustus_jobs(self, mode: str = "TM") -> bool:
        """Step 3: Create Augustus job list using createAugustusJoblist.pl."""
        logger.info(f"Step 3: Creating Augustus jobs for {mode} mode...")
        
        try:
            chr_lst_file = self.temp_dir / "chr.lst"
            jobs_lst_file = self.temp_dir / f"jobs_{mode}.lst"
            
            # Debug: Check if chr.lst file exists
            if not chr_lst_file.exists():
                logger.error(f"chr.lst file not found: {chr_lst_file}")
                return False
            
            # Debug: Check if jobs directory exists
            if not self.jobs_dir.exists():
                logger.info(f"Creating jobs directory: {self.jobs_dir}")
                self.jobs_dir.mkdir(parents=True, exist_ok=True)
            else:
                logger.info(f"Jobs directory exists: {self.jobs_dir}")
            
            # Choose configuration file based on mode
            cfg_file = self.args.tm_cfg if mode == "TM" else self.args.tmr_cfg
            
            # Debug: Check if config file exists
            if not os.path.exists(cfg_file):
                logger.error(f"Config file not found: {cfg_file}")
                return False
            
            logger.info(f"Using config file: {cfg_file}")
            
            # Build Augustus command with parameters that createAugustusJoblist.pl will handle
            # Note: createAugustusJoblist.pl automatically adds:
            # - {seqfile} (sequence file path)
            # - --predictionStart and --predictionEnd (based on chunking)
            # - --hintsfile (from chr.lst file)
            # - --outfile and --errfile (output and error files)
            aug_cmd = [
                "augustus",
                "--extrinsicCfgFile={}".format(cfg_file),
                "--UTR={}".format(int(self.args.utr)),
                "--alternatives-from-evidence=0",
                "--species={}".format(self.args.augustus_species),
                "--allow_hinted_splicesites=atac",
                "--protein=0",
                "--softmasking=1",
                "--/augustus/verbosity=0"
            ]
            
            aug_call = " ".join(aug_cmd)
            logger.info(f"Augustus command template: {aug_call}")
            
            # Create job list with genome-specific prefix to avoid conflicts
            create_jobs_cmd = [
                "createAugustusJoblist.pl",
                f"--sequences={chr_lst_file}",
                "--wrap=#",
                "--overlap=100000",
                "--chunksize=1100000",
                f"--outputdir={self.jobs_dir}",
                f"--joblist={jobs_lst_file}",
                f"--jobprefix={self.args.genome}_aug_{mode}_",
                f"--command={aug_call}"
            ]
            
            logger.info(f"Running: {' '.join(create_jobs_cmd)}")
            # createAugustusJoblist.pl --wrap writes job scripts into CWD (basenames
            # in the joblist). Run from temp_dir so SLURM/SGE (cd temp_dir; bash $JOB)
            # and local execution find the same files. --outputdir is only for .gff/.err.
            result = subprocess.run(
                create_jobs_cmd,
                capture_output=True,
                text=True,
                check=True,
                cwd=str(self.temp_dir),
            )
            
            # Debug: Log createAugustusJoblist.pl output
            logger.info(f"createAugustusJoblist.pl stdout: {result.stdout}")
            if result.stderr:
                logger.info(f"createAugustusJoblist.pl stderr: {result.stderr}")
            
            # Debug: Check if job list file was created
            if jobs_lst_file.exists():
                jobs_lst_size = os.path.getsize(jobs_lst_file)
                logger.info(f"Created job list with {jobs_lst_size} bytes")
                # Log content of job list
                with open(jobs_lst_file, 'r') as f:
                    job_files = f.readlines()
                    logger.info(f"Job list contains {len(job_files)} jobs:")
                    for i, job_file in enumerate(job_files, 1):
                        logger.info(f"  Job {i}: {job_file.strip()}")
            else:
                logger.error("Job list file was not created!")
                return False
            
            # Job scripts land in temp_dir (CWD), not jobs_dir
            job_files_created = list(self.temp_dir.glob(f"{self.args.genome}_aug_{mode}_*"))
            logger.info(f"Created {len(job_files_created)} job files in {self.temp_dir}")
            for job_file in job_files_created:
                file_size = os.path.getsize(job_file)
                logger.info(f"  {job_file.name}: {file_size} bytes")
            if not job_files_created:
                logger.error(
                    f"createAugustusJoblist.pl wrote jobs_{mode}.lst but no job scripts "
                    f"under {self.temp_dir}. Refusing to submit empty Augustus array."
                )
                return False
            
            logger.info(f"Created job list: {jobs_lst_file}")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating Augustus jobs: {e}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"Unexpected error in create_augustus_jobs: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            return False
    
    def create_initial_setup_slurm_script(self) -> str:
        """Create SLURM script for initial setup (genome splitting and chr.lst creation)."""
        logger.info("Creating SLURM script for initial setup (genome splitting and chr.lst)...")
        
        setup_mem = (
            getattr(self.args, 'slurm_setup_mem', None)
            or getattr(self.args, 'slurm_jobs_mem', '16G')
        )
        log_out, log_err = self._scheduler.job_log_paths(self.temp_dir, "augustus_setup")
        header = self._scheduler.header(
            job_name="augustus_setup",
            cpus=1,
            mem=setup_mem,
            walltime=getattr(self.args, 'slurm_setup_time', '04:00:00'),
            log_out=log_out,
            log_err=log_err,
            partition=getattr(self.args, 'slurm_partition', ''),
            queue=getattr(self.args, 'slurm_partition', ''),
        )
        slurm_script = self._wrap_cluster_header(header) + rf"""
# Enable strict error handling
set -e          # Exit immediately if a command exits with a non-zero status
set -u          # Treat unset variables as an error
set -o pipefail # Return value of a pipeline is the value of the last command to exit with a non-zero status

echo "Starting initial setup..."
echo "Working directory: {self.temp_dir}"
echo "Split directory: {self.split_dir}"

# Step 1: Split genome by chromosome
echo "Step 1: Splitting genome by chromosome..."
splitMfasta.pl {self.args.genome_fasta} --outputpath={self.split_dir}

# Rename files to use chromosome names
for split_file in {self.split_dir}/*.split.*; do
    if [ -f "$split_file" ]; then
        # Extract chromosome name from FASTA header
        chrom_name=$(head -n1 "$split_file" | sed 's/^>//' | awk '{{print $1}}')
        new_name="{self.split_dir}/${{chrom_name}}.fa"
        mv "$split_file" "$new_name"
        echo "Renamed $(basename $split_file) to $(basename $new_name)"
    fi
done

# Generate ACGT content summary for each split chromosome
summary_out="{self.temp_dir}/summary.out"
for split_file in {self.split_dir}/*.fa; do
    if [ -f "$split_file" ]; then
        echo "Running summarizeACGTcontent.pl on $(basename $split_file)"
        summarizeACGTcontent.pl "$split_file" >> "$summary_out"
    fi
done

echo "Genome splitting completed successfully"

# Step 2: Create hints location file (initial version without hints files)
echo "Step 2: Creating chromosome list..."
chr_lst_file="{self.temp_dir}/chr.lst"
chrom_list="{self.temp_dir}/chromosomes.txt"

# Parse summary.out to create chr.lst and chromosome list
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
echo "Initial setup completed successfully!"
"""

        slurm_script_path = self.temp_dir / "augustus_setup.slurm"
        with open(slurm_script_path, 'w') as f:
            f.write(slurm_script)
        
        logger.info(f"Created initial setup SLURM script: {slurm_script_path}")
        return str(slurm_script_path)
    
    def create_hints_generation_slurm_script(self, num_chromosomes: int, dependency_job_id: str = None) -> str:
        """Create cluster array script for parallel chromosome hints generation.

        Already 1-based; the body keeps using SLURM_ARRAY_TASK_ID for legacy
        log readability but reads from a backend-agnostic ${{TASK_ID}} variable
        so this works on SGE too.
        """
        logger.info(f"Creating cluster array script for hints generation ({num_chromosomes} chromosomes)...")

        dependency = self._scheduler.depends_on_job_id(dependency_job_id) if dependency_job_id else None
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        task_var = self._scheduler.task_id_env()

        log_out, log_err = self._scheduler.array_log_paths(self.temp_dir, "augustus_hints")
        header = self._scheduler.header(
            job_name="augustus_hints",
            cpus=1,
            mem=getattr(self.args, 'slurm_hints_mem', '64G'),
            walltime=getattr(self.args, 'slurm_hints_time', '04:00:00'),
            log_out=log_out,
            log_err=log_err,
            partition=getattr(self.args, 'slurm_partition', ''),
            queue=getattr(self.args, 'slurm_partition', ''),
            array=(1, num_chromosomes),
            max_concurrent=getattr(self.args, 'slurm_hints_concurrency', 10),
            dependency=dependency,
        )
        slurm_script = self._wrap_cluster_header(header) + rf"""
# Enable strict error handling
set -euo pipefail

TASK_ID="${{{task_var}}}"
CHROM=$(sed -n "${{TASK_ID}}p" {self.temp_dir}/chromosomes.txt)

echo "Processing chromosome: $CHROM"
echo "Task ID: $TASK_ID"

# Set up Python environment (PYTHONPATH may be unset; guard for `set -u`)
export PYTHONPATH="{project_root}:${{PYTHONPATH:-}}"

# Create hints file for this chromosome using Python
python3 -c "
import sys
import os
sys.path.insert(0, '{project_root}')
import tools.tm2hints
import tools.transcripts
import tools.psl
import tools.nameConversions
from tools.hintsDatabaseInterface import get_rnaseq_hints, reflect_hints_db

chrom = '$CHROM'
print(f'Generating hints for chromosome: {{chrom}}')

# Load required data
print('Loading alignment dictionaries...')
ref_psl_dict = tools.psl.get_alignment_dict('{self.args.ref_psl}')
tm_psl_dict = tools.psl.get_alignment_dict('{self.args.filtered_tm_psl}')
ref_tx_dict = tools.transcripts.get_gene_pred_dict('{self.args.annotation_gp}')
tx_dict = tools.transcripts.get_gene_pred_dict('{self.args.coding_gp}')

print(f'Loaded {{len(tx_dict)}} transcripts')

# augMP when --miniprot_hints_gff was passed (matches local create_chromosome_hints_files)
is_augmp_mode = {self.args.miniprot_hints_gff is not None}

# Load miniprot hints ONLY for augMP mode
mp_hints = None
if is_augmp_mode:
    print('Loading miniprot hints for augMP mode...')
    with open('{self.args.miniprot_hints_gff}') as mpf:
        mp_hints_all = mpf.read()
    
    # Filter miniprot hints for this chromosome
    mp_hints_chrom = []
    for line in mp_hints_all.split('\\n'):
        if line.strip() and not line.startswith('#'):
            parts = line.split('\\t')
            if len(parts) > 0 and parts[0] == chrom:
                mp_hints_chrom.append(line)
    mp_hints = '\\n'.join(mp_hints_chrom)

# Load TMR hints if available (NOT for augMP mode - only for augTMR)
tmr_hints_data = None
tmr_gtf = '{self.args.augustus_tmr_gtf if self.args.augustus_tmr_gtf else ""}'
hints_db = '{self.args.augustus_hints_db if self.args.augustus_hints_db else ""}'
if tmr_gtf and hints_db and not is_augmp_mode:
    try:
        print('Loading TMR hints database...')
        hints_db_file = hints_db
        speciesnames, seqnames, hints, featuretypes, session = reflect_hints_db(hints_db_file)
        tmr_hints_data = (speciesnames, seqnames, hints, featuretypes, session)
    except Exception as e:
        print(f'Could not load TMR hints database: {{e}}')

# Get transcripts for this chromosome
chrom_tx_ids = [tx_id for tx_id, tx in tx_dict.items() if tx.chromosome == chrom]
print(f'Found {{len(chrom_tx_ids)}} transcripts for chromosome {{chrom}}')

# Create hints file for this chromosome
hints_file = '{self.temp_dir}/' + chrom + '_hints.gff'
print(f'Writing hints to: {{hints_file}}')

with open(hints_file, 'w') as outf:
    # Process transcripts for this chromosome
    for tx_id in chrom_tx_ids:
        tx = tx_dict[tx_id]
        ref_tx = ref_tx_dict.get(tools.nameConversions.remove_alignment_number(tx_id))
        tm_psl = tm_psl_dict.get(tx_id)
        ref_psl = ref_psl_dict.get(tools.nameConversions.remove_alignment_number(tx_id))
        
        if all([ref_tx, tm_psl, ref_psl]):
            # Generate TM hints
            tm_hints = tools.tm2hints.tm_to_hints(tx, tm_psl, ref_psl)
            outf.write(tm_hints + '\\n')
    
    # Add miniprot hints for this chromosome
    if mp_hints:
        outf.write(mp_hints + '\\n')
    
    # Add TMR hints if available
    if tmr_hints_data:
        print('Adding TMR hints...')
        speciesnames, seqnames, hints, featuretypes, session = tmr_hints_data
        # Get RNA-seq hints for this chromosome
        rnaseq_hints = get_rnaseq_hints(
            '{self.args.genome}', chrom, 0, 999999999,  # Large range
            speciesnames, seqnames, hints, featuretypes, session
        )
        if rnaseq_hints:
            outf.write(rnaseq_hints + '\\n')
        session.close()

print(f'Successfully created hints file for chromosome {{chrom}}')
"

echo "Completed hints generation for chromosome: $CHROM"
"""

        slurm_script_path = self.temp_dir / "augustus_hints_array.slurm"
        with open(slurm_script_path, 'w') as f:
            f.write(slurm_script)
        
        logger.info(f"Created hints generation SLURM array script: {slurm_script_path}")
        return str(slurm_script_path)
    
    def create_joblist_generation_slurm_script(self, dependency_job_id: str = None) -> str:
        """Create cluster job script for Augustus joblist generation."""
        logger.info("Creating cluster script for Augustus job list generation...")

        dependency = self._scheduler.depends_on_job_id(dependency_job_id) if dependency_job_id else None

        setup_mem = (
            getattr(self.args, 'slurm_setup_mem', None)
            or getattr(self.args, 'slurm_jobs_mem', '16G')
        )
        log_out, log_err = self._scheduler.job_log_paths(self.temp_dir, "augustus_joblist")
        header = self._scheduler.header(
            job_name="augustus_joblist",
            cpus=1,
            mem=setup_mem,
            walltime=getattr(self.args, 'slurm_setup_time', '04:00:00'),
            log_out=log_out,
            log_err=log_err,
            partition=getattr(self.args, 'slurm_partition', ''),
            queue=getattr(self.args, 'slurm_partition', ''),
            dependency=dependency,
        )
        slurm_script = self._wrap_cluster_header(header) + rf"""
# Enable strict error handling
set -euo pipefail

echo "Creating Augustus job lists..."

# Set unique temporary directory to avoid conflicts between parallel runs
export TMPDIR="{self.temp_dir}/tmp_$$"
mkdir -p "$TMPDIR"

chr_lst_file="{self.temp_dir}/chr.lst"

# Debug: Verify chr.lst content
echo "Contents of chr.lst:"
cat "$chr_lst_file"
echo "---"

# TM mode
echo "Creating TM job list..."
cfg_file="{self.args.tm_cfg}"
aug_cmd="augustus --extrinsicCfgFile=${{cfg_file}} --UTR={int(self.args.utr)} --alternatives-from-evidence=0 --species={self.args.augustus_species} --allow_hinted_splicesites=atac --protein=0 --softmasking=1 --/augustus/verbosity=0"

# Change to the temp directory to ensure createAugustusJoblist.pl uses local paths
cd "{self.temp_dir}"

createAugustusJoblist.pl --sequences="$chr_lst_file" --wrap="#" --overlap=100000 --chunksize=1100000 --outputdir="{self.jobs_dir}" --joblist="{self.temp_dir}/jobs_TM.lst" --jobprefix="{self.args.genome}_aug_TM_" --command="$aug_cmd"

echo "Created TM job list: {self.temp_dir}/jobs_TM.lst"

# Debug: Show first job file content and verify it exists
if [ -f "{self.temp_dir}/jobs_TM.lst" ]; then
    echo "Job list content:"
    cat "{self.temp_dir}/jobs_TM.lst"
    first_job=$(head -1 "{self.temp_dir}/jobs_TM.lst")
    echo "First TM job file path: $first_job"
    if [ -f "$first_job" ]; then
        echo "First TM job exists with size: $(wc -c < "$first_job") bytes"
        echo "First TM job content (first 20 lines):"
        head -20 "$first_job"
    else
        echo "ERROR: First TM job file does not exist: $first_job"
        echo "Checking jobs directory:"
        ls -la "{self.jobs_dir}/"
    fi
else
    echo "ERROR: jobs_TM.lst was not created"
fi

# TMR mode (if requested)
if [[ -n "{self.args.augustus_tmr_gtf}" && -n "{self.args.tmr_cfg}" ]]; then
    echo "Creating TMR job list..."
    tmr_cfg_file="{self.args.tmr_cfg}"
    tmr_aug_cmd="augustus --extrinsicCfgFile=${{tmr_cfg_file}} --UTR={int(self.args.utr)} --alternatives-from-evidence=0 --species={self.args.augustus_species} --allow_hinted_splicesites=atac --protein=0 --softmasking=1 --/augustus/verbosity=0"
    
    createAugustusJoblist.pl --sequences="$chr_lst_file" --wrap="#" --overlap=100000 --chunksize=1100000 --outputdir="{self.jobs_dir}" --joblist="{self.temp_dir}/jobs_TMR.lst" --jobprefix="{self.args.genome}_aug_TMR_" --command="$tmr_aug_cmd"
    
    echo "Created TMR job list: {self.temp_dir}/jobs_TMR.lst"
    
    # Debug TMR job list
    if [ -f "{self.temp_dir}/jobs_TMR.lst" ]; then
        echo "TMR job list content:"
        cat "{self.temp_dir}/jobs_TMR.lst"
        first_tmr_job=$(head -1 "{self.temp_dir}/jobs_TMR.lst")
        echo "First TMR job file path: $first_tmr_job"
        if [ -f "$first_tmr_job" ]; then
            echo "First TMR job exists with size: $(wc -c < "$first_tmr_job") bytes"
        else
            echo "ERROR: First TMR job file does not exist: $first_tmr_job"
        fi
    fi
fi

# Cleanup temporary directory
rm -rf "$TMPDIR"

echo "Job list generation completed successfully!"
"""

        slurm_script_path = self.temp_dir / "augustus_joblist.slurm"
        with open(slurm_script_path, 'w') as f:
            f.write(slurm_script)
        
        logger.info(f"Created job list generation SLURM script: {slurm_script_path}")
        return str(slurm_script_path)

    def create_slurm_script(self, mode: str = "TM", num_jobs: int = 30, dependency_job_id: str = None) -> str:
        """Step 4: Create cluster array job script for parallel execution."""
        logger.info(f"Step 4: Creating cluster script for {mode} mode...")

        dependency = self._scheduler.depends_on_job_id(dependency_job_id) if dependency_job_id else None
        task_var = self._scheduler.task_id_env()

        log_out, log_err = self._scheduler.array_log_paths(self.temp_dir, f"augustus_{mode}")
        header = self._scheduler.header(
            job_name=f"augustus_{mode}",
            cpus=1,
            mem=getattr(self.args, 'slurm_jobs_mem', '16G'),
            walltime=getattr(self.args, 'slurm_jobs_time', '01:00:00'),
            log_out=log_out,
            log_err=log_err,
            partition=getattr(self.args, 'slurm_jobs_partition', ''),
            queue=getattr(self.args, 'slurm_jobs_partition', ''),
            array=(1, num_jobs),
            max_concurrent=getattr(self.args, 'slurm_jobs_concurrency', 100),
            dependency=dependency,
        )
        slurm_script = self._wrap_cluster_header(header) + f"""
# Change to temp directory (job paths in the list are relative to this directory)
cd {self.temp_dir}

# Get the job file for this array task (1-based — works on both SLURM and SGE).
TASK_ID="${{{task_var}}}"
JOB_FILE=$(sed -n "${{TASK_ID}}p" jobs_{mode}.lst)

# Run the Augustus job
if [ -f "$JOB_FILE" ]; then
    echo "Running Augustus job: $JOB_FILE"
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
        
        slurm_script_path = self.temp_dir / f"augustus_{mode}.slurm"
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
        """Wait for a cluster job to complete and verify it succeeded.

        Two-stage approach (works for both SLURM and SGE):
          1. Poll ``job_present()`` until the job has drained from the queue.
             This catches the case where the job hangs in PENDING with
             ``DependencyNeverSatisfied`` (we must time out) by treating
             "still pending" as still present.
          2. Ask the scheduler to validate the final state via
             ``verify_completed()`` (SLURM: ``sacct``; SGE: no-op since
             ``qacct`` is non-portable; Local: no-op).

        Returns False on failure so the caller can abort the dependency chain
        instead of pushing another job that will instantly land in
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

    def _should_run_jobs_locally(self) -> bool:
        """True when Augustus chunk jobs must run in-process (not via sbatch/qsub)."""
        mode = getattr(self.args, "execution_mode", "auto")
        if mode == "local" or getattr(self.args, "no_slurm_jobs", False):
            return True
        # Snakemake local path sets --no_slurm_preprocessing; keep jobs local too
        # so auto→sge detection cannot qsub orphan arrays without PATH/env.
        if not getattr(self.args, "use_slurm_preprocessing", True):
            return True
        return False

    def run_local_jobs(self, jobs_file: str, num_cpus: int = None) -> bool:
        """Run Augustus job scripts locally (same semantics as SLURM array body)."""
        import multiprocessing

        logger.info("Running Augustus jobs locally...")
        try:
            list_dir = os.path.dirname(os.path.abspath(jobs_file))
            with open(jobs_file, 'r') as f:
                job_entries = [line.strip() for line in f if line.strip()]

            jobs = []
            for entry in job_entries:
                if os.path.isabs(entry) and os.path.isfile(entry):
                    jobs.append(entry)
                    continue
                candidates = [
                    os.path.join(list_dir, entry),
                    os.path.join(str(self.temp_dir), entry),
                    os.path.join(str(self.jobs_dir), entry),
                    os.path.join(str(self.temp_dir), os.path.basename(entry)),
                    os.path.join(str(self.jobs_dir), os.path.basename(entry)),
                ]
                resolved = next((p for p in candidates if os.path.isfile(p)), None)
                if resolved is None:
                    logger.error(
                        f"Job script not found for list entry {entry!r}; tried: {candidates}"
                    )
                    return False
                jobs.append(resolved)

            if num_cpus is None:
                num_cpus = getattr(self.args, "num_cpus", None) or multiprocessing.cpu_count()
            num_cpus = max(1, int(num_cpus))
            logger.info(f"Total jobs to run: {len(jobs)}")
            logger.info(f"Using {num_cpus} CPUs for parallel execution")

            with multiprocessing.Pool(processes=num_cpus) as pool:
                results = pool.map(_run_augustus_job_script, jobs)

            failed = [(job, err) for ok, job, err in results if not ok]
            if failed:
                logger.error(f"{len(failed)} jobs failed:")
                for job, err in failed[:5]:
                    logger.error(f"Job: {job}\nError: {err}")
                return False
            logger.info(f"All {len(jobs)} jobs completed successfully")
            return True
        except Exception as e:
            logger.error(f"Error running local jobs: {e}")
            return False

    def run_slurm_jobs(self, slurm_script_path: str) -> bool:
        """Submit and wait for SLURM jobs to complete and verify success."""
        logger.info("Submitting SLURM jobs...")
        
        try:
            # Submit the job
            job_id = self._scheduler.submit(slurm_script_path)
            logger.info(f"Submitted {self._scheduler.name} job array: {job_id}")

            logger.info(f"Waiting for {self._scheduler.name} jobs to complete...")
            while True:
                if not self._scheduler.job_present(job_id):
                    logger.info(f"All {self._scheduler.name} jobs no longer in queue")
                    break
                time.sleep(60)

            # Per-task failure detection: SLURM reports via sacct (which
            # includes DependencyNeverSatisfied → never-COMPLETED parent
            # state). SGE skips this (qacct flavours diverge); callers that
            # need per-task SGE accuracy should wire sentinel files.
            result = self._scheduler.verify_completed(job_id)
            if not result.ok:
                logger.error(f"{self._scheduler.name} job array {job_id} failed: {result.detail}")
                return False
            logger.info(f"{self._scheduler.name} job array {job_id} completed")
            return True

        except (subprocess.CalledProcessError, RuntimeError) as e:
            logger.error(f"Error running cluster jobs: {e}")
            return False
    
    def save_intermediate_outputs(self, mode: str = "TM") -> bool:
        """Save intermediate Augustus outputs for debugging."""
        logger.info(f"Saving intermediate Augustus outputs for {mode} mode...")
        
        try:
            # Create intermediate outputs directory
            intermediate_dir = self.temp_dir / f"intermediate_{mode}_outputs"
            intermediate_dir.mkdir(exist_ok=True)
            logger.info(f"Created intermediate outputs directory: {intermediate_dir}")
            
            # Save all GFF files from jobs directory
            gff_files = list(self.jobs_dir.glob('*.gff'))
            logger.info(f"Found {len(gff_files)} GFF files to save")
            
            for gff_file in gff_files:
                if gff_file.exists() and os.path.getsize(gff_file) > 0:
                    # Copy to intermediate directory with descriptive name
                    intermediate_name = f"{mode}_{gff_file.name}"
                    intermediate_path = intermediate_dir / intermediate_name
                    shutil.copy2(gff_file, intermediate_path)
                    logger.info(f"Saved intermediate output: {intermediate_path}")
                    
                    # Log first few lines of the GFF file
                    with open(gff_file, 'r') as f:
                        lines = f.readlines()[:10]
                        logger.info(f"Preview of {gff_file.name}: {len(lines)} lines")
                        for i, line in enumerate(lines, 1):
                            logger.info(f"  Line {i}: {line.strip()}")
                else:
                    logger.warning(f"GFF file is empty or doesn't exist: {gff_file}")
            
            # Save all error files
            err_files = list(self.jobs_dir.glob('*.err'))
            logger.info(f"Found {len(err_files)} error files to save")
            
            for err_file in err_files:
                if err_file.exists():
                    intermediate_name = f"{mode}_{err_file.name}"
                    intermediate_path = intermediate_dir / intermediate_name
                    shutil.copy2(err_file, intermediate_path)
                    logger.info(f"Saved error file: {intermediate_path}")
                    
                    # Log error file content
                    with open(err_file, 'r') as f:
                        content = f.read().strip()
                        if content:
                            logger.info(f"Error content in {err_file.name}: {content}")
                        else:
                            logger.info(f"Error file {err_file.name} is empty")
            
            # Save job list and job files
            jobs_lst_file = self.temp_dir / f"jobs_{mode}.lst"
            if jobs_lst_file.exists():
                intermediate_jobs_lst = intermediate_dir / f"{mode}_jobs.lst"
                shutil.copy2(jobs_lst_file, intermediate_jobs_lst)
                logger.info(f"Saved job list: {intermediate_jobs_lst}")
            
            # Save actual job files
            job_files = list(self.temp_dir.glob(f"{self.args.genome}_aug_{mode}_*"))
            if not job_files:
                job_files = list(self.jobs_dir.glob(f"{self.args.genome}_aug_{mode}_*"))
            logger.info(f"Found {len(job_files)} job files to save")
            
            for job_file in job_files:
                if job_file.is_file():  # Only save files, not directories
                    intermediate_name = f"{mode}_{job_file.name}"
                    intermediate_path = intermediate_dir / intermediate_name
                    shutil.copy2(job_file, intermediate_path)
                    logger.info(f"Saved job file: {intermediate_path}")
            
            logger.info(f"Intermediate outputs saved to: {intermediate_dir}")
            return True
            
        except Exception as e:
            logger.error(f"Error saving intermediate outputs: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            return False

    def merge_augustus_output(self, mode: str = "TM") -> str:
        """Step 5: Merge Augustus output using join_aug_pred.pl to remove duplicates."""
        logger.info(f"Step 5: Merging Augustus output for {mode} mode using join_aug_pred.pl...")
        
        try:
            jobs_lst_file = self.temp_dir / f"jobs_{mode}.lst"
            merged_gff = self.temp_dir / f"augustus_{mode}_merged.gff"
            
            # Debug: List all files in jobs directory
            all_job_files = list(self.jobs_dir.glob('*'))
            logger.info(f"All files in jobs directory: {[f.name for f in all_job_files]}")
            
            # Prefer globbing the produced GFFs in the jobs directory (more robust than parsing job scripts)
            gff_files = sorted([str(p) for p in (self.jobs_dir).glob('*.gff')])
            logger.info(f"Found {len(gff_files)} GFF files: {[os.path.basename(f) for f in gff_files]}")

            # Debug: Check each GFF file
            valid_gff_files = []
            for gff_file in gff_files:
                if os.path.exists(gff_file):
                    file_size = os.path.getsize(gff_file)
                    logger.info(f"GFF file {os.path.basename(gff_file)}: {file_size} bytes")
                    if file_size > 0:
                        valid_gff_files.append(gff_file)
                        # Log first few lines
                        with open(gff_file, 'r') as f:
                            lines = f.readlines()[:5]
                            logger.info(f"  Preview: {len(lines)} lines")
                            for i, line in enumerate(lines, 1):
                                logger.info(f"    Line {i}: {line.strip()}")
                    else:
                        logger.warning(f"GFF file is empty: {gff_file}")
                else:
                    logger.warning(f"GFF file doesn't exist: {gff_file}")
            
            gff_files = valid_gff_files
            logger.info(f"Valid GFF files for merging: {len(gff_files)}")

            # Fallback: parse job scripts if no GFFs found yet
            if not gff_files and jobs_lst_file.exists():
                logger.info("No GFF files found, trying to parse job scripts...")
                with open(jobs_lst_file, 'r') as f:
                    job_files = [line.strip() for line in f if line.strip()]
                logger.info(f"Job files from list: {job_files}")
                
                for job_file_name in job_files:
                    # Job scripts are written to temp_dir (createAugustusJoblist.pl CWD);
                    # fall back to jobs_dir for older layouts.
                    if os.path.isabs(job_file_name):
                        job_file = job_file_name
                    else:
                        candidates = [
                            str(self.temp_dir / job_file_name),
                            str(self.jobs_dir / job_file_name),
                        ]
                        job_file = next((p for p in candidates if os.path.isfile(p)), candidates[0])
                    
                    logger.info(f"Parsing job file: {job_file_name} -> {job_file}")
                    try:
                        if os.path.exists(job_file):
                            with open(job_file, 'r') as jf:
                                for line_num, line in enumerate(jf, 1):
                                    if "--outfile=" in line:
                                        parts = line.split("--outfile=")
                                        if len(parts) > 1:
                                            outfile = parts[1].split()[0]
                                            logger.info(f"Found output file in job: {outfile}")
                                            if os.path.exists(outfile):
                                                file_size = os.path.getsize(outfile)
                                                logger.info(f"Output file exists: {outfile} ({file_size} bytes)")
                                                gff_files.append(outfile)
                                            else:
                                                logger.warning(f"Output file doesn't exist: {outfile}")
                                        break
                        else:
                            logger.warning(f"Job file doesn't exist: {job_file}")
                    except Exception as e:
                        logger.error(f"Error parsing job file {job_file}: {e}")
                        continue

            if not gff_files:
                logger.error("No valid GFF files found to merge")
                return None

            # Use join_aug_pred.pl to properly merge overlapping predictions
            # This script removes duplicates from overlapping chunks
            logger.info(f"Using join_aug_pred.pl to merge {len(gff_files)} GFF files and remove duplicates...")
            
            # First, concatenate all GFF files (join_aug_pred.pl reads from stdin)
            concatenated_gff = self.temp_dir / f"augustus_{mode}_concatenated.gff"
            logger.info(f"Concatenating {len(gff_files)} GFF files...")
            
            with open(concatenated_gff, 'w') as outf:
                for gff_file in gff_files:
                    logger.info(f"  Adding {os.path.basename(gff_file)}")
                    with open(gff_file, 'r') as inf:
                        outf.write(inf.read())
            
            logger.info(f"Created concatenated GFF: {concatenated_gff}")
            
            # Run join_aug_pred.pl (reads from stdin, writes to stdout)
            join_cmd = ["join_aug_pred.pl"]
            
            logger.info(f"Running: {' '.join(join_cmd)} < {concatenated_gff} > {merged_gff}")
            
            with open(concatenated_gff, 'r') as inf, open(merged_gff, 'w') as outf:
                result = subprocess.run(join_cmd, stdin=inf, stdout=outf, stderr=subprocess.PIPE, text=True, check=True)
            
            if result.stderr:
                logger.info(f"join_aug_pred.pl stderr: {result.stderr}")
            
            # Debug: Check merged file
            if merged_gff.exists():
                merged_size = os.path.getsize(merged_gff)
                with open(merged_gff, 'r') as f:
                    total_lines = sum(1 for _ in f)
                logger.info(f"Created deduplicated merged Augustus output: {merged_gff} ({merged_size} bytes, {total_lines} lines)")
                
                # Log first few lines of merged file
                with open(merged_gff, 'r') as f:
                    lines = f.readlines()[:10]
                    logger.info(f"Merged file preview ({len(lines)} lines):")
                    for i, line in enumerate(lines, 1):
                        logger.info(f"  Line {i}: {line.strip()}")
            else:
                logger.error("Merged GFF file was not created!")
                return None
            
            logger.info("Successfully merged and deduplicated Augustus output")
            return str(merged_gff)
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running join_aug_pred.pl: {e}")
            logger.error(f"stderr: {e.stderr}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            return None
        except Exception as e:
            logger.error(f"Error merging Augustus output: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            return None
    
    def munge_augustus_output(self, aug_output, mode, tm_tx):
        """
        Process Augustus output and convert to GTF format.
        Modified to handle multiple overlapping Augustus transcripts by selecting the best match.
        
        :param aug_output: List of Augustus output lines
        :param mode: Prediction mode ('TM' or 'TMR')
        :param tm_tx: TransMap transcript object
        :return: List of GTF records or None if processing failed
        """
        # extract the transcript lines
        tx_entries = [x.split() for x in aug_output if "\ttranscript\t" in x]
        valid_txs = [x[-1] for x in tx_entries if tm_tx.interval.overlap(tools.intervals.ChromosomeInterval(x[0], x[3],
                                                                                                        x[4], x[6]))]
        
        if len(valid_tx_candidates) != 0:
            return None
        
        valid_tx = valid_txs[0]
        tx_id = 'aug{}-{}'.format(mode, tm_tx.name)
        tx_lines = [x.split('\t') for x in aug_output if valid_tx in x and not x.startswith('#')]
        features = {"exon", "CDS", "start_codon", "stop_codon", "tts", "tss"}
        gtf = []
        for chrom, source, feature, start, stop, score, strand, frame, attributes in tx_lines:
            if feature not in features:
                continue
            new_attributes = 'transcript_id "{}"; gene_id "{}";'.format(tx_id, tm_tx.name2)
            gtf.append([chrom, source, feature, start, stop, score, strand, frame, new_attributes])
        return gtf

    def process_final_output(self, merged_gff: str, mode: str = "TM") -> bool:
        """Step 6: Process final output using the EXACT logic from working augustus.py.
        
        The key insight is that we need to process each Augustus output chunk individually
        with its corresponding TransMap transcript, just like the working version does.
        """
        logger.info(f"Step 6: Processing final output for {mode} mode using working augustus.py logic...")

        try:
            import tools.transcripts
            from tools.intervals import ChromosomeInterval

            # Load coding GP transcripts
            tm_tx_dict = tools.transcripts.get_gene_pred_dict(self.args.coding_gp)
            logger.info(f"Loaded {len(tm_tx_dict)} TM transcripts from {self.args.coding_gp}")
            
            # Check if merged GFF file exists and has content
            if not os.path.exists(merged_gff):
                logger.error(f"Merged GFF file does not exist: {merged_gff}")
                return False
            
            gff_size = os.path.getsize(merged_gff)
            logger.info(f"Merged GFF file size: {gff_size} bytes")
            
            if gff_size == 0:
                logger.warning(f"Merged GFF file is empty: {merged_gff}")
                # Create empty output file
                output_gtf = self.args.augustus_tm_gtf if mode == "TM" else self.args.augustus_tmr_gtf
                os.makedirs(os.path.dirname(output_gtf), exist_ok=True)
                with open(output_gtf, 'w') as out_f:
                    pass  # Create empty file
                logger.info(f"Created empty GTF file: {output_gtf}")
                return True

            # Read the merged GFF file
            with open(merged_gff, 'r') as f:
                all_lines = f.readlines()
            
            logger.info(f"Read {len(all_lines)} lines from merged GFF file")
            
            # OPTIMIZATION: Process transcripts with streaming output
            # SLURM mode now writes directly to the output file, returns record count
            # Local mode returns records list for compatibility
            
            result = self.process_transcripts_parallel(all_lines, mode, tm_tx_dict)
            
            # Check if result is a count (SLURM streaming mode) or list of records (local mode)
            if isinstance(result, int):
                # SLURM streaming mode - file already written
                record_count = result
                output_gtf = self.args.augustus_tm_gtf if mode == "TM" else self.args.augustus_tmr_gtf
                logger.info(f"Final GTF file already created by streaming: {output_gtf}")
                
                if record_count == 0 and len(tm_tx_dict) > 0:
                    # The merged GFF was non-empty (checked above), so Augustus DID predict
                    # genes, yet the transcript-processing array produced zero results. This
                    # means every batch was missing/failed (e.g. the SLURM array died on
                    # startup). Do NOT let this be silently treated as a successful (empty)
                    # run — that writes the .gtf.done sentinel and makes Snakemake think
                    # augMP is complete. Fail loudly so the rule re-runs.
                    logger.error(
                        f"No records written to {output_gtf} despite {len(tm_tx_dict)} input "
                        f"transcripts and a non-empty merged Augustus GFF. The transcript-"
                        f"processing array produced no output (all batches missing/failed). "
                        f"Treating augMP as FAILED so it is not marked complete."
                    )
                    return False
            else:
                # Local mode - need to write records
                all_gtf_records = result
                output_gtf = self.args.augustus_tm_gtf if mode == "TM" else self.args.augustus_tmr_gtf
                os.makedirs(os.path.dirname(output_gtf), exist_ok=True)
                
                logger.info(f"Writing {len(all_gtf_records)} records to {output_gtf}...")
                write_start = time.time()
                
                # Use large buffer for faster writing (16MB buffer)
                with open(output_gtf, 'w', buffering=16*1024*1024) as out_f:
                    # Pre-join records in batches for faster writing
                    batch_size = 10000
                    for i in range(0, len(all_gtf_records), batch_size):
                        batch = all_gtf_records[i:i+batch_size]
                        # Join batch records into lines
                        lines = ['\t'.join(map(str, rec)) + '\n' for rec in batch]
                        # Write entire batch at once
                        out_f.write(''.join(lines))
                        
                        # Log progress every 1M records
                        if (i + batch_size) % 1000000 == 0:
                            elapsed = time.time() - write_start
                            records_written = min(i + batch_size, len(all_gtf_records))
                            rate = records_written / elapsed if elapsed > 0 else 0
                            logger.info(f"Written {records_written:,}/{len(all_gtf_records):,} records "
                                      f"({rate:.0f} records/sec)")
                
                write_elapsed = time.time() - write_start
                logger.info(f"Created final GTF file: {output_gtf} (records: {len(all_gtf_records):,}, "
                           f"time: {write_elapsed:.1f}s, rate: {len(all_gtf_records)/write_elapsed:.0f} records/sec)")
                
                # Add debugging information
                if len(all_gtf_records) == 0:
                    logger.warning(f"No records written to {output_gtf}")
                    logger.warning(f"Total TM transcripts processed: {len(tm_tx_dict)}")
                    
                    # Show some examples of what was found
                    if tm_tx_dict:
                        logger.warning("Example TM transcript chromosomes:")
                        chroms = set(tx.chromosome for tx in tm_tx_dict.values())
                        for chrom in list(chroms)[:3]:
                            logger.warning(f"  {chrom}")
                    
                    # Debug: Show what's in the Augustus output
                    logger.warning("Augustus output preview:")
                    for i, line in enumerate(all_lines[:20]):
                        if not line.startswith('#'):
                            logger.warning(f"  Line {i}: {line.strip()}")
            
            return True

        except Exception as e:
            logger.error(f"Error processing final output: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            return False
    
    def process_transcripts_parallel(self, all_lines, mode, tm_tx_dict):
        """
        OPTIMIZED: Prefer SLURM array processing; fall back to local multiprocessing.
        Returns: int (record count) if SLURM streaming mode, list of records if local mode.
        """
        # Prefer SLURM path when enabled
        use_slurm = getattr(self.args, 'use_slurm_transcripts', True)
        if use_slurm:
            logger.info("Using SLURM-based parallel processing for transcript optimization...")
            try:
                result = self.process_transcripts_slurm(all_lines, mode, tm_tx_dict)
                # SLURM now returns int (record count) and writes file directly
                if result is not None and result >= 0:
                    return result
                else:
                    logger.warning("SLURM transcript processing produced no results; falling back to local multiprocessing")
            except Exception as e:
                logger.error(f"SLURM transcript processing failed: {e}")
                logger.warning("Falling back to local multiprocessing")

        logger.info("Using optimized local multiprocessing for transcript processing...")
        return self.process_transcripts_local_parallel(all_lines, mode, tm_tx_dict)
    
    def process_transcripts_local_parallel(self, all_lines, mode, tm_tx_dict):
        """
        Fallback: Local multiprocessing for transcript processing.
        Used when SLURM optimization is not available.
        """
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        import time
        
        # Honour --num_cpus (Snakemake threads in local mode); fall back to a
        # capped cpu_count when the flag is unset (standalone runs).
        requested = getattr(self.args, "num_cpus", None)
        if requested is not None and int(requested) > 0:
            num_processes = int(requested)
        else:
            num_processes = min(mp.cpu_count(), 32)
        num_processes = max(1, num_processes)
        logger.info(f"Processing {len(tm_tx_dict)} transcripts using local multiprocessing with {num_processes} processes")
        
        # Split transcripts into batches for parallel processing
        transcript_items = list(tm_tx_dict.items())
        batch_size = max(1, len(transcript_items) // num_processes)
        transcript_batches = [
            transcript_items[i:i + batch_size] 
            for i in range(0, len(transcript_items), batch_size)
        ]
        
        logger.info(f"Created {len(transcript_batches)} batches of ~{batch_size} transcripts each")
        
        # Use threads to avoid pickling issues on some clusters
        from concurrent.futures import ThreadPoolExecutor

        def process_batch_fn(batch_data):
            all_lines_local, mode_local, tm_tx_batch = batch_data
            results_local = []
            for tm_tx_id, tm_tx in tm_tx_batch:
                try:
                    gtf_records = self.munge_augustus_output(all_lines_local, mode_local, tm_tx)
                    if gtf_records:
                        results_local.extend(gtf_records)
                except Exception as e:
                    logger.error(f"Error processing transcript {tm_tx_id}: {e}")
                    continue
            return results_local

        # Process batches in parallel using threads
        start_time = time.time()
        all_gtf_records = []
        
        with ThreadPoolExecutor(max_workers=num_processes) as executor:
            # Submit all batches
            future_to_batch = {
                executor.submit(process_batch_fn, (all_lines, mode, batch)): i
                for i, batch in enumerate(transcript_batches)
            }
            
            # Collect results as they complete
            completed_batches = 0
            for future in future_to_batch:
                try:
                    batch_results = future.result()
                    all_gtf_records.extend(batch_results)
                    completed_batches += 1
                    
                    if completed_batches % 10 == 0:
                        elapsed = time.time() - start_time
                        logger.info(f"Completed {completed_batches}/{len(transcript_batches)} batches "
                                  f"in {elapsed:.1f}s ({completed_batches/elapsed:.1f} batches/sec)")
                        
                except Exception as e:
                    batch_idx = future_to_batch[future]
                    logger.error(f"Error processing batch {batch_idx}: {e}")
        
        elapsed = time.time() - start_time
        logger.info(f"Local parallel transcript processing completed in {elapsed:.1f}s")
        logger.info(f"Total GTF records found: {len(all_gtf_records)}")
        
        return all_gtf_records

    def process_transcripts_slurm(self, all_lines, mode, tm_tx_dict):
        """
        OPTIMIZED: Run transcript processing as SLURM array with streaming results.
        Each task processes a batch and streams to a shared output file with file locking.
        Returns list of GTF records by reading the single consolidated output file.
        """
        import math
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        # Working dirs
        slurm_root = self.work_dir / f"slurm_transcript_temp_{self.args.genome}"
        jobs_dir = slurm_root / "jobs"
        results_dir = slurm_root / "results"
        logs_dir = slurm_root / "logs"
        for p in (slurm_root, jobs_dir, results_dir, logs_dir):
            p.mkdir(parents=True, exist_ok=True)

        # Persist merged GFF to a temp file for workers to read
        merged_for_slurm = self.temp_dir / f"augustus_{mode}_merged_for_slurm.gff"
        with open(merged_for_slurm, 'w') as f:
            f.writelines(all_lines)
        logger.info(f"Created temporary merged GFF file: {merged_for_slurm}")

        # OPTIMIZATION 1: Use larger batches to reduce number of files (1000 transcripts per task)
        transcript_items = list(tm_tx_dict.items())
        total = len(transcript_items)
        transcripts_per_task = max(1, int(getattr(self.args, 'transcripts_per_task', 1000)))
        batch_size = transcripts_per_task
        batches = [transcript_items[i:i + batch_size] for i in range(0, total, batch_size)]
        logger.info(f"Creating {len(batches)} transcript batches for SLURM processing (batch_size={batch_size})...")

        # Write batch files (store transcript IDs per line)
        for idx, batch in enumerate(batches):
            job_file = jobs_dir / f"batch_{idx:04d}.txt"
            with open(job_file, 'w') as jf:
                for tx_id, _ in batch:
                    jf.write(tx_id + "\n")

        # OPTIMIZATION 2: Use a single consolidated output file instead of one per batch
        consolidated_output = results_dir / "all_results.gtf"
        
        # Create cluster array script with streaming output. Indexing is
        # 1-based at the scheduler level (works on SGE natively); the body
        # shifts to 0-based IDX so the batch_NNNN.txt filename convention
        # is preserved.
        array_script = slurm_root / "transcript_processing.slurm"
        array_count = len(batches)
        max_concurrent = getattr(self.args, 'slurm_transcripts_concurrency', 100)
        python_exe = sys.executable or "python3"
        task_var = self._scheduler.task_id_env()
        log_out, log_err = self._scheduler.array_log_paths(logs_dir, "tx")
        header = self._scheduler.header(
            job_name="aug-tx-proc",
            cpus=1,
            mem=getattr(self.args, 'slurm_transcripts_mem', '64G'),
            walltime=getattr(self.args, 'slurm_transcripts_time', '05:00:00'),
            log_out=log_out,
            log_err=log_err,
            partition=getattr(self.args, 'slurm_transcripts_partition', ''),
            queue=getattr(self.args, 'slurm_transcripts_partition', ''),
            array=(1, array_count),
            max_concurrent=max_concurrent,
        )
        with open(array_script, 'w') as sf:
            sf.write(self._wrap_cluster_header(header) + f"""
set -eo pipefail
export PYTHONPATH="{project_root}:${{PYTHONPATH:-}}"

TASK_ID="${{{task_var}}}"
IDX=$((TASK_ID - 1))
JOB_FILE=$(printf "{jobs_dir}/batch_%04d.txt" "$IDX")
RES_FILE=$(printf "{results_dir}/results_%04d.gtf" "$IDX")

{python_exe} {os.path.abspath(__file__)} --worker_mode --worker_merged_gff {merged_for_slurm} --worker_coding_gp {os.path.abspath(self.args.coding_gp)} --worker_run_mode {mode} --worker_batch_file "$JOB_FILE" --worker_result_file "$RES_FILE"
""")

        try:
            job_id = self._scheduler.submit(str(array_script))
            logger.info(f"Submitting {self._scheduler.name} array job for transcript processing...")
            logger.info(f"Submitted {self._scheduler.name} array job: {job_id}")
        except (subprocess.CalledProcessError, RuntimeError) as e:
            logger.error(f"Error submitting transcript array job: {e}")
            return []

        # Wait for completion via the scheduler's backend-appropriate poll.
        #
        # This is only a BACKSTOP: the loop exits as soon as the array leaves the
        # queue. Size it from the ACTUAL workload so we never abandon an array that
        # is merely working through a low concurrency limit. The previous fixed 2h
        # backstop fired while thousands of tasks were still queued behind
        # max_concurrent; the controller then assembled a partial result AND the
        # caller deleted the per-task batch files, so the still-pending tasks died
        # with FileNotFoundError (the "augMP empty" incident).
        def _walltime_to_hours(s, default_h=5.0):
            try:
                parts = [int(x) for x in str(s).split(":")]
                while len(parts) < 3:
                    parts.insert(0, 0)
                h, m, sec = parts[-3], parts[-2], parts[-1]
                return max(0.1, h + m / 60.0 + sec / 3600.0)
            except Exception:
                return default_h
        per_task_h = _walltime_to_hours(getattr(self.args, 'slurm_transcripts_time', '05:00:00'))
        waves = (array_count + max_concurrent - 1) // max(1, max_concurrent)
        # Allow every wave up to its full per-task walltime, plus a 1h buffer,
        # capped so a genuinely stuck scheduler can't hang the controller forever.
        timeout_sec = int(min(max(2 * 3600, waves * per_task_h * 3600 + 3600), 48 * 3600))
        logger.info(
            f"Transcript-processing wait backstop: {timeout_sec / 3600:.1f}h "
            f"(array={array_count}, concurrency={max_concurrent}, waves={waves}, "
            f"per_task<={per_task_h:.1f}h)"
        )
        start = time.time()
        timed_out = False
        while True:
            if not self._scheduler.job_present(job_id):
                break
            if time.time() - start > timeout_sec:
                logger.warning(
                    f"Timeout waiting for {self._scheduler.name} job {job_id} "
                    f"after {timeout_sec / 3600:.1f}h — tasks may still be queued/running."
                )
                timed_out = True
                break
            time.sleep(30)
        
        # Per-task validation is delegated to the downstream result-file
        # inspection (which is portable across SLURM and SGE). sacct gives a
        # richer answer on SLURM but doesn't have a portable SGE analog, so
        # we drop it here and rely on the existing "missing/empty output
        # file" detection below.
        logger.info(f"Cluster transcript processing job {job_id} returned; validating outputs below.")

        # OPTIMIZATION 3: Stream results directly to output file instead of loading into memory
        logger.info("Streaming results from worker files directly to output...")
        from concurrent.futures import ThreadPoolExecutor
        import threading
        
        # Get output path from args
        output_gtf = self.args.augustus_tm_gtf if mode == "TM" else self.args.augustus_tmr_gtf
        os.makedirs(os.path.dirname(output_gtf), exist_ok=True)
        
        missing = 0
        missing_lock = threading.Lock()
        total_records = 0
        records_lock = threading.Lock()
        write_start = time.time()
        
        # Open output file once for streaming writes
        with open(output_gtf, 'w', buffering=16*1024*1024) as out_f:  # 16MB buffer
            def stream_file(idx):
                nonlocal missing, total_records
                rf = results_dir / f"results_{idx:04d}.gtf"
                if not rf.exists() or os.path.getsize(rf) == 0:
                    logger.warning(f"Missing result file: {rf}")
                    with missing_lock:
                        missing += 1
                    return 0
                
                # Stream file contents directly to output (no parsing into memory)
                local_count = 0
                with open(rf, 'r', buffering=1024*1024) as f:  # 1MB buffer
                    for line in f:
                        out_f.write(line)
                        local_count += 1
                
                # Update total count
                with records_lock:
                    total_records += local_count
                    
                    # Log progress every 1M records
                    if total_records % 1000000 < local_count:
                        elapsed = time.time() - write_start
                        rate = total_records / elapsed if elapsed > 0 else 0
                        logger.info(f"Written {total_records:,} records ({rate:.0f} records/sec)")
                
                return local_count
            
            # Use ThreadPoolExecutor for parallel streaming I/O
            # Process files sequentially to maintain deterministic order
            for idx in range(array_count):
                try:
                    count = stream_file(idx)
                    
                    # Log progress every 100 files
                    if (idx + 1) % 100 == 0:
                        elapsed = time.time() - write_start
                        rate = total_records / elapsed if elapsed > 0 else 0
                        logger.info(f"Streamed {idx + 1}/{array_count} files, {total_records:,} records so far ({rate:.0f} records/sec)...")
                except Exception as e:
                    logger.error(f"Error streaming batch {idx}: {e}")

        write_elapsed = time.time() - write_start
        logger.info(f"Completed batches: {array_count - missing}/{array_count}")
        logger.info(f"Total GTF records: {total_records:,}")
        logger.info(f"Created final GTF file: {output_gtf} (time: {write_elapsed:.1f}s, rate: {total_records/write_elapsed:.0f} records/sec)")
        
        if missing > 0:
            # Print a few SLURM stderr/stdout lines to help debug missing batches.
            try:
                err_files = sorted([p for p in logs_dir.glob('tx_*_*.err') if p.stat().st_size > 0])
                out_files = sorted([p for p in logs_dir.glob('tx_*_*.out')])
                for coll, name in ((err_files[:3], 'err'), (out_files[:3], 'out')):
                    for p in coll:
                        try:
                            with open(p, 'r') as lf:
                                head = ''.join(lf.readlines()[:20]).strip()
                                if head:
                                    logger.warning(f"SLURM {name} sample {p.name}:\n{head}")
                        except Exception:
                            pass
            except Exception:
                pass

        # If we stopped waiting on the backstop while batches were still missing,
        # the array was abandoned mid-flight (tasks still queued behind the
        # concurrency limit). Do NOT return a silently-truncated result as success:
        # the caller would write the .gtf.done sentinel and delete the batch files,
        # and the still-pending tasks would die with FileNotFoundError. Fail loudly
        # so the whole augMP step re-runs cleanly.
        if timed_out and missing > 0:
            raise RuntimeError(
                f"Transcript-processing array {job_id} did not finish within the "
                f"{timeout_sec / 3600:.1f}h backstop: {missing}/{array_count} batches "
                f"still missing. Increase 'transcripts_concurrency' (currently "
                f"{max_concurrent}) and/or 'transcripts_time' for this rule so the "
                f"array can complete. Failing so augMP is not marked complete with a "
                f"truncated gene set."
            )

        # Return total count instead of records (to signal success/failure)
        return total_records
    
    def cleanup_temp_files(self):
        """Clean up temporary files while preserving SLURM logs and important debugging files."""
        logger.info("Cleaning up temporary files...")
        
        try:
            if self.temp_dir.exists():
                # Preserve important files for debugging
                preserve_patterns = [
                    '*.out',     # SLURM stdout
                    '*.err',     # SLURM stderr
                    '*.log',     # Log files
                    '*.slurm',   # SLURM scripts
                    'jobs_*.lst', # Job list files
                    'chr.lst',   # Chromosome list
                    'summary.out' # Summary file
                ]
                
                # Delete only intermediate processing files
                files_to_delete = []
                for item in self.temp_dir.rglob('*'):
                    if item.is_file():
                        # Check if file matches any preserve pattern
                        should_preserve = any(item.match(pattern) for pattern in preserve_patterns)
                        if not should_preserve:
                            files_to_delete.append(item)
                
                # Delete intermediate files
                for file_path in files_to_delete:
                    try:
                        file_path.unlink()
                        logger.debug(f"Deleted: {file_path}")
                    except Exception as e:
                        logger.warning(f"Could not delete {file_path}: {e}")
                
                # Remove empty subdirectories (but keep temp_dir itself with preserved files)
                for dirpath in sorted(self.temp_dir.rglob('*'), reverse=True):
                    if dirpath.is_dir() and not any(dirpath.iterdir()):
                        try:
                            dirpath.rmdir()
                            logger.debug(f"Removed empty directory: {dirpath}")
                        except Exception as e:
                            logger.warning(f"Could not remove directory {dirpath}: {e}")
                
                logger.info(f"Cleanup completed. Important files preserved in: {self.temp_dir}")
            else:
                logger.info("Temp directory does not exist, nothing to clean up")
        except Exception as e:
            logger.warning(f"Error during cleanup: {e}")
    
    def _coding_gp_empty(self) -> bool:
        gp = self.args.coding_gp
        if not gp or not os.path.exists(gp) or os.path.getsize(gp) == 0:
            return True
        with open(gp) as f:
            return not any(line.strip() for line in f)

    def _write_empty_gtf_outputs(self) -> bool:
        """Succeed with empty GTFs when there are no transcripts to annotate."""
        try:
            for path in (self.args.augustus_tm_gtf, getattr(self.args, "augustus_tmr_gtf", None)):
                if not path:
                    continue
                os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
                Path(path).touch()
                logger.info(f"Wrote empty Augustus output (no input transcripts): {path}")
            return True
        except Exception as e:
            logger.error(f"Failed to write empty Augustus outputs: {e}")
            return False

    def run_pipeline(self) -> bool:
        """Run the complete parallel Augustus pipeline."""
        logger.info("Starting parallel Augustus pipeline...")
        
        try:
            # Validate inputs
            if not self.validate_inputs():
                return False

            # Empty coding GP is common for sparse pairwise / ancestor modes — do
            # not run Augustus against missing hints files; emit empty GTFs.
            if self._coding_gp_empty():
                logger.warning(
                    f"No transcripts in coding_gp ({self.args.coding_gp}); "
                    "writing empty Augustus outputs and skipping prediction."
                )
                return self._write_empty_gtf_outputs()
            
            # Step 1: Run preprocessing (either on SLURM or locally)
            joblist_job_id = None  # Will be set if using SLURM preprocessing
            
            if not self.args.use_slurm_preprocessing:
                logger.info("Step 1: Running preprocessing steps locally...")
                if not self.split_genome_by_chromosome():
                    return False
                
                if not self.create_hints_location_file():
                    return False
                
                if not self.create_chromosome_hints_files():
                    return False
                
                if not self.create_augustus_jobs("TM"):
                    return False
                
                if self.args.augustus_tmr_gtf:
                    if not self.create_augustus_jobs("TMR"):
                        return False
            else:
                logger.info("Step 1: Running preprocessing steps on SLURM (3-stage parallel pipeline)...")
                
                # Stage 1: Initial setup (genome splitting + chr.lst creation)
                logger.info("Stage 1a: Genome splitting and chromosome list creation...")
                setup_script = self.create_initial_setup_slurm_script()
                setup_job_id = self.submit_slurm_job(setup_script, "initial setup")
                
                if not setup_job_id:
                    logger.error("Failed to submit initial setup job")
                    return False
                
                # Wait for setup to complete
                if not self.wait_for_slurm_job(setup_job_id, "initial setup", check_interval=30):
                    logger.error("Initial setup job failed")
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
                    slurm_out_files = list(self.temp_dir.glob("augustus_setup_*.out"))
                    slurm_err_files = list(self.temp_dir.glob("augustus_setup_*.err"))
                    
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
                hints_job_id = self.submit_slurm_job(hints_script, "hints generation array")
                
                if not hints_job_id:
                    logger.error("Failed to submit hints generation array job")
                    return False
                
                # Wait for all hints jobs to complete
                if not self.wait_for_slurm_job(hints_job_id, "hints generation array", check_interval=30):
                    logger.error("Hints generation array job failed")
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
                    # SLURM: augustus_hints_%A_%a.err; SGE: augustus_hints.$JOB_ID.$TASK_ID.err
                    slurm_err_files = sorted(self.temp_dir.glob("augustus_hints*.err"))
                    if slurm_err_files:
                        logger.error(f"\nCheck cluster error files in: {self.temp_dir}")
                        for err_file in slurm_err_files[:3]:
                            logger.error(f"  - {err_file.name}")
                    else:
                        logger.error(
                            f"\nNo augustus_hints*.err files found in {self.temp_dir}; "
                            "inspect qstat/qacct for the hints array job."
                        )
                    return False
                
                logger.info(f"Successfully created hints files for all {num_chromosomes} chromosomes")
                
                # Stage 3: Job list generation
                logger.info("Stage 3: Creating Augustus job lists...")
                joblist_script = self.create_joblist_generation_slurm_script(hints_job_id)
                joblist_job_id = self.submit_slurm_job(joblist_script, "joblist generation")
                
                if not joblist_job_id:
                    logger.error("Failed to submit joblist generation job")
                    return False
                
                # Wait for joblist generation to complete
                if not self.wait_for_slurm_job(joblist_job_id, "joblist generation", check_interval=30):
                    logger.error("Joblist generation job failed")
                    return False
                
                # Verify joblist outputs
                logger.info("Verifying joblist outputs...")
                expected_files = [self.temp_dir / "jobs_TM.lst"]
                if self.args.augustus_tmr_gtf:
                    expected_files.append(self.temp_dir / "jobs_TMR.lst")
                
                missing_files = [f for f in expected_files if not f.exists()]
                if missing_files:
                    logger.error("Joblist generation completed but expected files are missing:")
                    for f in missing_files:
                        logger.error(f"  - {f}")
                    
                    # Check for SLURM output files
                    slurm_out_files = list(self.temp_dir.glob("augustus_joblist_*.out"))
                    slurm_err_files = list(self.temp_dir.glob("augustus_joblist_*.err"))
                    
                    if slurm_out_files:
                        logger.error(f"\nCheck SLURM output file: {slurm_out_files[0]}")
                    if slurm_err_files:
                        logger.error(f"Check SLURM error file: {slurm_err_files[0]}")
                    
                    return False
                
                logger.info("All preprocessing stages completed successfully!")
            
            # Step 2: Run TM mode Augustus jobs
            logger.info("Step 2: Running TM mode Augustus jobs...")
            
            # Count number of TM jobs
            jobs_file = self.temp_dir / "jobs_TM.lst"
            if not jobs_file.exists():
                logger.error(f"TM jobs file not found: {jobs_file}")
                return False
                
            with open(jobs_file, 'r') as f:
                num_tm_jobs = len([line for line in f if line.strip()])
            
            # Use the last preprocessing job ID as dependency (joblist_job_id if SLURM, None if local)
            last_preprocessing_job_id = joblist_job_id if self.args.use_slurm_preprocessing else None
            if self._should_run_jobs_locally():
                logger.info(f"Running {num_tm_jobs} TM Augustus jobs locally...")
                if not self.run_local_jobs(str(jobs_file), getattr(self.args, "num_cpus", None)):
                    return False
            else:
                tm_slurm_script = self.create_slurm_script("TM", num_tm_jobs, last_preprocessing_job_id)
                if not self.run_slurm_jobs(tm_slurm_script):
                    return False
            
            # Save intermediate outputs for debugging
            if self.args.save_intermediate:
                logger.info("Saving intermediate TM outputs for debugging...")
                self.save_intermediate_outputs("TM")
            else:
                logger.info("Skipping intermediate TM output saving (disabled)")
            
            merged_gff = self.merge_augustus_output("TM")
            if not merged_gff:
                return False
            
            if not self.process_final_output(merged_gff, "TM"):
                return False
            
            # Step 3: Run TMR mode if requested
            if self.args.augustus_tmr_gtf:
                logger.info("Step 3: Running TMR mode Augustus jobs...")
                
                # Count number of TMR jobs
                tmr_jobs_file = self.temp_dir / "jobs_TMR.lst"
                if not tmr_jobs_file.exists():
                    logger.error(f"TMR jobs file not found: {tmr_jobs_file}")
                    return False
                    
                with open(tmr_jobs_file, 'r') as f:
                    num_tmr_jobs = len([line for line in f if line.strip()])
                
                if self._should_run_jobs_locally():
                    logger.info(f"Running {num_tmr_jobs} TMR Augustus jobs locally...")
                    if not self.run_local_jobs(str(tmr_jobs_file), getattr(self.args, "num_cpus", None)):
                        return False
                else:
                    tmr_slurm_script = self.create_slurm_script("TMR", num_tmr_jobs, last_preprocessing_job_id)
                    if not self.run_slurm_jobs(tmr_slurm_script):
                        return False
                
                # Save intermediate outputs for debugging
                if self.args.save_intermediate:
                    logger.info("Saving intermediate TMR outputs for debugging...")
                    self.save_intermediate_outputs("TMR")
                else:
                    logger.info("Skipping intermediate TMR output saving (disabled)")
                
                merged_gff = self.merge_augustus_output("TMR")
                if not merged_gff:
                    return False
                
                if not self.process_final_output(merged_gff, "TMR"):
                    return False
            
            logger.info("Parallel Augustus pipeline completed successfully!")
            return True
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            return False
        
        finally:
            if not self.args.keep_temp:
                self.cleanup_temp_files()


def main():
    """Main entry point."""
    # Check if we're in worker mode BEFORE parsing args to avoid required argument issues
    is_worker_mode = '--worker_mode' in sys.argv
    
    parser = argparse.ArgumentParser(description="Parallel Augustus Pipeline for CAT")
    
    # Input Files - only required if not in worker mode
    parser.add_argument("--genome_fasta", required=not is_worker_mode, 
                       help="Genome FASTA file for Augustus gene prediction.")
    parser.add_argument("--coding_gp", required=not is_worker_mode, 
                       help="GenePred file containing coding transcripts from TransMap.")
    parser.add_argument("--filtered_tm_psl", required=not is_worker_mode, 
                       help="Filtered TransMap PSL alignment file.")
    parser.add_argument("--ref_psl", required=not is_worker_mode, 
                       help="Reference genome PSL alignment file.")
    parser.add_argument("--annotation_gp", required=not is_worker_mode, 
                       help="Reference annotation in GenePred format.")
    parser.add_argument("--tm_cfg", required=not is_worker_mode, 
                       help="Augustus configuration file for TM (TransMap) mode.")
    parser.add_argument("--miniprot_hints_gff", required=False,
                       help="GFF file containing protein alignment hints from Miniprot (only needed for augMP mode).")
    
    # Augustus Parameters - only required if not in worker mode
    parser.add_argument("--genome", required=not is_worker_mode, 
                       help="Genome name identifier for retrieving RNA-seq hints from database.")
    parser.add_argument("--augustus_species", required=not is_worker_mode, 
                       help="Species parameter for Augustus (e.g., 'human', 'mouse', 'fly').")
    parser.add_argument("--utr", type=int, required=not is_worker_mode, choices=[0, 1], 
                       help="UTR prediction parameter for Augustus (0=no UTRs, 1=predict UTRs).")
    
    # Output Files - only required if not in worker mode
    parser.add_argument("--augustus_tm_gtf", required=not is_worker_mode, 
                       help="Output path for Augustus TM mode predictions in GTF format.")
    
    # TMR Mode Parameters (optional)
    parser.add_argument("--augustus_tmr_gtf", 
                       help="Output path for Augustus TMR mode predictions in GTF format. "
                            "If provided, TMR mode is activated using RNA-seq evidence.")
    parser.add_argument("--augustus_hints_db", 
                       help="SQLite database file containing RNA-seq hints. Required for TMR mode.")
    parser.add_argument("--tmr_cfg", 
                       help="Augustus configuration file for TMR (TransMap + RNA-seq) mode. Required for TMR mode.")
    
    # Pipeline options - only work_dir required if not in worker mode
    parser.add_argument("--work_dir", required=not is_worker_mode,
                       help="Working directory for temporary files.")
    parser.add_argument("--keep_temp", action="store_true", default=True,
                       help="Keep temporary files for debugging. Default: True")
    parser.add_argument("--use_slurm_preprocessing", action="store_true", default=True,
                       help="Use SLURM for preprocessing steps (genome splitting, hints creation, job list generation). Default: True")
    parser.add_argument("--no_slurm_preprocessing", action="store_true",
                       help="Disable SLURM preprocessing and run preprocessing steps locally.")
    parser.add_argument("--save_intermediate", action="store_true", default=True,
                       help="Save intermediate Augustus outputs for debugging. Default: True")
    parser.add_argument("--no_save_intermediate", action="store_true",
                       help="Disable saving intermediate Augustus outputs.")
    # SLURM transcript processing controls
    parser.add_argument("--use_slurm_transcripts", action="store_true", default=True,
                       help="Use SLURM array for transcript processing (default: enabled).")
    parser.add_argument("--no_slurm_transcripts", action="store_true",
                       help="Disable SLURM transcript processing and run transcripts locally.")
    parser.add_argument("--no_slurm_jobs", action="store_true",
                       help="Run Augustus chunk job scripts locally instead of a cluster array.")
    parser.add_argument("--transcripts_per_task", type=int, default=1000,
                       help="Number of transcripts to process per SLURM task (default: 1000, optimized for faster collection).")
    parser.add_argument("--array_concurrency", type=int, default=500,
                       help="Max concurrent SLURM array tasks (default: 500).")
    # SLURM resource configuration (applied when using SLURM mode)
    parser.add_argument("--slurm_partition", default="",
                       help="SLURM partition for preprocessing steps (setup/hints/joblist). Empty = cluster default.")
    parser.add_argument("--slurm_jobs_partition", default="",
                       help="SLURM partition for Augustus execution array jobs. Empty = cluster default.")
    parser.add_argument("--slurm_transcripts_partition", default="",
                       help="SLURM partition for transcript processing array jobs. Empty = cluster default.")
    parser.add_argument("--slurm_hints_mem", default="64G",
                       help="Memory per hints generation SLURM task. Default: 64G.")
    parser.add_argument("--slurm_jobs_mem", default="16G",
                       help="Memory per Augustus execution SLURM task. Default: 16G.")
    parser.add_argument("--slurm_transcripts_mem", default="64G",
                       help="Memory per transcript processing SLURM task. Default: 64G.")
    parser.add_argument("--slurm_setup_mem", default="",
                       help="Memory for setup/joblist jobs. Empty = use --slurm_jobs_mem.")
    parser.add_argument("--slurm_setup_time", default="04:00:00",
                       help="Time limit for setup and joblist SLURM jobs. Default: 04:00:00.")
    parser.add_argument("--slurm_hints_time", default="04:00:00",
                       help="Time limit for hints generation SLURM array jobs. Default: 04:00:00.")
    parser.add_argument("--slurm_jobs_time", default="01:00:00",
                       help="Time limit for Augustus execution SLURM array jobs. Default: 01:00:00.")
    parser.add_argument("--slurm_transcripts_time", default="05:00:00",
                       help="Time limit for transcript processing SLURM array jobs. Default: 05:00:00.")
    parser.add_argument("--slurm_hints_concurrency", type=int, default=10,
                       help="Max concurrent hints generation array tasks. Default: 10.")
    parser.add_argument("--slurm_jobs_concurrency", type=int, default=100,
                       help="Max concurrent Augustus execution array tasks. Default: 100.")
    parser.add_argument("--slurm_transcripts_concurrency", type=int, default=100,
                       help="Max concurrent transcript processing array tasks. Default: 100. "
                            "Low values here (with a large miniprot set) can leave thousands of "
                            "tasks queued and cause the controller's wait backstop to fire early.")
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
    parser.add_argument("--num_cpus", type=int, default=None,
                       help="Max local worker threads (Snakemake threads in local mode). "
                            "Defaults to min(cpu_count, 32) when unset.")
    # Worker mode (internal; used by SLURM array)
    # Worker mode (internal; use unique flag names to avoid conflicts)
    parser.add_argument("--worker_mode", action="store_true",
                       help="Internal: run transcript worker mode (called by SLURM array).")
    parser.add_argument("--worker_merged_gff", help="Merged GFF input for worker mode")
    parser.add_argument("--worker_coding_gp", help="Coding GP path for worker mode")
    parser.add_argument("--worker_run_mode", help="Mode for worker (TM/TMR)")
    parser.add_argument("--worker_batch_file", help="Batch file with transcript IDs for worker")
    parser.add_argument("--worker_result_file", help="Worker output GTF path")
    
    args = parser.parse_args()
    
    # Handle SLURM preprocessing options
    if args.no_slurm_preprocessing:
        args.use_slurm_preprocessing = False

    # Handle SLURM transcript options
    if args.no_slurm_transcripts:
        args.use_slurm_transcripts = False

    # Resolve auto → slurm/sge/local so local Snakemake runs don't qsub via auto-detect.
    from cat.scheduler import resolve_execution_mode
    args.execution_mode = resolve_execution_mode(args.execution_mode)
    
    # Handle intermediate output saving options
    if args.no_save_intermediate:
        args.save_intermediate = False
    
    # Validate TMR mode requirements
    if args.augustus_tmr_gtf and (not args.augustus_hints_db or not args.tmr_cfg):
        parser.error("--augustus_hints_db and --tmr_cfg are required when --augustus_tmr_gtf is specified.")
    
    # Worker mode: process a batch and exit
    if args.worker_mode:
        import tools.transcripts
        import tools.intervals
        # Read lines
        with open(args.worker_merged_gff, 'r') as f:
            aug_lines = f.readlines()
        tm_tx_dict = tools.transcripts.get_gene_pred_dict(args.worker_coding_gp)
        # Munge per transcript
        def munge_local(aug_lines, mode, tm_tx):
            tx_entries = [x.split() for x in aug_lines if "\ttranscript\t" in x]
            # Store entries with overlap size for selection
            valid_tx_candidates = []
            for x in tx_entries:
                try:
                    aug_interval = tools.intervals.ChromosomeInterval(x[0], int(x[3]), int(x[4]), x[6])
                    if tm_tx.interval.overlap(aug_interval):
                        # Calculate overlap size
                        intersection = tm_tx.interval.intersection(aug_interval)
                        overlap_size = (intersection.stop - intersection.start) if intersection else 0
                        valid_tx_candidates.append((x[-1], overlap_size, int(x[3]), int(x[4])))
                except (ValueError, IndexError):
                    continue
            
            if len(valid_tx_candidates) == 0:
                return []
            
            # Select the best match: highest overlap size, then shorter transcript (more complete prediction)
            # Sort by: -overlap_size (descending), then transcript_length (ascending)
            valid_tx_candidates.sort(key=lambda x: (-x[1], x[3] - x[2]))
            valid_tx = valid_tx_candidates[0][0]
            features = {"exon", "CDS", "start_codon", "stop_codon", "tts", "tss"}
            out = []
            for line in aug_lines:
                if line.startswith('#'):
                    continue
                if valid_tx in line:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 9:
                        continue
                    chrom, source, feature, start, stop, score, strand, frame, attributes = parts
                    if feature not in features:
                        continue
                    new_attr = f'transcript_id "aug{mode}-{tm_tx.name}"; gene_id "{tm_tx.name2}";'
                    out.append([chrom, source, feature, start, stop, score, strand, frame, new_attr])
            return out
        all_out = []
        with open(args.worker_batch_file, 'r') as bf:
            for tx_id in (ln.strip() for ln in bf if ln.strip()):
                tx = tm_tx_dict.get(tx_id)
                if tx is None:
                    continue
                all_out.extend(munge_local(aug_lines, args.worker_run_mode, tx))
        os.makedirs(os.path.dirname(args.worker_result_file), exist_ok=True)
        with open(args.worker_result_file, 'w') as outf:
            for rec in all_out:
                outf.write("\t".join(map(str, rec)) + "\n")
        sys.exit(0)

    # Run pipeline
    pipeline = ParallelAugustus(args)
    success = pipeline.run_pipeline()
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
