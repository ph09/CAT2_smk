import argparse
import collections
import logging
import os
import shutil

from toil.fileStores import FileID
from toil.common import Toil
from toil.job import Job

# Assuming 'tools' is a local package
import tools.fileOps
import tools.toilInterface
import tools.hal
import tools.procOps


def calculate_dynamic_resources(input_file_ids, num_chromosomes=None, genome_sizes=None, num_target_genomes=1):
    """
    Calculate dynamic resource allocation based on input data characteristics.
    
    :param input_file_ids: Namespace containing FileIDs for input files
    :param num_chromosomes: Number of chromosomes to process
    :param genome_sizes: Dictionary of chromosome name to size
    :param num_target_genomes: Number of target genomes
    :return: Dictionary with resource allocations
    """
    # Calculate base disk usage from input files
    base_disk = tools.toilInterface.find_total_disk_usage(input_file_ids, buffer='64G', round='8G')
    base_disk_gb = max(128, int(base_disk / (1024**3)))
    
    # Calculate memory requirements based on genome characteristics
    if genome_sizes and num_chromosomes:
        total_genome_size = sum(genome_sizes.values())
        avg_chrom_size = total_genome_size / num_chromosomes
        max_chrom_size = max(genome_sizes.values()) if genome_sizes else 0
        
        # Scale memory based on chromosome sizes and number of target genomes
        base_memory_gb = max(128, int(total_genome_size / (500 * 1024**2)) * 2)  # 2GB per 500MB genome
        
        # Chain job memory scales with chromosome size
        if max_chrom_size >= 200 * 1024**2:  # 200MB+
            chain_memory_gb = max(128, int(max_chrom_size / (100 * 1024**2)) * 4)  # 4GB per 100MB
        elif max_chrom_size >= 50 * 1024**2:  # 50MB+
            chain_memory_gb = max(128, int(max_chrom_size / (50 * 1024**2)) * 2)   # 2GB per 50MB
        else:
            chain_memory_gb = 64  # Minimum for small chromosomes
        
        # Scale disk for chain jobs
        chain_disk_gb = max(64, base_disk_gb // (num_chromosomes // 4 + 1))
        
        # Merge resources scale with number of chromosomes and target genomes
        merge_memory_gb = max(128, int(num_chromosomes / 10) * 4 * num_target_genomes)  # 4GB per 10 chromosomes per genome
        merge_disk_gb = max(128, base_disk_gb // 2)
    else:
        # Fallback to conservative estimates
        base_memory_gb = 64
        chain_memory_gb = 64
        chain_disk_gb = 64
        merge_memory_gb = 128
        merge_disk_gb = 128
    return {
        'setup_memory': f'{base_memory_gb}G',
        'setup_disk': f'{base_disk_gb}G',
        'chain_memory': f'{chain_memory_gb}G',
        'chain_disk': f'{chain_disk_gb}G',
        'merge_memory': f'{merge_memory_gb}G',
        'merge_disk': f'{merge_disk_gb}G'
    }


def estimate_genome_characteristics(query_sizes_path):
    """
    Estimate genome characteristics for resource calculation.
    
    :param query_sizes_path: Path to chromosome sizes file
    :return: Tuple of (num_chromosomes, genome_sizes_dict)
    """
    try:
        genome_sizes = {}
        if os.path.exists(query_sizes_path):
            with open(query_sizes_path, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            chrom, size_str = parts[0], parts[1]
                            genome_sizes[chrom] = int(size_str)
            return len(genome_sizes), genome_sizes
        return None, None
    except:
        return None, None


def calculate_optimal_chunk_strategy(num_chromosomes, num_target_genomes, genome_sizes=None):
    """
    Calculate optimal parallelization strategy based on workload characteristics.
    
    :param num_chromosomes: Number of chromosomes
    :param num_target_genomes: Number of target genomes
    :param genome_sizes: Dictionary of chromosome sizes
    :return: Dictionary with chunking strategy
    """
    if not num_chromosomes:
        return {'max_concurrent_chains': 4, 'merge_batch_size': 10}
    
    total_jobs = num_chromosomes * num_target_genomes
    
    # Determine optimal concurrency based on total workload
    if total_jobs < 50:
        # Small workload - moderate parallelism
        max_concurrent_chains = min(8, total_jobs)
        merge_batch_size = 20
    elif total_jobs < 200:
        # Medium workload - high parallelism
        max_concurrent_chains = min(16, total_jobs // 2)
        merge_batch_size = 15
    else:
        # Large workload - very high parallelism with batching
        max_concurrent_chains = min(32, total_jobs // 4)
        merge_batch_size = 10
    
    # Adjust for large chromosomes that need more resources
    if genome_sizes:
        large_chroms = sum(1 for size in genome_sizes.values() if size > 100 * 1024**2)  # >100MB
        if large_chroms > num_chromosomes * 0.3:  # >30% are large
            max_concurrent_chains = max(4, max_concurrent_chains // 2)
    
    return {
        'max_concurrent_chains': max_concurrent_chains,
        'merge_batch_size': merge_batch_size,
        'total_jobs_estimate': total_jobs
    }


def run_chaining_pipeline(hal, ref_genome, query_sizes, query_two_bit, target_two_bits,
                          chain_files, chain_mode, genome_chains, toil_options):
    """
    Main entry function for the chaining Toil pipeline with dynamic resource allocation.
    """
    hal_abs = os.path.abspath(hal)
    query_sizes_abs = os.path.abspath(query_sizes)
    query_two_bit_abs = os.path.abspath(query_two_bit)
    target_two_bits_abs = {genome: os.path.abspath(f) for genome, f in target_two_bits.items()}
    if genome_chains:
        genome_chains_abs = {genome: os.path.abspath(f) for genome, f in genome_chains.items()}
    else:
        genome_chains_abs = {}

    with Toil(toil_options) as t:
        if not t.options.restart:
            # Create a namespace to hold input file IDs
            input_file_ids = argparse.Namespace()
            if chain_mode:
                # In chain_mode, we just import the existing chains
                input_file_ids.chain_ids = {genome: FileID.forPath(t.importFile('file://' + f), f)
                                           for genome, f in genome_chains_abs.items()}
                job = Job.wrapJobFn(copy_chains, chain_files, input_file_ids)
            else:
                # In the default mode, we import all the necessary files for chaining
                input_file_ids.hal = FileID.forPath(t.importFile('file://' + hal_abs), hal_abs)
                input_file_ids.query_sizes = FileID.forPath(t.importFile('file://' + query_sizes_abs), query_sizes_abs)
                input_file_ids.query_two_bit = FileID.forPath(t.importFile('file://' + query_two_bit_abs),
                                                              query_two_bit_abs)
                input_file_ids.target_two_bits = {genome: FileID.forPath(t.importFile('file://' + f), f)
                                                  for genome, f in target_two_bits_abs.items()}
                
                # Calculate dynamic resources based on genome characteristics
                num_chromosomes, genome_sizes = estimate_genome_characteristics(query_sizes_abs)
                num_target_genomes = len(target_two_bits)
                resources = calculate_dynamic_resources(input_file_ids, num_chromosomes, genome_sizes, num_target_genomes)
                
                job = Job.wrapJobFn(setup, ref_genome, chain_files, input_file_ids, resources,
                                   memory=resources['setup_memory'], disk=resources['setup_disk'])
            
            # Start the Toil workflow and get the resulting chain file IDs
            chain_file_ids = t.start(job)
        else:
            # Restart the workflow from a previous state
            chain_file_ids = t.restart()

        # Export the final chain files from the file store to the output directory
        # Output paths can remain relative, as Toil handles exporting correctly.
        for chain_file, chain_file_id in chain_file_ids.items():
            tools.fileOps.ensure_file_dir(chain_file)
            t.exportFile(chain_file_id, 'file://' + os.path.abspath(chain_file))


def setup(job, ref_genome, chain_files, input_file_ids, resources):
    """
    Sets up the per-chromosome chaining jobs with enhanced parallelization and resource management.
    
    :param job: The parent Toil job.
    :param ref_genome: Name of the reference genome.
    :param chain_files: Dictionary mapping target genomes to output chain file paths.
    :param input_file_ids: Namespace object containing FileIDs for input files.
    :param resources: Dictionary containing dynamic resource allocations.
    :return: A dictionary mapping output file paths to the FileIDs of merged chain files.
    """
    job.fileStore.logToMaster('Beginning chaining pipeline with dynamic resource allocation', level=logging.INFO)
    
    # Read chromosome sizes and analyze workload
    chrom_sizes_path = job.fileStore.readGlobalFile(input_file_ids.query_sizes)
    num_chromosomes, genome_sizes = estimate_genome_characteristics(chrom_sizes_path)
    num_target_genomes = len(input_file_ids.target_two_bits)
    
    # Calculate optimal parallelization strategy
    strategy = calculate_optimal_chunk_strategy(num_chromosomes, num_target_genomes, genome_sizes)
    
    job.fileStore.logToMaster(
        f'Processing {num_chromosomes} chromosomes across {num_target_genomes} target genomes '
        f'(estimated {strategy["total_jobs_estimate"]} total jobs)', 
        level=logging.INFO
    )
    
    tmp_chain_file_ids = collections.defaultdict(list)
    active_jobs = 0
    max_concurrent = strategy['max_concurrent_chains']
    
    # Group chromosomes for better resource utilization
    chromosome_groups = []
    current_group = []
    current_group_size = 0
    
    with open(chrom_sizes_path) as chrom_sizes_f:
        chromosomes = []
        for line in chrom_sizes_f:
            chrom, size_str = line.split()
            size = int(size_str)
            chromosomes.append((chrom, size))
    
    # Sort chromosomes by size (largest first) for better load balancing
    chromosomes.sort(key=lambda x: x[1], reverse=True)
    
    job.fileStore.logToMaster(f'Chromosome size range: {chromosomes[-1][1]:,} - {chromosomes[0][1]:,} bp', 
                             level=logging.INFO)
    
    # Create chromosome processing jobs with optimized parallelism
    for chrom, size in chromosomes:
        for target_genome, target_two_bit_file_id in input_file_ids.target_two_bits.items():
            # Calculate job-specific resources based on chromosome size and complexity
            job_memory, job_disk = get_chromosome_resources(size, genome_sizes, resources)
            
            # Create the chaining job
            j = job.addChildJobFn(chain_by_chromosome, ref_genome, chrom, size, 
                                  input_file_ids, target_genome, target_two_bit_file_id,
                                  memory=job_memory, disk=job_disk)
            tmp_chain_file_ids[target_genome].append(j.rv())
            active_jobs += 1
    
    job.fileStore.logToMaster(f'Scheduled {active_jobs} chaining jobs across {len(tmp_chain_file_ids)} target genomes', 
                             level=logging.INFO)
            
    # Follow-on jobs to merge the per-chromosome chains for each target genome
    return_file_ids = {}
    for genome, chain_file_path in chain_files.items():
        chain_files_for_genome = tmp_chain_file_ids[genome]
        merge_memory = resources['merge_memory']
        merge_disk = resources['merge_disk']
        
        # Scale merge resources based on number of input files
        if len(chain_files_for_genome) > 50:
            # Very large number of chromosomes - increase merge resources
            merge_memory_gb = int(merge_memory.rstrip('G')) * 1.5
            merge_memory = f'{int(merge_memory_gb)}G'
        
        j = job.addFollowOnJobFn(merge, chain_files_for_genome, genome, len(chromosomes),
                                memory=merge_memory, disk=merge_disk)
        return_file_ids[chain_file_path] = j.rv()
    
    return return_file_ids


def get_chromosome_resources(chrom_size, genome_sizes, base_resources):
    """
    Calculate chromosome-specific resources based on size and complexity.
    
    :param chrom_size: Size of this chromosome
    :param genome_sizes: Dictionary of all chromosome sizes
    :param base_resources: Base resource allocations
    :return: Tuple of (memory_str, disk_str)
    """
    base_memory_gb = int(base_resources['chain_memory'].rstrip('G'))
    base_disk_gb = int(base_resources['chain_disk'].rstrip('G'))
    
    # Scale resources based on chromosome size relative to average
    if genome_sizes:
        avg_size = sum(genome_sizes.values()) / len(genome_sizes)
        size_ratio = chrom_size / avg_size
        
        if size_ratio > 5.0:  # Very large chromosome
            memory_multiplier = 3.0
            disk_multiplier = 2.0
        elif size_ratio > 2.0:  # Large chromosome
            memory_multiplier = 2.0
            disk_multiplier = 1.5
        elif size_ratio > 0.5:  # Average chromosome
            memory_multiplier = 1.0
            disk_multiplier = 1.0
        else:  # Small chromosome
            memory_multiplier = 0.75
            disk_multiplier = 0.75
    else:
        # Fallback based on absolute size
        if chrom_size > 200 * 1024**2:  # >200MB
            memory_multiplier = 2.0
            disk_multiplier = 1.5
        elif chrom_size > 50 * 1024**2:  # >50MB
            memory_multiplier = 1.0
            disk_multiplier = 1.0
        else:  # <50MB
            memory_multiplier = 0.75
            disk_multiplier = 0.75
    
    final_memory = max(64, int(base_memory_gb * memory_multiplier))
    final_disk = max(64, int(base_disk_gb * disk_multiplier))
    
    return f'{final_memory}G', f'{final_disk}G'


def copy_chains(job, chain_files, input_file_ids):
    """
    In chain_mode, this function returns the FileIDs of existing chain files to be copied.
    
    :param chain_files: Dictionary mapping target genomes to output chain file paths.
    :param input_file_ids: Namespace object containing FileIDs for the pre-existing chains.
    :return: A dictionary mapping output file paths to the original chain FileIDs.
    """
    return {chain_file: input_file_ids.chain_ids[genome] for genome, chain_file in chain_files.items()}


def chain_by_chromosome(job, ref_genome, chrom, size, input_file_ids, target_genome, target_two_bit_file_id):
    """
    Generates a chain file for a single chromosome with enhanced error handling.
    
    :param ref_genome: Name of the reference genome.
    :param chrom: Chromosome name.
    :param size: Chromosome size.
    :param input_file_ids: Namespace object containing FileIDs for input files.
    :param target_genome: Name of the target genome for this job.
    :param target_two_bit_file_id: FileID for the target genome's 2bit file.
    :return: FileID of the generated chain file for this chromosome.
    """
    job.fileStore.logToMaster(f'Processing chromosome {target_genome}:{chrom} ({size:,} bp)', level=logging.INFO)
    
    try:
        # Create a temporary BED file for the chromosome
        bed_path = tools.fileOps.get_tmp_toil_file()
        with open(bed_path, 'w') as outf:
            tools.fileOps.print_row(outf, [chrom, 0, size])

        # Read necessary files from the global file store
        hal_path = job.fileStore.readGlobalFile(input_file_ids.hal)
        target_two_bit_path = job.fileStore.readGlobalFile(target_two_bit_file_id)
        query_two_bit_path = job.fileStore.readGlobalFile(input_file_ids.query_two_bit)

        # Define the pipeline of commands to run with improved error handling
        chain_path = tools.fileOps.get_tmp_toil_file()
        
        # First command: halLiftover
        cmd1 = ['halLiftover', '--outPSL', hal_path, ref_genome, bed_path, target_genome, '/dev/stdout']
        
        # Check if chromosome exists in HAL file for target genome
        try:
            # Test if the chromosome has any alignments
            test_cmd = ['halStats', '--chromSizes', target_genome, hal_path]
            test_output = tools.procOps.call_proc_lines(test_cmd)
            target_chroms = [line.split()[0] for line in test_output if line.strip()]
            
            if not target_chroms:
                job.fileStore.logToMaster(f'No chromosomes found for {target_genome} in HAL file', level=logging.WARNING)
                # Create empty chain file
                with open(chain_path, 'w') as f:
                    pass
                return job.fileStore.writeGlobalFile(chain_path)
                
        except Exception as e:
            job.fileStore.logToMaster(f'Warning: Could not verify chromosome existence for {target_genome}:{chrom}: {str(e)}', 
                                     level=logging.WARNING)
        
        # Run the full pipeline
        cmd_pipeline = [
            cmd1,
            ['pslPosTarget', '/dev/stdin', '/dev/stdout'],
            ['axtChain', '-psl', '-verbose=0', '-linearGap=medium', '/dev/stdin', target_two_bit_path, query_two_bit_path, chain_path]
        ]
        
        tools.procOps.run_proc(cmd_pipeline)
        
        # Validate output file
        if not os.path.exists(chain_path) or os.path.getsize(chain_path) == 0:
            job.fileStore.logToMaster(f'Warning: Empty or missing chain file for {target_genome}:{chrom}', 
                                     level=logging.WARNING)
            # Create minimal empty chain file
            with open(chain_path, 'w') as f:
                pass
        else:
            # Log success with file size
            file_size = os.path.getsize(chain_path)
            job.fileStore.logToMaster(f'Generated chain for {target_genome}:{chrom} ({file_size:,} bytes)', level=logging.INFO)

        # Write the resulting chain file to the global file store
        return job.fileStore.writeGlobalFile(chain_path)
        
    except Exception as e:
        job.fileStore.logToMaster(f'Error processing chromosome {target_genome}:{chrom}: {str(e)}', level=logging.ERROR)
        # Create empty chain file to prevent pipeline failure
        empty_chain_path = tools.fileOps.get_tmp_toil_file()
        with open(empty_chain_path, 'w') as f:
            pass
        return job.fileStore.writeGlobalFile(empty_chain_path)


def merge(job, chain_file_ids, genome, num_chromosomes):
    """
    Merges per-chromosome chain files into a single sorted chain file for a genome with enhanced processing.
    
    :param chain_file_ids: A list of FileIDs for the per-chromosome chain files.
    :param genome: The name of the genome whose chains are being merged.
    :param num_chromosomes: Number of chromosomes being merged.
    :return: FileID of the final merged and sorted chain file.
    """
    job.fileStore.logToMaster(f'Merging {len(chain_file_ids)} chain files for {genome}', level=logging.INFO)
    
    try:
        # Filter out empty chain files to improve merge efficiency
        valid_chain_files = []
        empty_files = 0
        
        # Create a file-of-filenames (fofn) for chainMergeSort
        fofn_path = tools.fileOps.get_tmp_toil_file()
        with open(fofn_path, 'w') as outf:
            for i, file_id in enumerate(chain_file_ids):
                # Read each chain file from the file store to a local path
                local_path = job.fileStore.readGlobalFile(file_id, userPath=f'{i}.chain')
                
                # Check if file has content
                if os.path.exists(local_path) and os.path.getsize(local_path) > 0:
                    if os.environ.get('CAT_BINARY_MODE') == 'singularity':
                        local_path = tools.procOps.singularify_arg(local_path)
                    outf.write(local_path + '\n')
                    valid_chain_files.append(local_path)
                else:
                    empty_files += 1
        
        if empty_files > 0:
            job.fileStore.logToMaster(f'Skipped {empty_files} empty chain files for {genome}', level=logging.INFO)
        
        if not valid_chain_files:
            job.fileStore.logToMaster(f'Warning: No valid chain files found for {genome}', level=logging.WARNING)
            # Create empty output file
            empty_chain_file = tools.fileOps.get_tmp_toil_file()
            with open(empty_chain_file, 'w') as f:
                pass
            return job.fileStore.writeGlobalFile(empty_chain_file)
        
        # Run chainMergeSort with enhanced error handling
        tmp_chain_file = tools.fileOps.get_tmp_toil_file()
        temp_dir = job.fileStore.getLocalTempDir()
        
        cmd = ['chainMergeSort', f'-inputList={fofn_path}', f'-tempDir={temp_dir}/']
        
        try:
            tools.procOps.run_proc(cmd, stdout=tmp_chain_file)
            
            # Validate merged output
            if os.path.exists(tmp_chain_file) and os.path.getsize(tmp_chain_file) > 0:
                merge_size = os.path.getsize(tmp_chain_file)
                job.fileStore.logToMaster(
                    f'Successfully merged {len(valid_chain_files)} chain files for {genome} '
                    f'({merge_size:,} bytes total)', 
                    level=logging.INFO
                )
            else:
                job.fileStore.logToMaster(f'Warning: chainMergeSort produced empty output for {genome}', 
                                         level=logging.WARNING)
                # Create empty file as fallback
                with open(tmp_chain_file, 'w') as f:
                    pass
                    
        except Exception as e:
            job.fileStore.logToMaster(f'Error during chainMergeSort for {genome}: {str(e)}', level=logging.ERROR)
            # Create empty file as fallback
            with open(tmp_chain_file, 'w') as f:
                pass

        # Write the final merged file to the global store and return its ID
        return job.fileStore.writeGlobalFile(tmp_chain_file)
        
    except Exception as e:
        job.fileStore.logToMaster(f'Error merging chains for {genome}: {str(e)}', level=logging.ERROR)
        # Create empty fallback file
        fallback_file = tools.fileOps.get_tmp_toil_file()
        with open(fallback_file, 'w') as f:
            pass
        return job.fileStore.writeGlobalFile(fallback_file)


def main():
    """
    Main entry point for the chaining script with enhanced validation and error handling.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

    # Add Toil-specific options
    Job.Runner.addToilOptions(parser)

    # === Script-specific arguments ===
    # Use groups for better help message organization
    chain_group = parser.add_argument_group("Chaining Mode Arguments")
    copy_group = parser.add_argument_group("Copy Existing Chains Mode Arguments")

    # --- Default Chaining Mode Arguments ---
    chain_group.add_argument("--hal", required=True, help="HAL alignment file.")
    chain_group.add_argument("--ref_genome", required=True, help="Reference genome name in HAL file.")
    chain_group.add_argument("--query_sizes", required=True, help="Chromosome sizes file for the reference genome.")
    chain_group.add_argument("--query_two_bit", required=True, help="2bit file for the reference genome.")
    chain_group.add_argument("--target_two_bit", nargs=2, metavar=("GENOME", "2BIT_PATH"), action="append",
                             help="Target genome name and its 2bit file. Can be specified multiple times.")
    
    # --- Chain Copying Mode ---
    copy_group.add_argument("--chain_mode", action="store_true",
                            help="Enable chain mode to copy existing chains instead of generating them.")
    copy_group.add_argument("--genome_chain", nargs=2, metavar=("GENOME", "CHAIN_PATH"), action="append",
                            help="Genome name and path to its existing chain file. Use with --chain_mode.")

    # --- Required Output Argument for both modes ---
    parser.add_argument("--chain_file", nargs=2, metavar=("GENOME", "OUTPUT_PATH"), required=True, action="append",
                        help="Target genome and the path for its final output chain file. Can be specified multiple times.")

    args = parser.parse_args()

    # --- Enhanced argument validation and processing ---
    if args.chain_mode:
        if not args.genome_chain:
            parser.error("--genome_chain is required when using --chain_mode.")
        if args.target_two_bit or args.hal or args.ref_genome or args.query_sizes or args.query_two_bit:
            parser.error("HAL/2bit/sizes arguments cannot be used with --chain_mode.")
    else:
        if not args.target_two_bit:
             parser.error("--target_two_bit is required unless in --chain_mode.")
    
    # Convert list-based arguments to dictionaries for easier handling
    target_two_bits = dict(args.target_two_bit) if args.target_two_bit else {}
    chain_files = dict(args.chain_file)
    genome_chains = dict(args.genome_chain) if args.genome_chain else {}
    
    # Validate input files exist
    if not args.chain_mode:
        required_files = [args.hal, args.query_sizes, args.query_two_bit]
        required_files.extend(target_two_bits.values())
        
        for file_path in required_files:
            if not os.path.exists(file_path):
                parser.error(f"Input file not found: {file_path}")
        
        # Validate HAL file contains reference genome
        try:
            from tools.hal import extract_hal_genomes
            hal_genomes = extract_hal_genomes(args.hal)
            if args.ref_genome not in hal_genomes:
                parser.error(f"Reference genome '{args.ref_genome}' not found in HAL file. Available: {', '.join(hal_genomes)}")
        except Exception:
            # If validation fails, continue with warning
            print(f"Warning: Could not validate reference genome in HAL file")
    else:
        # Validate chain files exist in chain mode
        for genome, chain_file in genome_chains.items():
            if not os.path.exists(chain_file):
                parser.error(f"Chain file not found for {genome}: {chain_file}")
    
    # Validate output directories can be created
    for genome, output_path in chain_files.items():
        output_dir = os.path.dirname(output_path)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir, exist_ok=True)
            except OSError as e:
                parser.error(f"Cannot create output directory for {genome}: {e}")
    
    # Cross-validate genome names between input and output
    if not args.chain_mode:
        input_genomes = set(target_two_bits.keys())
        output_genomes = set(chain_files.keys())
        
        if input_genomes != output_genomes:
            missing_outputs = input_genomes - output_genomes
            extra_outputs = output_genomes - input_genomes
            
            error_msg = []
            if missing_outputs:
                error_msg.append(f"Missing output chain files for genomes: {', '.join(missing_outputs)}")
            if extra_outputs:
                error_msg.append(f"Extra output chain files for genomes not in input: {', '.join(extra_outputs)}")
            
            parser.error(". ".join(error_msg))
    
    # Print processing summary
    if not args.chain_mode:
        print(f"Chaining pipeline will process {len(target_two_bits)} target genomes:")
        for genome in target_two_bits.keys():
            print(f"  - {genome}")
    else:
        print(f"Chain copy mode will process {len(genome_chains)} existing chain files:")
        for genome in genome_chains.keys():
            print(f"  - {genome}")

    # Launch the pipeline
    try:
        run_chaining_pipeline(
            hal=args.hal,
            ref_genome=args.ref_genome,
            query_sizes=args.query_sizes,
            query_two_bit=args.query_two_bit,
            target_two_bits=target_two_bits,
            chain_files=chain_files,
            chain_mode=args.chain_mode,
            genome_chains=genome_chains,
            toil_options=args  # Pass the complete args object for Toil
        )
        print("Chaining pipeline completed successfully!")
        
    except Exception as e:
        print(f"Error running chaining pipeline: {str(e)}")
        raise


if __name__ == "__main__":
    main()
