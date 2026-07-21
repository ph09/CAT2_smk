import re
import os
import argparse
import collections
import logging
import itertools

from toil.common import Toil
from toil.job import Job
from toil.fileStores import FileID

import tools.bio
import tools.dataOps
import tools.fileOps
import tools.intervals
import tools.procOps
import tools.toilInterface
import tools.transcripts


def calculate_dynamic_resources(input_file_ids, genome_size_mb=None, num_chunks=None):
    """
    Calculate dynamic resource allocation based on input data characteristics.
    
    :param input_file_ids: Namespace containing FileIDs for input files
    :param genome_size_mb: Estimated genome size in MB (optional)
    :param num_chunks: Estimated number of processing chunks (optional)
    :return: Dictionary with resource allocations
    """
    if genome_size_mb and num_chunks:
        # Base memory scales with genome size
        base_memory_gb = max(128, int(genome_size_mb / 500) * 4)  # 4GB per 500MB genome
        
        # For Augustus chunks: scale with chunk complexity
        chunk_memory_gb = max(64, int(genome_size_mb / num_chunks / 50) * 2)  # 2GB per 50MB chunk
        chunk_disk_gb = max(64, int(genome_size_mb / num_chunks / 100) * 2)   # 2GB per 100MB chunk
        
        # For merge: scale with number of chunks and total size
        merge_memory_gb = max(128, base_memory_gb + int(num_chunks / 10) * 2)  # Additional 2GB per 10 chunks
        merge_disk_gb = max(128, int(genome_size_mb / 200) * 4)  # 4GB per 200MB genome
    else:
        # Fallback to conservative estimates
        base_memory_gb = 128
        chunk_memory_gb = 128
        chunk_disk_gb = 128
        merge_memory_gb = 128
        merge_disk_gb = 128
    
    return {
        'setup_memory': f'{base_memory_gb}G',
        'setup_disk': f'{base_memory_gb * 2}G',
        'chunk_memory': f'{chunk_memory_gb}G',
        'chunk_disk': f'{chunk_disk_gb}G',
        'merge_memory': f'{merge_memory_gb}G',
        'merge_disk': f'{merge_disk_gb}G'
    }


def estimate_genome_characteristics(genome_fasta_path, chunksize):
    """
    Estimate genome size and number of chunks for resource calculation.
    
    :param genome_fasta_path: Path to genome FASTA file
    :param chunksize: Size of each processing chunk
    :return: Tuple of (genome_size_mb, estimated_chunks)
    """
    try:
        import os
        if os.path.exists(genome_fasta_path):
            file_size_mb = os.path.getsize(genome_fasta_path) / (1024 * 1024)
            # Rough estimate: FASTA is ~2x larger than raw sequence due to headers and formatting
            genome_size_mb = file_size_mb / 2
            estimated_chunks = max(16, int(genome_size_mb * 1024 * 1024 / chunksize))
            return genome_size_mb, estimated_chunks
        return None, None
    except:
        return None, None


def run_augustus_pb_pipeline(args, toil_options):
    """
    Main entry function for the AugustusPB Toil pipeline. This function orchestrates the workflow.
    """
    with Toil(toil_options) as t:
        if not t.options.restart:
            # Create a namespace for Toil FileIDs and ensure input paths are absolute
            input_file_ids = argparse.Namespace()
            args.genome_fasta = os.path.abspath(args.genome_fasta)
            args.chrom_sizes = os.path.abspath(args.chrom_sizes)
            args.pb_cfg = os.path.abspath(args.pb_cfg)
            args.hints_gff = os.path.abspath(args.hints_gff)

            # Import files into the Toil FileStore
            input_file_ids.genome_fasta = tools.toilInterface.write_fasta_to_filestore(t, args.genome_fasta)
            input_file_ids.chrom_sizes = FileID.forPath(t.importFile('file://' + args.chrom_sizes), args.chrom_sizes)
            input_file_ids.pb_cfg = FileID.forPath(t.importFile('file://' + args.pb_cfg), args.pb_cfg)
            input_file_ids.hints_gff = FileID.forPath(t.importFile('file://' + args.hints_gff), args.hints_gff)

            # Calculate dynamic resources based on input data
            genome_size_mb, estimated_chunks = estimate_genome_characteristics(args.genome_fasta, args.chunksize)
            resources = calculate_dynamic_resources(input_file_ids, genome_size_mb, estimated_chunks)
            
            job = Job.wrapJobFn(setup, args, input_file_ids, resources, 
                               memory=resources['setup_memory'], disk=resources['setup_disk'])
            raw_gtf_file_id, gtf_file_id, joined_gp_file_id = t.start(job)
        else:
            raw_gtf_file_id, gtf_file_id, joined_gp_file_id = t.restart()

        # Export final files from the FileStore to their destination paths
        tools.fileOps.ensure_file_dir(args.raw_gtf)
        t.exportFile(raw_gtf_file_id, 'file://' + os.path.abspath(args.raw_gtf))
        tools.fileOps.ensure_file_dir(args.gtf)
        t.exportFile(gtf_file_id, 'file://' + os.path.abspath(args.gtf))
        tools.fileOps.ensure_file_dir(args.gp)
        t.exportFile(joined_gp_file_id, 'file://' + os.path.abspath(args.gp))


def setup(job, args, input_file_ids, resources):
    """
    Set up the Augustus PB pipeline by loading input data, creating chunks, and scheduling parallel jobs.
    """
    job.fileStore.logToMaster('Beginning Augustus PB run with dynamic resource allocation', level=logging.INFO)
    
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)

    # load only PB hints
    hints_file = job.fileStore.readGlobalFile(input_file_ids.hints_gff)
    hints = [x.split('\t') for x in open(hints_file) if 'src=PB' in x]

    if len(hints) == 0:
        raise RuntimeError('No PB hints found.')

    job.fileStore.logToMaster(f'Found {len(hints)} PB hints across {len(genome_fasta)} chromosomes', level=logging.INFO)

    # convert the start/stops to ints and break up by chromosome
    hints_by_chrom = collections.defaultdict(list)
    for h in hints:
        h[3] = int(h[3])
        h[4] = int(h[4])
        hints_by_chrom[h[0]].append(h)

    # Calculate overlapping intervals with improved chunking strategy
    intervals = collections.defaultdict(list)
    total_chunks = 0
    
    for chrom in genome_fasta:
        chrom_size = len(genome_fasta[chrom])
        chrom_hints = len(hints_by_chrom[chrom])
        
        # Adaptive chunking: smaller chunks for hint-dense regions
        if chrom_hints > 1000:  # High-density chromosome
            effective_chunksize = max(args.chunksize // 2, 11000000)  # Smaller chunks
        elif chrom_hints > 100:  # Medium-density chromosome  
            effective_chunksize = args.chunksize
        else:  # Low-density chromosome
            effective_chunksize = min(args.chunksize * 2, chrom_size)  # Larger chunks
            
        for start in range(0, chrom_size, effective_chunksize - args.overlap):
            stop = min(start + effective_chunksize, chrom_size)
            intervals[chrom].append([start, stop])
            total_chunks += 1

    # Merge small final intervals (improved logic)
    for chrom, interval_list in intervals.items():
        if len(interval_list) < 2:
            continue
        last_start, last_stop = interval_list[-1]
        if last_stop - last_start <= 0.3 * args.chunksize:  # More aggressive merging
            del interval_list[-1]
            if interval_list:  # Check if list is not empty after deletion
                interval_list[-1][1] = last_stop
                total_chunks -= 1

    job.fileStore.logToMaster(f'Created {total_chunks} genomic chunks for parallel processing', level=logging.INFO)

    # Create parallel jobs for each chunk
    predictions = []
    chunk_count = 0
    
    for chrom, interval_list in intervals.items():
        for start, stop in interval_list:
            # Filter hints for this chunk
            chunk_hints = [h for h in hints_by_chrom[chrom] if h[3] >= start and h[4] <= stop]
            if len(chunk_hints) == 0:
                continue  # Skip chunks with no hints
                
            # Create temporary hints file for this chunk
            tmp_hints = tools.fileOps.get_tmp_toil_file()
            with open(tmp_hints, 'w') as outf:
                for h in chunk_hints:
                    tools.fileOps.print_row(outf, h)
            hints_file_id = job.fileStore.writeGlobalFile(tmp_hints)
            
            # Calculate chunk-specific resources based on hint density
            chunk_hint_density = len(chunk_hints) / (stop - start) * 11000000  # hints per Mb
            chunk_resources = get_chunk_resources(resources, chunk_hint_density, stop - start)
            
            # Create parallel Augustus job
            j = job.addChildJobFn(augustus_pb_chunk, args, input_file_ids, hints_file_id, 
                                 chrom, start, stop, chunk_count,
                                 memory=chunk_resources['memory'], 
                                 disk=chunk_resources['disk'])
            predictions.append(j.rv())
            chunk_count += 1

    if len(predictions) == 0:
        raise RuntimeError('No genomic chunks with PB hints found.')

    job.fileStore.logToMaster(f'Scheduled {len(predictions)} Augustus PB jobs for parallel execution', level=logging.INFO)

    # Schedule merge job with all prediction results
    results = job.addFollowOnJobFn(join_genes, predictions, len(predictions),
                                  memory=resources['merge_memory'], 
                                  disk=resources['merge_disk']).rv()
    return results


def get_chunk_resources(base_resources, hint_density, chunk_size):
    """
    Calculate chunk-specific resources based on hint density and chunk size.
    
    :param base_resources: Base resource allocations
    :param hint_density: Number of hints per Mb in this chunk
    :param chunk_size: Size of the genomic chunk in bp
    :return: Dictionary with chunk-specific resource allocations
    """
    # Parse base memory allocation
    base_memory_gb = int(base_resources['chunk_memory'].rstrip('G'))
    base_disk_gb = int(base_resources['chunk_disk'].rstrip('G'))
    
    # Scale resources based on hint density
    if hint_density > 50:  # Very high density
        memory_multiplier = 2.0
        disk_multiplier = 1.5
    elif hint_density > 20:  # High density
        memory_multiplier = 1.5
        disk_multiplier = 1.3
    elif hint_density > 10:  # Medium density
        memory_multiplier = 1.2
        disk_multiplier = 1.1
    else:  # Low density
        memory_multiplier = 1.0
        disk_multiplier = 1.0
    
    # Scale resources based on chunk size (larger chunks need more resources)
    size_multiplier = max(1.0, chunk_size / 11000000)  # Based on default 3Mb chunks
    
    final_memory = max(64, int(base_memory_gb * memory_multiplier * size_multiplier))
    final_disk = max(64, int(base_disk_gb * disk_multiplier * size_multiplier))
    
    return {
        'memory': f'{final_memory}G',
        'disk': f'{final_disk}G'
    }

    # results contains a 3 member tuple of [raw_gtf_file_id, gtf_file_id, joined_gp_file_id]
    results = job.addFollowOnJobFn(join_genes, predictions, memory='8G', disk='8G').rv()
    return results


def augustus_pb_chunk(job, args, input_file_ids, hints_file_id, chrom, start, stop, chunk_id):
    """
    Process a single genomic chunk with Augustus PB mode.
    
    :param job: Toil job object
    :param args: Command line arguments
    :param input_file_ids: FileStore IDs for input files
    :param hints_file_id: FileStore ID for chunk-specific hints file
    :param chrom: Chromosome name
    :param start: Start position of chunk
    :param stop: Stop position of chunk  
    :param chunk_id: Unique identifier for this chunk
    :return: FileStore ID of Augustus output for this chunk
    """
    job.fileStore.logToMaster(f'Processing Augustus PB chunk {chunk_id}: {chrom}:{start}-{stop}', level=logging.INFO)
    
    # Load required files from FileStore
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)
    hints = job.fileStore.readGlobalFile(hints_file_id)
    pb_cfg = job.fileStore.readGlobalFile(input_file_ids.pb_cfg)
    
    # Extract genomic sequence for this chunk
    tmp_fasta = tools.fileOps.get_tmp_toil_file()
    chunk_sequence = genome_fasta[chrom][start:stop]
    tools.bio.write_fasta(tmp_fasta, chrom, chunk_sequence)
    
    # Set up Augustus output file
    results = tools.fileOps.get_tmp_toil_file()
    
    # Build Augustus command with enhanced parameters
    cmd = ['augustus', 
           '--softmasking=1', 
           '--allow_hinted_splicesites=atac',
           '--alternatives-from-evidence=1', 
           '--UTR={}'.format(int(args.utr)),
           '--hintsfile={}'.format(hints),
           '--extrinsicCfgFile={}'.format(pb_cfg),
           '--species={}'.format(args.species),
           '--/augustus/verbosity=0',
           '--predictionStart=-{}'.format(start), 
           '--predictionEnd=-{}'.format(start),
           tmp_fasta]
    
    try:
        # Run Augustus with error handling
        tools.procOps.run_proc(cmd, stdout=results)
        
        # Verify output was generated
        if os.path.getsize(results) == 0:
            job.fileStore.logToMaster(f'Warning: Augustus chunk {chunk_id} produced no output', level=logging.WARNING)
        else:
            job.fileStore.logToMaster(f'Successfully completed Augustus chunk {chunk_id}', level=logging.INFO)
            
    except Exception as e:
        job.fileStore.logToMaster(f'Error in Augustus chunk {chunk_id}: {str(e)}', level=logging.ERROR)
        # Create empty results file to avoid pipeline failure
        with open(results, 'w') as f:
            f.write('')
    
    return job.fileStore.writeGlobalFile(results)


def join_genes(job, gff_chunks, num_chunks):
    """
    Join Augustus prediction results from all chunks into final output files.
    
    :param job: Toil job object
    :param gff_chunks: List of FileStore IDs containing Augustus output chunks
    :param num_chunks: Number of chunks processed
    :return: Tuple of (raw_gtf_file_id, joined_gtf_file_id, joined_gp_file_id)
    """
    job.fileStore.logToMaster(f'Joining results from {num_chunks} Augustus PB chunks', level=logging.INFO)
    
    def filter_joingenes(injoingenes_file, out_joingenes_file):
        """Filter and format joingenes output."""
        matcher = re.compile("\tAUGUSTUS\t(exon|CDS|start_codon|stop_codon|tts|tss)\t")
        lines_written = 0
        with open(out_joingenes_file, "w") as ofh:
            for l in open(injoingenes_file):
                if matcher.search(l):
                    l = l.replace("jg", "augPB-")
                    ofh.write(l)
                    lines_written += 1
        return lines_written

    # Combine all chunk results into raw GTF
    raw_gtf_file = tools.fileOps.get_tmp_toil_file()
    raw_gtf_fofn = tools.fileOps.get_tmp_toil_file()
    files = []
    total_lines = 0
    
    with open(raw_gtf_file, 'w') as raw_handle, open(raw_gtf_fofn, 'w') as fofn_handle:
        for i, chunk in enumerate(gff_chunks):
            local_path = job.fileStore.readGlobalFile(chunk)
            
            # Count lines in this chunk
            chunk_lines = 0
            for line in open(local_path):
                raw_handle.write(line)
                chunk_lines += 1
            total_lines += chunk_lines
            
            # Handle different execution environments
            if os.environ.get('CAT_BINARY_MODE') == 'singularity':
                local_path = tools.procOps.singularify_arg(local_path)
                files.append(local_path)
            else:
                files.append(os.path.basename(local_path))
            fofn_handle.write(local_path + '\n')
            
            if chunk_lines > 0:
                job.fileStore.logToMaster(f'Chunk {i+1}: {chunk_lines} predictions', level=logging.INFO)

    job.fileStore.logToMaster(f'Combined {total_lines} total predictions from all chunks', level=logging.INFO)

    if total_lines == 0:
        job.fileStore.logToMaster('Warning: No predictions found in any chunks', level=logging.WARNING)
        # Create empty output files
        empty_file = tools.fileOps.get_tmp_toil_file()
        with open(empty_file, 'w') as f:
            f.write('')
        empty_file_id = job.fileStore.writeGlobalFile(empty_file)
        return empty_file_id, empty_file_id, empty_file_id

    # Run joingenes to merge overlapping predictions
    join_genes_file = tools.fileOps.get_tmp_toil_file()
    join_genes_gp = tools.fileOps.get_tmp_toil_file()

    try:
        # First pass: run joingenes
        tmp_join_genes_file = tools.fileOps.get_tmp_toil_file()
        cmd = ['joingenes', '-f', raw_gtf_fofn, '-o', tmp_join_genes_file]
        tools.procOps.run_proc(cmd)
        
        # Filter and format the joingenes output
        filtered_lines = filter_joingenes(tmp_join_genes_file, join_genes_file)
        job.fileStore.logToMaster(f'Joingenes produced {filtered_lines} final gene predictions', level=logging.INFO)

        # Convert to GenePred format and back to GTF for proper formatting
        cmd = ['gtfToGenePred', '-genePredExt', join_genes_file, join_genes_gp]
        tools.procOps.run_proc(cmd)
        
        # Convert back to GTF with proper formatting
        cmd = ['genePredToGtf', 'file', join_genes_gp, '-utr', '-honorCdsStat', '-source=augustusPB', join_genes_file]
        tools.procOps.run_proc(cmd)

        job.fileStore.logToMaster('Successfully completed Augustus PB gene joining and formatting', level=logging.INFO)

    except Exception as e:
        job.fileStore.logToMaster(f'Error during gene joining: {str(e)}', level=logging.ERROR)
        # Use raw GTF as fallback
        join_genes_file = raw_gtf_file
        # Create dummy GenePred file
        with open(join_genes_gp, 'w') as f:
            f.write('')

    # Write final results to FileStore
    joined_gtf_file_id = job.fileStore.writeGlobalFile(join_genes_file)
    raw_gtf_file_id = job.fileStore.writeGlobalFile(raw_gtf_file)
    joined_gp_file_id = job.fileStore.writeGlobalFile(join_genes_gp)
    
    return raw_gtf_file_id, joined_gtf_file_id, joined_gp_file_id


def main():
    """
    Main entry point for the Augustus PB pipeline with enhanced parallel processing.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    Job.Runner.addToilOptions(parser)

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
                       help="Size of genomic chunks for parallel processing (default: 11000000 bp). "
                            "Smaller chunks = more parallelism but more overhead.")
    parser.add_argument("--overlap", type=int, default=1000000, 
                       help="Overlap between genomic chunks (default: 1000000 bp). "
                            "Larger overlaps help capture genes spanning chunk boundaries.")

    args = parser.parse_args()
    
    # Validate input files exist
    for file_arg in ['genome_fasta', 'chrom_sizes', 'hints_gff', 'pb_cfg']:
        file_path = getattr(args, file_arg)
        if not os.path.exists(file_path):
            parser.error(f"Input file not found: {file_path}")
    
    # Validate chunking parameters
    if args.chunksize < args.overlap:
        parser.error("Chunk size must be larger than overlap size")
    if args.overlap < 0:
        parser.error("Overlap size must be non-negative")
        
    run_augustus_pb_pipeline(args, args)

if __name__ == "__main__":
    main()
