import argparse
import collections
import itertools
import logging
import os
import sys

from toil.common import Toil
from toil.job import Job
from toil.fileStores import FileID
from toil.lib.humanize import human2bytes

import tools.bio
import tools.dataOps
import tools.fileOps
import tools.nameConversions
import tools.parasail_wrapper
import tools.procOps
import tools.sqlInterface
import tools.toilInterface
import tools.transcripts


def calculate_dynamic_resources(input_file_ids, num_transcripts=None, chunk_size=1000):
    """
    Calculate dynamic memory and disk requirements based on input data size and transcript count.
    
    :param input_file_ids: File IDs for input files
    :param num_transcripts: Estimated number of transcripts to process
    :param chunk_size: Number of transcripts per chunk
    :return: dict with memory and disk requirements
    """
    # Calculate base disk usage from input files
    base_disk = tools.toilInterface.find_total_disk_usage(input_file_ids, buffer='8G', round='8G')
    base_disk_gb = max(32, int(base_disk / (1024**3)))  # Convert to GB, minimum 2GB
    
    # Dynamic memory calculation based on transcript count
    if num_transcripts:
        # Base memory: 4GB + 100MB per 1000 transcripts
        base_memory_gb = 32 + max(1, int(num_transcripts / 1000) * 0.1)
        
        # For alignment chunks: scale with chunk size
        chunk_memory_gb = max(64, int(chunk_size / 500) * 2)  # 2GB per 500 transcripts
        chunk_disk_gb = max(16, base_disk_gb // 4)  # Quarter of total disk per chunk
        
        # For merge: scale with total output size
        merge_memory_gb = max(64, base_memory_gb * 2)
        merge_disk_gb = max(32, base_disk_gb)
    else:
        # Fallback to conservative estimates
        base_memory_gb = 128
        chunk_memory_gb = 128
        chunk_disk_gb = 128
        merge_memory_gb = 128
        merge_disk_gb = 128
    
    return {
        'setup_memory': f'{base_memory_gb}G',
        'chunk_memory': f'{chunk_memory_gb}G',
        'chunk_disk': f'{chunk_disk_gb}G',
        'merge_memory': f'{merge_memory_gb}G',
        'merge_disk': f'{merge_disk_gb}G'
    }


def estimate_transcript_count(transcript_modes, annotation_gp):
    """
    Estimate the total number of transcripts to be processed.
    
    :param transcript_modes: Dictionary of transcript modes and their files
    :param annotation_gp: Path to annotation genePred file
    :return: Estimated transcript count
    """
    try:
        # Quick estimate based on annotation file size and average transcript modes
        import os
        if os.path.exists(annotation_gp):
            file_size_mb = os.path.getsize(annotation_gp) / (1024 * 1024)
            # Rough estimate: ~100 transcripts per MB of genePred file
            base_count = int(file_size_mb * 100)
            # Multiply by number of modes
            total_count = base_count * len(transcript_modes)
            return total_count
        return None
    except:
        return None


def run_align_transcripts_pipeline(args, toil_options):
    """
    Main entry function for the transcript alignment Toil pipeline.
    """
    # Re-structure the mode files argument into a more usable dictionary
    # The new structure will be: {'transMap': {'gp': '/path', 'mRNA': '/path', 'CDS': '/path'}, ...}
    args.transcript_modes = collections.defaultdict(dict)
    for mode, gp_path, mrna_path, cds_path in args.mode_files:
        args.transcript_modes[mode]['gp'] = gp_path
        args.transcript_modes[mode]['mRNA'] = mrna_path
        args.transcript_modes[mode]['CDS'] = cds_path

    with Toil(toil_options) as t:
        if not t.options.restart:
            input_file_ids = argparse.Namespace()
            
            # Make paths absolute and import files into the Toil filestore
            input_file_ids.ref_genome_fasta = tools.toilInterface.write_fasta_to_filestore(t, os.path.abspath(args.ref_genome_fasta))
            input_file_ids.genome_fasta = tools.toilInterface.write_fasta_to_filestore(t, os.path.abspath(args.genome_fasta))
            input_file_ids.annotation_gp = t.importFile('file://' + os.path.abspath(args.annotation_gp))
            input_file_ids.ref_db = t.importFile('file://' + os.path.abspath(args.ref_db_path))
            
            input_file_ids.modes = {}
            for mode, paths in args.transcript_modes.items():
                input_file_ids.modes[mode] = t.importFile('file://' + os.path.abspath(paths['gp']))
            
            # Calculate dynamic resources based on input data
            transcript_count = estimate_transcript_count(args.transcript_modes, args.annotation_gp)
            resources = calculate_dynamic_resources(input_file_ids, transcript_count)
            
            job = Job.wrapJobFn(setup, args, input_file_ids, resources, memory=resources['setup_memory'])
            results_file_ids = t.start(job)
        else:
            results_file_ids = t.restart()

        # Export all resulting PSL files from the filestore
        for file_path, file_id in results_file_ids.items():
            tools.fileOps.ensure_file_dir(file_path)
            t.exportFile(file_id, 'file://' + os.path.abspath(file_path))


def setup(job, args, input_file_ids, resources):
    job.fileStore.logToMaster('Beginning Align Transcripts run on {}'.format(args.genome), level=logging.INFO)
    # load all fileStore files necessary
    annotation_gp = job.fileStore.readGlobalFile(input_file_ids.annotation_gp)
    ref_genome_db = job.fileStore.readGlobalFile(input_file_ids.ref_db)
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)
    ref_genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.ref_genome_fasta,
                                                                     prefix='ref_genome', upper=False)
    # load required reference data into memory
    tx_biotype_map = tools.sqlInterface.get_transcript_biotype_map(ref_genome_db)
    ref_transcript_dict = tools.transcripts.get_gene_pred_dict(annotation_gp)
    # will hold a mapping of output file paths to lists of Promise objects containing output
    results = collections.defaultdict(list)
    for tx_mode in ['transMap', 'transMap_pairwise', 'augTM', 'augTM_pairwise', 'augTMR', 'augTMR_pairwise', 'augMP', 'txTM']:
        if tx_mode not in args.transcript_modes:
            continue
        # output file paths
        mrna_path = args.transcript_modes[tx_mode]['mRNA']
        cds_path = args.transcript_modes[tx_mode]['CDS']
        # begin loading transcripts and sequences
        gp_path = job.fileStore.readGlobalFile(input_file_ids.modes[tx_mode])
        transcript_dict = tools.transcripts.get_gene_pred_dict(gp_path)
        # Filter protein-coding transcripts (including NMD and LoF variants)
        protein_coding_biotypes = {
            'protein_coding',
            'nonsense_mediated_decay',
            'non_stop_decay',
            'protein_coding_LoF',
            'protein_coding_CDS_not_defined'
        }
        filtered_transcripts = {}
        for aln_id, tx in transcript_dict.items():
            ref_id = tools.nameConversions.alignment_id_to_ref_transcript_id(aln_id)
            biotype = tx_biotype_map.get(ref_id)
            if biotype in protein_coding_biotypes or biotype is None:
                filtered_transcripts[aln_id] = tx
        # Fallback: if filtering removes everything but there are inputs, use all transcripts
        if not filtered_transcripts and transcript_dict:
            filtered_transcripts = transcript_dict
        transcript_dict = filtered_transcripts
        for aln_mode, out_path in zip(*[['mRNA', 'CDS'], [mrna_path, cds_path]]):
            seq_iter = get_alignment_sequences(transcript_dict, ref_transcript_dict, genome_fasta,
                                               ref_genome_fasta, aln_mode)
            for chunk in group_transcripts(seq_iter):
                j = job.addChildJobFn(run_aln_chunk, chunk, 
                                    memory=resources['chunk_memory'], 
                                    disk=resources['chunk_disk'])
                results[out_path].append(j.rv())

    if len(results) == 0:
        err_msg = 'Align Transcripts pipeline did not detect any input genePreds for {}'.format(args.genome)
        raise RuntimeError(err_msg)
    # convert the results Promises into resolved values
    return job.addFollowOnJobFn(merge, results, args, 
                               memory=resources['merge_memory'], 
                               disk=resources['merge_disk']).rv()


def get_alignment_sequences(transcript_dict, ref_transcript_dict, genome_fasta, ref_genome_fasta, mode):
    """Generator that yields tuples of (tx_id, tx_seq, ref_tx_id, ref_tx_seq) for alignment."""
    assert mode in ['mRNA', 'CDS']
    for tx_id, tx in transcript_dict.items():
        ref_tx_id = tools.nameConversions.alignment_id_to_ref_transcript_id(tx_id)
        # Some inputs (e.g., raw Augustus IDs like g1.t1) may not map to a reference transcript.
        # Skip those gracefully instead of raising a KeyError.
        ref_tx = ref_transcript_dict.get(ref_tx_id)
        if ref_tx is None:
            continue
        tx_seq = tx.get_mrna(genome_fasta) if mode == 'mRNA' else tx.get_cds(genome_fasta)
        ref_tx_seq = ref_tx.get_mrna(ref_genome_fasta) if mode == 'mRNA' else ref_tx.get_cds(ref_genome_fasta)
        # Parasail has issues with very short sequences
        if len(ref_tx_seq) > 20 and len(tx_seq) > 20:
            yield tx_id, tx_seq, ref_tx_id, ref_tx_seq


def run_aln_chunk(job, chunk):
    """
    Runs a chunk of sequences through parasail alignment.
    """
    results = []
    for tx_id, tx_seq, ref_tx_id, ref_tx_seq in chunk:
        # Perform nucleotide alignment
        p = tools.parasail_wrapper.aln_nucleotides(tx_seq, tx_id, ref_tx_seq, ref_tx_id)
        # Convert result to PSL format
        psl_str = '\t'.join(p.psl_string())
        results.append(psl_str)
    return results


def merge(job, results, genome):
    """
    Merges the PSL results from all alignment chunks into final output files.
    """
    job.fileStore.logToMaster(f'Merging alignment output for {genome}', level=logging.INFO)
    results_file_ids = {}
    for out_path, result_promises in results.items():
        tmp_results_file = tools.fileOps.get_tmp_toil_file()
        with open(tmp_results_file, 'w') as outf:
            # Chain together the lists of PSL strings from each chunk
            for psl_line in itertools.chain.from_iterable(result_promises):
                if psl_line:
                    outf.write(psl_line + '\n')
        results_file_ids[out_path] = job.fileStore.writeGlobalFile(tmp_results_file)
    return results_file_ids


def group_transcripts(tx_iter, num_bases=10**6, max_seqs=1000):
    """
    Greedily groups transcripts into bins of approximately num_bases, without exceeding max_seqs.
    This helps balance the workload for each parallel job.
    """
    try:
        tx_id, tx_seq, ref_tx_id, ref_tx_seq = next(tx_iter)
    except StopIteration:
        return # Handle empty iterator
    this_bin = [(tx_id, tx_seq, ref_tx_id, ref_tx_seq)]
    bin_base_count = len(tx_seq)
    num_seqs = 1
    for tx_id, tx_seq, ref_tx_id, ref_tx_seq in tx_iter:
        bin_base_count += len(tx_seq)
        num_seqs += 1
        if bin_base_count >= num_bases or num_seqs >= max_seqs:
            yield this_bin
            this_bin = []
            bin_base_count = 0
            num_seqs = 0
        this_bin.append((tx_id, tx_seq, ref_tx_id, ref_tx_seq))
    if this_bin:
        yield this_bin


def main():
    """Parse arguments and launch the alignment pipeline.

    Supports two top-level modes:

    - ``--mode toil``    – original Toil-based pipeline used by ``execution_mode: local``.
    - ``--mode cluster`` – submit array jobs via cat.scheduler (SLURM or SGE),
                          selected by ``--execution-mode`` (defaults to slurm).
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--mode", choices=['toil', 'cluster'], default='toil',
                        help="'toil' = single-machine via Toil; 'cluster' = SLURM or SGE array jobs.")

    parser.add_argument("--ref-genome-fasta", required=True, help="FASTA file for the reference genome.")
    parser.add_argument("--genome-fasta", required=True, help="FASTA file for the target genome.")
    parser.add_argument("--annotation-gp", required=True, help="GenePred annotation for the reference genome.")
    parser.add_argument("--ref-db-path", required=True, help="Path to the reference genome database.")
    parser.add_argument("--genome", required=True, help="Name of the target genome.")

    parser.add_argument("--mode-files", nargs=4, action='append', required=True,
                        metavar=("MODE", "INPUT_GP", "MRNA_PSL", "CDS_PSL"),
                        help="Define a transcript mode and its files. Provide 4 values: "
                             "the mode name (e.g., 'augTM'), the input genePred path, "
                             "the output mRNA PSL path, and the output CDS PSL path. "
                             "This argument can be specified multiple times for different modes.")

    # Cluster-mode arguments. Names mirror cat/align_transcripts_cluster.py so
    # we can forward them through unchanged.
    parser.add_argument("--execution-mode", choices=("auto", "slurm", "sge", "local"), default="auto")
    parser.add_argument("--partition", default="short")
    parser.add_argument("--exclude-nodes", default="")
    parser.add_argument("--module-load", default="")
    parser.add_argument("--sge-parallel-env", default="smp")
    parser.add_argument("--sge-memory-flag", default="h_vmem")
    parser.add_argument("--memory", type=int, default=16)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--time", default="01:00:00")
    parser.add_argument("--max-jobs", type=int, default=200)
    parser.add_argument("--timeout-hours", type=int, default=12)
    parser.add_argument("--chunk-size", type=int, default=500)
    parser.add_argument("--cleanup", action="store_true")

    if '--mode' in sys.argv and sys.argv[sys.argv.index('--mode') + 1] == 'cluster':
        args = parser.parse_args()
        from cat.align_transcripts_cluster import run_cluster_alignment_pipeline
        run_cluster_alignment_pipeline(args)
    else:
        # Toil mode adds Toil's own CLI options on top.
        Job.Runner.addToilOptions(parser)
        args = parser.parse_args()
        run_align_transcripts_pipeline(args, args)


if __name__ == "__main__":
    main()
