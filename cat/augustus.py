import argparse
import itertools
import os
import logging
import collections

from toil.fileStores import FileID
from toil.common import Toil
from toil.job import Job

import tools.bio
import tools.dataOps
import tools.fileOps
import tools.intervals
import tools.nameConversions
import tools.procOps
import tools.psl
import tools.misc
import tools.tm2hints
import tools.toilInterface
import tools.transcripts
from tools.hintsDatabaseInterface import reflect_hints_db, get_rnaseq_hints


def estimate_transcript_characteristics(coding_gp_path):
    """
    Estimate transcript count and complexity for resource calculation.
    
    :param coding_gp_path: Path to coding GenePred file
    :return: Estimated transcript count
    """
    try:
        if os.path.exists(coding_gp_path):
            # Count lines in GenePred file (each line = one transcript)
            with open(coding_gp_path, 'r') as f:
                transcript_count = sum(1 for line in f if line.strip())
            return transcript_count
        return None
    except:
        return None

def run_augustus_pipeline(args, toil_options):
    """
    Main entry function for the Augustus TM/TMR Toil pipeline with static resource allocation.
    """
    with Toil(toil_options) as t:
        if not t.options.restart:
            input_file_ids = argparse.Namespace()
            
            # Make paths absolute
            args.genome_fasta = os.path.abspath(args.genome_fasta)
            args.tm_cfg = os.path.abspath(args.tm_cfg)
            args.coding_gp = os.path.abspath(args.coding_gp)
            args.ref_psl = os.path.abspath(args.ref_psl)
            args.tm_psl = os.path.abspath(args.filtered_tm_psl)
            args.annotation_gp = os.path.abspath(args.annotation_gp)
            args.miniprot_hints_gff = os.path.abspath(args.miniprot_hints_gff)
            
            file_ids = [args.genome_fasta, args.coding_gp, args.ref_psl,
                        args.tm_psl, args.annotation_gp, args.tm_cfg, args.miniprot_hints_gff]
            
            # Handle TMR mode files
            if args.augustus_tmr_gtf:
                args.augustus_hints_db = os.path.abspath(args.augustus_hints_db)
                input_file_ids.augustus_hints_db = FileID.forPath(t.importFile('file://' + args.augustus_hints_db), args.augustus_hints_db)
                args.tmr_cfg = os.path.abspath(args.tmr_cfg)
                input_file_ids.tmr_cfg = FileID.forPath(t.importFile('file://' + args.tmr_cfg), args.tmr_cfg)
                file_ids.append(args.augustus_hints_db)
            
            # Import files to FileStore
            input_file_ids.genome_fasta = tools.toilInterface.write_fasta_to_filestore(t, args.genome_fasta)
            input_file_ids.tm_cfg = FileID.forPath(t.importFile('file://' + args.tm_cfg), args.tm_cfg)
            input_file_ids.coding_gp = FileID.forPath(t.importFile('file://' + args.coding_gp), args.coding_gp)
            input_file_ids.ref_psl = FileID.forPath(t.importFile('file://' + args.ref_psl), args.ref_psl)
            input_file_ids.tm_psl = FileID.forPath(t.importFile('file://' + args.tm_psl), args.tm_psl)
            input_file_ids.annotation_gp = FileID.forPath(t.importFile('file://' + args.annotation_gp), args.annotation_gp)
            input_file_ids.miniprot_hints_gff = FileID.forPath(t.importFile('file://' + args.miniprot_hints_gff), args.miniprot_hints_gff)

            # Calculate dynamic resources based on transcript count
            disk_usage = tools.toilInterface.find_total_disk_usage(file_ids)
            job = Job.wrapJobFn(setup, args, input_file_ids, disk_usage, disk=disk_usage)
            tm_file_id, tmr_file_id = t.start(job)
        else:
            tm_file_id, tmr_file_id = t.restart()
            
        # Export output files
        tools.fileOps.ensure_file_dir(args.augustus_tm_gtf)
        t.exportFile(tm_file_id, 'file://' + os.path.abspath(args.augustus_tm_gtf))
        if tmr_file_id is not None:
            tools.fileOps.ensure_file_dir(args.augustus_tmr_gtf)
            t.exportFile(tmr_file_id, 'file://' + os.path.abspath(args.augustus_tmr_gtf))


def setup(job, args, input_file_ids, disk_usage):
    """
    Set up the Augustus pipeline by loading data and creating optimally-sized chunks.
    """
    def start_jobs(mode, chunk_size, cfg_file_id):
        results = []
        for chunk in tools.dataOps.grouper(iter(tx_dict.items()), chunk_size):
            grouped_recs = {}
            for tx_id, tx in chunk:
                grouped_recs[tx_id] = [tx, ref_tx_dict[tools.nameConversions.remove_alignment_number(tx_id)], tm_psl_dict[tx_id], ref_psl_dict[tools.nameConversions.remove_alignment_number(tx_id)]]
            j = job.addChildJobFn(run_augustus_chunk, args, grouped_recs, input_file_ids, mode, cfg_file_id, disk=disk_usage)
            results.append(j.rv())
            
        return results
    
    # Load required input files
    ref_psl = job.fileStore.readGlobalFile(input_file_ids.ref_psl)
    tm_psl = job.fileStore.readGlobalFile(input_file_ids.tm_psl)
    annotation_gp = job.fileStore.readGlobalFile(input_file_ids.annotation_gp)
    coding_gp = job.fileStore.readGlobalFile(input_file_ids.coding_gp)
    
    # Create lookup dictionaries
    ref_psl_dict = tools.psl.get_alignment_dict(ref_psl)
    tm_psl_dict = tools.psl.get_alignment_dict(tm_psl)
    ref_tx_dict = tools.transcripts.get_gene_pred_dict(annotation_gp)
    tx_dict = tools.transcripts.get_gene_pred_dict(coding_gp)
    
    # Start TM mode jobs
    tm_results = start_jobs('TM', 25, input_file_ids.tm_cfg)
    
    # Start TMR mode jobs if requested
    if args.augustus_tmr_gtf:
        tmr_results = start_jobs('TMR', 15, input_file_ids.tmr_cfg)
    else:
        tmr_results = None
    return job.addFollowOnJobFn(merge, tm_results, tmr_results).rv()


def run_augustus_chunk(job, args, grouped_recs, input_file_ids, mode, cfg_file_id, padding=20000):
    """
    Process a chunk of transcripts with Augustus in specified mode.
    """
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)
    cfg_file = job.fileStore.readGlobalFile(cfg_file_id)
    miniprot_gff = job.fileStore.readGlobalFile(input_file_ids.miniprot_hints_gff)
    
    with open(miniprot_gff) as mpf:
        mp_hints = mpf.read()
    
    if args.augustus_tmr_gtf:
        try:
            hints_db_file = job.fileStore.readGlobalFile(input_file_ids.augustus_hints_db)
            speciesnames, seqnames, hints, featuretypes, session = reflect_hints_db(hints_db_file)
        except Exception as e:
            job.fileStore.logToMaster(f'Error loading hints database for TMR mode: {str(e)}', level=logging.ERROR)
            mode = 'TM'

    results = []
    for tm_tx, ref_tx, tm_psl, ref_psl in grouped_recs.values():
        if len(tm_tx) > 3 * 10 ** 6:
            continue
            
        chromosome = tm_tx.chromosome
        start = max(tm_tx.start - padding, 0)
        stop = min(tm_tx.stop + padding, len(genome_fasta[chromosome]))
            
        tm_hints = tools.tm2hints.tm_to_hints(tm_tx, tm_psl, ref_psl)
        hint_components = [tm_hints.strip(), mp_hints.strip()]
        if args.augustus_tmr_gtf:
            rnaseq_hints = get_rnaseq_hints(args.genome, chromosome, start, stop, 
                                                    speciesnames, seqnames, hints, featuretypes, session)
            hint_components.append(rnaseq_hints.strip())
            hint = '\n'.join(h for h in hint_components if h)
        else:
            hint = '\n'.join(h for h in hint_components if h)

            
        transcript = run_augustus(hint, genome_fasta, tm_tx, cfg_file, start, stop, 
                                     args.augustus_species, mode, args.utr)
        if transcript is not None:
            results.extend(transcript)
            
    if args.augustus_tmr_gtf:
        session.close()
    return results



def run_augustus(hint, fasta, tm_tx, cfg_file, start, stop, species, mode, utr):
    """
    Run Augustus gene prediction with enhanced error handling.
    
    :param hint: Hint string for Augustus
    :param fasta: Genome FASTA dictionary
    :param tm_tx: TransMap transcript object
    :param cfg_file: Augustus configuration file path
    :param start: Start position (0-based)
    :param stop: Stop position (0-based)
    :param species: Augustus species parameter
    :param mode: Prediction mode ('TM' or 'TMR')
    :param utr: UTR prediction flag
    :param tx_id: Transcript identifier for logging
    :return: List of GTF records or None if failed
    """
    # Create temporary FASTA file for this transcript region
    tmp_fasta = tools.fileOps.get_tmp_toil_file()
    tools.bio.write_fasta(tmp_fasta, tm_tx.chromosome, fasta[tm_tx.chromosome][start:stop])
        
    # Create temporary hints file
    hints_out = tools.fileOps.get_tmp_toil_file()
    with open(hints_out, 'w') as outf:
        outf.write(hint)

    # Build Augustus command
    cmd = ['augustus', tmp_fasta, 
           '--predictionStart=-{}'.format(start), 
           '--predictionEnd=-{}'.format(start),
           '--extrinsicCfgFile={}'.format(cfg_file), 
           '--hintsfile={}'.format(hints_out), 
           '--UTR={}'.format(int(utr)),
           '--alternatives-from-evidence=0', 
           '--species={}'.format(species), 
           '--allow_hinted_splicesites=atac',
           '--protein=0', 
           '--softmasking=1', 
           '--/augustus/verbosity=0']
    aug_output = tools.procOps.call_proc_lines(cmd)
    transcript = munge_augustus_output(aug_output, mode, tm_tx)
    return transcript


def merge(job, tm_results, tmr_results):
    """
    Merge Augustus prediction results from all chunks.
    
    :param job: Toil job object
    :param tm_results: List of TM mode results from chunks
    :param tmr_results: List of TMR mode results from chunks (or None)
    :return: Tuple of (tm_file_id, tmr_file_id)
    """
    tmp_results_file = tools.fileOps.get_tmp_toil_file()
    tools.fileOps.print_rows(tmp_results_file, list(itertools.chain.from_iterable(tm_results)))
    tm_results_file_id = job.fileStore.writeGlobalFile(tmp_results_file)
    if tmr_results is not None:
        tmp_tmr_results_file = tools.fileOps.get_tmp_toil_file()
        tools.fileOps.print_rows(tmp_tmr_results_file, list(itertools.chain.from_iterable(tmr_results)))
        tmr_results_file_id = job.fileStore.writeGlobalFile(tmp_tmr_results_file)
    else:
        tmr_results_file_id = None
    return tm_results_file_id, tmr_results_file_id


def munge_augustus_output(aug_output, mode, tm_tx):
    """
    Process Augustus output and convert to GTF format with enhanced error handling.
    
    :param aug_output: List of Augustus output lines
    :param mode: Prediction mode ('TM' or 'TMR')
    :param tm_tx: TransMap transcript object
    :return: List of GTF records or None if processing failed
    """
    tx_entries = [x.split() for x in aug_output if "\ttranscript\t" in x]
    valid_txs = [x[-1] for x in tx_entries if tm_tx.interval.overlap(tools.intervals.ChromosomeInterval(x[0], x[3],
                                                                                                        x[4], x[6]))]
    if len(valid_txs) != 1:
        return None  # Expect exactly one valid transcript
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


def main():
    """
    Main entry point for the Augustus TM/TMR pipeline with enhanced parallel processing.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    Job.Runner.addToilOptions(parser)

    # Input Files
    parser.add_argument("--genome_fasta", required=True, 
                       help="Genome FASTA file for Augustus gene prediction.")
    parser.add_argument("--coding_gp", required=True, 
                       help="GenePred file containing coding transcripts from TransMap.")
    parser.add_argument("--filtered_tm_psl", required=True, 
                       help="Filtered TransMap PSL alignment file.")
    parser.add_argument("--ref_psl", required=True, 
                       help="Reference genome PSL alignment file.")
    parser.add_argument("--annotation_gp", required=True, 
                       help="Reference annotation in GenePred format.")
    parser.add_argument("--tm_cfg", required=True, 
                       help="Augustus configuration file for TM (TransMap) mode.")
    parser.add_argument("--miniprot_hints_gff", required=True,
                       help="GFF file containing protein alignment hints from Miniprot.")
    
    # Augustus Parameters
    parser.add_argument("--genome", required=True, 
                       help="Genome name identifier for retrieving RNA-seq hints from database.")
    parser.add_argument("--augustus_species", required=True, 
                       help="Species parameter for Augustus (e.g., 'human', 'mouse', 'fly').")
    parser.add_argument("--utr", type=int, required=True, choices=[0, 1], 
                       help="UTR prediction parameter for Augustus (0=no UTRs, 1=predict UTRs).")
    
    # Output Files
    parser.add_argument("--augustus_tm_gtf", required=True, 
                       help="Output path for Augustus TM mode predictions in GTF format.")
    
    # TMR Mode Parameters (optional)
    parser.add_argument("--augustus_tmr_gtf", 
                       help="Output path for Augustus TMR mode predictions in GTF format. "
                            "If provided, TMR mode is activated using RNA-seq evidence.")
    parser.add_argument("--augustus_hints_db", 
                       help="SQLite database file containing RNA-seq hints. Required for TMR mode.")
    parser.add_argument("--tmr_cfg", 
                       help="Augustus configuration file for TMR (TransMap + RNA-seq) mode. Required for TMR mode.")

    args = parser.parse_args()
    
    # Validate TMR mode requirements
    if args.augustus_tmr_gtf and (not args.augustus_hints_db or not args.tmr_cfg):
        parser.error("--augustus_hints_db and --tmr_cfg are required when --augustus_tmr_gtf is specified.")
    
    # Validate input files exist
    required_files = [
        args.genome_fasta, args.coding_gp, args.filtered_tm_psl, 
        args.ref_psl, args.annotation_gp, args.tm_cfg, args.miniprot_hints_gff
    ]
    
    if args.augustus_tmr_gtf:
        required_files.extend([args.augustus_hints_db, args.tmr_cfg])
    
    for file_path in required_files:
        if not os.path.exists(file_path):
            parser.error(f"Input file not found: {file_path}")
    
    # Validate output directories
    for output_path in [args.augustus_tm_gtf, args.augustus_tmr_gtf]:
        if output_path:
            output_dir = os.path.dirname(output_path)
            if output_dir and not os.path.exists(output_dir):
                try:
                    os.makedirs(output_dir, exist_ok=True)
                except OSError as e:
                    parser.error(f"Cannot create output directory {output_dir}: {e}")
    
    run_augustus_pipeline(args, args)


if __name__ == "__main__":
    main()
