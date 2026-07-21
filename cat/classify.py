import argparse
import bisect
import collections
import sqlite3
import pickle
import pandas as pd

import tools.bio
import tools.fileOps
import tools.intervals
import tools.mathOps
import tools.nameConversions
import tools.psl
import tools.sqlInterface
import tools.transcripts

# Distance allowed between intron locations to be considered equivalent
FUZZ_DISTANCE = 7


def run_classification(args):
    ref_tx_dict = tools.transcripts.get_gene_pred_dict(args.annotation_gp)
    tx_biotype_map = tools.sqlInterface.get_transcript_biotype_map(args.ref_db_path)
    seq_dict = tools.bio.get_sequence_dict(args.fasta)
    transcript_modes = collections.defaultdict(dict)
    for mode, gp_path, mrna_path, cds_path in args.mode_files:
        transcript_modes[mode]['gp'] = gp_path
        transcript_modes[mode]['mRNA'] = mrna_path
        transcript_modes[mode]['CDS'] = cds_path
    results = []
    for tx_mode, path_dict in transcript_modes.items():
        tx_dict = tools.transcripts.get_gene_pred_dict(path_dict['gp'])
        for aln_mode in ['CDS', 'mRNA']:
            psl_path = path_dict.get(aln_mode)
            psl_iter = list(tools.psl.psl_iterator(psl_path))
            mc_df = metrics_classify(aln_mode, ref_tx_dict, tx_dict, tx_biotype_map, psl_iter, seq_dict)
            ec_df = evaluation_classify(aln_mode, ref_tx_dict, tx_dict, tx_biotype_map, psl_iter, seq_dict)
            metrics_tbl_name = tools.sqlInterface.tables[aln_mode][tx_mode]['metrics'].__tablename__
            eval_tbl_name = tools.sqlInterface.tables[aln_mode][tx_mode]['evaluation'].__tablename__
            results.append((metrics_tbl_name, mc_df))
            results.append((eval_tbl_name, ec_df))
    
    with open(args.resolved_df, 'wb') as f:
        pickle.dump(results, f)
            
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation-gp", required=True, help="Reference annotation genePred file.")
    parser.add_argument("--ref-db-path", required=True, help="Path to the reference genome database.")
    parser.add_argument("--fasta", required=True, help="Target genome FASTA file.")
    parser.add_argument("--db-path", required=True, help="Path to the SQLite database to write results to.")
    parser.add_argument("--resolved-df", required=True, help="Path to save the resolved DataFrame.")
    parser.add_argument("--mode-files", nargs=4, action='append', required=True,
                        metavar=("MODE", "INPUT_GP", "MRNA_PSL", "CDS_PSL"),
                        help="Define a transcript mode and its files. This argument can be specified multiple times.")

    args = parser.parse_args()
    results = run_classification(args)
    return results
    


def metrics_classify(aln_mode, ref_tx_dict, tx_dict, tx_biotype_map, psl_iter, seq_dict):
    r = []
    for ref_tx, tx, psl, biotype in tx_iter(psl_iter, ref_tx_dict, tx_dict, tx_biotype_map):
        original_intron_vector = calculate_original_intron_vector(ref_tx, tx, psl, aln_mode)
        adj_start, adj_stop = find_adj_start_stop(tx, seq_dict)
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'AlnCoverage', 100 * psl.target_coverage])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'AlnIdentity', 100 * psl.identity])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'AlnGoodness', 100 * (1 - psl.badness)])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'PercentUnknownBases', psl.percent_n])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'OriginalIntrons', original_intron_vector])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'ValidStart', tools.transcripts.has_start_codon(seq_dict, tx)])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'ValidStop', tools.transcripts.has_stop_codon(seq_dict, tx)])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'ProperOrf', tx.cds_size % 3 == 0])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'AdjStart', adj_start])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'AdjStop', adj_stop])
    columns = ['GeneId', 'TranscriptId', 'AlignmentId', 'classifier', 'value']
    df = pd.DataFrame(r, columns=columns)
    df = df.sort_values(columns)
    df = df.set_index('AlignmentId')
    assert len(r) == len(df)
    return df


def evaluation_classify(aln_mode, ref_tx_dict, tx_dict, tx_biotype_map, psl_iter, seq_dict):
    """
    Calculates the evaluation metrics on this transcript_chunk
    :return: DataFrame
    """
    r = []
    for ref_tx, tx, psl, biotype in tx_iter(psl_iter, ref_tx_dict, tx_dict, tx_biotype_map):
        r.extend(find_indels(tx, psl, aln_mode))
        if biotype == 'protein_coding':
            line = in_frame_stop(tx, seq_dict)
            if line is not None:
                r.append(line)
    columns = ['AlignmentId', 'chromosome', 'start', 'stop', 'name', 'score', 'strand', 'thickStart',
               'thickStop', 'rgb', 'blockCount', 'blockSizes', 'blockStarts']
    df = pd.DataFrame(r, columns=columns)
    df = df.sort_values(columns)
    df = df.set_index('AlignmentId')
    assert len(r) == len(df)
    return df

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
        """Converts either an insertion or a deletion into a output transcript"""
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
        if q_offset != 0:  # query insertion -> insertion in target sequence
            left_pos = q_start - q_offset
            right_pos = q_start
            row = parse_indel(left_pos, right_pos, coordinate_fn, tx, q_offset, 'Insertion')
            if row is not None:
                r.append(row)
        if t_offset != 0:  # target insertion -> insertion in reference sequence
            left_pos = right_pos = q_start
            row = parse_indel(left_pos, right_pos, coordinate_fn, tx, t_offset, 'Deletion')
            if row is not None:
                r.append(row)
        q_pos = q_start
        t_pos = t_start
    return r

def tx_iter(psl_iter, ref_tx_dict, tx_dict, tx_biotype_map):
    """
    yields tuples of (GenePredTranscript <reference> , GenePredTranscript <target>, PslRow, biotype
    """
    for psl in psl_iter:
        # this psl is target-referenced
        ref_tx = ref_tx_dict[psl.t_name]
        tx = tx_dict[psl.q_name]
        biotype = tx_biotype_map[psl.t_name]
        yield ref_tx, tx, psl, biotype


def convert_cds_frames(ref_tx, tx, aln_mode):
    if aln_mode == 'CDS':
        if ref_tx.offset != 0:
            ref_tx = convert_cds_frame(ref_tx)
        if tx.offset != 0:
            tx = convert_cds_frame(tx)
    return ref_tx, tx


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
    # remove None which means this transcript is protein_coding and that exon is entirely non-coding
    return [x for x in introns if x is not None]


def get_exon_intervals(tx, aln_mode):
    if aln_mode == 'CDS':
        tx = convert_cds_frame(tx)
    exons = {}
    for exon in tx.exon_intervals:
        start = tx.chromosome_coordinate_to_mrna(exon.start)
        stop = tx.chromosome_coordinate_to_mrna(exon.stop - 1)  # zero based, half open
        if tx.strand == '-':
            start, stop = stop, start
        i = tools.intervals.ChromosomeInterval(None, start, stop + 1, '.')
        exons[exon] = i
    return exons

if __name__ == "__main__":
    main()
