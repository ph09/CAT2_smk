#!/usr/bin/env python3
"""
Generate consensus gene set bigBed track files.
"""

import argparse
import os
import sys
import subprocess
import tempfile
import pandas as pd
from collections import namedtuple

GenePred = namedtuple('GenePred', [
    'name', 'chromosome', 'strand', 'start', 'stop', 'cds_start', 'cds_end',
    'block_count', 'exon_starts', 'exon_ends', 'score', 'name2',
    'cds_start_stat', 'cds_end_stat', 'exon_frames', 'thick_start', 'thick_stop'
])


def find_rgb(info):
    """Determine RGB color based on transcript biotype.
    Protein coding predicted only by augPB uses same color as unknown_likely_coding.
    """
    if info.transcript_biotype == 'unknown_likely_coding' or info.alignment_mode == 'augPB':
        return '135,76,212'
    elif info.transcript_biotype == 'protein_coding' or info.transcript_biotype == 'mRNA':
        return '76,85,212'
    elif "pseudogene" in info.transcript_biotype:
        return '255,51,255'
    return '85,212,76'


def parse_genepred_line(line):
    """Parse a genePred format line."""
    fields = line.strip().split('\t')
    if len(fields) < 10:
        return None
    
    name = fields[0]
    chromosome = fields[1]
    strand = fields[2]
    tx_start = int(fields[3])
    tx_end = int(fields[4])
    cds_start = int(fields[5])
    cds_end = int(fields[6])
    exon_count = int(fields[7])
    exon_starts = [int(x) for x in fields[8].rstrip(',').split(',')]
    exon_ends = [int(x) for x in fields[9].rstrip(',').split(',')]
    
    score = int(fields[10]) if len(fields) > 10 else 0
    name2 = fields[11] if len(fields) > 11 else name
    cds_start_stat = fields[12] if len(fields) > 12 else 'none'
    cds_end_stat = fields[13] if len(fields) > 13 else 'none'
    exon_frames = [int(x) for x in fields[14].rstrip(',').split(',')] if len(fields) > 14 else [-1] * exon_count
    
    return GenePred(
        name=name, chromosome=chromosome, strand=strand,
        start=tx_start, stop=tx_end, cds_start=cds_start, cds_end=cds_end,
        block_count=exon_count, exon_starts=exon_starts, exon_ends=exon_ends,
        score=score, name2=name2, cds_start_stat=cds_start_stat,
        cds_end_stat=cds_end_stat, exon_frames=exon_frames,
        thick_start=cds_start, thick_stop=cds_end
    )


def create_bed_info_gp(tx):
    """Create BED-style block information from genePred."""
    block_starts = ','.join(str(start - tx.start) for start in tx.exon_starts) + ','
    block_sizes = ','.join(str(end - start) for start, end in zip(tx.exon_starts, tx.exon_ends)) + ','
    exon_frames = ','.join(str(f) for f in tx.exon_frames) + ','
    return block_starts, block_sizes, exon_frames


def construct_consensus_gp_as(has_rna, has_pb):
    """Dynamically generate an autosql file for consensus."""
    consensus_gp_as = '''table bigCat
"bigCat gene models"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name"
    uint score;         "Score (0-1000)"
    char[1] strand;     "+ or - for strand"
    uint thickStart;    "Start of where display should be thick (start codon)"
    uint thickEnd;      "End of where display should be thick (stop codon)"
    uint reserved;       "RGB value (use R,G,B string in input file)"
    int blockCount;     "Number of blocks"
    int[blockCount] blockSizes; "Comma separated list of block sizes"
    int[blockCount] chromStarts; "Start positions relative to chromStart"
    string name2;       "Gene name"
    string cdsStartStat; "Status of CDS start annotation"
    string cdsEndStat;   "Status of CDS end annotation"
    int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame for exon"
    string txId; "Transcript ID"
    string type;        "Transcript type"
    string geneName;    "Gene ID"
    string geneType;    "Gene type"
    string sourceGene;    "Source gene ID"
    string sourceTranscript;    "Source transcript ID"
    string alignmentId;  "Alignment ID"
    lstring alternativeSourceTranscripts;    "Alternative source transcripts"
    string frameshift;  "Frameshifted relative to source?"
    lstring exonAnnotationSupport;   "Exon support in reference annotation"
    lstring intronAnnotationSupport;   "Intron support in reference annotation"
    string transcriptClass;    "Transcript class"
    string transcriptModes;    "Transcript mode(s)"
    string validStart;         "Valid start codon"
    string validStop;          "Valid stop codon"
    string properOrf;           "Proper multiple of 3 ORF"
'''
    if has_rna:
        consensus_gp_as += '    lstring intronRnaSupport;   "RNA intron support"\n'
        consensus_gp_as += '    lstring exonRnaSupport;  "RNA exon support"\n'
    if has_pb:
        consensus_gp_as += '    string pbIsoformSupported;   "Is this transcript supported by IsoSeq?"'
    consensus_gp_as += '\n)\n'
    return consensus_gp_as


def generate_consensus_template(genome, path):
    """Generate trackDb entry."""
    return f"""track consensus
shortLabel CAT Consensus
longLabel CAT Consensus gene set
type bigBed 12 +
bigDataUrl {path}
visibility pack
itemRgb on
searchIndex name
"""


def run_command(cmd, stderr_redirect=None):
    """Run a shell command and handle errors."""
    try:
        if stderr_redirect == '/dev/null':
            subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)
        else:
            subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(cmd)}", file=sys.stderr)
        raise


def main():
    parser = argparse.ArgumentParser(
        description='Generate consensus gene set bigBed track files',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--consensus-gp', required=True,
                        help='Input consensus genePred file')
    parser.add_argument('--consensus-gp-info', required=True,
                        help='Input consensus genePred info TSV file')
    parser.add_argument('--chrom-sizes', required=True,
                        help='Chromosome sizes file')
    parser.add_argument('--genome', required=True,
                        help='Genome name/identifier')
    parser.add_argument('--output-bigbed', required=True,
                        help='Output bigBed file path')
    parser.add_argument('--output-trackdb', required=True,
                        help='Output trackDb file path')
    parser.add_argument('--has-rnaseq', action='store_true',
                        help='Include RNA-seq support columns')
    parser.add_argument('--has-pacbio', action='store_true',
                        help='Include PacBio isoform support column')
    
    args = parser.parse_args()
    
    # Read consensus info
    print("Reading consensus gene prediction info...")
    consensus_gp_info = pd.read_csv(
        args.consensus_gp_info,
        sep='\t',
        header=0,
        na_filter=False
    ).set_index('transcript_id')
    
    # Auto-detect if we have RNA-seq or PacBio data if not explicitly specified
    has_rnaseq = args.has_rnaseq or 'intron_rna_support' in consensus_gp_info.columns
    has_pb = args.has_pacbio or 'pacbio_isoform_supported' in consensus_gp_info.columns
    
    # Create temporary files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_gp, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.as', delete=False) as as_file:
        
        tmp_gp_path = tmp_gp.name
        as_file_path = as_file.name
        
        try:
            # Generate bed file
            print("Generating BED12+ file...")
            with open(args.consensus_gp, 'r') as inf:
                for line in inf:
                    if line.startswith('#'):
                        continue
                    
                    tx = parse_genepred_line(line)
                    if tx is None:
                        continue
                    
                    if tx.name not in consensus_gp_info.index:
                        print(f"Warning: transcript {tx.name} not found in info file, skipping", 
                              file=sys.stderr)
                        continue
                    
                    info = consensus_gp_info.loc[tx.name]
                    block_starts, block_sizes, exon_frames = create_bed_info_gp(tx)
                    tx_name = info.source_transcript_name if info.source_transcript_name != 'N/A' else tx.name
                    
                    row = [
                        tx.chromosome, tx.start, tx.stop, tx_name, tx.score, tx.strand,
                        tx.thick_start, tx.thick_stop, find_rgb(info), tx.block_count,
                        block_sizes, block_starts,
                        info.source_gene_common_name, tx.cds_start_stat, tx.cds_end_stat, exon_frames,
                        tx.name, info.transcript_biotype, tx.name2, info.gene_biotype, info.source_gene,
                        info.source_transcript, info.alignment_id, info.alternative_source_transcripts,
                        info.frameshift, info.exon_annotation_support,
                        info.intron_annotation_support, info.transcript_class, info.alignment_mode,
                        info.valid_start, info.valid_stop, info.proper_orf
                    ]
                    
                    if has_rnaseq:
                        row.extend([info.intron_rna_support, info.exon_rna_support])
                    if has_pb:
                        row.append(info.pacbio_isoform_supported)
                    
                    tmp_gp.write('\t'.join(str(x) for x in row) + '\n')
            
            # Generate AutoSQL file
            print("Generating AutoSQL file...")
            as_str = construct_consensus_gp_as(has_rnaseq, has_pb)
            as_file.write(as_str)
        
        finally:
            tmp_gp.close()
            as_file.close()
        
        try:
            # Sort bed file
            print("Sorting BED file...")
            run_command(['bedSort', tmp_gp_path, tmp_gp_path])
            
            # Convert to bigBed
            print("Converting to bigBed...")
            cmd = [
                'bedToBigBed',
                '-extraIndex=name,name2,txId,geneName,sourceGene,sourceTranscript,alignmentId',
                '-type=bed12+23',
                '-tab',
                f'-as={as_file_path}',
                tmp_gp_path,
                args.chrom_sizes,
                args.output_bigbed
            ]
            run_command(cmd, stderr_redirect='/dev/null')
            
            print(f"Successfully created bigBed file: {args.output_bigbed}")
            
            # Generate trackDb
            print("Generating trackDb file...")
            with open(args.output_trackdb, 'w') as outf:
                trackdb_content = generate_consensus_template(
                    args.genome,
                    os.path.basename(args.output_bigbed)
                )
                outf.write(trackdb_content)
            
            print(f"Successfully created trackDb file: {args.output_trackdb}")
            
        finally:
            # Cleanup temporary files
            if os.path.exists(tmp_gp_path):
                os.unlink(tmp_gp_path)
            if os.path.exists(as_file_path):
                os.unlink(as_file_path)
    
    print("Done!")


if __name__ == '__main__':
    main()
