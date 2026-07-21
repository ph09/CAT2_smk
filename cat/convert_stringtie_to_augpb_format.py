#!/usr/bin/env python3
"""
Convert StringTie GTF format to a clean GTF with strg- prefixed names,
compatible with the augPB processing pipeline for consensus building.
"""

import sys
import re

def convert_name(name, name_type='gene'):
    """
    Convert StringTie naming to strg naming.
    
    Examples:
        MSTRG.1 -> strg-1 (gene)
        MSTRG.1.1 -> strg-1.t1 (transcript)
        MSTRG.2.2 -> strg-2.t2 (transcript)
    """
    if name_type == 'gene':
        match = re.match(r'M?STRG\.(\d+)', name)
        if match:
            return f"strg-{match.group(1)}"
    elif name_type == 'transcript':
        match = re.match(r'M?STRG\.(\d+)\.(\d+)', name)
        if match:
            return f"strg-{match.group(1)}.t{match.group(2)}"
    return name

def parse_attributes(attr_string):
    """Parse GTF attributes into a dictionary."""
    attrs = {}
    pattern = r'(\w+)\s+"([^"]+)"'
    for match in re.finditer(pattern, attr_string):
        key, value = match.groups()
        attrs[key] = value
    return attrs

def format_attributes(attrs, feature_type):
    """
    Format attributes to match augPB style.
    
    For transcript features:
        gene_id "X"; transcript_id "Y";  gene_name "Z";
    
    For exon features:
        gene_id "X"; transcript_id "Y"; exon_number "N"; exon_id "Y.N"; gene_name "Z";
    """
    parts = []
    
    if 'gene_id' in attrs:
        parts.append(f'gene_id "{attrs["gene_id"]}"')
    
    if 'transcript_id' in attrs:
        parts.append(f'transcript_id "{attrs["transcript_id"]}"')
    
    if feature_type == 'transcript':
        if 'gene_name' in attrs:
            parts.append(f' gene_name "{attrs["gene_name"]}"')
        return '; '.join(parts) + ';'
    else:
        if 'exon_number' in attrs:
            parts.append(f'exon_number "{attrs["exon_number"]}"')
        
        if 'exon_id' in attrs:
            parts.append(f'exon_id "{attrs["exon_id"]}"')
        elif 'transcript_id' in attrs and 'exon_number' in attrs:
            exon_id = f"{attrs['transcript_id']}.{attrs['exon_number']}"
            parts.append(f'exon_id "{exon_id}"')
        
        if 'gene_name' in attrs:
            parts.append(f'gene_name "{attrs["gene_name"]}"')
        
        return '; '.join(parts) + ';'

def convert_gtf_line(line):
    """Convert a single GTF line from StringTie format to strg format."""
    if line.startswith('#'):
        return None
    
    fields = line.rstrip('\n').split('\t')
    if len(fields) != 9:
        return None
    
    chrom, source, feature, start, end, score, strand, frame, attributes = fields
    
    if strand == '.':
        strand = '+'
    
    source = 'StringTie'
    
    attrs = parse_attributes(attributes)
    
    if 'gene_id' in attrs:
        attrs['gene_id'] = convert_name(attrs['gene_id'], 'gene')
    
    if 'transcript_id' in attrs:
        attrs['transcript_id'] = convert_name(attrs['transcript_id'], 'transcript')
    
    if 'gene_name' not in attrs and 'gene_id' in attrs:
        attrs['gene_name'] = attrs['gene_id']
    else:
        attrs['gene_name'] = convert_name(attrs['gene_name'], 'gene')
    
    if feature == 'exon' and 'exon_id' not in attrs:
        if 'transcript_id' in attrs and 'exon_number' in attrs:
            attrs['exon_id'] = f"{attrs['transcript_id']}.{attrs['exon_number']}"
    
    new_attributes = format_attributes(attrs, feature)
    
    # Reconstruct line
    new_line = '\t'.join([chrom, source, feature, start, end, score, strand, frame, new_attributes])
    return new_line

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_stringtie.gtf> <output_strg_format.gtf>", file=sys.stderr)
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            converted_line = convert_gtf_line(line)
            if converted_line:
                outfile.write(converted_line + '\n')
    
    print(f"Conversion complete: {input_file} -> {output_file}", file=sys.stderr)

if __name__ == '__main__':
    main()
