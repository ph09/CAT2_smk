#!/usr/bin/env python3
"""
Fix GFF3 escape sequences by properly escaping special characters in attribute values.
In GFF3 format:
- % must be escaped as %25
- ; (semicolon) in values must be escaped as %3B
- = (equals) in values must be escaped as %3D
- , (comma) in values should be escaped as %2C (though often tolerated)
- \t (tab) in values must be escaped as %09
- \n (newline) in values must be escaped as %0A

Also fixes duplicate transcript IDs by appending suffixes to make them unique.
"""

import sys
import re
from collections import defaultdict

def escape_gff3_value(value):
    """Escape special characters in GFF3 attribute values."""
    if not value:
        return value
    
    # First, escape % that is not already part of an escape sequence
    value = re.sub(r'%(?![0-9A-Fa-f]{2})', '%25', value)
    
    # Then escape other special characters that should be URL-encoded
    # But be careful not to double-escape
    # Semicolon is critical - it separates attributes
    value = value.replace(';', '%3B')
    # Equals is also critical - it separates key from value
    value = value.replace('=', '%3D')
    # Tab and newline should also be escaped
    value = value.replace('\t', '%09')
    value = value.replace('\n', '%0A')
    # Note: We don't escape commas as they're commonly used unescaped in GFF3
    # (e.g., in Dbxref attributes) and many parsers tolerate them
    
    return value

def parse_and_fix_attributes(attr_string):
    """Properly parse GFF3 attributes and fix escape sequences."""
    if not attr_string:
        return attr_string
    
    # GFF3 attributes are key=value pairs separated by semicolons
    # We need to handle cases where semicolons appear unescaped in values
    # Strategy: use regex to find all key=value patterns, where value extends
    # until the next ;key= pattern or end of string
    
    result_parts = []
    
    # Pattern: find key=value where:
    # - key is non-empty and doesn't contain = or ;
    # - value extends until next ;key= or end of string
    # We'll match: (key)=(value) where value can contain unescaped semicolons
    pattern = r'([^=;]+)=([^;]*(?:;(?![^=;]+=)[^;]*)*)'
    
    last_end = 0
    for match in re.finditer(pattern, attr_string):
        # Handle any text before this match (shouldn't happen in valid GFF3, but be safe)
        if match.start() > last_end:
            # There's text we didn't match - preserve it
            result_parts.append(attr_string[last_end:match.start()])
        
        key = match.group(1)
        value = match.group(2)
        
        # Escape special characters in the value
        escaped_value = escape_gff3_value(value)
        result_parts.append(f"{key}={escaped_value}")
        
        last_end = match.end()
    
    # Handle any remaining text
    if last_end < len(attr_string):
        remaining = attr_string[last_end:]
        # If it starts with ;, it's just a trailing separator, skip it
        if not remaining.startswith(';'):
            result_parts.append(remaining)
    
    return ';'.join(result_parts) if result_parts else attr_string

def fix_gff3_escapes(input_file, output_file):
    """Fix escape sequences in a GFF3 file and filter out duplicate transcript IDs."""
    seen_transcript_ids = set()
    duplicate_ids_to_skip = set()
    duplicate_count = 0
    
    with open(input_file, 'r') as inf, open(output_file, 'w') as outf:
        for line in inf:
            if line.startswith('#') or not line.strip():
                # Write comments and empty lines as-is
                outf.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) != 9:
                # Write malformed lines as-is
                outf.write(line)
                continue
            
            feature_type = fields[2]
            attributes = fields[8]
            
            # Extract ID and Parent from attributes
            transcript_id = None
            parent_id = None
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    if key == 'ID':
                        transcript_id = value
                    elif key == 'Parent':
                        parent_id = value
            
            # Check for duplicate transcript IDs
            if feature_type == 'transcript':
                if transcript_id in seen_transcript_ids:
                    # Skip this duplicate transcript
                    duplicate_ids_to_skip.add(transcript_id)
                    duplicate_count += 1
                    print(f"  Skipping duplicate transcript ID: {transcript_id}")
                    continue
                else:
                    seen_transcript_ids.add(transcript_id)
            
            # Skip child features (exon, CDS, etc.) of duplicate transcripts
            if feature_type in ['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'start_codon', 'stop_codon']:
                if parent_id in duplicate_ids_to_skip:
                    # Skip this child feature
                    continue
            
            # Fix the attributes column (9th column, index 8)
            fields[8] = parse_and_fix_attributes(fields[8])
            outf.write('\t'.join(fields) + '\n')
    
    if duplicate_count > 0:
        print(f"  Filtered out {duplicate_count} duplicate transcript IDs")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: fix_gff3_escape.py <input.gff3> <output.gff3>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    print(f"Fixing escape sequences in {input_file}...")
    fix_gff3_escapes(input_file, output_file)
    print(f"Output written to {output_file}")

