#!/usr/bin/env python3
"""
Fix augMP genePred files
"""
import argparse
import sqlite3
import os
import sys
import tempfile
import shutil


def build_tx_to_gene_map(ref_db):
    """Build transcript ID -> gene ID mapping from the reference annotation DB."""
    conn = sqlite3.connect(ref_db)
    rows = conn.execute("SELECT TranscriptId, GeneId FROM annotation").fetchall()
    conn.close()
    return {tx_id: gene_id for tx_id, gene_id in rows}


def fix_genepred(gp_path, tx_to_gene):
    """Fix a single augMP genePred file in-place."""
    if not os.path.exists(gp_path) or os.path.getsize(gp_path) == 0:
        print(f"Skipping empty or missing file: {gp_path}")
        return

    fd, tmp_path = tempfile.mkstemp(dir=os.path.dirname(gp_path))
    os.close(fd)
    fixed = 0
    total = 0

    with open(gp_path) as inf, open(tmp_path, 'w') as outf:
        for line in inf:
            line = line.rstrip('\n')
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 12:
                outf.write(line + '\n')
                continue
            total += 1

            name = cols[0]
            if name.startswith('augTM-'):
                cols[0] = 'augMP-' + name[6:]

            tx_id = cols[11]
            if tx_id in tx_to_gene:
                cols[11] = tx_to_gene[tx_id]
                fixed += 1

            outf.write('\t'.join(cols) + '\n')

    shutil.move(tmp_path, gp_path)
    print(f"Processed {gp_path}: {fixed}/{total} name2 fields mapped to gene IDs")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ref-db', required=True, help='Reference annotation database')
    parser.add_argument('--augmp-files', nargs='+', required=True, help='augMP genePred file(s)')
    args = parser.parse_args()

    tx_to_gene = build_tx_to_gene_map(args.ref_db)
    print(f"Loaded {len(tx_to_gene)} transcript->gene mappings from {args.ref_db}")

    for gp_path in args.augmp_files:
        fix_genepred(gp_path, tx_to_gene)


if __name__ == '__main__':
    main()
