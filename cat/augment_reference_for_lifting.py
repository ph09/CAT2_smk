#!/usr/bin/env python3
"""
Augment CAT2 reference files so that pseudogenes / single-feature ncRNA
"gene" entries in the input GFF3 — features that gff3ToGenePred and gffread
skip because they have no transcript / exon children — get represented as
synthetic single-exon transcripts in the reference GenePred, transcript FASTA,
and gp_attrs.
"""

import argparse
import json
import subprocess
import sqlite3
import sys
import time
from pathlib import Path


# ---------------------------------------------------------------------------
# GFF3 sqlite (gffutils) query
# ---------------------------------------------------------------------------

def _load_existing_gp_ids(gp_path):
    ids = set()
    with open(gp_path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith('#'):
                ids.add(s.split('\t', 1)[0])
    return ids


def _iter_gene_only_features(gff3_db, existing_ids):
    """
    Yield (gene_id, chrom, start_1based, end_inclusive, strand, biotype, name) for
    every `gene` (or pseudogene-family) feature whose ID is not already a
    transcript in existing_ids and which has no transcript-type children.

    Uses raw SQL against the gffutils sqlite schema (much faster than
    instantiating gffutils.FeatureDB for millions of rows).
    """
    conn = sqlite3.connect(gff3_db)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    # gffutils tables: features(id, seqid, start, end, strand, featuretype, attributes JSON ...)
    # relations(parent, child, level)
    # find all gene-like features
    cur.execute(
        "SELECT id, seqid, start, end, strand, featuretype, attributes "
        "FROM features WHERE featuretype IN ('gene','pseudogene','ncRNA_gene')"
    )

    rows = cur.fetchall()
    print(f"  scanning {len(rows):,} gene-level features", flush=True)

    # children lookup
    has_children = set()
    cur.execute("SELECT DISTINCT parent FROM relations WHERE level = 1")
    for (parent,) in cur.fetchall():
        has_children.add(parent)

    skipped_have_tx = skipped_existing = skipped_tx_id = 0
    for r in rows:
        gid = r['id']
        if gid in existing_ids:
            skipped_existing += 1
            continue
        if gid in has_children:
            skipped_have_tx += 1
            continue
        attrs = {}
        try:
            attrs = json.loads(r['attributes']) if r['attributes'] else {}
        except (json.JSONDecodeError, TypeError):
            attrs = {}

        def _first(key, default=''):
            val = attrs.get(key, default)
            if isinstance(val, list):
                return val[0] if val else default
            return val or default

        biotype = (_first('gene_biotype') or _first('biotype') or r['featuretype'])
        name = (_first('gene_name') or _first('Name') or gid)
        transcript_id = _first('transcript_id') or name
        gene_id = _first('gene_id') or name
        if transcript_id in existing_ids or name in existing_ids:
            skipped_tx_id += 1
            continue

        yield {
            'id': transcript_id,
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': r['seqid'],
            'start': int(r['start']),
            'end': int(r['end']),
            'strand': r['strand'] or '+',
            'biotype': biotype,
            'name': name,
        }

    print(f"  {skipped_existing:,} skipped (already in GenePred), "
          f"{skipped_have_tx:,} skipped (have transcript children), "
          f"{skipped_tx_id:,} skipped (transcript_id already in GenePred)",
          flush=True)
    conn.close()


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def _gp_record(feat):
    """Build a 15-field genePredExt row for a single-exon synthetic transcript."""
    start0 = feat['start'] - 1            # gffutils stores 1-based inclusive starts
    end_excl = feat['end']
    return '\t'.join([
        feat['id'],          # name
        feat['chrom'],
        feat['strand'],
        str(start0),
        str(end_excl),
        str(end_excl),       # cdsStart = cdsEnd → non-coding
        str(end_excl),
        '1',
        f"{start0},",
        f"{end_excl},",
        '0',
        feat['name'],
        'none', 'none',
        '-1,',
    ]) + '\n'


def _faidx_seq(genome_fa, chrom, start_1, end_1):
    """Return upper-case sequence (single string) for chrom:start_1-end_1 (1-based inclusive)."""
    region = f"{chrom}:{start_1}-{end_1}"
    proc = subprocess.run(
        ["samtools", "faidx", genome_fa, region],
        capture_output=True, text=True, check=True,
    )
    lines = proc.stdout.splitlines()
    if not lines:
        return ''
    return ''.join(lines[1:]).upper()


def _revcomp(seq):
    table = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(table)[::-1]


def _wrap_fasta(seq, width=60):
    out = []
    for i in range(0, len(seq), width):
        out.append(seq[i:i + width])
    return '\n'.join(out) + '\n'


# ---------------------------------------------------------------------------
# Bulk extraction (one samtools faidx call) — much faster than per-feature
# ---------------------------------------------------------------------------

def _extract_bulk(genome_fa, regions):
    """
    Use samtools faidx -r to extract many regions in one call.
    regions: list of (key, chrom, start_1, end_1).
    Returns dict: key -> sequence (uppercase, forward-strand of target).
    """
    if not regions:
        return {}
    gfa = Path(genome_fa)
    list_file = gfa.with_suffix(gfa.suffix + '.augregions.tmp')
    with open(list_file, 'w') as fh:
        for _, chrom, s, e in regions:
            fh.write(f"{chrom}:{s}-{e}\n")
    try:
        proc = subprocess.run(
            ["samtools", "faidx", str(genome_fa), "-r", str(list_file)],
            capture_output=True, text=True, check=True,
        )
    finally:
        list_file.unlink(missing_ok=True)

    seqs = {}
    current_region = None
    current = []
    for line in proc.stdout.splitlines():
        if line.startswith('>'):
            if current_region is not None:
                seqs[current_region] = ''.join(current).upper()
            current_region = line[1:].split()[0]
            current = []
        else:
            current.append(line)
    if current_region is not None:
        seqs[current_region] = ''.join(current).upper()

    out = {}
    for key, chrom, s, e in regions:
        region_key = f"{chrom}:{s}-{e}"
        out[key] = seqs.get(region_key, '')
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gff3-db', required=True)
    ap.add_argument('--ref-gp', required=True)
    ap.add_argument('--ref-fa', required=True)
    ap.add_argument('--genome-fa', required=True)
    ap.add_argument('--gp-attrs', required=True)
    ap.add_argument('--min-len', type=int, default=15,
                    help='Skip features whose genomic span is shorter than this many bp.')
    ap.add_argument('--max-len', type=int, default=2_000_000,
                    help='Skip features whose genomic span is longer than this (safety).')
    args = ap.parse_args()

    t0 = time.time()
    print(f"[{time.time()-t0:6.1f}s] Loading existing GenePred ids: {args.ref_gp}", flush=True)
    existing = _load_existing_gp_ids(args.ref_gp)
    print(f"[{time.time()-t0:6.1f}s] {len(existing):,} ids already present", flush=True)

    print(f"[{time.time()-t0:6.1f}s] Scanning GFF3 db for gene-only features", flush=True)
    feats = []
    for f in _iter_gene_only_features(args.gff3_db, existing):
        span = f['end'] - f['start'] + 1
        if span < args.min_len or span > args.max_len:
            continue
        feats.append(f)
    print(f"[{time.time()-t0:6.1f}s] {len(feats):,} synthetic transcripts to add", flush=True)
    if not feats:
        return

    print(f"[{time.time()-t0:6.1f}s] Extracting genomic sequences via samtools faidx", flush=True)
    regions = [(f['id'], f['chrom'], f['start'], f['end']) for f in feats]
    BATCH = 50000
    seqs = {}
    for i in range(0, len(regions), BATCH):
        chunk = regions[i:i + BATCH]
        seqs.update(_extract_bulk(args.genome_fa, chunk))
        print(f"[{time.time()-t0:6.1f}s]   extracted {min(i+BATCH, len(regions)):,}/{len(regions):,}",
              flush=True)

    print(f"[{time.time()-t0:6.1f}s] Appending to reference files", flush=True)
    n_gp = n_fa = n_attrs = 0
    with open(args.ref_gp, 'a') as gp_out, \
         open(args.ref_fa, 'a') as fa_out, \
         open(args.gp_attrs, 'a') as at_out:
        for f in feats:
            seq = seqs.get(f['id'], '')
            if not seq:
                continue
            if f['strand'] == '-':
                seq = _revcomp(seq)
            gp_out.write(_gp_record(f))
            n_gp += 1
            fa_out.write(f">{f['id']}\n")
            fa_out.write(_wrap_fasta(seq))
            n_fa += 1
            at_out.write(f"{f['id']}\tgene_biotype\t{f['biotype']}\n")
            at_out.write(f"{f['id']}\tgene_name\t{f['name']}\n")
            at_out.write(f"{f['id']}\tgene_id\t{f['gene_id']}\n")
            at_out.write(f"{f['id']}\ttranscript_id\t{f['transcript_id']}\n")
            at_out.write(f"{f['id']}\ttranscript_name\t{f['name']}\n")
            at_out.write(f"{f['id']}\ttranscript_biotype\t{f['biotype']}\n")
            n_attrs += 6

    print(f"[{time.time()-t0:6.1f}s] Re-indexing transcript FASTA", flush=True)
    fai = args.ref_fa + '.fai'
    Path(fai).unlink(missing_ok=True)
    subprocess.run(["samtools", "faidx", args.ref_fa], check=True)

    print(f"[{time.time()-t0:6.1f}s] Done. Added {n_gp:,} GenePred rows, "
          f"{n_fa:,} FASTA records, {n_attrs:,} gp_attrs lines.", flush=True)


if __name__ == '__main__':
    main()
