#!/usr/bin/env python3
"""
Prepare reference GFF3 files for CAT / gff3ToGenePred / transMap.

Fixes applied (in order):
  1. Remove tRNA gene loci (transMapPslToGenePred cannot lift TRN* genePreds).
  2. Remove orphan gene shells with no transcript/exon/CDS children.  gff3ToGenePred
     still emits genePreds from these, but they break transMapPslToGenePred.
  3. Normalize transcript feature types: mRNA / lncRNA -> transcript.
  4. Insert a transcript layer for genes whose exons are parented directly to gene
     (common for NCBI Gnomon pseudogenes on sex chromosomes).
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict

TRANSCRIPT_TYPES = frozenset({
    'mRNA', 'lncRNA', 'transcript', 'ncRNA', 'primary_transcript',
})
TRANSCRIPT_TYPES_TO_NORMALIZE = frozenset({'mRNA', 'lncRNA'})
CHILD_TYPES_TO_REPARENT = frozenset({
    'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'UTR',
    'start_codon', 'stop_codon',
})


def parse_attrs(attr_str: str) -> dict[str, str]:
    out = {}
    for part in attr_str.strip().split(';'):
        if not part or '=' not in part:
            continue
        k, v = part.split('=', 1)
        out[k] = v
    return out


def format_attrs(attrs: dict[str, str], extra: list[tuple[str, str]] | None = None) -> str:
    keys_done = set()
    parts = []
    for k, v in (extra or []):
        parts.append(f'{k}={v}')
        keys_done.add(k)
    for k, v in attrs.items():
        if k in keys_done:
            continue
        parts.append(f'{k}={v}')
    return ';'.join(parts)


def rewrite_parent(attr_str: str, gene_parent: str, rna_parent: str) -> str:
    attrs = parse_attrs(attr_str)
    parents = [p.strip() for p in attrs.get('Parent', '').split(',') if p.strip()]
    new_parents = [rna_parent if p == gene_parent else p for p in parents]
    if not new_parents:
        new_parents = [rna_parent]
    attrs['Parent'] = ','.join(new_parents)
    return format_attrs(attrs)


def pick_transcript_id(gene_attrs: dict[str, str], gene_id_key: str) -> str:
    for key in ('transcript_id', 'Name', 'gene', 'gene_id', 'ID'):
        val = gene_attrs.get(key, '')
        if not val:
            continue
        if key == 'ID' and val.startswith('gene-'):
            val = val[5:]
        return val
    return gene_id_key[5:] if gene_id_key.startswith('gene-') else gene_id_key


def build_transcript_line(gene_parts: list[str], gene_attrs: dict[str, str], gene_id_key: str,
                          tx_start: int, tx_end: int) -> str:
    chrom, source, _gene, _start, _end, score, strand, phase = gene_parts[:8]
    biotype = gene_attrs.get('gene_biotype', gene_attrs.get('gene_type', 'pseudogene'))
    tx_id = pick_transcript_id(gene_attrs, gene_id_key)
    rna_id = f'rna-{tx_id}'
    tx_biotype = gene_attrs.get('transcript_biotype', biotype)
    extra = [
        ('ID', rna_id),
        ('Parent', gene_id_key),
        ('transcript_id', tx_id),
        ('gene_id', gene_attrs.get('gene_id', tx_id)),
        ('gene_biotype', biotype),
        ('transcript_biotype', tx_biotype),
    ]
    if 'Name' in gene_attrs:
        extra.append(('Name', gene_attrs['Name']))
    if 'gene' in gene_attrs:
        extra.append(('gene', gene_attrs['gene']))
    if 'gbkey' in gene_attrs:
        extra.append(('gbkey', 'mRNA'))
    attr_str = format_attrs(gene_attrs, extra=extra)
    return '\t'.join([
        chrom, source, 'transcript', str(tx_start), str(tx_end), score, strand, phase, attr_str,
    ]) + '\n'


def remove_trna_loci(lines: list[str]) -> tuple[list[str], int]:
    """Drop tRNA genes and all child features."""
    trna_genes = set()
    trna_transcripts = set()

    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            continue
        attrs = parse_attrs(parts[8])
        feat = parts[2]
        if feat == 'gene' and attrs.get('gene_biotype') == 'tRNA':
            gid = attrs.get('ID', '')
            if gid:
                trna_genes.add(gid)
        elif feat == 'tRNA':
            for par in attrs.get('Parent', '').split(','):
                par = par.strip()
                if par in trna_genes:
                    tid = attrs.get('ID', '')
                    if tid:
                        trna_transcripts.add(tid)

    drop_parents = trna_genes | trna_transcripts
    if not drop_parents:
        return lines, 0

    out = []
    removed = 0
    for line in lines:
        if line.startswith('#') or not line.strip():
            out.append(line)
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            out.append(line)
            continue
        attrs = parse_attrs(parts[8])
        fid = attrs.get('ID', '')
        if fid and fid in drop_parents:
            removed += 1
            continue
        if parts[2] == 'gene' and attrs.get('gene_biotype') == 'tRNA':
            removed += 1
            continue
        if parts[2] == 'tRNA':
            removed += 1
            continue
        parents = [p.strip() for p in attrs.get('Parent', '').split(',') if p.strip()]
        if any(p in drop_parents for p in parents):
            removed += 1
            continue
        out.append(line)
    return out, removed


def _gene_ids(lines: list[str]) -> set[str]:
    gene_ids = set()
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        gid = parse_attrs(parts[8]).get('ID', '')
        if gid:
            gene_ids.add(gid)
    return gene_ids


def _gene_has_feature_children(lines: list[str]) -> set[str]:
    """Return gene IDs that have at least one non-gene child feature."""
    gene_ids = _gene_ids(lines)
    has_child = set()
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9 or parts[2] == 'gene':
            continue
        parents = [p.strip() for p in parse_attrs(parts[8]).get('Parent', '').split(',') if p.strip()]
        for par in parents:
            if par in gene_ids:
                has_child.add(par)
    return has_child


def remove_orphan_genes(lines: list[str]) -> tuple[list[str], int]:
    """Drop gene rows with no child features, plus any dangling descendants.

    NCBI annotation sometimes lists loci as gene-only shells with no exons or
    transcripts.  gff3ToGenePred still emits a genePred, but transMapPslToGenePred
    can abort on the resulting PSL rows.
    """
    has_child = _gene_has_feature_children(lines)
    drop_genes = set()
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        gid = parse_attrs(parts[8]).get('ID', '')
        if gid and gid not in has_child:
            drop_genes.add(gid)

    if not drop_genes:
        return lines, 0

    # Collect transcript/ncRNA IDs parented to dropped genes so exons/CDS are removed too.
    drop_parents = set(drop_genes)
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            continue
        attrs = parse_attrs(parts[8])
        parents = [p.strip() for p in attrs.get('Parent', '').split(',') if p.strip()]
        if any(p in drop_genes for p in parents):
            fid = attrs.get('ID', '')
            if fid:
                drop_parents.add(fid)

    out = []
    removed = 0
    for line in lines:
        if line.startswith('#') or not line.strip():
            out.append(line)
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            out.append(line)
            continue
        attrs = parse_attrs(parts[8])
        fid = attrs.get('ID', '')
        if fid and fid in drop_parents:
            removed += 1
            continue
        parents = [p.strip() for p in attrs.get('Parent', '').split(',') if p.strip()]
        if any(p in drop_parents for p in parents):
            removed += 1
            continue
        out.append(line)
    return out, removed


def normalize_transcript_features(lines: list[str]) -> tuple[list[str], int]:
    out = []
    changed = 0
    for line in lines:
        if line.startswith('#') or not line.strip():
            out.append(line)
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) >= 9 and parts[2] in TRANSCRIPT_TYPES_TO_NORMALIZE:
            parts[2] = 'transcript'
            changed += 1
            out.append('\t'.join(parts) + '\n')
        else:
            out.append(line)
    return out, changed


def analyze_gff3(lines: list[str]):
    gene_rows = {}
    has_transcript_child = set()
    direct_children = defaultdict(list)

    for idx, line in enumerate(lines):
        if line.startswith('#') or not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            continue
        feat = parts[2]
        attrs = parse_attrs(parts[8])
        if feat == 'gene':
            gid = attrs.get('ID', '')
            if gid:
                gene_rows[gid] = (idx, parts, attrs)
            continue
        parents = [p.strip() for p in attrs.get('Parent', '').split(',') if p.strip()]
        for par in parents:
            if not par.startswith('gene-') or par not in gene_rows:
                continue
            if feat in TRANSCRIPT_TYPES:
                has_transcript_child.add(par)
            elif feat in CHILD_TYPES_TO_REPARENT:
                start, end = int(parts[3]), int(parts[4])
                direct_children[par].append((idx, start, end))

    needs_fix = {}
    for gid, (_gidx, gparts, gattrs) in gene_rows.items():
        if gid in has_transcript_child:
            continue
        kids = direct_children.get(gid)
        if not kids:
            continue
        starts = [x[1] for x in kids]
        ends = [x[2] for x in kids]
        needs_fix[gid] = {
            'gene_parts': gparts,
            'gene_attrs': gattrs,
            'tx_start': min(starts),
            'tx_end': max(ends),
            'child_indices': [x[0] for x in kids],
        }
    return needs_fix


def fix_direct_exon_parents(lines: list[str], needs_fix: dict) -> tuple[list[str], int]:
    if not needs_fix:
        return lines, 0

    child_reparent = {}
    for gid, info in needs_fix.items():
        rna_id = f'rna-{pick_transcript_id(info["gene_attrs"], gid)}'
        for cidx in info['child_indices']:
            child_reparent[cidx] = (gid, rna_id)

    out = []
    n_tx = 0
    for idx, line in enumerate(lines):
        if line.startswith('#') or not line.strip():
            out.append(line)
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            out.append(line)
            continue
        attrs = parse_attrs(parts[8])
        gid = attrs.get('ID', '')
        if parts[2] == 'gene' and gid in needs_fix:
            out.append(line)
            info = needs_fix[gid]
            out.append(build_transcript_line(
                info['gene_parts'], info['gene_attrs'], gid,
                info['tx_start'], info['tx_end'],
            ))
            n_tx += 1
            continue
        if idx in child_reparent:
            gene_par, rna_par = child_reparent[idx]
            parts[8] = rewrite_parent(parts[8], gene_par, rna_par)
            out.append('\t'.join(parts) + '\n')
            continue
        out.append(line)
    return out, n_tx


def prepare_cat_reference_gff3(lines: list[str], also_pc: bool = False) -> tuple[list[str], dict]:
    stats = {}
    lines, stats['trna_features_removed'] = remove_trna_loci(lines)
    lines, removed = remove_orphan_genes(lines)
    stats['orphan_genes_removed'] = removed
    lines, stats['transcript_types_normalized'] = normalize_transcript_features(lines)
    needs_fix = analyze_gff3(lines)
    if not also_pc:
        needs_fix = {
            gid: info for gid, info in needs_fix.items()
            if info['gene_attrs'].get('gene_biotype') != 'protein_coding'
        }
    lines, stats['transcript_layers_inserted'] = fix_direct_exon_parents(lines, needs_fix)
    return lines, stats


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input_gff3', help='Input GFF3 path')
    ap.add_argument('output_gff3', help='Output GFF3 path')
    ap.add_argument(
        '--also-pc',
        action='store_true',
        help='Also fix protein_coding genes with direct exons (default: non-PC only).',
    )
    ap.add_argument(
        '--skip-trna-removal',
        action='store_true',
        help='Do not remove tRNA loci (not recommended for CAT transMap).',
    )
    ap.add_argument(
        '--skip-orphan-removal',
        action='store_true',
        help='Do not remove gene-only loci with no transcript/exon/CDS children.',
    )
    ap.add_argument(
        '--skip-transcript-normalize',
        action='store_true',
        help='Do not rewrite mRNA/lncRNA feature types to transcript.',
    )
    args = ap.parse_args()

    with open(args.input_gff3) as fh:
        lines = fh.readlines()

    stats = {}
    working = lines
    if not args.skip_trna_removal:
        working, stats['trna_features_removed'] = remove_trna_loci(working)
    else:
        stats['trna_features_removed'] = 0
    if not args.skip_orphan_removal:
        working, stats['orphan_genes_removed'] = remove_orphan_genes(working)
    else:
        stats['orphan_genes_removed'] = 0
    if not args.skip_transcript_normalize:
        working, stats['transcript_types_normalized'] = normalize_transcript_features(working)
    else:
        stats['transcript_types_normalized'] = 0
    needs_fix = analyze_gff3(working)
    if not args.also_pc:
        needs_fix = {
            gid: info for gid, info in needs_fix.items()
            if info['gene_attrs'].get('gene_biotype') != 'protein_coding'
        }
    out_lines, stats['transcript_layers_inserted'] = fix_direct_exon_parents(working, needs_fix)

    with open(args.output_gff3, 'w') as fh:
        fh.writelines(out_lines)

    print(
        f'Wrote {args.output_gff3}: '
        + ', '.join(f'{k}={v}' for k, v in stats.items()),
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
