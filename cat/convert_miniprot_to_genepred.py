#!/usr/bin/env python3
"""
Convert miniprot PAF output to GenePred (and real PSL) format.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from typing import Iterable, List, Tuple

CIGAR_OP_RE = re.compile(r'(\d+)([MIDNUVFG])')


# ---------------------------------------------------------------------------
# PAF / CIGAR parsing
# ---------------------------------------------------------------------------

def _parse_tags(tag_fields: Iterable[str]) -> dict:
    """Parse PAF tag fields of the form ``KEY:TYPE:VALUE`` into ``{KEY: VALUE}``.

    Integer-typed tags (``i``) are coerced to int; everything else stays as str.
    """
    out = {}
    for tag in tag_fields:
        if not tag or ':' not in tag:
            continue
        parts = tag.split(':', 2)
        if len(parts) != 3:
            continue
        key, typ, val = parts
        if typ == 'i':
            try:
                out[key] = int(val)
            except ValueError:
                out[key] = 0
        else:
            out[key] = val
    return out


def parse_cigar_to_exons(
    cg: str, t_start: int,
) -> Tuple[List[Tuple[int, int]], int, int, int]:
    """Walk a miniprot ``cg:Z:`` CIGAR and return:

    * ``exons``      -- list of ``(genomic_start, genomic_end)`` half-open intervals.
    * ``matched_aa`` -- sum of ``M`` operations (true aa-aa matches in BLOSUM space).
    * ``aligned_aa`` -- ``M + I + U + V`` (aa positions consumed by alignment, the
                       denominator for "fraction-of-query aligned" identity).
    * ``frameshifts``-- count of ``F`` / ``G`` ops in the alignment (each one is a
                       1-bp insertion or deletion in the target).

    Exon boundaries are opened/closed at ``N`` / ``U`` / ``V`` intron operations.
    For simplicity, the 1–2 bp of a split codon spanning a U/V intron is folded
    into the intron span; Augustus refines exon boundaries downstream anyway.
    """
    exons: List[Tuple[int, int]] = []
    cur_start = t_start
    t_pos = t_start
    matched_aa = 0
    aligned_aa = 0
    frameshifts = 0

    for n_str, op in CIGAR_OP_RE.findall(cg):
        n = int(n_str)
        if op == 'M':
            t_pos += 3 * n
            matched_aa += n
            aligned_aa += n
        elif op == 'I':
            aligned_aa += n
        elif op == 'D':
            t_pos += n
        elif op in ('N', 'U', 'V'):
            if t_pos > cur_start:
                exons.append((cur_start, t_pos))
            t_pos += n
            cur_start = t_pos
            if op in ('U', 'V'):
                aligned_aa += 1
        elif op == 'F':
            t_pos += 1
            frameshifts += 1
        elif op == 'G':
            t_pos -= 1
            frameshifts += 1

    if t_pos > cur_start:
        exons.append((cur_start, t_pos))
    return exons, matched_aa, aligned_aa, frameshifts


def parse_paf_row(line: str) -> dict | None:
    """Parse a single miniprot PAF row into a record dict.

    Returns ``None`` for rows that lack required fields or whose strand /
    target name indicate an unmapped read.
    """
    f = line.rstrip('\n').split('\t')
    if len(f) < 12:
        return None
    try:
        rec = {
            'q_name': f[0], 'q_len': int(f[1]),
            'q_start': int(f[2]), 'q_end': int(f[3]),
            'strand': f[4],
            't_name': f[5], 't_len': int(f[6]),
            't_start': int(f[7]), 't_end': int(f[8]),
            'paf_match': int(f[9]),
            'paf_aln_len': int(f[10]),
            'mapq': int(f[11]),
        }
    except (ValueError, IndexError):
        return None
    if rec['strand'] not in ('+', '-') or rec['t_name'] == '*':
        return None
    rec['tags'] = _parse_tags(f[12:])
    return rec


# ---------------------------------------------------------------------------
# Record construction
# ---------------------------------------------------------------------------

def _exon_frames(exons: List[Tuple[int, int]], strand: str) -> List[int]:
    """genePredExt-style per-exon reading frame.

    Frames count cumulative coding bases up to the start of each exon, in the
    direction of transcription. For a fully-coding miniprot record every exon is
    coding, so cumulative bases roll over modulo 3.
    """
    n = len(exons)
    frames = [-1] * n
    cds_bases = 0
    order = range(n) if strand == '+' else range(n - 1, -1, -1)
    for i in order:
        es, ee = exons[i]
        frames[i] = cds_bases % 3
        cds_bases += ee - es
    return frames


def build_gp_line(
    name: str, rec: dict, exons: List[Tuple[int, int]],
) -> str:
    chrom = rec['t_name']
    strand = rec['strand']
    tx_start = exons[0][0]
    tx_end = exons[-1][1]
    cds_start = tx_start
    cds_end = tx_end
    exon_count = len(exons)
    exon_starts = ','.join(str(e[0]) for e in exons) + ','
    exon_ends = ','.join(str(e[1]) for e in exons) + ','
    score = int(rec['tags'].get('AS', 0))

    # cdsStartStat: complete if miniprot's start-codon tag (st:i:) is 1.
    cds_start_stat = 'cmpl' if int(rec['tags'].get('st', 0)) == 1 else 'incmpl'
    # cdsEndStat: complete iff alignment reaches the C-terminus of the protein.
    cds_end_stat = 'cmpl' if rec['q_end'] == rec['q_len'] else 'incmpl'

    frames = _exon_frames(exons, strand)

    return '\t'.join([
        name, chrom, strand,
        str(tx_start), str(tx_end),
        str(cds_start), str(cds_end),
        str(exon_count),
        exon_starts, exon_ends,
        str(score),
        rec['q_name'],
        cds_start_stat, cds_end_stat,
        ','.join(str(fr) for fr in frames) + ',',
    ])


def build_psl_line(
    name: str, rec: dict,
    exons: List[Tuple[int, int]],
    matched_aa: int, aligned_aa: int,
    frameshifts: int,
) -> str:
    """Build a PSL row in protein-residue query space.

    The PSL is consumed by ``store_psl_metrics.py``, which reads only:
        col0  matches
        col1  misMatches
        col2  repMatches
        col9  qName
        col10 qSize
        col18 blockSizes  (comma-terminated)

    and derives:
        AlnIdentity = matches / (matches + misMatches + repMatches)
        AlnCoverage = sum(blockSizes) / qSize

    We keep the rest of the standard PSL columns sensible so the file also
    validates against PSL-aware tooling. Block sizes are reported in aa to
    keep coverage in protein-residue space; qSize = protein length.
    """
    matches = matched_aa
    mismatches = max(0, aligned_aa - matched_aa)
    rep_matches = 0
    ns = 0
    q_inserts = 0
    q_insert_bases = 0
    t_inserts = max(0, len(exons) - 1)
    t_insert_bases = (exons[-1][1] - exons[0][0]) - sum(e[1] - e[0] for e in exons)

    block_sizes_aa = [max(1, (e[1] - e[0]) // 3) for e in exons]
    block_sizes = ','.join(str(b) for b in block_sizes_aa) + ','

    # qStarts in aa, distributed along blocks. Begin at q_start.
    q_starts_list = []
    cum = rec['q_start']
    for b in block_sizes_aa:
        q_starts_list.append(cum)
        cum += b
    q_starts = ','.join(str(x) for x in q_starts_list) + ','

    t_starts = ','.join(str(e[0]) for e in exons) + ','

    return '\t'.join([
        str(matches), str(mismatches), str(rep_matches), str(ns),
        str(q_inserts), str(q_insert_bases),
        str(t_inserts), str(t_insert_bases),
        rec['strand'],
        name, str(rec['q_len']), str(rec['q_start']), str(rec['q_end']),
        rec['t_name'], str(rec['t_len']),
        str(exons[0][0]), str(exons[-1][1]),
        str(len(exons)),
        block_sizes, q_starts, t_starts,
    ])


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------

def convert(
    paf_path: str, gp_path: str, psl_path: str | None,
    min_coverage: float, min_identity: float, min_mapq: int, min_score: int,
):
    n_rows = 0
    n_dropped_struct = 0
    n_dropped_cov = 0
    n_dropped_id = 0
    n_dropped_mapq = 0
    n_dropped_score = 0
    n_dropped_no_exons = 0
    n_written = 0

    copy_counter: dict[str, int] = defaultdict(int)

    fpsl = open(psl_path, 'w') if psl_path else None
    try:
        with open(paf_path) as fin, open(gp_path, 'w') as fgp:
            for line in fin:
                if not line.strip() or line.startswith('#'):
                    continue
                n_rows += 1
                rec = parse_paf_row(line)
                if rec is None:
                    n_dropped_struct += 1
                    continue

                if rec['mapq'] < min_mapq:
                    n_dropped_mapq += 1
                    continue
                if int(rec['tags'].get('AS', 0)) < min_score:
                    n_dropped_score += 1
                    continue

                cg = rec['tags'].get('cg', '')
                if not cg:
                    n_dropped_struct += 1
                    continue

                exons, matched_aa, aligned_aa, frameshifts = parse_cigar_to_exons(
                    cg, rec['t_start']
                )
                if not exons:
                    n_dropped_no_exons += 1
                    continue

                # Protein-space coverage / identity (Liftoff-style, kept for parity with the external Liftoff tool's metric definitions).
                coverage = aligned_aa / rec['q_len'] if rec['q_len'] else 0.0
                identity = matched_aa / aligned_aa if aligned_aa else 0.0
                if coverage < min_coverage:
                    n_dropped_cov += 1
                    continue
                if identity < min_identity:
                    n_dropped_id += 1
                    continue

                copy_counter[rec['q_name']] += 1
                idx = copy_counter[rec['q_name']]
                name = rec['q_name'] if idx == 1 else f"{rec['q_name']}_{idx}"

                fgp.write(build_gp_line(name, rec, exons) + '\n')
                if fpsl is not None:
                    fpsl.write(build_psl_line(
                        name, rec, exons, matched_aa, aligned_aa, frameshifts,
                    ) + '\n')
                n_written += 1
    finally:
        if fpsl is not None:
            fpsl.close()

    print(
        f"convert_miniprot_to_genepred: read {n_rows:,} PAF rows; "
        f"wrote {n_written:,} records "
        f"(dropped struct={n_dropped_struct:,} no-exons={n_dropped_no_exons:,} "
        f"cov={n_dropped_cov:,} id={n_dropped_id:,} "
        f"mapq={n_dropped_mapq:,} score={n_dropped_score:,})",
        file=sys.stderr,
    )
    return n_written


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Convert miniprot PAF output to GenePred and (optionally) a real "
            "PSL. Parses the cg:Z: CIGAR for proper exon structure; emits one "
            "record per PAF row with _2, _3, ... copy suffixes for paralogs; "
            "filters by protein-space coverage / identity / mapq / score."
        ),
    )
    ap.add_argument('paf', help='miniprot PAF input')
    ap.add_argument('gp',  help='Output genePredExt')
    ap.add_argument('--psl', default=None,
                    help='Optional real PSL output for downstream metric '
                         'derivation (store_psl_metrics.py).')
    ap.add_argument('--min-coverage', type=float, default=0.0,
                    help='Drop records whose (aligned_aa / q_len) < this. '
                         'Default 0.0 = keep everything (real metrics are '
                         'preserved in the PSL output and downstream consensus '
                         'filtering can act on them honestly). Raise to e.g. '
                         '0.50 to filter at converter time.')
    ap.add_argument('--min-identity', type=float, default=0.0,
                    help='Drop records whose (matched_aa / aligned_aa) < this. '
                         'Default 0.0 = keep everything (same rationale as '
                         '--min-coverage).')
    ap.add_argument('--min-mapq', type=int, default=0,
                    help='Drop records whose mapping quality < this '
                         '(default 0, i.e. no filter).')
    ap.add_argument('--min-score', type=int, default=0,
                    help='Drop records whose AS:i: alignment score < this '
                         '(default 0, i.e. no filter).')
    args = ap.parse_args()

    convert(
        args.paf, args.gp, args.psl,
        min_coverage=args.min_coverage,
        min_identity=args.min_identity,
        min_mapq=args.min_mapq,
        min_score=args.min_score,
    )


if __name__ == '__main__':
    main()
