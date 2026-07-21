#!/usr/bin/env python3
"""
Transcript-level minimap2 mapping — CAT2 ``txTM`` mode
"""

import bisect
import collections
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Shell helper
# ---------------------------------------------------------------------------

def _run(cmd, log=None, check=True):
    kwargs = dict(shell=True, executable='/bin/bash')
    if log:
        with open(log, 'a') as lf:
            result = subprocess.run(cmd, stdout=lf, stderr=lf, **kwargs)
    else:
        result = subprocess.run(cmd, **kwargs)
    if check and result.returncode != 0:
        sys.exit(f"Command failed (exit {result.returncode}): {cmd}")
    return result.returncode


def _tool_path(tool_name: str) -> str:
    """
    Return an executable path for an external dependency.

    Prefer PATH, but fall back to CAT2's vendored UCSC binaries in standalones/.
    This is important on SLURM nodes where user PATH can be minimal.
    """
    p = shutil.which(tool_name)
    if p:
        return p
    local = Path(__file__).resolve().parents[1] / "standalones" / tool_name
    if local.exists():
        return str(local)
    return tool_name


# ---------------------------------------------------------------------------
# CNV copy enumeration
# ---------------------------------------------------------------------------

def _psl_alignment_score(psl_fields):
    """Sort key: higher is better (identity, coverage, raw matches)."""
    try:
        matches    = int(psl_fields[0])
        mismatches = int(psl_fields[1])
        rep_match  = int(psl_fields[2])
        q_size     = int(psl_fields[10])
        block_sizes = [int(x) for x in psl_fields[18].rstrip(',').split(',') if x]
    except Exception:
        return (0.0, 0.0, 0)
    aligned = matches + mismatches + rep_match
    ident = (matches / aligned) if aligned > 0 else 0.0
    cov = (sum(block_sizes) / q_size) if q_size > 0 else 0.0
    return (ident, cov, matches)


def enumerate_copies(in_psl, out_psl):
    """
    Rename secondary alignments of the same transcript as <tx_id>-2, <tx_id>-3 …
    so every row has a unique query name (PSL column 9).
    """
    rows = []
    with open(in_psl) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#') or not line.strip() or len(fields) < 21:
                rows.append(('comment', line))
                continue
            if len(fields) not in (21, 23):
                continue
            rows.append(('aln', fields))

    # Rank alignments per transcript so the best keeps the original ID.
    per_q = collections.defaultdict(list)
    ordered = []
    for kind, data in rows:
        if kind == 'comment':
            ordered.append((kind, data))
        else:
            qname = data[9]
            per_q[qname].append(data)

    query_seen = collections.defaultdict(int)
    with open(out_psl, 'w') as fh:
        # Emit comments in original order, then alignments grouped by query.
        for kind, data in ordered:
            if kind == 'comment':
                fh.write(data)
        for qname, alns in per_q.items():
            alns_sorted = sorted(alns, key=_psl_alignment_score, reverse=True)
            for fields in alns_sorted:
                query_seen[qname] += 1
                if query_seen[qname] > 1:
                    fields = list(fields)
                    # Liftoff-style CNV suffix: _2, _3, ...
                    fields[9] = f"{qname}_{query_seen[qname]}"
                fh.write('\t'.join(fields) + '\n')

    total = sum(query_seen.values())
    multi = sum(1 for v in query_seen.values() if v > 1)
    print(f"    Copy enumeration: {total} alignments, "
          f"{len(query_seen)} unique transcripts, "
          f"{multi} with multiple copies")


# ---------------------------------------------------------------------------
# PSL filtering
# ---------------------------------------------------------------------------

def filter_psl(in_psl, out_psl, min_coverage=0.80, min_identity=0.80):
    """
    Keep alignments above absolute coverage and identity thresholds.
    Retains all CNV copies that pass the filter.

    Coverage = sum(blockSizes) / qSize
    Identity = matches / (matches + misMatches + repMatches)
    """
    kept = total = 0
    with open(in_psl) as fin, open(out_psl, 'w') as fout:
        for line in fin:
            if line.startswith('#') or not line.strip():
                fout.write(line)
                continue
            fields = line.split('\t')
            if len(fields) < 21:
                fout.write(line)
                continue
            if len(fields) not in (21, 23):
                continue
            try:
                matches    = int(fields[0])
                mismatches = int(fields[1])
                rep_match  = int(fields[2])
                q_size     = int(fields[10])
                block_sizes = [int(x) for x in fields[18].rstrip(',').split(',') if x]
            except ValueError:
                continue
            total += 1
            if q_size == 0:
                continue
            coverage = sum(block_sizes) / q_size
            aligned  = matches + mismatches + rep_match
            identity = matches / aligned if aligned > 0 else 0.0
            if coverage >= min_coverage and identity >= min_identity:
                fout.write(line)
                kept += 1
    print(f"    PSL filter: kept {kept}/{total} alignments "
          f"(cov≥{min_coverage:.0%}, id≥{min_identity:.0%})")


# Copy suffix convention:
# - Liftoff-style CNV copies are recognized downstream by an "_<int>" suffix.
# - Using "_" (not "-") also lets consensus code detect these as CNV txTM
#   placements and apply appropriate fragment handling.
_COPY_SUFFIX_PAT = re.compile(r'_\d+$')


def renumber_filtered_psl_copy_suffixes(psl_path):
    """
    Copy suffixes (-2, -3, …) are assigned before cov/id/copy filtering. If the
    best-ranked alignment is filtered out, the surviving hit can still be named
    TX-2 even though it is now the only (or best) mapping — transcript IDs no
    longer match txTM / RefSeq names and look 'lost' in GP diffs.

    Renumber each base transcript's *kept* rows so the highest-scoring alignment
    uses the bare transcript ID; extras become base-2, base-3, … (same scheme
    as enumerate_copies).
    """
    with open(psl_path) as fh:
        lines = fh.readlines()

    preamble = []
    groups = collections.defaultdict(list)
    for line in lines:
        if not line.strip():
            preamble.append(line)
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 21:
            preamble.append(line)
            continue
        try:
            int(fields[0])
        except ValueError:
            preamble.append(line)
            continue
        groups[_COPY_SUFFIX_PAT.sub('', fields[9])].append(fields)

    n_renamed = 0
    out_align = []
    for base in sorted(groups.keys()):
        alns = sorted(groups[base], key=_psl_alignment_score, reverse=True)
        for i, fields in enumerate(alns, start=1):
            old = fields[9]
            fields[9] = base if i == 1 else f"{base}_{i}"
            if old != fields[9]:
                n_renamed += 1
            out_align.append('\t'.join(fields) + '\n')

    with open(psl_path, 'w') as fh:
        for line in preamble:
            fh.write(line)
        fh.writelines(out_align)

    print(f"    PSL renumber: {len(out_align):,} rows, {n_renamed:,} query names adjusted "
          f"(best kept hit → bare transcript ID)")


def _build_ref_strand(ref_gp_path):
    """Return dict: transcript_name -> '+'/'-' from reference genePred."""
    strand = {}
    with open(ref_gp_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 3:
                continue
            strand[f[0]] = f[2]
    return strand


# ---------------------------------------------------------------------------
# Two-pass helpers
# ---------------------------------------------------------------------------

def _get_psl_query_names(psl_path):
    names = set()
    with open(psl_path) as f:
        for line in f:
            fields = line.split('\t')
            if len(fields) >= 21:
                try:
                    names.add(fields[9])
                except IndexError:
                    pass
    return names


def _extract_fasta_sequences(fa_path, names, out_fa):
    if not names:
        open(out_fa, 'w').close()
        return 0
    names_file = str(out_fa) + '.names'
    with open(names_file, 'w') as f:
        for n in sorted(names):
            f.write(n + '\n')
    _run(f"samtools faidx {fa_path} -r {names_file} > {out_fa}")
    return len(names)


def _merge_psls(psl_a, psl_b, out_psl):
    header_written = False
    with open(out_psl, 'w') as fout:
        for src in (psl_a, psl_b):
            with open(src) as fin:
                for line in fin:
                    if len(line.split('\t')) < 21:
                        if not header_written:
                            fout.write(line)
                        continue
                    header_written = True
                    fout.write(line)


# ---------------------------------------------------------------------------
# Reference GenePred parsing helpers
# ---------------------------------------------------------------------------

def _build_ref_exon_tx_space(ref_gp_path):
    """
    Parse reference GenePred → dict: name → list of (tx_start, tx_end) per exon,
    in 5'→3' transcript space (0-based, exclusive end).
    """
    result = {}
    with open(ref_gp_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            if len(f) < 10:
                continue
            try:
                name   = f[0]
                strand = f[2]
                ex_s   = [int(x) for x in f[8].rstrip(',').split(',') if x]
                ex_e   = [int(x) for x in f[9].rstrip(',').split(',') if x]
            except (ValueError, IndexError):
                continue
            exon_pairs = list(zip(ex_s, ex_e))
            if strand == '-':
                exon_pairs = list(reversed(exon_pairs))
            boundaries = []
            tx_pos = 0
            for es, ee in exon_pairs:
                length = ee - es
                boundaries.append((tx_pos, tx_pos + length))
                tx_pos += length
            result[name] = boundaries
    return result


def _build_ref_cds_tx_space(ref_gp_path):
    """
    Parse reference GenePred → dict: name → (cds_tx_start, cds_tx_end) or None.
    0-based, 5'→3' transcript space, exclusive end.
    """
    result = {}
    with open(ref_gp_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            if len(f) < 10:
                continue
            try:
                name   = f[0]
                strand = f[2]
                cds_s  = int(f[5])
                cds_e  = int(f[6])
                ex_s   = [int(x) for x in f[8].rstrip(',').split(',') if x]
                ex_e   = [int(x) for x in f[9].rstrip(',').split(',') if x]
            except (ValueError, IndexError):
                continue
            if cds_s >= cds_e:
                result[name] = None
                continue
            tx_start = tx_end = None
            tx_pos = 0
            if strand == '+':
                for es, ee in zip(ex_s, ex_e):
                    if tx_start is None and es <= cds_s < ee:
                        tx_start = tx_pos + (cds_s - es)
                    if tx_end is None and es < cds_e <= ee:
                        tx_end = tx_pos + (cds_e - es)
                    tx_pos += ee - es
            else:
                for es, ee in reversed(list(zip(ex_s, ex_e))):
                    if tx_start is None and es < cds_e <= ee:
                        tx_start = tx_pos + (ee - cds_e)
                    if tx_end is None and es <= cds_s < ee:
                        tx_end = tx_pos + (ee - cds_s)
                    tx_pos += ee - es
            if tx_start is not None and tx_end is not None and tx_start < tx_end:
                result[name] = (tx_start, tx_end)
            else:
                result[name] = None
    return result


def _build_ref_name2(ref_gp_path):
    """Return dict: transcript_name → gene_name (name2, field 11 of genePredExt)."""
    name2 = {}
    with open(ref_gp_path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 12:
                name2[f[0]] = f[11]
    return name2


def _build_ref_intron_transcript_order(ref_gp_path):
    """
    For each reference transcript, intron lengths between consecutive exons in
    transcript 5'→3' order (same indexing as _build_ref_exon_tx_space).

    Used to reject minimap2 chimeras that splice a 5' block into a distant locus
    with a target intron far larger than the reference allows between the same
    exon pair (e.g. XM_017021679.3 picking up ACOT2 sequence before DNAL1).
    """
    introns = {}
    with open(ref_gp_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            if len(f) < 10:
                continue
            try:
                name = f[0]
                strand = f[2]
                ex_s = [int(x) for x in f[8].rstrip(',').split(',') if x]
                ex_e = [int(x) for x in f[9].rstrip(',').split(',') if x]
            except (ValueError, IndexError):
                continue
            n = len(ex_s)
            if n < 2:
                introns[name] = []
                continue
            if strand == '+':
                tx_introns = [ex_s[i + 1] - ex_e[i] for i in range(n - 1)]
            else:
                # Minus strand: transcript 5' exon is the largest genomic block.
                tx_introns = []
                for k in range(n - 1):
                    gi = n - 1 - k
                    gj = n - 2 - k
                    tx_introns.append(ex_s[gi] - ex_e[gj])
            introns[name] = tx_introns
    return introns


def _trim_exon_target_chimeric(exon_target, introns_tx, ratio, pad_bp, floor_bp):
    """
    Split exon_target into runs separated by (a) non-consecutive ref exon index or
    (b) target intron >> reference intron between those exons. Keep the run with
    the largest total exonic span on the target.

    exon_target: ref_exon_index → [t_start, t_end]
    introns_tx: list length n_exons-1, intron after transcript exon i.
    """
    if not exon_target or ratio is None or float(ratio) <= 0:
        return exon_target
    if not introns_tx:
        return exon_target

    items = sorted(exon_target.items(), key=lambda kv: kv[1][0])
    runs = []
    current = [items[0]]

    def allowed_gap_after_exon(ei):
        if ei >= len(introns_tx):
            return 10**9
        ri = introns_tx[ei]
        return max(int(float(ri) * float(ratio)) + int(pad_bp), int(floor_bp))

    for i in range(len(items) - 1):
        (ei_a, (sa, ea)) = items[i]
        (ei_b, (sb, eb)) = items[i + 1]
        tg = sb - ea
        bad = ei_b != ei_a + 1 or tg > allowed_gap_after_exon(ei_a)
        if bad:
            runs.append(current)
            current = [items[i + 1]]
        else:
            current.append(items[i + 1])
    runs.append(current)

    def run_score(run):
        return sum(b - a for _, (a, b) in run)

    best = max(runs, key=run_score)
    return dict(best)


def write_gp_attrs(out_gp_path, ref_db_path, out_attrs_path):
    """
    Write a gp_attrs file for the transcripts in out_gp_path by looking up
    biotype information from the reference annotation database.

    The gp_attrs format (tab-separated, 3 columns) is:
        transcript_id  attribute_name  attribute_value

    The consensus module reads gene_biotype (and gene_name) from this file
    to assign biotypes to txTM/transcript_map transcripts that don't match
    entries in the reference DB by alignment ID (because transcript_map uses
    a different naming convention for multi-mapped copies).
    """
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    import tools.sqlInterface as sql

    # Load reference annotation: TranscriptId → GeneBiotype, GeneName
    ref_df = sql.load_annotation(ref_db_path)
    # Build lookup by base transcript ID (strip _N CNV suffix)
    tx_to_gene_biotype = {}
    tx_to_gene_name = {}
    for _, row in ref_df.iterrows():
        tx_id = row['TranscriptId']
        base = re.sub(r'_\d+$', '', tx_id)
        tx_to_gene_biotype[tx_id]  = row.get('GeneBiotype', '')
        tx_to_gene_biotype[base]   = row.get('GeneBiotype', '')
        tx_to_gene_name[tx_id]     = row.get('GeneName', '')
        tx_to_gene_name[base]      = row.get('GeneName', '')

    # Read transcript IDs from the output genePred (column 0 = name)
    written = 0
    with open(out_gp_path) as gp_in, open(out_attrs_path, 'w') as af:
        for line in gp_in:
            f = line.rstrip('\n').split('\t')
            if not f or not f[0]:
                continue
            tx_id = f[0]
            # Strip CNV suffix (_2, _3 …) added by enumerate_copies
            base_tx = re.sub(r'_\d+$', '', tx_id)
            # Also strip _N txTM CNV suffix
            base_tx2 = re.sub(r'_\d+$', '', base_tx)

            biotype = (tx_to_gene_biotype.get(tx_id)
                       or tx_to_gene_biotype.get(base_tx)
                       or tx_to_gene_biotype.get(base_tx2)
                       or '')
            gene_name = (tx_to_gene_name.get(tx_id)
                         or tx_to_gene_name.get(base_tx)
                         or tx_to_gene_name.get(base_tx2)
                         or '')

            if biotype:
                af.write(f"{tx_id}\tgene_biotype\t{biotype}\n")
                written += 1
            if gene_name:
                af.write(f"{tx_id}\tgene_name\t{gene_name}\n")

    print(f"    gp_attrs: wrote biotypes for {written:,} transcripts → {out_attrs_path}")


# ---------------------------------------------------------------------------
# CDS coordinate mapping
# ---------------------------------------------------------------------------

def _map_tx_pos_to_target(psl_fields, tx_pos):
    """
    Map a 0-based 5'→3' transcript position (inclusive) to target genome
    position (inclusive, forward-strand). Returns None if tx_pos falls in
    an unaligned region.
    """
    strand  = psl_fields[8][0]
    qsize   = int(psl_fields[10])
    bsizes  = [int(x) for x in psl_fields[18].rstrip(',').split(',') if x]
    qstarts = [int(x) for x in psl_fields[19].rstrip(',').split(',') if x]
    tstarts = [int(x) for x in psl_fields[20].rstrip(',').split(',') if x]

    if strand == '+':
        for qs, ts, bs in zip(qstarts, tstarts, bsizes):
            if qs <= tx_pos < qs + bs:
                return ts + (tx_pos - qs)
    else:
        for qs_rc, ts, bs in zip(qstarts, tstarts, bsizes):
            tx_blk_start = qsize - qs_rc - bs
            tx_blk_end   = qsize - qs_rc
            if tx_blk_start <= tx_pos < tx_blk_end:
                offset = tx_pos - tx_blk_start
                return ts + bs - 1 - offset
    return None


def _parse_psl_blocks(psl_fields):
    """Pre-parse a PSL row's block arrays once. Returns a dict with all
    coordinate arrays and (for '-' strand) per-block tx-space ranges sorted
    by tx_start to enable binary search."""
    strand  = psl_fields[8][0]
    qsize   = int(psl_fields[10])
    bsizes  = [int(x) for x in psl_fields[18].rstrip(',').split(',') if x]
    qstarts = [int(x) for x in psl_fields[19].rstrip(',').split(',') if x]
    tstarts = [int(x) for x in psl_fields[20].rstrip(',').split(',') if x]
    if strand == '+':
        tx_blk_start = qstarts
        tx_blk_end   = [qs + bs for qs, bs in zip(qstarts, bsizes)]
    else:
        tx_blk_start = [qsize - qs - bs for qs, bs in zip(qstarts, bsizes)]
        tx_blk_end   = [qsize - qs for qs in qstarts]
    order = sorted(range(len(bsizes)), key=lambda i: tx_blk_start[i])
    tx_starts_sorted = [tx_blk_start[i] for i in order]
    tx_ends_sorted   = [tx_blk_end[i]   for i in order]
    t_starts_sorted  = [tstarts[i]      for i in order]
    bsizes_sorted    = [bsizes[i]       for i in order]
    return {
        'strand': strand, 'qsize': qsize,
        'tx_s': tx_starts_sorted, 'tx_e': tx_ends_sorted,
        't_s':  t_starts_sorted,  'bs':   bsizes_sorted,
    }


def _map_tx_pos_via_parsed(parsed, tx_pos):
    """Fast O(log n_blocks) lookup using pre-parsed/sorted block arrays."""
    tx_s = parsed['tx_s']
    if not tx_s:
        return None
    i = bisect.bisect_right(tx_s, tx_pos) - 1
    if i < 0:
        return None
    if not (parsed['tx_s'][i] <= tx_pos < parsed['tx_e'][i]):
        return None
    bs = parsed['bs'][i]
    t_start = parsed['t_s'][i]
    if parsed['strand'] == '+':
        return t_start + (tx_pos - parsed['tx_s'][i])
    else:
        offset = tx_pos - parsed['tx_s'][i]
        return t_start + bs - 1 - offset


def _map_tx_pos_tolerant_parsed(parsed, tx_pos, snap_max=0):
    """Same as the tolerant version but takes pre-parsed blocks and uses
    binary search. Snap walks neighboring aligned positions rather than
    iterating by 1 bp at a time, so cost is O(log n_blocks) regardless of
    snap_max."""
    base = _map_tx_pos_via_parsed(parsed, tx_pos)
    if base is not None or snap_max <= 0:
        return base
    tx_s = parsed['tx_s']
    tx_e = parsed['tx_e']
    if not tx_s:
        return None
    # Find the closest aligned position within snap_max on either side.
    # Right side: the block whose start is >= tx_pos
    j = bisect.bisect_left(tx_s, tx_pos)
    right_pos = tx_s[j] if j < len(tx_s) else None
    # Left side: the block whose end is <= tx_pos
    i = bisect.bisect_right(tx_e, tx_pos) - 1
    left_pos  = tx_e[i] - 1 if i >= 0 else None
    candidates = []
    if right_pos is not None and abs(right_pos - tx_pos) <= snap_max:
        candidates.append((abs(right_pos - tx_pos), right_pos))
    if left_pos is not None and abs(left_pos - tx_pos) <= snap_max:
        candidates.append((abs(left_pos - tx_pos), left_pos))
    if not candidates:
        return None
    candidates.sort()
    return _map_tx_pos_via_parsed(parsed, candidates[0][1])


def _map_tx_pos_to_target_tolerant(psl_fields, tx_pos, snap_max=0):
    """
    Like _map_tx_pos_to_target, but if tx_pos falls in an unaligned region
    of the query, snap to the nearest mapped query position within snap_max
    bp and return the target coordinate of that position.

    Used by gene-region projection so that an exon endpoint that falls
    inside a small insertion/deletion in the gene-region alignment does
    not kill the entire projection. Quality is still gated by the
    overall PSL cov/id threshold, so a modest snap (e.g. 50 bp) does not
    introduce noise — it only extends the reach of already-good
    alignments into their unaligned fringes.

    snap_max=0 reproduces the strict behavior exactly.
    """
    # Backwards-compatible wrapper: parse + delegate. Callers in tight loops
    # should use _parse_psl_blocks + _map_tx_pos_tolerant_parsed directly.
    parsed = _parse_psl_blocks(psl_fields)
    return _map_tx_pos_tolerant_parsed(parsed, tx_pos, snap_max)


# ---------------------------------------------------------------------------
# Exon frame computation
# ---------------------------------------------------------------------------

def _compute_exon_frames(exon_starts, exon_ends, cds_start, cds_end, strand):
    """
    Compute per-exon reading frames (genePredExt field 14).
    Returns list of ints: -1 for non-CDS exons, 0/1/2 for CDS exons.
    """
    n = len(exon_starts)
    if cds_start >= cds_end:
        return [-1] * n
    frames   = [-1] * n
    cds_bases = 0
    order = range(n) if strand == '+' else range(n - 1, -1, -1)
    for i in order:
        es, ee = exon_starts[i], exon_ends[i]
        seg_s = max(es, cds_start)
        seg_e = min(ee, cds_end)
        if seg_e <= seg_s:
            continue
        frames[i] = cds_bases % 3
        cds_bases += seg_e - seg_s
    return frames


def _clip_cds_range_to_exon_hull(cds_s, cds_e, ex_starts, ex_ends):
    """
    Intersect half-open genomic CDS [cds_s, cds_e) with genePred exons.

    After chimeric exon trimming, PSL-mapped CDS endpoints can fall in introns
    relative to the written exon list; clip to exonic sequence so thickStart /
    thickEnd are valid for GFF/CDS mapping.

    Returns (thick_start, thick_end) with thick_end exclusive (genePred cdsEnd).
    Returns (None, None) if there is no exonic overlap or total exonic CDS < 3 bp.
    """
    if cds_s >= cds_e or len(ex_starts) != len(ex_ends):
        return None, None
    segs = []
    for es, ee in zip(ex_starts, ex_ends):
        lo = max(cds_s, es)
        hi = min(cds_e, ee)
        if lo < hi:
            segs.append((lo, hi))
    if not segs:
        return None, None
    exonic = sum(hi - lo for lo, hi in segs)
    if exonic < 3:
        return None, None
    return segs[0][0], segs[-1][1]


# ---------------------------------------------------------------------------
# PSL → genePredExt  (reference-guided exon merging)
# ---------------------------------------------------------------------------

_COPY_PAT = re.compile(r'_\d+$')


def _merged_exon_cap(ref_n_exons, is_copy, *, factor, extra, abs_cap):
    """
    Max merged exon count allowed vs reference exon count.

    Applied **after** reference-guided block merging (not on raw minimap2 block
    count, which is typically 10–50× ref exons for good spliced alignments).
    """
    if ref_n_exons is None or ref_n_exons <= 0 or factor is None or float(factor) <= 0:
        return None
    f = float(factor)
    e = int(extra)
    cap = int(max(ref_n_exons * f + e, ref_n_exons + e))
    return min(int(abs_cap), cap)


def _psl_to_genepred(
    in_psl,
    out_gp,
    min_intron=50,
    ref_gp_path=None,
    chimeric_intron_ratio=8.0,
    chimeric_intron_pad_bp=8000,
    chimeric_intron_floor_bp=35000,
    merged_exon_max_factor=3.5,
    merged_exon_max_extra=15,
    merged_exon_max_abs=200,
    copy_merged_exon_max_factor=2.5,
    copy_merged_exon_max_extra=8,
    copy_merged_exon_max_abs=80,
):
    """
    Pure-Python PSL → 15-field genePredExt conversion.

    Exon merging strategy (in priority order):
    1. Reference-guided (when ref_gp_path is provided):
       Each PSL block is mapped to 5'→3' transcript space and assigned to
       the reference exon it overlaps most. Target coords within each exon
       group are merged, producing exactly the reference exon count for
       well-aligned transcripts — regardless of within-exon D-operation
       false splits or species-specific micro-insertions.
       Optional chimeric trim: if consecutive merged exons are separated on the
       target by far more than the reference intron between those exons, keep
       only the contiguous run with the largest total exonic span (reduces
       minimap2 over-merging into neighboring loci).
    2. Gap-based fallback (no ref, or transcript absent from ref):
       Adjacent blocks with target gap < min_intron bp are merged.

    Transcript txStart/txEnd are always the hull of merged exons (not raw PSL
    bounds), so stray minimap2 micro-blocks outside merged exons cannot inflate
    the genePred span.

    cdsStart/cdsEnd are set to txStart/txEnd as placeholders; corrected by
    _add_cds_to_genepred.
    """
    ref_exon_tx = {}
    ref_introns_tx = {}
    if ref_gp_path and Path(ref_gp_path).exists():
        ref_exon_tx = _build_ref_exon_tx_space(ref_gp_path)
        ref_introns_tx = _build_ref_intron_transcript_order(ref_gp_path)
        print(f"    Reference-guided merging for {len(ref_exon_tx):,} transcripts")

    written = skipped_exon_cap = 0
    with open(in_psl) as fin, open(out_gp, 'w') as fout:
        for line in fin:
            f = line.rstrip('\n').split('\t')
            if len(f) < 21:
                continue
            try:
                int(f[0])
            except ValueError:
                continue
            try:
                strand      = f[8][0]
                name        = f[9]
                chrom       = f[13]
                t_start     = int(f[15])
                t_end       = int(f[16])
                q_size      = int(f[10])
                block_sizes = [int(x) for x in f[18].rstrip(',').split(',') if x]
                q_starts    = [int(x) for x in f[19].rstrip(',').split(',') if x]
                t_starts    = [int(x) for x in f[20].rstrip(',').split(',') if x]
            except (ValueError, IndexError):
                continue
            if not block_sizes:
                continue

            base = _COPY_PAT.sub('', name)
            is_copy = bool(_COPY_PAT.search(name))
            ref_exons = ref_exon_tx.get(base)
            ref_n_exons = len(ref_exons) if ref_exons else None

            if ref_exons:
                # Reference-guided: group PSL blocks by reference exon overlap
                exon_target = {}
                for qs, ts, bs in zip(q_starts, t_starts, block_sizes):
                    if strand == '+':
                        tx_blk_s = qs
                        tx_blk_e = qs + bs
                    else:
                        tx_blk_s = q_size - qs - bs
                        tx_blk_e = q_size - qs
                    best_ei = best_ov = None
                    for ei, (ex_ts, ex_te) in enumerate(ref_exons):
                        ov = min(tx_blk_e, ex_te) - max(tx_blk_s, ex_ts)
                        if best_ov is None or ov > best_ov:
                            best_ov = ov
                            best_ei = ei
                    if best_ei is not None and best_ov > 0:
                        if best_ei not in exon_target:
                            exon_target[best_ei] = [ts, ts + bs]
                        else:
                            exon_target[best_ei][0] = min(exon_target[best_ei][0], ts)
                            exon_target[best_ei][1] = max(exon_target[best_ei][1], ts + bs)

                if exon_target:
                    introns_tx = ref_introns_tx.get(base, [])
                    exon_target = _trim_exon_target_chimeric(
                        exon_target,
                        introns_tx,
                        chimeric_intron_ratio,
                        chimeric_intron_pad_bp,
                        chimeric_intron_floor_bp,
                    )
                    sorted_groups = sorted(exon_target.values(), key=lambda x: x[0])
                    exon_starts = [g[0] for g in sorted_groups]
                    exon_ends   = [g[1] for g in sorted_groups]
                else:
                    ref_exons = None  # degenerate — fall through to gap-based

            if not ref_exons:
                # Gap-based fallback
                exon_starts = [t_starts[0]]
                exon_ends   = [t_starts[0] + block_sizes[0]]
                for i in range(1, len(t_starts)):
                    gap = t_starts[i] - exon_ends[-1]
                    if gap < min_intron:
                        exon_ends[-1] = t_starts[i] + block_sizes[i]
                    else:
                        exon_starts.append(t_starts[i])
                        exon_ends.append(t_starts[i] + block_sizes[i])

            if not exon_starts:
                continue

            if ref_n_exons is not None:
                if is_copy:
                    cap = _merged_exon_cap(
                        ref_n_exons, True,
                        factor=copy_merged_exon_max_factor,
                        extra=copy_merged_exon_max_extra,
                        abs_cap=copy_merged_exon_max_abs,
                    )
                else:
                    cap = _merged_exon_cap(
                        ref_n_exons, False,
                        factor=merged_exon_max_factor,
                        extra=merged_exon_max_extra,
                        abs_cap=merged_exon_max_abs,
                    )
                if cap is not None and len(exon_starts) > cap:
                    skipped_exon_cap += 1
                    continue

            # Hull of merged exons — never use raw PSL tStart/tEnd (they can include
            # blocks that were not assigned to any reference exon).
            gp_tx_s = min(exon_starts)
            gp_tx_e = max(exon_ends)

            exon_frames = ','.join(['-1'] * len(exon_starts)) + ','
            fout.write('\t'.join([
                name, chrom, strand,
                str(gp_tx_s), str(gp_tx_e),
                str(gp_tx_s), str(gp_tx_e),   # cdsStart/cdsEnd placeholders
                str(len(exon_starts)),
                ','.join(str(x) for x in exon_starts) + ',',
                ','.join(str(x) for x in exon_ends) + ',',
                '0',       # score
                name,      # name2 placeholder (filled by _add_cds_to_genepred)
                'unk', 'unk',
                exon_frames,
            ]) + '\n')
            written += 1
    msg = f"    psl_to_genepred: {written:,} records written"
    if skipped_exon_cap:
        msg += f"; {skipped_exon_cap:,} dropped (merged exon count >> reference)"
    print(msg)


# ---------------------------------------------------------------------------
# CDS transfer
# ---------------------------------------------------------------------------

def _add_cds_to_genepred(in_gp, psl_path, ref_gp_path, out_gp):
    """
    Post-process _psl_to_genepred output (15-field genePredExt) to add:
      - cdsStart / cdsEnd (mapped from reference CDS via PSL blocks)
      - name2 (gene name)
      - cdsStartStat / cdsEndStat
      - exonFrames

    Failed CDS transfer for ref-coding transcripts falls back to non-coding
    placeholders (divergent / bad mappings should be removed earlier in PSL
    filtering, not here).
    """
    ref_cds   = _build_ref_cds_tx_space(ref_gp_path)
    ref_name2 = _build_ref_name2(ref_gp_path)

    psl_index = {}
    with open(psl_path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 21:
                psl_index[f[9]] = f

    transferred = fallback = 0
    with open(in_gp) as fin, open(out_gp, 'w') as fout:
        for line in fin:
            raw = line.rstrip('\n')
            if not raw or raw.startswith('#'):
                fout.write(raw + '\n')
                continue
            gp = raw.split('\t')
            if len(gp) < 15:
                fout.write(raw + '\n')
                continue

            name   = gp[0]
            strand = gp[2]
            tx_end = gp[4]
            ex_starts = [int(x) for x in gp[8].rstrip(',').split(',') if x]
            ex_ends   = [int(x) for x in gp[9].rstrip(',').split(',') if x]

            # Strip minimap2 CNV copy suffix (_2, _3, …) for fallback to primary transcript.
            base = _COPY_PAT.sub('', name)
            # Prefer exact ref transcript id (e.g. NM_*.2_1) for CDS/name2; fall back to base id.
            if name in ref_name2:
                gp[11] = ref_name2[name]
            elif base in ref_name2:
                gp[11] = ref_name2[base]
            else:
                gp[11] = base

            if name in ref_cds:
                cds_info = ref_cds[name]
            else:
                cds_info = ref_cds.get(base)
            psl_rec  = psl_index.get(name)
            ok = False
            if cds_info is not None and psl_rec is not None:
                cds_tx_s, cds_tx_e = cds_info
                pos_s = _map_tx_pos_to_target(psl_rec, cds_tx_s)
                pos_e = _map_tx_pos_to_target(psl_rec, cds_tx_e - 1)
                if pos_s is not None and pos_e is not None:
                    cds_s = min(pos_s, pos_e)
                    cds_e = max(pos_s, pos_e) + 1
                    cl_s, cl_e = _clip_cds_range_to_exon_hull(
                        cds_s, cds_e, ex_starts, ex_ends)
                    if cl_s is not None:
                        cds_s, cds_e = cl_s, cl_e
                        gp[5] = str(cds_s)
                        gp[6] = str(cds_e)
                        gp[12] = 'cmpl'
                        gp[13] = 'cmpl'
                        frames = _compute_exon_frames(
                            ex_starts, ex_ends, cds_s, cds_e, strand)
                        gp[14] = ','.join(str(fr) for fr in frames) + ','
                        transferred += 1
                        ok = True

            if not ok:
                gp[5] = tx_end
                gp[6] = tx_end
                gp[12] = 'none'
                gp[13] = 'none'
                gp[14] = ','.join(['-1'] * len(ex_starts)) + ','
                fallback += 1

            fout.write('\t'.join(gp) + '\n')

    total = transferred + fallback
    pct = (transferred / total) if total else 0.0
    print(f"    CDS transfer: {transferred:,}/{total:,} ({pct:.1%}); "
          f"{fallback:,} set non-coding")


def _psl_matches_by_qname(psl_path):
    """Map PSL query name (col 9) -> matches (col 0). Last row wins if duplicates."""
    out = {}
    with open(psl_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            flds = line.rstrip('\n').split('\t')
            if len(flds) < 21:
                continue
            try:
                out[flds[9]] = int(flds[0])
            except ValueError:
                continue
    return out


def _filter_gp_drop_weak_gene_loci(
    gp_path,
    psl_path,
    *,
    min_transcript_frac: float,
    min_match_frac: float,
    keep_if_match_frac_ge: Optional[float],
    escape_min_transcript_frac: float,
):
    """
    Drop genePred rows whose (gene, chromosome) locus is weak relative to the
    same gene's strongest locus elsewhere, using only:

      - transcript counts per chromosome (isoform / mapping structure), and
      - sum of PSL ``matches`` for those transcripts (alignment strength).

    Let nr = N(c)/Nmax, sr = S(c)/Smax (floats). A chromosome ``c`` is **kept** if:

      (nr >= min_transcript_frac and sr >= min_match_frac)
      OR (keep_if_match_frac_ge is not None and sr >= keep_if_match_frac_ge
          and nr >= escape_min_transcript_frac)

    The second clause rescues real CNVR second copies that have fewer isoforms
    on the duplicate but nearly as much total alignment mass as the primary
    locus. Paralog-style hits (e.g. LYPLA1 on X) stay weak in **both** nr and sr.

    Genes with a single chromosome are unchanged.

    Returns (n_rows_dropped, set of kept transcript names col0).
    """
    matches = _psl_matches_by_qname(psl_path)

    with open(gp_path, encoding='utf-8', errors='replace') as fh:
        raw_lines = fh.readlines()

    entries = []
    for line in raw_lines:
        s = line.rstrip('\n')
        if not s or s.startswith('#'):
            entries.append(('meta', line))
            continue
        flds = s.split('\t')
        if len(flds) < 12:
            entries.append(('meta', line))
            continue
        entries.append(('gp', line, flds))

    # gene -> chrom -> list of (transcript_id, matches)
    by_gene_chrom = collections.defaultdict(lambda: collections.defaultdict(list))
    for item in entries:
        if item[0] != 'gp':
            continue
        _, _line, flds = item
        gene = flds[11]
        chrom = flds[1]
        tx = flds[0]
        m = matches.get(tx, 0)
        by_gene_chrom[gene][chrom].append((tx, m))

    # gene -> chrom -> (N, S)
    ns = {}
    for gene, chmap in by_gene_chrom.items():
        ns[gene] = {}
        for chrom, rows in chmap.items():
            n = len(rows)
            ssum = sum(m for _tx, m in rows)
            ns[gene][chrom] = (n, ssum)

    dropped = 0
    kept_tx = set()
    out_lines = []
    tf = float(min_transcript_frac)
    mf = float(min_match_frac)
    esc = float(keep_if_match_frac_ge) if keep_if_match_frac_ge is not None else None
    emn = float(escape_min_transcript_frac)

    for item in entries:
        if item[0] == 'meta':
            out_lines.append(item[1])
            continue
        _, line, flds = item
        gene = flds[11]
        chrom = flds[1]
        gmap = ns.get(gene, {})
        if len(gmap) <= 1:
            out_lines.append(line)
            kept_tx.add(flds[0])
            continue
        nmax = max(v[0] for v in gmap.values())
        smax = max(v[1] for v in gmap.values())
        if nmax <= 0 or smax <= 0:
            out_lines.append(line)
            kept_tx.add(flds[0])
            continue
        nc, sc = gmap[chrom]
        nr = nc / float(nmax)
        sr = sc / float(smax) if smax > 0 else 0.0
        keep_main = (nr >= tf and sr >= mf)
        keep_esc = (
            esc is not None
            and sr >= esc
            and nr >= emn
        )
        if not (keep_main or keep_esc):
            dropped += 1
            continue
        out_lines.append(line)
        kept_tx.add(flds[0])

    with open(gp_path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.writelines(out_lines)
    return dropped, kept_tx


def _filter_psl_to_transcripts(psl_path, keep_tx: set):
    """Drop alignment rows whose query (col 9) is not in keep_tx; keep headers."""
    kept_lines = []
    n_drop = 0
    with open(psl_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                kept_lines.append(line)
                continue
            flds = line.rstrip('\n').split('\t')
            if len(flds) < 21:
                kept_lines.append(line)
                continue
            try:
                int(flds[0])
            except ValueError:
                kept_lines.append(line)
                continue
            if flds[9] in keep_tx:
                kept_lines.append(line)
            else:
                n_drop += 1
    with open(psl_path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.writelines(kept_lines)
    return n_drop


# ---------------------------------------------------------------------------
# Gene-region rescue (Liftoff-style)
#
# For transcripts that don't survive transcript-level mapping (no minimap2 hit,
# coverage < threshold on best hit, or dropped by later filters), align the
# **gene region** (txStart - flank .. txEnd + flank from the reference genome)
# to the target instead. The gene region carries intron sequence as anchor mass,
# so a transcript whose exons alone only cover 30% can still be placed via a
# contiguous gene-region alignment that the 50% transcript-level filter rejects.
#
# Output: one GenePred record per accepted alignment, with the original exon
# structure projected from reference-genome coords into target-genome coords
# via the PSL block alignment.
# ---------------------------------------------------------------------------


def _ref_gp_coords_for_ids(ref_gp_path, wanted_ids):
    """Parse reference GenePred for IDs in wanted_ids.

    Returns dict: tx_id -> dict(chrom, strand, tx_start, tx_end,
                                exon_starts, exon_ends, cds_start, cds_end,
                                gene_name).
    Coordinates are 0-based half-open in reference genome space.
    """
    out = {}
    with open(ref_gp_path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            f = s.split('\t')
            if len(f) < 10:
                continue
            tx = f[0]
            if tx not in wanted_ids:
                continue
            try:
                out[tx] = dict(
                    chrom       = f[1],
                    strand      = f[2],
                    tx_start    = int(f[3]),
                    tx_end      = int(f[4]),
                    cds_start   = int(f[5]),
                    cds_end     = int(f[6]),
                    exon_starts = [int(x) for x in f[8].rstrip(',').split(',') if x],
                    exon_ends   = [int(x) for x in f[9].rstrip(',').split(',') if x],
                    gene_name   = f[11] if len(f) > 11 else tx,
                )
            except (ValueError, IndexError):
                continue
    return out


def _build_gene_region_fasta(ref_coords, ref_genome_fa, out_fa,
                              flank_bp=2000, chrom_sizes=None):
    """Extract gene-region genomic sequence for each tx_id in ref_coords.

    Writes FASTA records keyed by tx_id (one record per tx). Returns dict:
        tx_id -> region_start_0  (0-based ref-genome coord where the FASTA seq
                                  starts; needed later to convert ref-genome
                                  positions into query coordinates).
    """
    if not ref_coords:
        return {}

    # Build samtools faidx region list (1-based inclusive). samtools clamps to
    # chrom end automatically, but we explicitly clamp to 0 on the low side.
    regions = []   # list of (tx_id, chrom, start_1, end_1, region_start_0)
    for tx, rec in ref_coords.items():
        chrom = rec['chrom']
        r_start_0 = max(0, rec['tx_start'] - flank_bp)
        r_end_0   = rec['tx_end'] + flank_bp
        if chrom_sizes and chrom in chrom_sizes:
            r_end_0 = min(r_end_0, chrom_sizes[chrom])
        if r_end_0 <= r_start_0:
            continue
        regions.append((tx, chrom, r_start_0 + 1, r_end_0, r_start_0))

    if not regions:
        return {}

    out_fa_path = Path(out_fa)
    list_file = out_fa_path.with_suffix(out_fa_path.suffix + '.regions.tmp')
    with open(list_file, 'w') as fh:
        for _, chrom, s_1, e_1, _ in regions:
            fh.write(f"{chrom}:{s_1}-{e_1}\n")

    try:
        proc = subprocess.run(
            ["samtools", "faidx", str(ref_genome_fa), "-r", str(list_file)],
            capture_output=True, text=True, check=True,
        )
    finally:
        list_file.unlink(missing_ok=True)

    # Parse samtools output into a {region_str: seq} map
    region_to_seq = {}
    cur_key = None
    cur_buf = []
    for line in proc.stdout.splitlines():
        if line.startswith('>'):
            if cur_key is not None:
                region_to_seq[cur_key] = ''.join(cur_buf).upper()
            cur_key = line[1:].split()[0]
            cur_buf = []
        else:
            cur_buf.append(line)
    if cur_key is not None:
        region_to_seq[cur_key] = ''.join(cur_buf).upper()

    offsets = {}
    with open(out_fa, 'w') as fa_out:
        for tx, chrom, s_1, e_1, r_start_0 in regions:
            seq = region_to_seq.get(f"{chrom}:{s_1}-{e_1}")
            if not seq:
                continue
            fa_out.write(f">{tx}\n")
            for i in range(0, len(seq), 60):
                fa_out.write(seq[i:i + 60] + '\n')
            offsets[tx] = r_start_0
    return offsets


def _merge_chained_psl_records(in_psl, out_psl,
                               target_cluster_gap_bp=1_000_000,
                               query_overlap_tolerance_bp=500):
    """Merge PSL records that came from the same minimap2 chained alignment
    but were split into primary + supplementary BAM records.

    Why this exists: minimap2 chains large genomic alignments (e.g. a 200 kb
    BRAF gene region against the syntenic macaque locus) into ONE chain, but
    when that chain spans gaps larger than what a single CIGAR can encode
    cleanly, the SAM/BAM spec forces it to emit one primary plus N
    supplementary records (same query, same target chromosome, same strand,
    non-overlapping query ranges, with hard-clips marking the offsets).
    `bamToPsl` then turns each into a separate PSL row, and downstream
    coverage/identity filters reject each tiny piece (e.g. 1 kb / 200 kb
    query = 0.5% coverage) even though the union spans ~97% of the query.

    Merge rule:
      1. Group rows by (qname, tname, strand).
      2. Sort each group by tStart and cluster into target-loci such that
         consecutive rows within `target_cluster_gap_bp` join the same
         cluster (separates real chained alignments at the syntenic locus
         from off-target secondaries elsewhere on the same chrom).
      3. Within each cluster, sort by qStart; if query ranges are
         non-overlapping (or overlap by at most `query_overlap_tolerance_bp`
         which is the typical hard-clip seam slop), concatenate the blocks
         and sum counts. Otherwise keep records separate (they're true
         duplicate alignments to the same locus).
    """
    header_lines = []
    rows = []
    with open(in_psl) as fh:
        for line in fh:
            if not line.strip():
                header_lines.append(line); continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 21:
                header_lines.append(line); continue
            try:
                int(f[0])
            except ValueError:
                header_lines.append(line); continue
            rows.append({'qname': f[9], 'tname': f[13],
                         'strand': f[8], 'f': f})

    groups = collections.defaultdict(list)
    for r in rows:
        groups[(r['qname'], r['tname'], r['strand'])].append(r['f'])

    n_in = len(rows)
    n_clusters_merged = 0
    out_records = []

    def _merge_cluster(srt):
        """Try to merge srt (already sorted by qStart). Returns one row if
        merge ok, else returns srt unchanged."""
        if len(srt) == 1:
            return [srt[0]]
        prev_qe = -1
        for f in srt:
            qs = int(f[11]); qe = int(f[12])
            if qs < prev_qe - query_overlap_tolerance_bp:
                return srt
            if qs < prev_qe:
                pass  # tolerable seam overlap
            prev_qe = max(prev_qe, qe)
        m = list(srt[0])
        ints_sum = [0, 1, 2, 3, 4, 5, 6, 7]
        for ci in ints_sum:
            m[ci] = str(sum(int(f[ci]) for f in srt))
        m[11] = str(min(int(f[11]) for f in srt))
        m[12] = str(max(int(f[12]) for f in srt))
        m[15] = str(min(int(f[15]) for f in srt))
        m[16] = str(max(int(f[16]) for f in srt))
        block_sizes = []
        q_starts = []
        t_starts = []
        for f in srt:
            block_sizes.extend(x for x in f[18].rstrip(',').split(',') if x)
            q_starts.extend(  x for x in f[19].rstrip(',').split(',') if x)
            t_starts.extend(  x for x in f[20].rstrip(',').split(',') if x)
        m[17] = str(len(block_sizes))
        m[18] = ','.join(block_sizes) + ','
        m[19] = ','.join(q_starts) + ','
        m[20] = ','.join(t_starts) + ','
        return [m]

    emitted = set()
    for r in rows:
        k = (r['qname'], r['tname'], r['strand'])
        if k in emitted:
            continue
        emitted.add(k)
        grp = groups[k]
        if len(grp) == 1:
            out_records.append(grp[0])
            continue
        # Cluster by target proximity: walk sorted-by-tStart, start a new
        # cluster whenever the gap to the previous tEnd exceeds the cutoff.
        by_t = sorted(grp, key=lambda f: int(f[15]))
        clusters = []
        cur = [by_t[0]]
        cur_t_end = int(by_t[0][16])
        for f in by_t[1:]:
            ts = int(f[15])
            if ts - cur_t_end > target_cluster_gap_bp:
                clusters.append(cur)
                cur = [f]
                cur_t_end = int(f[16])
            else:
                cur.append(f)
                cur_t_end = max(cur_t_end, int(f[16]))
        clusters.append(cur)
        for cl in clusters:
            srt = sorted(cl, key=lambda f: int(f[11]))
            merged = _merge_cluster(srt)
            if len(merged) == 1 and len(srt) > 1:
                n_clusters_merged += 1
            out_records.extend(merged)

    with open(out_psl, 'w') as fout:
        for hl in header_lines:
            fout.write(hl)
        for f in out_records:
            fout.write('\t'.join(f) + '\n')

    print(f"    chained-PSL merge: {n_in} rows -> {len(out_records)} rows "
          f"({n_clusters_merged} target-clusters merged)")


def _project_gene_region_psl_to_gp(
    psl_path, out_gp, ref_coords, region_offsets,
    min_coverage=0.50, min_identity=0.50, max_span_ratio=3.60,
    snap_max_bp=0,
):
    """Project each gene-region PSL alignment back into a GenePred record.

    For each PSL row whose query is a gene-region of <tx_id>:
      1. Filter by cov/id (computed against the gene region, which is what
         minimap2 actually aligned — same yardstick Liftoff uses).
      2. For each reference exon (in reference-genome coords), compute its
         position in the query sequence (= ref_coord - region_offset).
      3. Map each query position into the target via the PSL blocks.
      4. Emit a GenePred row whose exons are the projected target coords,
         preserving the reference exon count and CDS structure when possible.
    """
    written = drop_covid = drop_proj = drop_span = 0

    with open(psl_path) as fin, open(out_gp, 'w') as fout:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 21:
                continue
            try:
                int(f[0])
            except ValueError:
                continue

            qname = f[9]
            base = _COPY_PAT.sub('', qname)
            rec = ref_coords.get(base)
            off = region_offsets.get(base)
            if rec is None or off is None:
                continue

            try:
                matches    = int(f[0])
                mismatches = int(f[1])
                rep_match  = int(f[2])
                q_size     = int(f[10])
                t_start    = int(f[15])
                t_end      = int(f[16])
                block_sizes = [int(x) for x in f[18].rstrip(',').split(',') if x]
            except (ValueError, IndexError):
                continue
            if q_size == 0:
                continue

            coverage = sum(block_sizes) / q_size
            denom = matches + mismatches + rep_match
            identity = matches / denom if denom > 0 else 0.0
            if coverage < min_coverage or identity < min_identity:
                drop_covid += 1
                continue

            ref_span = rec['tx_end'] - rec['tx_start']
            t_span   = t_end - t_start
            if ref_span > 0 and t_span > max(int(ref_span * max_span_ratio),
                                              ref_span + 28000):
                drop_span += 1
                continue

            # Project each reference exon endpoint to target.
            # For highly fragmented gene-region alignments (e.g. a 200 kb gene
            # whose alignment has small unaligned gaps at intron boundaries),
            # a single exon falling in an unaligned region used to kill the
            # whole transcript. We now keep the transcript as long as we
            # successfully project at least `min_proj_frac` of its exons, and
            # silently skip exons that fall in gaps.
            min_proj_frac = 0.50
            n_exons_total = len(rec['exon_starts'])
            min_proj_required = max(1, int(n_exons_total * min_proj_frac))

            # Cheap pre-filter: how many exons fall inside the [qStart,qEnd]
            # of this PSL alignment? If fewer than min_proj_required do,
            # there's no point invoking the expensive per-position projector.
            # This avoids spending O(snap_max_bp * n_blocks) per exon on
            # secondary alignments that only cover 5 kb of a 200 kb gene.
            q_start = int(f[11])
            q_end   = int(f[12])
            n_exons_in_window = 0
            for es_ref, ee_ref in zip(rec['exon_starts'], rec['exon_ends']):
                q_s = es_ref - off
                q_e = ee_ref - off - 1
                if 0 <= q_s < q_size and 0 <= q_e < q_size:
                    if q_s >= q_start and q_e < q_end:
                        n_exons_in_window += 1
            if n_exons_in_window < min_proj_required:
                drop_proj += 1
                continue

            exon_targets = []
            n_skipped = 0
            max_skipped = n_exons_total - min_proj_required
            parsed = _parse_psl_blocks(f)
            for es_ref, ee_ref in zip(rec['exon_starts'], rec['exon_ends']):
                if n_skipped > max_skipped:
                    break
                q_s = es_ref - off
                q_e = ee_ref - off - 1
                if q_s < 0 or q_e < 0 or q_s >= q_size or q_e >= q_size:
                    n_skipped += 1
                    continue
                t_s = _map_tx_pos_tolerant_parsed(parsed, q_s, snap_max_bp)
                t_e = _map_tx_pos_tolerant_parsed(parsed, q_e, snap_max_bp)
                if t_s is None or t_e is None:
                    n_skipped += 1
                    continue
                lo = min(t_s, t_e)
                hi = max(t_s, t_e) + 1
                exon_targets.append((lo, hi))
            if len(exon_targets) < min_proj_required:
                drop_proj += 1
                continue

            exon_targets.sort(key=lambda x: x[0])
            # Sanitize: snap-tolerance can produce exons that overlap on the
            # target (an exon's end snaps a few bp into the next exon's
            # region). That makes downstream intron intervals invert
            # (start > stop), which trips ChromosomeInterval's assertion in
            # consensus_runner. Walk the sorted exons and either drop or
            # fuse overlaps so introns are always valid.
            sanitized = []
            for lo, hi in exon_targets:
                if hi <= lo:
                    continue
                if sanitized and lo < sanitized[-1][1]:
                    # Overlap with previous exon: extend previous to cover
                    # this one. This merges adjacent fragments at the cost
                    # of losing an intron, which is the correct call here
                    # since the alignment didn't see a clear gap anyway.
                    prev_lo, prev_hi = sanitized[-1]
                    sanitized[-1] = (prev_lo, max(prev_hi, hi))
                else:
                    sanitized.append((lo, hi))
            exon_targets = sanitized
            if len(exon_targets) < min_proj_required:
                drop_proj += 1
                continue
            ts_list = [x[0] for x in exon_targets]
            te_list = [x[1] for x in exon_targets]

            # PSL query strand: '+' means query forward against target,
            # '-' means query reverse-complemented. Output GenePred strand =
            # ref strand XOR PSL strand.
            psl_strand = f[8][0]
            ref_strand = rec['strand']
            out_strand = ref_strand if psl_strand == '+' else (
                '+' if ref_strand == '-' else '-')

            gp_tx_s = ts_list[0]
            gp_tx_e = te_list[-1]

            # CDS projection
            cds_s_gp = gp_tx_e
            cds_e_gp = gp_tx_e
            cds_stat = 'none'
            exon_frames = [-1] * len(ts_list)
            if rec['cds_end'] > rec['cds_start']:
                q_cs = rec['cds_start'] - off
                q_ce = rec['cds_end'] - off - 1
                if 0 <= q_cs < q_size and 0 <= q_ce < q_size:
                    t_cs = _map_tx_pos_to_target(f, q_cs)
                    t_ce = _map_tx_pos_to_target(f, q_ce)
                    if t_cs is not None and t_ce is not None:
                        cds_s = min(t_cs, t_ce)
                        cds_e = max(t_cs, t_ce) + 1
                        cl_s, cl_e = _clip_cds_range_to_exon_hull(
                            cds_s, cds_e, ts_list, te_list)
                        if cl_s is not None:
                            cds_s_gp, cds_e_gp = cl_s, cl_e
                            cds_stat = 'cmpl'
                            exon_frames = _compute_exon_frames(
                                ts_list, te_list, cds_s_gp, cds_e_gp, out_strand)

            fout.write('\t'.join([
                qname, f[13], out_strand,
                str(gp_tx_s), str(gp_tx_e),
                str(cds_s_gp), str(cds_e_gp),
                str(len(ts_list)),
                ','.join(str(x) for x in ts_list) + ',',
                ','.join(str(x) for x in te_list) + ',',
                '0',
                rec['gene_name'],
                cds_stat, cds_stat,
                ','.join(str(fr) for fr in exon_frames) + ',',
            ]) + '\n')
            written += 1

    print(f"    gene-region projection: wrote {written:,} records "
          f"(dropped cov/id={drop_covid:,} span={drop_span:,} proj-failed={drop_proj:,})")
    return written


def _read_chrom_sizes(genome_fa):
    """Build {chrom: length} dict from genome_fa.fai (created by samtools faidx)."""
    fai = Path(str(genome_fa) + '.fai')
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(genome_fa)], check=True)
    sizes = {}
    with open(fai) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2:
                try:
                    sizes[p[0]] = int(p[1])
                except ValueError:
                    pass
    return sizes


def _run_gene_region_rescue(
    target_fa, ref_genome_fa, ref_gp_path,
    missing_ids,
    out_dir, tmp_dir, threads, log_path,
    min_coverage=0.50, min_identity=0.50,
    flank_bp=2000, max_secondary=10, secondary_ratio=0.5,
    minimap2_preset='asm10',
    extra_minimap2='',
    max_span_ratio=3.60,
    snap_max_bp=0,
    tag='rescue',
):
    """Run Liftoff-style gene-region rescue for missing_ids.

    Returns (rescue_gp_path, rescue_psl_path) Path objects (rescue_gp_path may
    point to an empty file if no rescue records were produced).
    """
    out_dir = Path(out_dir)
    tmp_dir = Path(tmp_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # All rescue artifacts live in *tmp_dir* (which the caller already
    # places under a per-genome subdirectory). Writing the .psl / .gp into
    # *out_dir* would race when multiple genomes run in parallel, since
    # *out_dir* (e.g. ``work_dir/txTM``) is shared across genomes.
    rescue_fa      = tmp_dir / f'{tag}_gene_regions.fa'
    rescue_bam     = tmp_dir / f'{tag}.bam'
    rescue_raw_psl = tmp_dir / f'{tag}_raw.psl'
    rescue_psl     = tmp_dir / f'{tag}.psl'
    rescue_gp      = tmp_dir / f'{tag}.gp'

    if not missing_ids:
        rescue_gp.write_text('')
        return rescue_gp, rescue_psl

    print(f"\n  Gene-region rescue: collecting coords for {len(missing_ids):,} ids")
    ref_coords = _ref_gp_coords_for_ids(ref_gp_path, missing_ids)
    print(f"    found {len(ref_coords):,}/{len(missing_ids):,} in reference GP")
    if not ref_coords:
        rescue_gp.write_text('')
        return rescue_gp, rescue_psl

    chrom_sizes = _read_chrom_sizes(ref_genome_fa)
    print(f"  Gene-region rescue: extracting gene-region FASTA (flank={flank_bp:,} bp)")
    offsets = _build_gene_region_fasta(
        ref_coords, ref_genome_fa, str(rescue_fa),
        flank_bp=flank_bp, chrom_sizes=chrom_sizes,
    )
    print(f"    wrote {len(offsets):,} gene-region sequences")
    if not offsets:
        rescue_gp.write_text('')
        return rescue_gp, rescue_psl

    print(f"  Gene-region rescue: minimap2 -x {minimap2_preset} → bam → psl")
    # -r 1000,500000  : broaden chain bandwidth so fragmented alignments across
    #                   large introns/repeats are bridged into one PSL record.
    # -z 2000,2000    : relax z-drop so local low-score regions (divergent
    #                   intron stretches, repeats) don't terminate the chain.
    # These complement -ax {asm10/asm20} for long-gene genomic rescue where the
    # default minimap2 chaining params produce dozens of short fragments rather
    # than one full-length alignment of the gene region.
    _run(
        f"minimap2 -ax {minimap2_preset} --MD --eqx "
        f"-r 1000,500000 -z 2000,2000 "
        f"-p {secondary_ratio} -N {max_secondary} "
        f"{extra_minimap2} -t {threads} "
        f"{target_fa} {rescue_fa} "
        f"| samtools sort {_samtools_sort_flags(threads)} -o {rescue_bam} "
        f"&& samtools index {rescue_bam}",
        log=log_path,
    )
    bam_to_psl = _tool_path('bamToPsl')
    _run(f"{bam_to_psl} {rescue_bam} {rescue_raw_psl}", log=log_path)

    # Merge primary + supplementary BAM records (now separate PSL rows) back
    # into single PSL rows per chained alignment. Without this, downstream
    # cov/id filters reject each tiny chunk of a long chained alignment that
    # SAM had to split across multiple records.
    rescue_raw_psl_merged = tmp_dir / f'{tag}_raw_merged.psl'
    _merge_chained_psl_records(str(rescue_raw_psl), str(rescue_raw_psl_merged))

    enumerate_copies(str(rescue_raw_psl_merged), str(rescue_psl))

    _project_gene_region_psl_to_gp(
        str(rescue_psl), str(rescue_gp), ref_coords, offsets,
        min_coverage=min_coverage, min_identity=min_identity,
        max_span_ratio=max_span_ratio,
        snap_max_bp=snap_max_bp,
    )

    return rescue_gp, rescue_psl


def _slurm_mem_gb() -> float | None:
    """Return SLURM memory request in GiB when running under sbatch/srun."""
    raw = os.environ.get("SLURM_MEM_PER_NODE") or os.environ.get("SLURM_MEM_PER_CPU")
    if not raw:
        return None
    try:
        mb = float(raw)
    except ValueError:
        return None
    if mb <= 0:
        return None
    # SLURM_MEM_PER_CPU is per-core; multiply when only that is set.
    if os.environ.get("SLURM_MEM_PER_NODE"):
        return mb / 1024.0
    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "1") or "1")
    return (mb * max(1, cpus)) / 1024.0


def _effective_minimap2_threads(requested_threads: int) -> int:
    """Cap minimap2 -t so minimap2|samtools sort fits in the SLURM memory request.

    minimap2 and samtools sort run concurrently in a pipe; peak RSS is roughly
    (minimap2_threads * per_thread_gb) + (sort_threads * sort_mem_per_thread).
    """
    requested = max(1, int(requested_threads))
    mem_gb = _slurm_mem_gb()
    if mem_gb is None:
        return requested

    sort_threads = max(1, min(8, requested))
    sort_mem_gb = 0.768
    sort_gb = sort_threads * sort_mem_gb
    overhead_gb = 20.0
    per_thread_gb = 2.0
    budget_gb = mem_gb * 0.9 - sort_gb - overhead_gb
    if budget_gb <= 0:
        return max(1, min(requested, 4))
    cap = max(1, int(budget_gb / per_thread_gb))
    return max(1, min(requested, cap))


def _samtools_sort_flags(threads):
    """Cap sort parallelism; keep per-thread -m modest (pipe runs with minimap2)."""
    sort_threads = max(1, min(8, int(threads)))
    return f"-@ {sort_threads} -m 768M"


# ---------------------------------------------------------------------------
# Alignment pass
# ---------------------------------------------------------------------------

def _align_and_to_psl(target_fa, ref_tx_fa, bam, raw_psl,
                      threads, secondary_ratio, max_secondary,
                      extra_flags, log_path, label,
                      k=14, w=None):
    """Run minimap2 → samtools sort + index → bamToPsl for one pass.

    k/w override minimap2's minimizer k-mer size and window. Defaults match
    minimap2's ``-x splice`` defaults; pass2 uses smaller k/w to catch short
    ncRNA (tRNA/miRNA) that pass1's ``-k14`` cannot seed.
    """
    w_flag = f"-w {int(w)} " if w is not None else ""
    mm_threads = _effective_minimap2_threads(threads)
    if mm_threads < int(threads):
        print(
            f"  [{label}] capping minimap2 -t {threads} -> {mm_threads} "
            f"(SLURM memory budget)",
            flush=True,
        )
    _run(
        # Liftoff uses minimap2 with many secondaries (-N 50) and permissive
        # secondary ratio (-p 0.5), plus --eqx and --end-bonus 5.
        # We keep splice-aware transcript→genome alignment (-x splice -uf).
        #
        # -G 1500000 raises max intron length from minimap2's splice default
        # (200kb) to 1.5Mb. Many primate mRNAs (BRAF intron 1 ≈88kb is fine,
        # but ADK has multi-100kb introns, GPC5 ~1Mb, etc.) get fragmented at
        # the 200kb cap and the consensus then sees only short partial tx,
        # which fail the downstream min_pc_len_ratio_vs_reference filter.
        f"minimap2 -ax splice -uf -G 1500000 -k{int(k)} {w_flag}--MD --eqx --end-bonus 5 "
        f"-p {secondary_ratio} -N {max_secondary} "
        f"{extra_flags} -t {mm_threads} "
        f"{target_fa} {ref_tx_fa} "
        f"| samtools sort {_samtools_sort_flags(threads)} -o {bam} "
        f"&& samtools index {bam}",
        log=log_path,
    )
    print(f"  [{label}] bamToPsl")
    bam_to_psl = _tool_path("bamToPsl")
    _run(f"{bam_to_psl} {bam} {raw_psl}", log=log_path)


# ---------------------------------------------------------------------------
# Main pipeline entry point
# ---------------------------------------------------------------------------

def run_transcript_map(work_dir, ref_genome, genome, threads, out_dir, tmp_dir,
                       log_path, min_coverage, min_identity, secondary_ratio,
                       copy_min_identity=0.80,
                       max_span_ratio=3.60,
                       max_span_extra_bp=28000,
                       copy_min_coverage=0.78,
                       max_blocks_factor=0.0,
                       max_blocks_extra=28,
                       max_blocks_abs=220,
                       primary_max_blocks_factor=0.0,
                       primary_max_blocks_extra=40,
                       primary_max_blocks_abs=350,
                       merged_exon_max_factor=3.5,
                       merged_exon_max_extra=15,
                       merged_exon_max_abs=200,
                       copy_merged_exon_max_factor=2.5,
                       copy_merged_exon_max_extra=8,
                       copy_merged_exon_max_abs=80,
                       min_target_span_ratio=None,
                       min_ref_span_for_span_ratio_bp=15000,
                       end_block_min_bp=14,
                       end_block_max_gap_bp=2600,
                       filter_weak_gene_loci=True,
                       gene_locus_min_transcript_frac=0.45,
                       gene_locus_min_match_frac=0.38,
                       gene_locus_keep_if_match_frac_ge=0.80,
                       gene_locus_escape_min_transcript_frac=0.15,
                       chimeric_intron_ratio=8.0,
                       chimeric_intron_pad_bp=8000,
                       chimeric_intron_floor_bp=35000,
                       gene_region_rescue=True,
                       gene_region_flank_bp=2000,
                       gene_region_minimap2_preset='asm10',
                       gene_region_min_coverage=0.50,
                       gene_region_min_identity=0.50,
                       gene_region_max_secondary=10,
                       gene_region_secondary_ratio=0.5,
                       gene_region_snap_max_bp=50,
                       gene_region_rescue_deep=True,
                       gene_region_deep_flank_bp=5000,
                       gene_region_deep_minimap2_preset='asm20',
                       gene_region_deep_min_coverage=0.50,
                       gene_region_deep_min_identity=0.50,
                       gene_region_deep_max_secondary=50,
                       gene_region_deep_secondary_ratio=0.5,
                       gene_region_deep_extra_minimap2='--end-bonus 5',
                       gene_region_deep_snap_max_bp=2000,
                       extra_minimap2='', two_pass=False, ref_db_path=None):
    """
    Run the full transcript-map pipeline for one target genome.

    Returns (output_gp_path, filtered_psl_path) as Path objects.
    """
    wd = Path(work_dir)
    od = Path(out_dir)
    od.mkdir(parents=True, exist_ok=True)

    ref_tx_fa = wd / f"reference/{ref_genome}.fa"
    ref_gp    = wd / f"reference/{ref_genome}.gp"
    target_fa = wd / f"genome_files/{genome}.fa"

    for p, label in [
        (ref_tx_fa, "reference transcript FASTA"),
        (ref_gp,    "reference GenePred"),
        (target_fa, "target genome FASTA"),
    ]:
        if not p.exists():
            sys.exit(f"Missing required file ({label}): {p}")

    tmp = Path(tmp_dir) / genome
    tmp.mkdir(parents=True, exist_ok=True)

    bam1         = tmp / f"{genome}_pass1.bam"
    raw_psl1     = tmp / f"{genome}_pass1_raw.psl"
    copies_psl   = tmp / f"{genome}_copies.psl"
    filtered_psl = od  / f"{genome}_filtered.psl"
    output_gp    = od  / f"{genome}.gp"

    ref_strand = _build_ref_strand(ref_gp)
    # Reference transcript genomic spans (txEnd - txStart), used to drop
    # chimeric/minimap2 over-extended spliced alignments that jump into nearby genes.
    # Liftoff is conservative here; without this, a few transcripts can spuriously
    # absorb distant micro-blocks and inflate gene spans, triggering consensus
    # spurious-overlap removals (e.g. HINFP→ABCG4/NLRX1/NHERF4 neighborhood).
    ref_span = {}
    ref_exons = {}
    with open(ref_gp) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 5:
                continue
            try:
                ref_span[f[0]] = int(f[4]) - int(f[3])
                # genePred: name chrom strand txStart txEnd cdsStart cdsEnd exonCount ...
                if len(f) > 7:
                    ref_exons[f[0]] = int(f[7])
            except ValueError:
                continue

    def _fails_end_block_sanity(t_start, t_end, t_starts, block_sizes):
        """
        Reject alignments whose genomic span is extended by tiny far-away blocks.

        This is a common failure mode of spliced aligners with permissive secondary
        reporting: a mostly-correct alignment gets a few 1–10 bp blocks far away,
        inflating txStart/txEnd and causing downstream spurious-overlap removals.

        We implement a simple Liftoff-like conservatism:
        - If the first block is tiny and separated from the next by a large gap → reject
        - If the last block is tiny and separated from the previous by a large gap → reject
        """
        if not t_starts or not block_sizes or len(t_starts) != len(block_sizes):
            return False
        if len(t_starts) < 2:
            return False
        # Blocks in PSL are in increasing target coordinate order for '+' strand; bamToPsl
        # emits them in target order. To be safe, sort by target start.
        blocks = sorted([(ts, ts + bs, bs) for ts, bs in zip(t_starts, block_sizes)], key=lambda x: x[0])
        (s0, e0, b0) = blocks[0]
        (s1, e1, b1) = blocks[1]
        (sp, ep, bp) = blocks[-2]
        (sl, el, bl) = blocks[-1]
        if b0 < end_block_min_bp and (s1 - e0) > end_block_max_gap_bp:
            return True
        if bl < end_block_min_bp and (sl - ep) > end_block_max_gap_bp:
            return True
        return False

    total_steps = 7 if two_pass else 5

    print(f"\n[{genome}] Step 1/{total_steps}: minimap2 (splice-aware, CNV-aware)")
    _align_and_to_psl(target_fa, ref_tx_fa, bam1, raw_psl1,
                      threads, secondary_ratio, 50,
                      f"{extra_minimap2}",
                      log_path, "pass1")

    if two_pass:
        pass1_aligned  = _get_psl_query_names(raw_psl1)
        all_tx_names   = [
            line[1:].split()[0]
            for line in open(ref_tx_fa)
            if line.startswith('>')
        ]
        unaligned_names = [n for n in all_tx_names if n not in pass1_aligned]
        print(f"\n[{genome}] Step 2/{total_steps}: extract unaligned transcripts "
              f"({len(unaligned_names):,} of {len(all_tx_names):,})")

        unaligned_fa = tmp / f"{genome}_unaligned.fa"
        _extract_fasta_sequences(ref_tx_fa, unaligned_names, unaligned_fa)

        bam2     = tmp / f"{genome}_pass2.bam"
        raw_psl2 = tmp / f"{genome}_pass2_raw.psl"
        print(f"[{genome}] Step 3/{total_steps}: minimap2 pass 2 "
              f"(permissive, smaller k for short ncRNA)")
        # Smaller k/w lets minimap2 seed alignments for short ncRNA (tRNA, miRNA,
        # snRNA) whose 22–80 bp sequences cannot generate any k=14 minimizer.
        # Liftoff masks this by aligning the whole gene-region (incl. introns)
        # so the pre-miRNA's ~80 bp is embedded in 1–2 kb of flanking sequence.
        _align_and_to_psl(target_fa, unaligned_fa, bam2, raw_psl2,
                          threads, 0.1, 500,
                          f"--end-bonus 10 {extra_minimap2}",
                          log_path, "pass2",
                          k=8, w=5)
        merged_psl = tmp / f"{genome}_merged.psl"
        print(f"[{genome}] Step 4/{total_steps}: merge pass 1 + pass 2 PSLs")
        _merge_psls(raw_psl1, raw_psl2, merged_psl)
        input_for_copies = merged_psl
    else:
        input_for_copies = raw_psl1

    step = 5 if two_pass else 2
    print(f"[{genome}] Step {step}/{total_steps}: enumerate CNV copies")
    enumerate_copies(input_for_copies, copies_psl)

    step += 1
    print(f"[{genome}] Step {step}/{total_steps}: filter by coverage / identity")
    # Enforce reference strand consistency (reduces chrom/strand flips vs txTM),
    # and apply txTM-like copy threshold: only keep extra copies if they meet
    # a stricter identity (copy_min_identity, analogous to txTM -sc).
    kept = total = 0
    copy_kept = copy_total = 0
    strand_mismatch = 0
    span_filtered = 0
    span_short_filtered = 0
    endblock_filtered = 0
    blockcount_filtered = 0
    primary_block_filtered = 0
    copy_cov_filtered = 0
    use_copy_block_cap = (
        max_blocks_factor is not None and float(max_blocks_factor) > 0.0
    )
    use_primary_block_cap = (
        primary_max_blocks_factor is not None
        and float(primary_max_blocks_factor) > 0.0
    )
    with open(copies_psl) as fin, open(filtered_psl, 'w') as fout:
        for line in fin:
            if line.startswith('#') or not line.strip():
                fout.write(line)
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 21:
                continue
            if len(fields) not in (21, 23):
                continue
            qname = fields[9]
            # Copies are suffixed as _2, _3, ... (see enumerate_copies and _COPY_SUFFIX_PAT).
            # We must strip that to look up reference strand/span by original transcript ID.
            base = _COPY_SUFFIX_PAT.sub('', qname)
            # Strand in PSL is query-vs-target; constrain to match reference transcript strand.
            psl_strand = fields[8][0]
            ref_s = ref_strand.get(base)
            # Prefer reference strand consistency (txTM-like), but do not
            # hard-drop mismatches: for some loci minimap2 may only return
            # opposite-strand mappings. Hard-dropping can remove entire genes.
            if ref_s and psl_strand != ref_s:
                strand_mismatch += 1
            try:
                matches    = int(fields[0])
                mismatches = int(fields[1])
                rep_match  = int(fields[2])
                q_size     = int(fields[10])
                t_start    = int(fields[15])
                t_end      = int(fields[16])
                block_sizes = [int(x) for x in fields[18].rstrip(',').split(',') if x]
                t_starts    = [int(x) for x in fields[20].rstrip(',').split(',') if x]
            except ValueError:
                continue
            if q_size == 0:
                continue
            total += 1
            coverage = sum(block_sizes) / q_size
            aligned  = matches + mismatches + rep_match
            identity = matches / aligned if aligned > 0 else 0.0

            is_copy = bool(_COPY_SUFFIX_PAT.search(qname))
            nblocks = len(block_sizes)
            rex = ref_exons.get(base)

            # Block-count sanity: spliced alignments can contain hundreds of tiny blocks
            # (often 1–5 bp) that inflate downstream exon structure. Liftoff is more
            # conservative; reject mappings whose blockCount is wildly larger than the
            # reference exon count for this transcript.
            # CNV copies use a tight cap; primary alignments use a looser cap so normal
            # minimap2 fragmentation still passes, but pathological divergence is dropped
            # regardless of coding status.
            if rex is not None and use_copy_block_cap and is_copy:
                max_allowed_blocks = min(
                    int(max_blocks_abs),
                    int(max(float(rex) * float(max_blocks_factor), float(rex) + float(max_blocks_extra))),
                )
                if nblocks > max_allowed_blocks:
                    blockcount_filtered += 1
                    continue
            elif rex is not None and use_primary_block_cap and not is_copy:
                max_primary = min(
                    int(primary_max_blocks_abs),
                    int(max(
                        float(rex) * float(primary_max_blocks_factor),
                        float(rex) + float(primary_max_blocks_extra),
                    )),
                )
                if nblocks > max_primary:
                    primary_block_filtered += 1
                    continue

            # Span sanity: if this transcript exists in the reference annotation, reject
            # target mappings whose genomic span is wildly inflated versus the reference
            # transcript span. This catches chimeric splice alignments that hop into a
            # neighboring gene locus but still satisfy cov/id due to many micro-blocks.
            ref_s = ref_span.get(base)
            t_span = t_end - t_start
            if ref_s is not None:
                max_allowed = max(int(ref_s * float(max_span_ratio)), int(ref_s + max_span_extra_bp))
                if t_span > max_allowed:
                    span_filtered += 1
                    continue
                if (
                    min_target_span_ratio is not None
                    and float(min_target_span_ratio) > 0.0
                    and ref_s >= int(min_ref_span_for_span_ratio_bp)
                    and t_span < int(float(ref_s) * float(min_target_span_ratio))
                ):
                    span_short_filtered += 1
                    continue

            # Tiny far-away end blocks: drop alignments whose span is driven by micro-blocks.
            if _fails_end_block_sanity(t_start, t_end, t_starts, block_sizes):
                endblock_filtered += 1
                continue

            # txTM-like base thresholds
            if not (coverage >= min_coverage and identity >= min_identity):
                continue
            if is_copy:
                copy_total += 1
                if identity < copy_min_identity:
                    continue
                # Liftoff extra copies (-copies) default requires near-perfect identity
                # AND full exon/CDS coverage. Imitate that conservatism: don't keep
                # multi-locus copies unless they are close to full-length.
                if coverage < float(copy_min_coverage):
                    copy_cov_filtered += 1
                    continue
                copy_kept += 1

            fout.write('\t'.join(fields) + '\n')
            kept += 1

    print(f"    PSL filter: kept {kept}/{total} alignments "
          f"(cov≥{min_coverage:.0%}, id≥{min_identity:.0%}); "
          f"strand-mismatch dropped {strand_mismatch:,}; "
          f"span-filtered {span_filtered:,}; "
          f"span-too-short-vs-ref {span_short_filtered:,}; "
          f"end-block-filtered {endblock_filtered:,}; "
          f"copy-blockcount-filtered {blockcount_filtered:,}; "
          f"primary-blockcount-filtered {primary_block_filtered:,}")
    if copy_total:
        print(f"    Copy filter: kept {copy_kept:,}/{copy_total:,} copies "
              f"(id≥{copy_min_identity:.0%}, cov≥{copy_min_coverage:.0%}; "
              f"dropped_low_cov={copy_cov_filtered:,})")

    renumber_filtered_psl_copy_suffixes(str(filtered_psl))

    step += 1
    print(f"[{genome}] Step {step}/{total_steps}: PSL→GenePred + CDS transfer")
    tmp_gp = tmp / f"{genome}_no_cds.gp"
    _psl_to_genepred(
        str(filtered_psl),
        str(tmp_gp),
        ref_gp_path=str(ref_gp),
        chimeric_intron_ratio=chimeric_intron_ratio,
        chimeric_intron_pad_bp=chimeric_intron_pad_bp,
        chimeric_intron_floor_bp=chimeric_intron_floor_bp,
        merged_exon_max_factor=merged_exon_max_factor,
        merged_exon_max_extra=merged_exon_max_extra,
        merged_exon_max_abs=merged_exon_max_abs,
        copy_merged_exon_max_factor=copy_merged_exon_max_factor,
        copy_merged_exon_max_extra=copy_merged_exon_max_extra,
        copy_merged_exon_max_abs=copy_merged_exon_max_abs,
    )
    _add_cds_to_genepred(str(tmp_gp), str(filtered_psl), str(ref_gp), str(output_gp))

    gp_tx = set()
    with open(output_gp, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            gp_tx.add(s.split('\t', 1)[0])
    n_psl_sync = _filter_psl_to_transcripts(str(filtered_psl), gp_tx)
    if n_psl_sync:
        print(f"    PSL sync: dropped {n_psl_sync:,} rows not written to genePred "
              f"({len(gp_tx):,} transcripts in GP)")

    # Drop weak extra loci per gene using only PSL match totals and per-chromosome
    # transcript counts (no chromosome-type rules).
    if filter_weak_gene_loci:
        print(f"[{genome}] Filter weak multi-locus genes (transcript count + PSL matches)")
        n_g, keep_tx = _filter_gp_drop_weak_gene_loci(
            str(output_gp),
            str(filtered_psl),
            min_transcript_frac=float(gene_locus_min_transcript_frac),
            min_match_frac=float(gene_locus_min_match_frac),
            keep_if_match_frac_ge=gene_locus_keep_if_match_frac_ge,
            escape_min_transcript_frac=float(gene_locus_escape_min_transcript_frac),
        )
        n_p = _filter_psl_to_transcripts(str(filtered_psl), keep_tx)
        print(f"    Gene-locus filter: dropped {n_g:,} genePred rows, "
              f"{n_p:,} PSL rows (kept {len(keep_tx):,} transcripts)")

    # ------------------------------------------------------------------
    # Liftoff-style gene-region rescue.
    # Reference transcripts whose transcript-level alignment failed (no PSL
    # row in filtered_psl, or dropped during PSL → GP conversion / weak-locus
    # filter) get one more chance: align their gene region (txStart-flank ..
    # txEnd+flank from the reference genome) to the target with minimap2
    # asm-style, then project original exons through the alignment.
    # ------------------------------------------------------------------
    if gene_region_rescue:
        ref_genome_fa = wd / f"genome_files/{ref_genome}.fa"
        if not ref_genome_fa.exists():
            print(f"  [gene-region rescue] skip: missing reference genome FASTA "
                  f"{ref_genome_fa}")
        else:
            kept_base_ids = {_COPY_PAT.sub('', t) for t in
                             _read_gp_first_col(str(output_gp))}
            all_ref_tx = _read_gp_first_col(str(ref_gp))
            missing = all_ref_tx - kept_base_ids
            print(f"\n[{genome}] Gene-region rescue: "
                  f"{len(missing):,} of {len(all_ref_tx):,} reference transcripts "
                  f"missing after transcript-level mapping")
            if missing:
                rescue_gp, rescue_psl = _run_gene_region_rescue(
                    str(target_fa), str(ref_genome_fa), str(ref_gp),
                    missing, str(od), str(tmp), threads, log_path,
                    min_coverage=gene_region_min_coverage,
                    min_identity=gene_region_min_identity,
                    flank_bp=gene_region_flank_bp,
                    max_secondary=gene_region_max_secondary,
                    secondary_ratio=gene_region_secondary_ratio,
                    minimap2_preset=gene_region_minimap2_preset,
                    max_span_ratio=max_span_ratio,
                    snap_max_bp=gene_region_snap_max_bp,
                    tag='rescue',
                )
                _append_gp(str(rescue_gp), str(output_gp),
                           skip_tx=kept_base_ids)
                # Also fold the rescue PSL into the main filtered PSL so that
                # downstream consumers (store_psl_metrics, evaluator) see metrics
                # for these transcripts too.
                _append_psl(str(rescue_psl), str(filtered_psl),
                            only_tx=_read_gp_first_col(str(output_gp)))

            # ----------------------------------------------------------
            # Deep rescue pass: for transcripts still missing after the
            # first rescue, re-run gene-region alignment with a wider
            # divergence preset (asm20), larger flank, more secondary
            # alignments and Liftoff's --end-bonus 5. Same 50% cov/id
            # threshold, so this only adds hits that we would have
            # accepted today had minimap2 found them in pass 1.
            # ----------------------------------------------------------
            if gene_region_rescue_deep:
                kept_base_ids = {_COPY_PAT.sub('', t) for t in
                                 _read_gp_first_col(str(output_gp))}
                missing2 = all_ref_tx - kept_base_ids
                print(f"\n[{genome}] Deep gene-region rescue: "
                      f"{len(missing2):,} of {len(all_ref_tx):,} reference "
                      f"transcripts still missing after first rescue")
                if missing2:
                    rescue2_gp, rescue2_psl = _run_gene_region_rescue(
                        str(target_fa), str(ref_genome_fa), str(ref_gp),
                        missing2, str(od), str(tmp), threads, log_path,
                        min_coverage=gene_region_deep_min_coverage,
                        min_identity=gene_region_deep_min_identity,
                        flank_bp=gene_region_deep_flank_bp,
                        max_secondary=gene_region_deep_max_secondary,
                        secondary_ratio=gene_region_deep_secondary_ratio,
                        minimap2_preset=gene_region_deep_minimap2_preset,
                        extra_minimap2=gene_region_deep_extra_minimap2,
                        max_span_ratio=max_span_ratio,
                        snap_max_bp=gene_region_deep_snap_max_bp,
                        tag='rescue_deep',
                    )
                    _append_gp(str(rescue2_gp), str(output_gp),
                               skip_tx=kept_base_ids)
                    _append_psl(str(rescue2_psl), str(filtered_psl),
                                only_tx=_read_gp_first_col(str(output_gp)))

    # Write gp_attrs so the consensus module can assign biotypes to these transcripts
    if ref_db_path is not None:
        attrs_path = od / f"{genome}_txTM.gp_attrs"
        print(f"\n[{genome}] Writing gp_attrs from reference DB...")
        write_gp_attrs(str(output_gp), str(ref_db_path), str(attrs_path))

    return output_gp, filtered_psl


def _read_gp_first_col(gp_path):
    """Return set of column-0 names from a genePred."""
    ids = set()
    with open(gp_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith('#'):
                ids.add(s.split('\t', 1)[0])
    return ids


def _append_gp(src_gp, dst_gp, skip_tx=None):
    """Append rows from src_gp to dst_gp, skipping any whose base name (after
    stripping _N copy suffix) is in skip_tx. Skips silently if src_gp is empty
    or doesn't exist.
    """
    src = Path(src_gp)
    if not src.exists() or src.stat().st_size == 0:
        return 0
    skip_tx = skip_tx or set()
    n = 0
    with open(src) as fin, open(dst_gp, 'a') as fout:
        for line in fin:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            base = _COPY_PAT.sub('', s.split('\t', 1)[0])
            if base in skip_tx:
                continue
            fout.write(line if line.endswith('\n') else line + '\n')
            n += 1
    if n:
        print(f"  appended {n:,} gene-region rescue rows → {dst_gp}")
    return n


def _append_psl(src_psl, dst_psl, only_tx=None):
    """Append rows from src_psl to dst_psl, keeping only rows whose query name
    (col 9) is in only_tx. Skips silently if src_psl doesn't exist or is empty.
    """
    src = Path(src_psl)
    if not src.exists() or src.stat().st_size == 0:
        return 0
    n = 0
    with open(src) as fin, open(dst_psl, 'a') as fout:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            flds = line.rstrip('\n').split('\t')
            if len(flds) < 21:
                continue
            try:
                int(flds[0])
            except ValueError:
                continue
            if only_tx is not None and flds[9] not in only_tx:
                continue
            fout.write(line if line.endswith('\n') else line + '\n')
            n += 1
    if n:
        print(f"  appended {n:,} gene-region rescue PSL rows → {dst_psl}")
    return n
