#!/usr/bin/env python3
"""Correct consensus gene structures and names from miniprot protein hits.

Projection (transMap / augTM) can fuse neighbouring paralogs, drop long first
introns, and carry a wrong ``source_gene_common_name``. Miniprot already mapped
the protein DB; this pass uses those hits to:

1. Confirm CAT names. A CAT name is corrected only when another peptide is
   clearly better *and* that symbol is not carried anywhere else in the
   genome, so corrections move names rather than adding copies.
2. Name unnamed coding models only when they are real protein-coding copies:
   intact ORF, near-full-length vs the matched protein, intron structure
   consistent with the parent, unambiguous best symbol. Retrocopies,
   fragments and broken ORFs stay unnamed and lose protein_coding.
3. Replace chimeric / truncated CDS with the owned full-protein miniprot locus.
4. Rescue reference genes absent from the genome whose protein hits sit in a
   clean syntenic gap.

Intended to run after ``generate_consensus`` and before ``annotate_novel_genes``.
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import math
import statistics
import shutil
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from cat.convert_miniprot_to_genepred import (
    _exon_frames,
    paf_identity,
    parse_cigar_to_exons,
    parse_paf_row,
)
from cat.peptide_qc import clearly_broken_peptide, recode_selenocysteine, protein_from_transcript
from tools.bio import get_sequence_dict, reverse_complement
from tools.transcripts import GenePredTranscript, gene_pred_iterator

logger = logging.getLogger("protein_refine")

GN_RE = re.compile(r"(?:^|\s)GN=([^\s;]+)")
COPY_SUFFIX_RE = re.compile(r"_\d+$")
NA = {"", "N/A", "NA", "None", "none", ".", "nan", "NaN"}

MIN_HIT_PID = 0.50
MIN_HIT_QCOV = 0.40
MIN_PRIMARY_PID = 0.78
MIN_NAMED_PID = 0.72
MIN_PRIMARY_QCOV = 0.50
MIN_ORTHOLOG_PID = 0.68
CAT_OVERRIDE_MARGIN = 0.08
MIN_COPY_PID = 0.80
MIN_COPY_QCOV = 0.85
MIN_COPY_AA_FRAC = 0.85
MAX_COPY_AA_FRAC = 1.25
MIN_COPY_EXON_FRAC = 0.70
MIN_INTRON_BP = 50
NAME_MARGIN = 0.02
DEMOTE_CALLS = {"broken_orf", "fragment", "chimeric", "retrocopy"}
RESCUE_CANDIDATES_PER_SYMBOL = 3
RESCUE_MAX_NAMED_OVERLAP = 0.10
FULL_PROTEIN_QCOV = 0.80
MIN_PROTEIN_SPAN = 2_000
# Nucleotide projection: a clean long transMap/augTM parse is the real locus.
# Protein modes (miniprot / augMP) often emit a compact isoform — use typical
# length for those, never let them overwrite a clean projection peptide.
PROJECTION_MODES = {"transMap", "transMap_pairwise", "augTM", "augTM_pairwise"}
PROTEIN_MODES = {"miniprot", "augMP"}


def _na(val) -> bool:
    if val is None:
        return True
    s = str(val).strip()
    return s in NA


def _as_float(val, default=0.0) -> float:
    try:
        if _na(val):
            return default
        return float(val)
    except (TypeError, ValueError):
        return default


def overlap_bp(a: tuple[int, int], b: tuple[int, int]) -> int:
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    return max(0, hi - lo + 1)


def _percentile(vals: list[int], p: float) -> float:
    s = sorted(vals)
    if not s:
        return 0.0
    k = (len(s) - 1) * p / 100.0
    lo = int(k)
    hi = min(lo + 1, len(s) - 1)
    if lo == hi:
        return float(s[lo])
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


def load_query_symbols(
    protein_fasta: str | None, gp_attrs: str | None
) -> tuple[dict[str, str], dict[str, int], set[str]]:
    """Map protein query IDs to gene symbols, plus typical UniProt length per symbol.

    Parent length comes from the protein DB (reviewed ``sp|`` median, else a
    stub/readthrough-aware typical of all GN= peptides), not from GENCODE tags.
    """
    out: dict[str, str] = {}
    if gp_attrs and Path(gp_attrs).exists():
        with open(gp_attrs) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                tx_id, attr, value = parts[0], parts[1], parts[2]
                if attr == "gene_name" and value and value not in NA:
                    out[tx_id] = value
                    out[tx_id.split(".")[0]] = value
        logger.info("  transcript→symbol from gp_attrs: %s", f"{len(out):,}")
    by_sp: dict[str, list[int]] = defaultdict(list)
    by_all: dict[str, list[int]] = defaultdict(list)
    seleno_symbols: set[str] = set()
    if protein_fasta and Path(protein_fasta).exists():
        n_gn = 0
        cur_sym = None
        cur_sp = False
        buf: list[str] = []

        def _flush() -> None:
            if not cur_sym or not buf:
                return
            seq = "".join(buf)
            n = len(seq)
            if n <= 0:
                return
            by_all[cur_sym].append(n)
            if cur_sp:
                by_sp[cur_sym].append(n)
            if "U" in seq.upper():
                seleno_symbols.add(cur_sym)

        with open(protein_fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    _flush()
                    header = line[1:].strip()
                    qid = header.split()[0]
                    gn = GN_RE.search(header)
                    cur_sym = gn.group(1) if gn else None
                    cur_sp = "|sp|" in qid or qid.startswith("sp|")
                    buf = []
                    if gn:
                        out[qid] = gn.group(1)
                        n_gn += 1
                    continue
                buf.append(line.strip())
            _flush()
        logger.info("  UniProt GN= symbols from protein FASTA: %s", f"{n_gn:,}")
    uniprot_aa: dict[str, int] = {}
    for sym in set(by_all) | set(by_sp):
        if by_sp.get(sym):
            uniprot_aa[sym] = int(statistics.median(by_sp[sym]))
        elif by_all.get(sym):
            uniprot_aa[sym] = typical_peptide_length(by_all[sym])
    logger.info("  UniProt typical lengths: %s symbols", f"{len(uniprot_aa):,}")
    if seleno_symbols:
        logger.info("  UniProt selenoproteins (U in peptide): %s symbols", f"{len(seleno_symbols):,}")
    return out, uniprot_aa, seleno_symbols


def symbol_for_query(qid: str, query_symbols: dict[str, str]) -> str | None:
    if qid in query_symbols:
        return query_symbols[qid]
    base = COPY_SUFFIX_RE.sub("", qid)
    if base in query_symbols:
        return query_symbols[base]
    nover = base.split(".")[0]
    if nover in query_symbols:
        return query_symbols[nover]
    gn = GN_RE.search(qid)
    if gn:
        return gn.group(1)
    return None


def load_ref_gene_order(ref_gp: str | None, query_symbols: dict[str, str]) -> dict[str, list[str]]:
    """Human chromosome → gene symbols in genomic order (one interval per symbol)."""
    by_chrom: dict[str, dict[str, tuple[int, int]]] = defaultdict(dict)
    if not ref_gp or not Path(ref_gp).exists():
        return {}
    for tx in gene_pred_iterator(ref_gp):
        sym = query_symbols.get(tx.name) or query_symbols.get(tx.name.split(".")[0])
        if not sym:
            sym = tx.name2 if tx.name2 and tx.name2 not in NA else None
        if not sym:
            continue
        prev = by_chrom[tx.chromosome].get(sym)
        if prev is None:
            by_chrom[tx.chromosome][sym] = (tx.start, tx.stop)
        else:
            by_chrom[tx.chromosome][sym] = (min(prev[0], tx.start), max(prev[1], tx.stop))
    order: dict[str, list[str]] = {}
    for chrom, genes in by_chrom.items():
        order[chrom] = [s for s, _ in sorted(genes.items(), key=lambda kv: kv[1][0])]
    logger.info("  reference gene order on %s chromosomes", len(order))
    return order


def real_exon_count(segs: list[tuple[int, int]]) -> int:
    """Exons after merging blocks split by gaps too short to be real introns."""
    segs = sorted(segs)
    if not segs:
        return 0
    n = 1
    for (_a0, a1), (b0, _b1) in zip(segs, segs[1:]):
        if b0 - a1 >= MIN_INTRON_BP:
            n += 1
    return n


def coding_exon_count(tx) -> int:
    segs = []
    for e in tx.exon_intervals:
        s, t = max(e.start, tx.thick_start), min(e.stop, tx.thick_stop)
        if t > s:
            segs.append((s, t))
    return real_exon_count(segs)


def typical_peptide_length(vals: list[int]) -> int:
    """Typical length from a pile of UniProt peptides (no reviewed subset).

    Drop fragments shorter than 40% of the 75th percentile, then drop a
    readthrough that is >1.5× the core median. Stops RPP21 / ALG10 / EBP
    picking a stub or a neighbour-fusion when there is no ``sp|`` entry.
    """
    if not vals:
        return 0
    if len(vals) == 1:
        return int(vals[0])
    q75 = _percentile(vals, 75)
    core = [v for v in vals if v >= 0.40 * q75] or list(vals)
    longest = max(core)
    median_aa = statistics.median(core)
    if longest > 1.5 * median_aa and median_aa > 0:
        return int(median_aa)
    return int(longest)


def typical_isoform(
    recs: list[tuple[int, int]], uniprot_aa: int | None = None
) -> tuple[int, int]:
    """Pick (cds_aa, coding_exons) for a gene.

    Prefer the reference isoform nearest the UniProt typical peptide length
    (reviewed median). That is annotation-tag agnostic: DDI2 stays ~400 aa
    because UniProt DDI2 is ~399, not because of a MANE tag, while GDI1 /
    NOTCH2 / ALG10 keep their full proteins instead of GENCODE stubs.

    Without a UniProt length, drop fragment isoforms shorter than 40% of the
    75th percentile, then apply the 1.5× median readthrough rule on that core.
    """
    if not recs:
        return (0, 0)
    recs = list(recs)
    if uniprot_aa and uniprot_aa > 0:
        return min(recs, key=lambda r: (abs(r[0] - uniprot_aa), -r[1]))
    recs.sort()
    if len(recs) == 1:
        return recs[0]
    q75 = _percentile([r[0] for r in recs], 75)
    core = [r for r in recs if r[0] >= 0.40 * q75] or recs
    longest = max(core, key=lambda r: r[0])
    median_aa = statistics.median(r[0] for r in core)
    if longest[0] > 1.5 * median_aa and median_aa > 0:
        return min(core, key=lambda r: (abs(r[0] - median_aa), -r[1]))
    return longest


def load_ref_structure(
    ref_gp: str | None,
    query_symbols: dict[str, str],
    uniprot_aa: dict[str, int] | None = None,
) -> dict[str, tuple[int, int]]:
    """Reference symbol → (cds_aa, coding_exons) of its typical coding isoform.

    Completeness is judged against the gene, not the matched peptide, which
    can itself be a short isoform or a UniProt fragment.
    """
    uniprot_aa = uniprot_aa or {}
    by_sym: dict[str, list[tuple[int, int]]] = defaultdict(list)
    if not ref_gp or not Path(ref_gp).exists():
        return {}
    for tx in gene_pred_iterator(ref_gp):
        if tx.cds_size <= 0:
            continue
        sym = query_symbols.get(tx.name) or query_symbols.get(tx.name.split(".")[0])
        if not sym:
            sym = tx.name2 if tx.name2 and tx.name2 not in NA else None
        if not sym:
            continue
        by_sym[sym].append((tx.cds_size // 3, coding_exon_count(tx)))
    out: dict[str, tuple[int, int]] = {}
    for sym, recs in by_sym.items():
        out[sym] = typical_isoform(recs, uniprot_aa.get(sym))
    return out


def parent_structure(sym: str, ref_struct: dict[str, tuple[int, int]]) -> tuple[int, int]:
    return ref_struct.get(sym, (0, 0))


def readthrough_parts(sym: str, ref_symbols: set[str]) -> tuple[str, str] | None:
    """Split A-B only when both sides are reference gene symbols."""
    if not sym or "-" not in sym:
        return None
    for i, ch in enumerate(sym):
        if ch != "-" or i == 0 or i == len(sym) - 1:
            continue
        left, right = sym[:i], sym[i + 1 :]
        if left in ref_symbols and right in ref_symbols:
            return (left, right)
    return None


def remap_readthrough_hit(
    h: Hit,
    ref_symbols: set[str],
    ref_struct: dict[str, tuple[int, int]],
    present: Counter | None = None,
) -> Hit:
    """Retarget a readthrough peptide to a missing constituent.

    If one constituent is already a separate named gene, this locus is not a
    readthrough (SPECC1L next to ADORA2A labeled SPECC1L-ADORA2A). Otherwise
    keep A-B only when the hit is a full A-B model; a short hit is remapped
    to a still-missing constituent it matches.
    """
    parts = readthrough_parts(h.symbol, ref_symbols)
    if not parts:
        return h
    left, right = parts
    left_p = bool(present and present.get(left, 0) > 0)
    right_p = bool(present and present.get(right, 0) > 0)
    if left_p ^ right_p:
        missing = right if left_p else left
        return Hit(
            h.chrom, h.start, h.end, h.strand, h.query, missing,
            h.pid, h.qcov, h.qlen, h.n_exons,
        )
    covered = int(h.qcov * h.qlen) if h.qlen else 0
    if is_complete_copy(h.n_exons, covered, h.symbol, ref_struct):
        return h
    for part in parts:
        if present is not None and present.get(part, 0) > 0:
            continue
        if is_complete_copy(h.n_exons, covered, part, ref_struct, require_exons=False):
            return Hit(
                h.chrom, h.start, h.end, h.strand, h.query, part, h.pid, h.qcov, h.qlen, h.n_exons
            )
    return h


def hits_conflict(a: Hit, b: Hit) -> bool:
    if a.chrom != b.chrom or a.strand != b.strand:
        return False
    ov = overlap_bp(a.span, b.span)
    if ov <= 50:
        return False
    la, lb = max(1, a.length), max(1, b.length)
    return ov / la >= 0.50 or ov / lb >= 0.50


def pick_unique_hits(
    by_sym: dict[str, list[Hit]], ref_struct: dict[str, tuple[int, int]]
) -> list[Hit]:
    """One locus per symbol. Highest-completeness hit wins; a losing symbol
    keeps its next non-overlapping hit. Single pass — no swap loop."""
    pairs: list[tuple] = []
    for sym, hs in by_sym.items():
        for h in hs:
            pairs.append((completeness(h, ref_struct), h.pid, sym, h))
    pairs.sort(key=lambda x: (x[0], x[1]), reverse=True)
    used_sym: set[str] = set()
    taken: list[Hit] = []
    for _sc, _pid, sym, h in pairs:
        if sym in used_sym:
            continue
        if any(hits_conflict(h, t) for t in taken):
            continue
        used_sym.add(sym)
        taken.append(h)
    return taken


def completeness(h, ref_struct: dict[str, tuple[int, int]]) -> float:
    """Identity × fraction of the reference gene's peptide covered by the hit."""
    p_aa = parent_structure(h.symbol, ref_struct)[0]
    if not p_aa or not h.qlen:
        return h.score
    return h.pid * min(1.0, h.qcov * h.qlen / p_aa)


def is_complete_copy(
    n_exons: int,
    aa: int,
    sym: str,
    ref_struct: dict[str, tuple[int, int]],
    require_exons: bool = True,
) -> bool:
    """Full-length vs the parent peptide.

    ``require_exons`` rejects a 1-exon alignment of a ≥3-coding-exon parent
    (retrocopy). It does not demand matching UTR-inflated exon counts: TLR9
    (1 CDS exon, 2-exon genePred) and miniprot's single CDS block both pass
    when the peptide itself is complete.
    """
    p_aa, p_ex = parent_structure(sym, ref_struct)
    if not p_aa:
        return False
    if not MIN_COPY_AA_FRAC * p_aa <= aa <= MAX_COPY_AA_FRAC * p_aa:
        return False
    if require_exons and p_ex >= 3 and n_exons < 2:
        return False
    return True


def extra_copy_complete(
    n_exons: int,
    aa: int,
    sym: str,
    ref_struct: dict[str, tuple[int, int]],
    pid: float | None = None,
    qcov: float | None = None,
) -> bool:
    """Intact extra copy of a gene that is already named elsewhere.

    Stricter than rescue: the parent must have splice structure (≥3 coding
    exons) and the hit must keep that intron pattern. A complete peptide is
    not enough — that is how D4Z4/DUX4, USP17, NUMTs and GAPDH retrocopies
    inflate copy number. Absent-gene rescue still uses
    ``require_exons=False`` so a thin miniprot block can recover EBP/TLR9.
    """
    if pid is not None and pid < MIN_COPY_PID:
        return False
    if qcov is not None and qcov < MIN_COPY_QCOV:
        return False
    p_aa, p_ex = parent_structure(sym, ref_struct)
    if not p_aa or p_ex < 3:
        return False
    if not MIN_COPY_AA_FRAC * p_aa <= aa <= MAX_COPY_AA_FRAC * p_aa:
        return False
    min_ex = max(2, math.ceil(MIN_COPY_EXON_FRAC * p_ex))
    if n_exons < min_ex:
        return False
    return True


def extra_copy_locus_ok(
    chrom: str,
    symbol: str,
    named_spans: dict[tuple[str, str], list[tuple[int, int]]],
) -> bool:
    """Extra copies expand a local array, they do not seed new chromosomes.

    TSPY units stay on chrY next to already-named TSPY; a GAPDH peptide on
    chr7 does not become a second GAPDH.
    """
    return bool(named_spans.get((symbol, chrom)))


def _span_dist(pos: int, span: tuple[int, int]) -> int:
    s, e = span
    if s <= pos <= e:
        return 0
    return s - pos if pos < s else pos - e


def nearest_named_symbol(
    chrom: str,
    pos: int,
    named_spans: dict[tuple[str, str], list[tuple[int, int]]],
) -> str | None:
    """Closest already-named gene on this chromosome."""
    best: tuple[int, str] | None = None
    for (sym, c), spans in named_spans.items():
        if c != chrom:
            continue
        for span in spans:
            d = _span_dist(pos, span)
            if best is None or d < best[0]:
                best = (d, sym)
    return best[1] if best else None


def pick_extra_from_cluster(
    cluster: list[Hit],
    named_spans: dict[tuple[str, str], list[tuple[int, int]]],
    ref_struct: dict[str, tuple[int, int]],
) -> Hit:
    """When AMY1 and AMY2 both fit, keep the paralog already named next door.

    Peptide identity cannot split recent tandem arrays. Prefer the nearest
    named symbol when it has a complete hit within NAME_MARGIN of the best.
    """
    best_sc = max(completeness(h, ref_struct) for h in cluster)
    near = nearest_named_symbol(cluster[0].chrom, cluster[0].start, named_spans)
    if near:
        local = [
            h
            for h in cluster
            if h.symbol == near and completeness(h, ref_struct) >= best_sc - NAME_MARGIN
        ]
        if local:
            cluster = local
    return max(cluster, key=lambda h: (completeness(h, ref_struct), h.pid, -h.length))


def resolve_tied_hit(
    chrom: str,
    pos: int,
    hits: list[Hit],
    named_spans: dict[tuple[str, str], list[tuple[int, int]]],
    ref_struct: dict[str, tuple[int, int]],
) -> Hit | None:
    """When AMY1/AMY2 or TSPY1/TSPY3 all fit, name for the nearest named copy."""
    complete = [
        h
        for h in hits
        if extra_copy_complete(
            h.n_exons, int(h.qcov * h.qlen), h.symbol, ref_struct, h.pid, h.qcov
        )
    ]
    if not complete:
        return None
    return pick_extra_from_cluster(complete, named_spans, ref_struct)


def is_fragment_model(g: dict, sym: str, ref_struct: dict[str, tuple[int, int]]) -> bool:
    p_aa, p_ex = parent_structure(sym, ref_struct)
    if not p_aa:
        return False
    aa = int(g.get("cds_size") or 0) // 3
    if aa < MIN_COPY_AA_FRAC * p_aa:
        return True
    if p_ex >= 3:
        min_ex = max(2, math.ceil(MIN_COPY_EXON_FRAC * p_ex))
        if int(g.get("coding_exons") or 0) < min_ex:
            return True
    return False


def demote_redundant_fragments(
    genes: list[dict],
    assign: dict[int, dict],
    calls: dict[int, str],
    present: Counter,
    ref_struct: dict[str, tuple[int, int]],
) -> int:
    """Un-name truncated extra projections of a gene that already has a full copy.

    transMap leaves 7-exon 340 aa AMY1B/AMY2A models in the gaps of the CHM13
    array after it has already placed the complete genes. Those are not extra
    copies.
    """

    def _sym(gi: int, g: dict) -> str:
        if gi in assign:
            return assign[gi]["symbol"]
        return g["source_name"]

    complete = set()
    for gi, g in enumerate(genes):
        if g.get("biotype") != "protein_coding":
            continue
        sym = _sym(gi, g)
        if sym and not is_fragment_model(g, sym, ref_struct):
            complete.add(sym)
    n = 0
    for gi, g in enumerate(genes):
        if g.get("biotype") != "protein_coding":
            continue
        sym = _sym(gi, g)
        if not sym or sym not in complete or not is_fragment_model(g, sym, ref_struct):
            continue
        calls[gi] = "fragment"
        if gi in assign:
            del assign[gi]
        if g["source_name"]:
            present[g["source_name"]] -= 1
            g["source_name"] = ""
        n += 1
    return n


def load_fasta_chroms(path: str | None, chroms: set[str]) -> dict[str, str]:
    if not path or not Path(path).exists() or not chroms:
        return {}
    wanted = set(chroms)
    out: dict[str, str] = {}
    cur, keep, buf = None, False, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and cur:
                    out[cur] = "".join(buf)
                cur = line[1:].split()[0]
                keep = cur in wanted
                buf = []
            elif keep:
                buf.append(line.strip())
        if keep and cur:
            out[cur] = "".join(buf)
    return out


def load_proteins(path: str | None, wanted: set[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path or not Path(path).exists():
        return out
    cur, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(buf)
                tid = line[1:].split()[0]
                cur = tid if tid in wanted else None
                buf = []
            elif cur is not None:
                buf.append(line.strip())
    if cur is not None:
        out[cur] = "".join(buf)
    return out


class Hit:
    __slots__ = ("chrom", "start", "end", "strand", "query", "symbol", "pid", "qcov", "qlen", "n_exons")

    def __init__(self, chrom, start, end, strand, query, symbol, pid, qcov, qlen=0, n_exons=0):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.query = query
        self.symbol = symbol
        self.pid = float(pid)
        self.qcov = float(qcov)
        self.qlen = int(qlen)
        self.n_exons = int(n_exons)

    @property
    def score(self) -> float:
        return self.pid * self.qcov

    @property
    def span(self) -> tuple[int, int]:
        return (self.start, self.end)

    @property
    def length(self) -> int:
        return self.end - self.start


def iter_paf_hits(paf_path: str, query_symbols: dict[str, str], with_exons: bool = False):
    """Yield Hit (and optional exon list) for aligned proteins that pass floors."""
    with open(paf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            rec = parse_paf_row(line)
            if rec is None:
                continue
            qlen = rec["q_len"] or 0
            if qlen <= 0:
                continue
            cg = rec["tags"].get("cg", "")
            if not cg:
                continue
            exons, _matched_aa, aligned_aa, _fs = parse_cigar_to_exons(
                cg, rec["t_start"], rec["t_end"], rec["strand"]
            )
            if not exons:
                continue
            qcov = aligned_aa / qlen
            pid = paf_identity(rec)
            if qcov < MIN_HIT_QCOV or pid < MIN_HIT_PID:
                continue
            sym = symbol_for_query(rec["q_name"], query_symbols)
            if not sym:
                continue
            h = Hit(
                rec["t_name"], rec["t_start"], rec["t_end"], rec["strand"],
                rec["q_name"], sym, pid, qcov, qlen, real_exon_count(exons),
            )
            if with_exons:
                yield h, exons
            else:
                yield h


def load_compact_hits(paf_path: str, query_symbols: dict[str, str]) -> list[Hit]:
    hits = list(iter_paf_hits(paf_path, query_symbols, with_exons=False))
    logger.info("  miniprot hits kept: %s", f"{len(hits):,}")
    return hits


def fetch_exons(
    paf_path: str,
    query_symbols: dict[str, str],
    needed: set[tuple],
) -> dict[tuple, list[tuple[int, int]]]:
    """needed keys: (chrom, start, end, query)."""
    out: dict[tuple, list[tuple[int, int]]] = {}
    if not needed:
        return out
    for h, exons in iter_paf_hits(paf_path, query_symbols, with_exons=True):
        key = (h.chrom, h.start, h.end, h.query)
        if key in needed and key not in out:
            out[key] = exons
            if len(out) == len(needed):
                break
    return out


def load_genes(gp_path: str, gp_info_path: str) -> tuple[list[dict], pd.DataFrame]:
    info = pd.read_csv(gp_info_path, sep="\t", low_memory=False)
    info_by_tx = {}
    if "transcript_id" in info.columns:
        for rec in info.to_dict("records"):
            info_by_tx[str(rec["transcript_id"])] = rec
    genes: list[dict] = []
    for tx in gene_pred_iterator(gp_path):
        rec = info_by_tx.get(tx.name, {})
        genes.append(
            {
                "tx": tx,
                "info": rec,
                "gene_id": str(rec.get("gene_id") or tx.name2 or tx.name),
                "transcript_id": tx.name,
                "chrom": tx.chromosome,
                "start": tx.start,
                "end": tx.stop,
                "strand": tx.strand,
                "biotype": str(rec.get("gene_biotype") or rec.get("transcript_biotype") or ""),
                "source_name": "" if _na(rec.get("source_gene_common_name")) else str(rec["source_gene_common_name"]),
                "transcript_class": str(rec.get("transcript_class") or ""),
                "alignment_mode": str(rec.get("alignment_mode") or ""),
                "cds_size": tx.cds_size,
                "coding_exons": coding_exon_count(tx),
            }
        )
    return genes, info


def alignment_mode_of(g: dict) -> str:
    return str(g.get("alignment_mode") or (g.get("info") or {}).get("alignment_mode") or "")


def is_projection_mode(g: dict) -> bool:
    mode = alignment_mode_of(g)
    return mode in PROJECTION_MODES or mode.startswith("transMap") or mode.startswith("augTM")


def is_protein_mode(g: dict) -> bool:
    mode = alignment_mode_of(g)
    return mode in PROTEIN_MODES or mode.startswith("augMP") or mode == "miniprot"


def pick_representatives(
    genes: list[dict],
    proteins: dict[str, str] | None = None,
    typical_aa: dict[str, int] | None = None,
    seleno_symbols: set[str] | None = None,
) -> list[int]:
    """Best transcript index per gene_id.

    A clean peptide beats a broken one (alignment frameshift flags are
    ignored). Among clean isoforms, a transMap/augTM model wins by longest
    peptide — those are spliced projections of the real locus. miniprot /
    augMP (and Augustus/StringTie extras) still pick the peptide nearest the
    UniProt typical length, because those modes invent compact or extended
    ORFs that are not the parent protein.
    """
    proteins = proteins or {}
    typical_aa = typical_aa or {}

    def peptide_len(g: dict, aa: str | None) -> int:
        if aa:
            seq = aa[:-1] if aa.endswith("*") else aa
            if seq and "*" not in seq:
                return len(seq)
        return int(g.get("cds_size") or 0) // 3

    def rank(g: dict) -> tuple:
        aa = proteins.get(g["transcript_id"])
        typ = typical_aa.get(g.get("source_name") or "") or 0
        plen = peptide_len(g, aa)
        clean = peptide_is_clean(aa, seleno=(g.get("source_name") or "") in (seleno_symbols or ()))
        proj = is_projection_mode(g)
        near_typical = -abs(plen - typ) if typ else 0
        return (
            not clearly_broken_peptide(aa),
            clean,
            clean and proj,
            plen if clean and proj else 0,
            near_typical if not proj else 0,
            near_typical,
            g["cds_size"],
            g["end"] - g["start"],
        )

    best: dict[str, int] = {}
    for i, g in enumerate(genes):
        prev = best.get(g["gene_id"])
        if prev is None or rank(g) > rank(genes[prev]):
            best[g["gene_id"]] = i
    return sorted(best.values(), key=lambda i: (genes[i]["chrom"], genes[i]["start"]))


def hits_by_chrom(hits: list[Hit]) -> dict[str, list[Hit]]:
    by = defaultdict(list)
    for h in hits:
        by[h.chrom].append(h)
    return by


class LocusIndex:
    """Hits on one chromosome, queried by overlap with a gene span."""

    def __init__(self, hits: list[Hit]):
        self.hits = sorted(hits, key=lambda h: h.start)
        self.starts = [h.start for h in self.hits]
        self.maxlen = max((h.length for h in self.hits), default=0)

    def overlapping(self, span: tuple[int, int], min_ov: int) -> list[Hit]:
        lo = bisect_left(self.starts, span[0] - self.maxlen)
        hi = bisect_right(self.starts, span[1])
        return [h for h in self.hits[lo:hi] if overlap_bp(span, h.span) >= min_ov]


def best_by_symbol(hits: list[Hit], ref_struct: dict[str, tuple[int, int]]) -> dict[str, Hit]:
    """Most complete hit per symbol; a short fragment query never beats the full protein."""
    by: dict[str, tuple[float, float, Hit]] = {}
    for h in hits:
        key = (completeness(h, ref_struct), h.pid)
        prev = by.get(h.symbol)
        if prev is None or key > prev[:2]:
            by[h.symbol] = (*key, h)
    return {s: v[2] for s, v in by.items()}


def orf_intact(g: dict, protein: str | None) -> bool:
    """True when the model CDS is a single clean ORF.

    A trailing ``*`` is optional: consensus FASTA for Augustus often omits the
    stop, which is why intact TSPY copies were labeled ``broken_orf``. Internal
    stops, a missing start, a frameshift flag, or a CDS not in frame still fail.
    """
    if g["cds_size"] <= 0 or g["cds_size"] % 3:
        return False
    if str(g.get("info", {}).get("frameshift", "")).lower() in {"true", "1"}:
        return False
    if protein is None:
        return True
    seq = protein[:-1] if protein.endswith("*") else protein
    seq = seq.replace("U", "A")
    if not seq or not seq.startswith("M"):
        return False
    return "*" not in seq


def peptide_is_clean(protein: str | None, *, seleno: bool = False) -> bool:
    """True when the translated peptide is one start-to-stop ORF.

    Alignment ``frameshift`` flags are ignored. ``U`` is Sec and is only
    accepted (along with leftover ``*`` below the broken floor) when this
    gene's parent protein is a selenoprotein. A premature TGA in anything
    else stays a stop and is not clean.
    """
    if not protein:
        return False
    seq = protein[:-1] if protein.endswith("*") else protein
    if not seq or not seq.startswith("M"):
        return False
    seq = seq.replace("U", "A")
    if seleno:
        return not clearly_broken_peptide(seq)
    return "*" not in seq


def peptide_aa_len(protein: str | None, cds_size: int, *, seleno: bool = False) -> int:
    """Usable peptide length. Sec counts as a residue only for selenoproteins."""
    if protein:
        seq = protein[:-1] if protein.endswith("*") else protein
        if not seq:
            return cds_size // 3 if cds_size else 0
        if peptide_is_clean(protein, seleno=seleno):
            return len(seq)
        parts = [p for p in seq.split("*") if p]
        return max((len(p) for p in parts), default=0)
    return cds_size // 3 if cds_size else 0


def _near_typical(aa: int, typ_aa: int) -> bool:
    return bool(typ_aa) and MIN_COPY_AA_FRAC * typ_aa <= aa <= MAX_COPY_AA_FRAC * typ_aa


def cds_needs_replace(
    g: dict,
    protein: str | None,
    symbol: str,
    ref_struct: dict[str, tuple[int, int]],
    old_span: int,
    new_span: int,
    new_aa: int | None = None,
    seleno_symbols: set[str] | None = None,
) -> str | None:
    """Whether miniprot exons should replace the current CDS.

    Genomic span alone is not a fusion signal. JAK1's real first intron is
    ~200 kb, so a compact miniprot block looks like a shrink and destroys a
    complete transMap protein. A clean peptide near the parent length is
    never replaced, including when an alignment frameshift flag is set.
    A clean peptide may change only when it is the wrong length versus the
    parent *and* the new CDS is closer to that length (fusion shrink or
    truncation expand). Broken peptides may still be replaced.
    """
    if new_span + 50 >= old_span and new_span - 50 <= old_span:
        return None
    typ_aa, _ = parent_structure(symbol, ref_struct)
    seleno = symbol in (seleno_symbols or ())
    old_aa = peptide_aa_len(protein, int(g.get("cds_size") or 0), seleno=seleno)
    if peptide_is_clean(protein, seleno=seleno):
        # A clean transMap/augTM peptide is the spliced locus. Miniprot's
        # compact block is often a short isoform of the same protein and
        # must not replace it, even when UniProt "typical" is that stub.
        if is_projection_mode(g) and new_span < old_span:
            return None
        if not typ_aa:
            return None
        if _near_typical(old_aa, typ_aa):
            return None
        if new_aa is not None:
            if abs(new_aa - typ_aa) >= abs(old_aa - typ_aa):
                return None
            if new_span < old_span and old_aa <= MAX_COPY_AA_FRAC * typ_aa:
                return None
            if new_span > old_span and old_aa >= MIN_COPY_AA_FRAC * typ_aa:
                return None
            if not _near_typical(new_aa, typ_aa):
                return None
            return "shrink" if new_span < old_span else "expand"
        if old_aa > MAX_COPY_AA_FRAC * typ_aa and new_span < old_span:
            return "shrink"
        if old_aa < MIN_COPY_AA_FRAC * typ_aa and new_span > old_span:
            return "expand"
        return None
    return "shrink" if new_span < old_span else "expand"


def isoform_worth_keeping(
    g: dict,
    protein: str | None,
    symbol: str,
    ref_struct: dict[str, tuple[int, int]],
) -> bool:
    """True if this isoform should survive after its representative is dropped.

    A superseded work transcript can be a wrecked or chimeric parse. Sibling
    isoforms of the same gene_id are kept when they still translate as a
    non-fusion protein. Fusion-length siblings of a split chimera are not.
    """
    if clearly_broken_peptide(protein):
        return False
    aa = peptide_aa_len(protein, int(g.get("cds_size") or 0))
    if aa < 50:
        return False
    typ_aa, _ = parent_structure(symbol, ref_struct)
    if typ_aa and aa > MAX_COPY_AA_FRAC * typ_aa:
        return False
    if protein and not peptide_is_clean(protein) and not clearly_broken_peptide(protein):
        return aa >= 50
    return True


def intact_sibling(
    work_gene: dict,
    siblings: list[dict],
    proteins: dict[str, str],
    symbol: str,
    ref_struct: dict[str, tuple[int, int]],
) -> dict | None:
    """Best clean non-fusion isoform of the same gene, excluding ``work_gene``."""
    typ_aa, _ = parent_structure(symbol, ref_struct)
    best = None
    best_key = None
    work_tid = work_gene.get("transcript_id")
    for g in siblings:
        if g.get("transcript_id") == work_tid:
            continue
        if g.get("gene_id") != work_gene.get("gene_id"):
            continue
        prot = proteins.get(g["transcript_id"])
        if not isoform_worth_keeping(g, prot, symbol, ref_struct):
            continue
        if prot and not peptide_is_clean(prot):
            continue
        plen = peptide_aa_len(prot, int(g.get("cds_size") or 0))
        key = (-abs(plen - typ_aa) if typ_aa else 0, plen)
        if best is None or key > best_key:
            best, best_key = g, key
    return best


def prune_superseded_isoforms(
    keep_tx: set[str],
    work: list[dict],
    superseded: set[int],
    all_genes: list[dict],
    proteins: dict[str, str],
    ref_struct: dict[str, tuple[int, int]],
) -> tuple[int, int]:
    """Drop superseded representatives; keep intact non-fusion siblings.

    Replacing a chimeric or still-broken work transcript must not delete the
    rest of the gene. Woolly JAK1 lost 29 intact 1154 aa isoforms because the
    representative was superseded and every ``gene_id`` member was discarded.
    Returns ``(n_dropped_siblings, n_kept_siblings)``.
    """
    by_gid: dict[str, list[dict]] = defaultdict(list)
    for g in all_genes:
        by_gid[g["gene_id"]].append(g)
    n_drop = n_keep = 0
    for wi in superseded:
        g = work[wi]
        keep_tx.discard(g["transcript_id"])
        sym = g.get("source_name") or ""
        for sib in by_gid[g["gene_id"]]:
            if sib["transcript_id"] == g["transcript_id"]:
                continue
            prot = proteins.get(sib["transcript_id"])
            if isoform_worth_keeping(sib, prot, sym or sib.get("source_name") or "", ref_struct):
                n_keep += 1
                continue
            keep_tx.discard(sib["transcript_id"])
            n_drop += 1
    return n_drop, n_keep


def coding_call(
    g: dict,
    by_sym: dict[str, Hit],
    ref_struct: dict[str, tuple[int, int]],
    protein: str | None,
) -> tuple[str, Hit | None]:
    """Decide whether an unnamed coding model is a real protein-coding copy.

    Returns (call, best_hit). Only ``real_copy`` gets a name; calls in
    DEMOTE_CALLS lose protein_coding; ``no_homology``, ``noncoding_parent``,
    ``low_identity`` and ``ambiguous`` are left as they are. A ``real_copy``
    overlapping a locus that already carries its symbol later becomes
    ``redundant`` in ``assign_names``.

    Intron loss is judged on the gene model. A thin miniprot hit (fragment
    peptide or a distant alignment that collapses exons) must not demote a
    model that itself has parent-like intron structure.
    """
    if not by_sym:
        return "no_homology", None
    ranked = sorted(
        by_sym.values(), key=lambda h: (completeness(h, ref_struct), h.pid), reverse=True
    )
    best = ranked[0]
    p_aa, p_ex = parent_structure(best.symbol, ref_struct)
    if not p_aa:
        return "noncoding_parent", best
    if not orf_intact(g, protein):
        # Augustus often drops a start or slips frame on tandem arrays while
        # miniprot still has a complete protein at the same locus (TSPY).
        # Require parent-like exons so a 1-exon GAPDH retrocopy stays broken.
        covered = int(best.qcov * best.qlen) if best.qlen else 0
        if extra_copy_complete(
            best.n_exons, covered, best.symbol, ref_struct, best.pid, best.qcov
        ):
            if len(ranked) > 1 and (
                completeness(best, ref_struct) - completeness(ranked[1], ref_struct) < NAME_MARGIN
            ):
                return "ambiguous", best
            return "real_copy", best
        return "broken_orf", best
    model_aa = g["cds_size"] // 3
    covered = best.qcov * best.qlen if best.qlen else best.qcov * p_aa
    if (
        best.qcov < MIN_COPY_QCOV
        or covered < MIN_COPY_AA_FRAC * p_aa
        or model_aa < MIN_COPY_AA_FRAC * p_aa
    ):
        return "fragment", best
    if model_aa > MAX_COPY_AA_FRAC * p_aa:
        return "chimeric", best
    min_ex = max(2, math.ceil(MIN_COPY_EXON_FRAC * p_ex))
    if p_ex >= 2 and g["coding_exons"] < min_ex:
        return "retrocopy", best
    if best.pid < MIN_COPY_PID:
        return "low_identity", best
    if len(ranked) > 1 and (
        completeness(best, ref_struct) - completeness(ranked[1], ref_struct) < NAME_MARGIN
    ):
        return "ambiguous", best
    return "real_copy", best


def _assignment(symbol: str, hit: Hit, label: str) -> dict:
    return {"symbol": symbol, "pid": hit.pid, "qcov": hit.qcov, "hit": hit, "label": label}


def should_correct_name(
    own: Hit | None, rival: Hit, gene: dict, ref_struct: dict[str, tuple[int, int]]
) -> bool:
    """Overwrite a CAT name only when the model is the rival gene, not the CAT gene.

    A complete own hit is never replaced (GDI1, FCGBP). A rival whose parent
    length does not fit the model (RPS6KA2 on a 1-exon RAMACL; ZNF519 on a
    20-exon FCGBP) is not a correction.
    """
    if rival.pid < MIN_PRIMARY_PID:
        return False
    own_pid = own.pid if own else 0.0
    if rival.pid < own_pid + CAT_OVERRIDE_MARGIN:
        return False
    if own and is_complete_copy(
        own.n_exons, int(own.qcov * own.qlen), own.symbol, ref_struct, require_exons=False
    ):
        return False
    p_aa = parent_structure(rival.symbol, ref_struct)[0]
    model_aa = int(gene.get("cds_size") or 0) // 3
    if not p_aa or not model_aa:
        return False
    if not MIN_COPY_AA_FRAC * p_aa <= model_aa <= MAX_COPY_AA_FRAC * p_aa:
        return False
    return is_complete_copy(
        rival.n_exons, int(rival.qcov * rival.qlen), rival.symbol, ref_struct, require_exons=False
    )


def assign_names(
    genes: list[dict],
    index_by_chrom: dict[str, LocusIndex],
    ref_symbols: set[str],
    ref_struct: dict[str, tuple[int, int]],
    proteins: dict[str, str],
) -> tuple[dict[int, dict], dict[int, str], Counter]:
    """Genome-wide naming: confirm CAT names, gate unnamed models.

    Returns (assign, calls, present). ``present`` counts loci per symbol after
    naming; corrections never target a symbol already present.
    """
    present: Counter = Counter(g["source_name"] for g in genes if g["source_name"])
    assign: dict[int, dict] = {}
    calls: dict[int, str] = {}
    corrections = []
    copies = []
    ambiguous: list[int] = []

    for gi, g in enumerate(genes):
        idx = index_by_chrom.get(g["chrom"])
        local = idx.overlapping((int(g["start"]), int(g["end"])), 100) if idx else []
        local = [
            h for h in local
            if h.symbol in ref_symbols and (not g["strand"] or h.strand == g["strand"])
        ]
        by_sym = best_by_symbol(local, ref_struct)
        src = g["source_name"]
        if src:
            own = by_sym.get(src)
            rivals = [
                h for s, h in by_sym.items()
                if s != src and s in ref_struct and h.qcov >= MIN_PRIMARY_QCOV
            ]
            rival = max(rivals, key=lambda h: (h.pid, h.qcov), default=None)
            parts = readthrough_parts(rival.symbol, ref_symbols) if rival else None
            if rival and parts and src in parts:
                if own and own.pid >= MIN_NAMED_PID:
                    assign[gi] = _assignment(src, own, "confirmed")
                continue
            if rival and should_correct_name(own, rival, g, ref_struct):
                corrections.append((rival.pid, rival.qcov, gi, rival, own))
            elif own and own.pid >= MIN_NAMED_PID:
                assign[gi] = _assignment(src, own, "confirmed")
            continue
        remapped = {}
        for h in by_sym.values():
            rh = remap_readthrough_hit(h, ref_symbols, ref_struct, present)
            remapped[rh.symbol] = rh
        call, hit = coding_call(g, remapped or by_sym, ref_struct, proteins.get(g["transcript_id"]))
        calls[gi] = call
        if call == "real_copy":
            copies.append((completeness(hit, ref_struct), hit.pid, gi, hit))
        elif call == "ambiguous":
            ambiguous.append(gi)

    corrections.sort(key=lambda x: (x[0], x[1]), reverse=True)
    for _pid, _qcov, gi, rival, own in corrections:
        src = genes[gi]["source_name"]
        if present[rival.symbol] > 0:
            if own and own.pid >= MIN_NAMED_PID:
                assign[gi] = _assignment(src, own, "confirmed")
            continue
        assign[gi] = _assignment(rival.symbol, rival, "corrected")
        present[rival.symbol] += 1
        present[src] -= 1

    named_spans: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for gi, g in enumerate(genes):
        sym = assign[gi]["symbol"] if gi in assign else g["source_name"]
        if sym:
            named_spans[(sym, g["chrom"])].append((int(g["start"]), int(g["end"])))

    for gi in ambiguous:
        g = genes[gi]
        if not orf_intact(g, proteins.get(g["transcript_id"])):
            continue
        idx = index_by_chrom.get(g["chrom"])
        local = idx.overlapping((int(g["start"]), int(g["end"])), 100) if idx else []
        local = [
            remap_readthrough_hit(h, ref_symbols, ref_struct, present)
            for h in local
            if h.symbol in ref_symbols and (not g["strand"] or h.strand == g["strand"])
        ]
        picked = resolve_tied_hit(g["chrom"], int(g["start"]), local, named_spans, ref_struct)
        if picked is None:
            continue
        calls[gi] = "real_copy"
        copies.append((completeness(picked, ref_struct), picked.pid, gi, picked))

    copies.sort(key=lambda x: (x[0], x[1]), reverse=True)
    for _sc, _pid, gi, hit in copies:
        g = genes[gi]
        span = (int(g["start"]), int(g["end"]))
        hit = remap_readthrough_hit(hit, ref_symbols, ref_struct, present)
        if present.get(hit.symbol, 0) > 0:
            idx = index_by_chrom.get(g["chrom"])
            local = idx.overlapping(span, 100) if idx else []
            alt = None
            for h in local:
                h = remap_readthrough_hit(h, ref_symbols, ref_struct, present)
                if present.get(h.symbol, 0) > 0:
                    continue
                if h.symbol not in ref_symbols:
                    continue
                if not is_complete_copy(h.n_exons, int(h.qcov * h.qlen), h.symbol, ref_struct):
                    continue
                alt = h
                break
            if alt:
                hit = alt
            elif any(overlap_bp(span, s) > 0 for s in named_spans[(hit.symbol, g["chrom"])]):
                calls[gi] = "redundant"
                continue
            else:
                near = nearest_named_symbol(g["chrom"], int(g["start"]), named_spans)
                if near and near != hit.symbol:
                    for h in local:
                        if h.symbol != near:
                            continue
                        if extra_copy_complete(
                            h.n_exons, int(h.qcov * h.qlen), h.symbol, ref_struct, h.pid, h.qcov
                        ):
                            hit = h
                            break
                if not extra_copy_complete(
                    hit.n_exons, int(hit.qcov * hit.qlen), hit.symbol, ref_struct, hit.pid, hit.qcov
                ) or not extra_copy_locus_ok(g["chrom"], hit.symbol, named_spans):
                    calls[gi] = "retrocopy"
                    continue
        key = (hit.symbol, g["chrom"])
        if any(overlap_bp(span, s) > 0 for s in named_spans[key]):
            calls[gi] = "redundant"
            continue
        named_spans[key].append(span)
        label = "named_copy" if present[hit.symbol] > 0 else "named_missing"
        assign[gi] = _assignment(hit.symbol, hit, label)
        present[hit.symbol] += 1

    return assign, calls, present


def two_complete_hits_inside(
    g: dict,
    idx: LocusIndex | None,
    ref_struct: dict[str, tuple[int, int]],
    ref_symbols: set[str],
) -> list[Hit]:
    """Non-overlapping complete protein hits whose CDS sits inside ``g``.

    Two such hits mean ``g`` is a chimera replacing two genes.
    """
    if idx is None:
        return []
    span = (int(g["start"]), int(g["end"]))
    segs = _coding_segments(g["tx"])
    local = []
    for h in idx.overlapping(span, 100):
        if h.symbol not in ref_symbols:
            continue
        if g["strand"] and h.strand != g["strand"]:
            continue
        if h.start < span[0] - 50 or h.end > span[1] + 50:
            continue
        if not is_complete_copy(
            h.n_exons, int(h.qcov * h.qlen), h.symbol, ref_struct, require_exons=False
        ):
            continue
        ov = sum(max(0, min(e, h.end) - max(s, h.start)) for s, e in segs)
        if ov < 0.30 * max(1, h.length):
            continue
        local.append(h)
    by_sym = best_by_symbol(local, ref_struct)
    hits = sorted(by_sym.values(), key=lambda h: h.start)
    picked: list[Hit] = []
    for h in hits:
        if any(overlap_bp(h.span, p.span) > 50 for p in picked):
            continue
        picked.append(h)
    return picked if len(picked) >= 2 else []


def split_chimeras(
    genes: list[dict],
    calls: dict[int, str],
    assign: dict[int, dict],
    present: Counter,
    index_by_chrom: dict[str, LocusIndex],
    ref_struct: dict[str, tuple[int, int]],
    ref_symbols: set[str],
) -> tuple[list[Hit], set[int]]:
    """Turn a fusion that covers two complete genes into two rescue candidates.

    The chimeric consensus model is dropped (superseded). Each contained
    complete hit is offered to ``place_rescues`` if that symbol is still absent.
    """
    extra: list[Hit] = []
    split: set[int] = set()
    for gi, g in enumerate(genes):
        oversized = (
            calls.get(gi) == "chimeric"
            or (g["cds_size"] // 3) >= 800 and gi not in assign
        )
        pair = two_complete_hits_inside(
            g, index_by_chrom.get(g["chrom"]), ref_struct, ref_symbols
        )
        if len(pair) < 2:
            continue
        if gi in assign:
            own = assign[gi]["symbol"]
            own_hit = assign[gi].get("hit")
            p_aa = parent_structure(own, ref_struct)[0]
            model_aa = g["cds_size"] // 3
            own_complete = False
            if own_hit is not None:
                own_complete = is_complete_copy(
                    own_hit.n_exons, int(own_hit.qcov * own_hit.qlen), own, ref_struct
                )
            gene_span = int(g["end"]) - int(g["start"])
            unit = max((h.length for h in pair), default=0)
            tandem = unit > 0 and gene_span > 2.5 * unit
            # A named CAT ortholog that matches its own parent is not a fusion
            # unless transMap glued several tandem units into one 80 kb model.
            if (
                g["source_name"] == own
                and own_complete
                and p_aa
                and model_aa <= MAX_COPY_AA_FRAC * p_aa
                and not tandem
            ):
                continue
            if not oversized and not any(h.symbol != own for h in pair) and not tandem:
                continue
        if gi in assign:
            present[assign[gi]["symbol"]] -= 1
            del assign[gi]
        if g["source_name"]:
            g["source_name"] = ""
        calls[gi] = "chimeric_split"
        split.add(gi)
        extra.extend(pair)
    return extra, split


def _is_rescue_complete(h: Hit, ref_struct: dict[str, tuple[int, int]]) -> bool:
    p_aa, _p_ex = parent_structure(h.symbol, ref_struct)
    if not p_aa or h.pid < MIN_COPY_PID or h.qcov < MIN_COPY_QCOV:
        return False
    covered = h.qcov * h.qlen
    if not MIN_COPY_AA_FRAC * p_aa <= covered <= MAX_COPY_AA_FRAC * p_aa:
        return False
    return is_complete_copy(
        h.n_exons, int(covered), h.symbol, ref_struct, require_exons=False
    )


def extra_copy_candidates(
    hits: list[Hit],
    present: Counter,
    named_spans: dict[tuple[str, str], list[tuple[int, int]]],
    ref_struct: dict[str, tuple[int, int]],
    ref_symbols: set[str] | None = None,
) -> list[Hit]:
    """Complete proteins at extra loci of genes already named elsewhere.

    transMap can only project GRCh38 copy number. T2T arrays (TSPY, RBMY, …)
    have additional intact copies whose Augustus models are often broken_orf.
    Same symbol may be placed at many non-overlapping loci on a chromosome
    that already carries that symbol. When two paralogs both fit (AMY1 vs
    AMY2), the copy is named for the nearest already-named gene. Processed
    pseudogenes (complete peptide, lost introns, other chromosomes) are not
    extras.
    """
    ref_symbols = ref_symbols or set(ref_struct)
    complete: list[Hit] = []
    for h in hits:
        if h.symbol not in ref_symbols or present.get(h.symbol, 0) <= 0:
            continue
        if not extra_copy_complete(
            h.n_exons, int(h.qcov * h.qlen), h.symbol, ref_struct, h.pid, h.qcov
        ):
            continue
        if not extra_copy_locus_ok(h.chrom, h.symbol, named_spans):
            continue
        already = False
        for (sym, chrom), spans in named_spans.items():
            if chrom != h.chrom:
                continue
            if any(overlap_bp(h.span, s) > 50 for s in spans):
                already = True
                break
        if already:
            continue
        complete.append(h)
    complete.sort(key=lambda h: (h.chrom, h.strand, h.start, -h.pid))
    clusters: list[list[Hit]] = []
    for h in complete:
        for c in clusters:
            if c[0].chrom == h.chrom and c[0].strand == h.strand and any(
                hits_conflict(h, x) for x in c
            ):
                c.append(h)
                break
        else:
            clusters.append([h])
    out = [pick_extra_from_cluster(c, named_spans, ref_struct) for c in clusters]
    out.sort(key=lambda h: (completeness(h, ref_struct), h.pid, -h.length), reverse=True)
    return out


def rescue_candidates(
    hits: list[Hit],
    present: Counter,
    ref_struct: dict[str, tuple[int, int]],
    ref_symbols: set[str] | None = None,
) -> list[Hit]:
    """Complete, intron-consistent hits to reference genes absent from the genome.

    Up to RESCUE_CANDIDATES_PER_SYMBOL per symbol, most complete first. Symbols
    that share a locus are resolved so each missing gene keeps its own interval
    when an alternative complete hit exists.
    """
    ref_symbols = ref_symbols or set(ref_struct)
    by_sym: dict[str, list[Hit]] = defaultdict(list)
    for h in hits:
        h = remap_readthrough_hit(h, ref_symbols, ref_struct, present)
        if present.get(h.symbol, 0) > 0:
            continue
        if not _is_rescue_complete(h, ref_struct):
            continue
        by_sym[h.symbol].append(h)
    for sym, hs in by_sym.items():
        hs.sort(key=lambda h: (completeness(h, ref_struct), h.pid, -h.length), reverse=True)
        by_sym[sym] = hs[:RESCUE_CANDIDATES_PER_SYMBOL]
    picked = pick_unique_hits(by_sym, ref_struct)
    # Keep the chosen hit plus unused alternates for place_rescues fallback.
    chosen_keys = {(h.symbol, h.chrom, h.start, h.end) for h in picked}
    out = list(picked)
    for hs in by_sym.values():
        for h in hs:
            key = (h.symbol, h.chrom, h.start, h.end)
            if key not in chosen_keys:
                out.append(h)
                chosen_keys.add(key)
    out.sort(key=lambda h: (h.symbol in {p.symbol for p in picked}, completeness(h, ref_struct), h.pid), reverse=True)
    return out


def _coding_segments(tx) -> list[tuple[int, int]]:
    segs = []
    for e in tx.exon_intervals:
        s, t = max(e.start, tx.thick_start), min(e.stop, tx.thick_stop)
        if t > s:
            segs.append((s, t))
    return segs


def _segment_overlap(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    return sum(max(0, min(x1, y1) - max(x0, y0)) for x0, x1 in a for y0, y1 in b)


def rescue_blocked_by_named(
    rescue_segs: list[tuple[int, int]],
    rescue_cds: int,
    named_segs: list[tuple[int, int]],
    named_cds: int,
) -> bool:
    """Block only when the rescue is the same gene or a fragment of the named gene.

    A small complete gene on the 5' end of a longer complete gene (EEF1AKMT4
    ⊂ ECE2) does not block the longer rescue. A small complete gene swallowed
    inside a chimeric larger model (EBP ⊂ TBC1D25 that starts on the EBP exon)
    is also allowed: overlap covers the rescue but is only a fraction of the
    named CDS. Same-size / same-interval overlap still blocks.
    """
    ov = _segment_overlap(rescue_segs, named_segs)
    if ov <= 0 or rescue_cds <= 0:
        return ov > 0
    if ov > 0.50 * rescue_cds:
        if named_cds >= 1.5 * rescue_cds and ov < 0.50 * named_cds:
            return False
        return True
    if named_cds > 0 and ov > 0.50 * named_cds and rescue_cds <= MAX_COPY_AA_FRAC * named_cds:
        return True
    return False


def _named_protein_coding(g: dict, gi: int, assign: dict[int, dict]) -> bool:
    if gi in assign:
        return True
    return bool(g["source_name"] and g["biotype"] == "protein_coding")


def nearest_named_gene(
    chrom: str, pos: int, genes: list[dict], assign: dict[int, dict]
) -> dict | None:
    best: tuple[int, dict] | None = None
    for gi, g in enumerate(genes):
        if g.get("chrom") != chrom or g.get("tx") is None:
            continue
        if not _named_protein_coding(g, gi, assign):
            continue
        d = _span_dist(pos, (int(g["start"]), int(g["end"])))
        if best is None or d < best[0]:
            best = (d, g)
    return best[1] if best else None


def _best_exon_window(
    query: str, chrom_seq: str, lo: int, hi: int, min_ident: float = 0.82
) -> tuple[int, int] | None:
    q = query.upper()
    n = len(q)
    if n < 12 or hi - lo < n:
        return None
    best: tuple[float, int, int] | None = None
    for i in range(lo, hi - n + 1):
        t = chrom_seq[i : i + n].upper()
        if len(t) < n:
            break
        ident = sum(a == b for a, b in zip(q, t)) / n
        if best is None or ident > best[0]:
            best = (ident, i, i + n)
            if ident >= 0.98:
                break
    if best and best[0] >= min_ident:
        return best[1], best[2]
    return None


def project_template_exons(
    template,
    chrom: str,
    strand: str,
    window: tuple[int, int],
    seq_dict: dict[str, str],
) -> list[tuple[int, int]] | None:
    """Copy the nearest named gene's coding exons onto a tandem extra locus.

    Miniprot often keeps a 10-exon AMY parse that is stop-ridden on an
    inverted unit; the named neighbour's 11–12 exon structure still translates.
    """
    tmpl_chrom = getattr(template, "chromosome", chrom)
    if chrom not in seq_dict or tmpl_chrom not in seq_dict:
        return None
    seq = seq_dict[chrom]
    tmpl_seq = seq_dict[tmpl_chrom]
    if not seq or not tmpl_seq:
        return None
    segs = []
    for e in template.exon_intervals:
        s, t = max(e.start, template.thick_start), min(e.stop, template.thick_stop)
        if t - s >= 12:
            segs.append((s, t))
    if len(segs) < 3:
        return None
    if strand != template.strand:
        segs = list(reversed(segs))
    lo, hi = max(0, window[0]), min(len(seq), window[1])
    found: list[tuple[int, int]] = []
    cursor = lo
    for s, t in segs:
        q = tmpl_seq[s:t]
        if strand != template.strand:
            q = reverse_complement(q)
        hit = _best_exon_window(q, seq, cursor, hi)
        if hit is None:
            continue
        found.append(hit)
        cursor = hit[1] + MIN_INTRON_BP
    if len(found) < max(3, math.ceil(MIN_COPY_EXON_FRAC * len(segs))):
        return None
    for (_a, b), (c, _d) in zip(found, found[1:]):
        if c < b:
            return None
    return found


def place_rescues(
    cands: list[Hit],
    exon_lookup: dict[tuple, list[tuple[int, int]]],
    genes: list[dict],
    assign: dict[int, dict],
    ref_struct: dict[str, tuple[int, int]],
    extra_cands: list[Hit] | None = None,
    seq_dict: dict[str, str] | None = None,
    seleno_symbols: set[str] | None = None,
) -> tuple[list[tuple[Hit, GenePredTranscript]], set[int]]:
    """Accept one rescue per absent symbol, then extra copies of present genes.

    A rescue is blocked when it is essentially the same CDS as a named
    protein-coding gene, not when it only shares a 5' exon with a neighbour
    or overlaps a projected GRCh38 pseudogene.
    """
    by_chrom: dict[str, list[int]] = defaultdict(list)
    for gi, g in enumerate(genes):
        by_chrom[g["chrom"]].append(gi)
    segs_cache: dict[int, list[tuple[int, int]]] = {}
    accepted: list[tuple[Hit, GenePredTranscript]] = []
    taken: dict[tuple[str, str], list[tuple[list[tuple[int, int]], int]]] = defaultdict(list)
    superseded: set[int] = set()
    used: set[str] = set()

    def _try_place(h: Hit, one_per_symbol: bool, extra: bool = False) -> None:
        if one_per_symbol and h.symbol in used:
            return
        exons = exon_lookup.get((h.chrom, h.start, h.end, h.query))
        if not exons:
            return
        tx = GenePredTranscript(
            gp_line_from_exons(f"MPR-{_safe(h.symbol)}-{h.start}", h.chrom, h.strand, exons, h.symbol).split("\t")
        )
        aa = tx.cds_size // 3
        n_ex = coding_exon_count(tx)
        if extra:
            if not extra_copy_complete(n_ex, aa, h.symbol, ref_struct, h.pid, h.qcov):
                return
            if seq_dict and tx.chromosome in seq_dict:
                try:
                    prot = protein_from_transcript(
                        tx, seq_dict, seleno=h.symbol in (seleno_symbols or ())
                    )
                except (KeyError, AssertionError, ValueError, RuntimeError, TypeError):
                    prot = None
                if prot is not None and not orf_intact(
                    {"cds_size": tx.cds_size, "info": {}}, prot
                ):
                    tmpl = nearest_named_gene(h.chrom, h.start, genes, assign)
                    rebuilt = None
                    if tmpl is not None:
                        tspan = int(tmpl["end"]) - int(tmpl["start"])
                        pad = max(h.length, tspan)
                        center = (h.start + h.end) // 2
                        new_exons = project_template_exons(
                            tmpl["tx"],
                            h.chrom,
                            h.strand,
                            (center - pad, center + pad),
                            seq_dict,
                        )
                        if new_exons:
                            rebuilt = GenePredTranscript(
                                gp_line_from_exons(
                                    f"MPR-{_safe(h.symbol)}-{h.start}",
                                    h.chrom,
                                    h.strand,
                                    new_exons,
                                    h.symbol,
                                ).split("\t")
                            )
                    if rebuilt is None:
                        return
                    try:
                        prot = protein_from_transcript(
                            rebuilt, seq_dict, seleno=h.symbol in (seleno_symbols or ())
                        )
                    except (KeyError, AssertionError, ValueError, RuntimeError, TypeError):
                        return
                    aa = rebuilt.cds_size // 3
                    n_ex = coding_exon_count(rebuilt)
                    if not extra_copy_complete(
                        n_ex, aa, h.symbol, ref_struct, h.pid, h.qcov
                    ) or not orf_intact({"cds_size": rebuilt.cds_size, "info": {}}, prot):
                        return
                    tx = rebuilt
        else:
            if not is_complete_copy(n_ex, aa, h.symbol, ref_struct, require_exons=False):
                return
            # A rescue is a new protein-coding gene, so its CDS must translate
            # without internal stops (Sec recoded). A missing start is allowed,
            # as in peptide_qc. Rejecting here lets the symbol's next candidate
            # hit be tried.
            if seq_dict and tx.chromosome in seq_dict:
                try:
                    prot = protein_from_transcript(
                        tx, seq_dict, seleno=h.symbol in (seleno_symbols or ())
                    )
                except (KeyError, AssertionError, ValueError, RuntimeError, TypeError):
                    prot = None
                if prot is not None and "*" in prot.rstrip("*"):
                    return
        segs = _coding_segments(tx)
        if any(
            rescue_blocked_by_named(segs, tx.cds_size, t_segs, t_cds)
            for t_segs, t_cds in taken[(h.chrom, h.strand)]
        ):
            return
        blocked = False
        unnamed = []
        for gi in by_chrom.get(h.chrom, []):
            g = genes[gi]
            if g["strand"] != h.strand or overlap_bp((int(g["start"]), int(g["end"])), h.span) == 0:
                continue
            if gi not in segs_cache:
                segs_cache[gi] = _coding_segments(g["tx"])
            ov = _segment_overlap(segs, segs_cache[gi])
            if ov == 0:
                continue
            if _named_protein_coding(g, gi, assign):
                named_hit = assign[gi].get("hit") if gi in assign else None
                if named_hit is not None and overlap_bp(named_hit.span, h.span) <= 50:
                    continue
                if rescue_blocked_by_named(segs, tx.cds_size, segs_cache[gi], int(g["cds_size"])):
                    blocked = True
                    break
            else:
                unnamed.append(gi)
        if blocked:
            return
        accepted.append((h, tx))
        superseded.update(unnamed)
        taken[(h.chrom, h.strand)].append((segs, tx.cds_size))
        used.add(h.symbol)

    for h in cands:
        _try_place(h, one_per_symbol=True, extra=False)
    for h in extra_cands or []:
        _try_place(h, one_per_symbol=False, extra=True)
    return accepted, superseded


def choose_protein_span(
    gene: dict, hit: Hit, idx: LocusIndex | None, ref_struct: dict[str, tuple[int, int]]
) -> Hit:
    """Owned full-protein hit: rank by completeness vs the reference gene, then span."""
    gs, ge = int(gene["start"]), int(gene["end"])
    nearby = idx.overlapping((gs - 40_000, ge + 40_000), 1) if idx else []
    cands = [h for h in nearby if h.symbol == hit.symbol and h.pid >= MIN_ORTHOLOG_PID]
    if not cands:
        cands = [hit]
    full = [h for h in cands if h.qcov >= FULL_PROTEIN_QCOV and h.length >= MIN_PROTEIN_SPAN]
    pool = full or cands
    return max(pool, key=lambda h: (completeness(h, ref_struct), h.length, h.pid))


def overlaps_other_gene(
    span: tuple[int, int],
    self_gi: int,
    genes: list[dict],
    same_chrom: list[int],
) -> bool:
    for gj in same_chrom:
        if gj == self_gi:
            continue
        g = genes[gj]
        if overlap_bp(span, (int(g["start"]), int(g["end"]))) >= 100:
            return True
    return False


def gp_line_from_exons(name: str, chrom: str, strand: str, exons: list[tuple[int, int]], name2: str) -> str:
    exons = sorted(exons)
    tx_start, tx_end = exons[0][0], exons[-1][1]
    frames = _exon_frames(exons, strand)
    return "\t".join(
        [
            name,
            chrom,
            strand,
            str(tx_start),
            str(tx_end),
            str(tx_start),
            str(tx_end),
            str(len(exons)),
            ",".join(str(e[0]) for e in exons) + ",",
            ",".join(str(e[1]) for e in exons) + ",",
            "0",
            name2,
            "cmpl",
            "cmpl",
            ",".join(str(fr) for fr in frames) + ",",
        ]
    )


def apply_structure(
    genes: list[dict],
    assign: dict[int, dict],
    index_by_chrom: dict[str, LocusIndex],
    exon_lookup: dict[tuple, list[tuple[int, int]]],
    ref_struct: dict[str, tuple[int, int]],
    proteins: dict[str, str] | None = None,
    seleno_symbols: set[str] | None = None,
) -> tuple[int, int]:
    proteins = proteins or {}
    n_shrink = n_expand = 0
    genes_by_chrom: dict[str, list[int]] = defaultdict(list)
    for gi, g in enumerate(genes):
        genes_by_chrom[g["chrom"]].append(gi)
    for gi, a in sorted(assign.items(), key=lambda kv: (genes[kv[0]]["chrom"], genes[kv[0]]["start"])):
        g = genes[gi]
        chosen = choose_protein_span(g, a["hit"], index_by_chrom.get(g["chrom"]), ref_struct)
        if chosen.qcov < FULL_PROTEIN_QCOV or chosen.length < MIN_PROTEIN_SPAN:
            continue
        if overlaps_other_gene(chosen.span, gi, genes, genes_by_chrom[g["chrom"]]):
            continue
        key = (chosen.chrom, chosen.start, chosen.end, chosen.query)
        exons = exon_lookup.get(key)
        if not exons:
            continue
        old_span = g["end"] - g["start"] + 1
        new_span = chosen.length
        line = gp_line_from_exons(
            g["transcript_id"], chosen.chrom, chosen.strand, exons, a["symbol"]
        )
        new_tx = GenePredTranscript(line.split("\t"))
        reason = cds_needs_replace(
            g,
            proteins.get(g["transcript_id"]),
            a["symbol"],
            ref_struct,
            old_span,
            new_span,
            new_aa=new_tx.cds_size // 3,
            seleno_symbols=seleno_symbols,
        )
        if reason is None:
            continue
        genes[gi]["tx"] = new_tx
        genes[gi]["start"] = new_tx.start
        genes[gi]["end"] = new_tx.stop
        genes[gi]["strand"] = new_tx.strand
        genes[gi]["cds_size"] = new_tx.cds_size
        a["hit"] = chosen
        if reason == "shrink":
            n_shrink += 1
        else:
            n_expand += 1
    return n_shrink, n_expand


def info_row_for_orphan(template: dict, gene_id: str, tx_id: str, hit: Hit, genome: str) -> dict:
    row = {k: ("N/A" if k not in ("gene_id", "transcript_id") else "") for k in template}
    row.update(template)
    row["gene_id"] = gene_id
    row["transcript_id"] = tx_id
    row["source_transcript"] = hit.query
    row["source_transcript_name"] = hit.query
    row["source_gene"] = hit.symbol
    row["source_gene_common_name"] = hit.symbol
    row["source_gene_biotype"] = "protein_coding"
    row["gene_biotype"] = "protein_coding"
    row["transcript_biotype"] = "protein_coding"
    row["alignment_id"] = f"miniprot-{hit.query}"
    row["alignment_mode"] = "miniprot"
    row["transcript_class"] = "ortholog"
    row["protein_only_novel"] = False
    row["augMP_recovered"] = True
    row["score"] = int(round(100 * hit.pid))
    row["transcript_score"] = hit.pid * hit.qcov
    row["valid_start"] = True
    row["valid_stop"] = True
    row["proper_orf"] = True
    row["frameshift"] = False
    row["original_gene_id"] = hit.query
    return row


def write_outputs(
    genes: list[dict],
    info: pd.DataFrame,
    assign: dict[str, dict],
    calls: dict[str, str],
    orphans: list[tuple[dict, Hit]],
    keep_tx: set[str],
    args,
    stats: dict,
    seleno_symbols: set[str] | None = None,
) -> None:
    from cat.consensus import write_consensus_fastas, write_consensus_gff3

    info_cols = list(info.columns)
    info_by_tx = {str(r["transcript_id"]): r for r in info.to_dict("records")}
    out_info_rows = []
    gene_dict: dict = defaultdict(dict)

    for g in genes:
        if g["transcript_id"] not in keep_tx:
            continue
        rec = dict(info_by_tx.get(g["transcript_id"], {}))
        rec["gene_id"] = g["gene_id"]
        rec["transcript_id"] = g["transcript_id"]
        a = assign.get(g["gene_id"])
        if a is not None:
            rec["protein_refine"] = a["label"]
            if a["label"] != "confirmed":
                rec["source_gene_common_name"] = a["symbol"]
                rec["source_gene"] = a["symbol"]
                rec["transcript_class"] = "ortholog"
                rec["gene_biotype"] = "protein_coding"
                rec["transcript_biotype"] = "protein_coding"
        else:
            rec["protein_refine"] = "unchanged"
        call = calls.get(g["gene_id"])
        rec["protein_refine_call"] = call or "N/A"
        if call in DEMOTE_CALLS:
            for col in ("gene_biotype", "transcript_biotype"):
                if rec.get(col) == "protein_coding":
                    rec[col] = "unknown_likely_coding"
            if call == "fragment":
                rec["source_gene_common_name"] = ""
        out_info_rows.append(rec)
        attrs = {k: rec[k] for k in rec}
        gene_dict[g["chrom"]].setdefault(g["gene_id"], []).append((g["tx"], attrs))

    for rec, hit in orphans:
        out_info_rows.append(rec)
        tx = rec["_tx"]
        attrs = {k: v for k, v in rec.items() if k != "_tx"}
        if not attrs.get("protein_refine") or attrs.get("protein_refine") in NA:
            attrs["protein_refine"] = "rescued"
        gene_dict[tx.chromosome].setdefault(rec["gene_id"], []).append((tx, attrs))

    gp_lines = []
    for chrom in sorted(gene_dict):
        for gene_id, tx_list in gene_dict[chrom].items():
            for tx_obj, attrs in tx_list:
                name2 = attrs.get("source_gene_common_name") or attrs.get("gene_id") or tx_obj.name2
                gp_lines.append("\t".join(tx_obj.get_gene_pred(name2=name2)))

    Path(args.output_gp).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_gp, "w") as fh:
        fh.write("\n".join(gp_lines) + ("\n" if gp_lines else ""))

    out_df = pd.DataFrame(out_info_rows)
    if "_tx" in out_df.columns:
        out_df = out_df.drop(columns=["_tx"])
    # Keep original column order, append new columns at the end.
    ordered = [c for c in info_cols if c in out_df.columns]
    extra = [c for c in out_df.columns if c not in ordered]
    out_df[ordered + extra].to_csv(args.output_gp_info, sep="\t", index=False)

    write_consensus_gff3(gene_dict, args.output_gff3)
    write_consensus_fastas(
        gene_dict, args.output_fasta, args.output_protein_fasta, args.fasta, seleno_symbols
    )

    with open(args.output_metrics_json, "w") as fh:
        json.dump(stats, fh, indent=2)
        fh.write("\n")


def copy_through(args) -> int:
    shutil.copy2(args.consensus_gp, args.output_gp)
    shutil.copy2(args.consensus_gp_info, args.output_gp_info)
    shutil.copy2(args.consensus_gff3, args.output_gff3)
    shutil.copy2(args.consensus_fasta, args.output_fasta)
    shutil.copy2(args.consensus_protein_fasta, args.output_protein_fasta)
    src_json = args.consensus_metrics_json
    if src_json and Path(src_json).exists():
        shutil.copy2(src_json, args.output_metrics_json)
    else:
        Path(args.output_metrics_json).write_text("{}\n")
    logger.info("copied consensus outputs unchanged")
    return 0


def run(args) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    if args.copy_through or not args.miniprot_paf or not Path(args.miniprot_paf).exists():
        return copy_through(args)

    logger.info("protein_refine %s", args.genome)
    query_symbols, uniprot_aa, seleno_symbols = load_query_symbols(args.protein_fasta, args.ref_gp_attrs)
    ref_order = load_ref_gene_order(args.ref_gp, query_symbols)
    ref_symbols = set()
    for genes in ref_order.values():
        ref_symbols.update(genes)
    if not ref_symbols:
        ref_symbols = set(query_symbols.values())
    logger.info("  reference symbols: %s", f"{len(ref_symbols):,}")

    hits = load_compact_hits(args.miniprot_paf, query_symbols)
    all_genes, info = load_genes(args.consensus_gp, args.consensus_gp_info)
    logger.info("  consensus transcripts: %s", f"{len(all_genes):,}")

    pc_genes = [
        g
        for g in all_genes
        if g["biotype"] == "protein_coding"
        or str(g["transcript_class"]).startswith("putative_novel")
        or g["cds_size"] > 0
    ]
    proteins = load_proteins(
        args.consensus_protein_fasta,
        {g["transcript_id"] for g in pc_genes},
    )
    work = [pc_genes[i] for i in pick_representatives(pc_genes, proteins, uniprot_aa, seleno_symbols)]
    chrom_hits = hits_by_chrom(hits)
    index_by_chrom = {c: LocusIndex(hs) for c, hs in chrom_hits.items()}
    ref_struct = load_ref_structure(args.ref_gp, query_symbols, uniprot_aa)

    assign, calls, present = assign_names(work, index_by_chrom, ref_symbols, ref_struct, proteins)
    n_frag = demote_redundant_fragments(work, assign, calls, present, ref_struct)
    split_hits, split = split_chimeras(
        work, calls, assign, present, index_by_chrom, ref_struct, ref_symbols
    )
    cands = rescue_candidates(hits, present, ref_struct, ref_symbols)
    split_extras: list[Hit] = []
    if split_hits:
        seen = {(h.symbol, h.chrom, h.start, h.end) for h in cands}
        for h in split_hits:
            key = (h.symbol, h.chrom, h.start, h.end)
            if key in seen:
                continue
            if present.get(h.symbol, 0) <= 0:
                cands.append(h)
                seen.add(key)
            else:
                split_extras.append(h)
        cands.sort(key=lambda h: (completeness(h, ref_struct), h.pid), reverse=True)
    named_spans: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for wi, g in enumerate(work):
        sym = assign[wi]["symbol"] if wi in assign else g["source_name"]
        if sym:
            named_spans[(sym, g["chrom"])].append((int(g["start"]), int(g["end"])))
    copies = split_extras + extra_copy_candidates(hits, present, named_spans, ref_struct, ref_symbols)
    labels = Counter(a["label"] for a in assign.values())
    call_counts = Counter(calls.values())
    logger.info("  names: %s", ", ".join(f"{k}={v}" for k, v in sorted(labels.items())))
    if n_frag:
        logger.info("  demoted %s redundant fragment copies", n_frag)
    logger.info(
        "  unnamed-model calls: %s",
        ", ".join(f"{k}={v}" for k, v in call_counts.most_common()),
    )
    logger.info(
        "  rescue candidates: %s hits for %s absent symbols; extra copies: %s",
        len(cands), len({h.symbol for h in cands}), len(copies),
    )

    needed = set()
    for wi, a in assign.items():
        chosen = choose_protein_span(work[wi], a["hit"], index_by_chrom.get(work[wi]["chrom"]), ref_struct)
        needed.add((chosen.chrom, chosen.start, chosen.end, chosen.query))
    for h in cands:
        needed.add((h.chrom, h.start, h.end, h.query))
    for h in copies:
        needed.add((h.chrom, h.start, h.end, h.query))
    exon_lookup = fetch_exons(args.miniprot_paf, query_symbols, needed)
    logger.info("  exon records fetched: %s / %s", len(exon_lookup), len(needed))

    n_shrink, n_expand = apply_structure(
        work, assign, index_by_chrom, exon_lookup, ref_struct, proteins, seleno_symbols
    )
    logger.info("  shrunk %s chimeric; expanded %s truncated", n_shrink, n_expand)

    seq_dict = get_sequence_dict(args.fasta) if args.fasta else {}
    superseded = set(split)
    n_drop_broken = 0
    n_promote = 0
    broken_only: set[str] = set()
    if seq_dict:
        for wi, g in enumerate(work):
            try:
                symbol = assign[wi]["symbol"] if wi in assign else g.get("source_name") or ""
                aa = protein_from_transcript(g["tx"], seq_dict, seleno=symbol in seleno_symbols)
            except (ValueError, KeyError, AssertionError, RuntimeError, TypeError):
                continue
            if not clearly_broken_peptide(aa):
                continue
            sibling = intact_sibling(g, all_genes, proteins, symbol, ref_struct)
            if sibling is not None:
                broken_only.add(g["transcript_id"])
                work[wi] = sibling
                n_promote += 1
                continue
            if wi in assign:
                present[assign[wi]["symbol"]] -= 1
                del assign[wi]
            g["source_name"] = ""
            superseded.add(wi)
            n_drop_broken += 1
        if n_drop_broken:
            logger.info("  dropped %s still-broken coding models for miniprot rescue", n_drop_broken)
        if n_promote:
            logger.info("  kept %s genes by promoting an intact sibling isoform", n_promote)

    if not seq_dict:
        seq_dict = load_fasta_chroms(args.fasta, {h.chrom for h in copies})
    rescues, rescued_super = place_rescues(
        cands, exon_lookup, work, assign, ref_struct, extra_cands=copies, seq_dict=seq_dict, seleno_symbols=seleno_symbols
    )
    superseded.update(rescued_super)
    n_extra_placed = sum(1 for h, _ in rescues if present.get(h.symbol, 0) > 0)
    logger.info(
        "  rescued %s genes (%s extra copies); superseded %s unnamed/chimeric models; "
        "split %s fusions",
        len(rescues), n_extra_placed, len(superseded), len(split),
    )

    # Broadcast names onto every transcript of an assigned gene.
    by_gid = defaultdict(list)
    for i, g in enumerate(all_genes):
        by_gid[g["gene_id"]].append(i)
    work_by_gid = {work[wi]["gene_id"]: (wi, a) for wi, a in assign.items()}
    keep_tx = {g["transcript_id"] for g in all_genes}
    n_drop_iso = 0
    for gid, (wi, a) in work_by_gid.items():
        refined = work[wi]
        rspan = (int(refined["start"]), int(refined["end"]))
        for gi in by_gid[gid]:
            g = all_genes[gi]
            if g["transcript_id"] == refined["transcript_id"]:
                all_genes[gi] = refined
                continue
            ov = overlap_bp((int(g["start"]), int(g["end"])), rspan)
            gl = max(1, int(g["end"]) - int(g["start"]) + 1)
            if ov / gl < 0.50:
                prot = proteins.get(g["transcript_id"])
                if isoform_worth_keeping(g, prot, a.get("symbol") or "", ref_struct):
                    continue
                keep_tx.discard(g["transcript_id"])
                n_drop_iso += 1

    n_drop_sib, n_keep_sib = prune_superseded_isoforms(
        keep_tx, work, superseded, all_genes, proteins, ref_struct
    )
    for tid in broken_only:
        keep_tx.discard(tid)
    if n_keep_sib:
        logger.info(
            "  kept %s intact sibling isoforms after superseding a broken/chimeric representative",
            n_keep_sib,
        )

    orphan_rows = []
    for h, tx in rescues:
        tx_id = tx.name
        gene_id = tx_id
        template = info.iloc[0].to_dict() if len(info) else {}
        rec = info_row_for_orphan(template, gene_id, tx_id, h, args.genome)
        rec["_tx"] = tx
        extra = present.get(h.symbol, 0) > 0
        rec["protein_refine"] = "named_copy" if extra else "rescued"
        rec["protein_refine_call"] = "named_copy" if extra else "rescued"
        orphan_rows.append((rec, h))

    stats = {
        "genome": args.genome,
        "names": dict(labels),
        "unnamed_model_calls": dict(call_counts),
        "n_rescued": len(orphan_rows),
        "n_extra_copies": n_extra_placed,
        "n_superseded": len(superseded),
        "n_shrink": n_shrink,
        "n_expand": n_expand,
        "n_drop_isoform": n_drop_iso,
        "n_drop_superseded_siblings": n_drop_sib,
        "n_keep_superseded_siblings": n_keep_sib,
        "n_promote_sibling": n_promote,
        "n_fragment_demoted": n_frag,
    }
    if args.consensus_metrics_json and Path(args.consensus_metrics_json).exists():
        try:
            prev = json.loads(Path(args.consensus_metrics_json).read_text())
            if isinstance(prev, dict):
                prev = dict(prev)
                prev["protein_refine"] = stats
                stats = prev
        except json.JSONDecodeError:
            pass

    assign_by_gid = {work[wi]["gene_id"]: a for wi, a in assign.items()}
    calls_by_gid = {work[wi]["gene_id"]: c for wi, c in calls.items()}
    write_outputs(
        all_genes, info, assign_by_gid, calls_by_gid, orphan_rows, keep_tx, args, stats, seleno_symbols
    )
    logger.info("  wrote %s", args.output_gp)
    return 0


def _safe(sym: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", sym)[:40]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--consensus-gp", required=True)
    ap.add_argument("--consensus-gp-info", required=True)
    ap.add_argument("--consensus-gff3", required=True)
    ap.add_argument("--consensus-fasta", required=True)
    ap.add_argument("--consensus-protein-fasta", required=True)
    ap.add_argument("--consensus-metrics-json", default="")
    ap.add_argument("--miniprot-paf", default="")
    ap.add_argument("--protein-fasta", default="")
    ap.add_argument("--ref-gp-attrs", default="")
    ap.add_argument("--ref-gp", default="")
    ap.add_argument("--fasta", required=True, help="Target genome FASTA")
    ap.add_argument("--output-gp", required=True)
    ap.add_argument("--output-gp-info", required=True)
    ap.add_argument("--output-gff3", required=True)
    ap.add_argument("--output-fasta", required=True)
    ap.add_argument("--output-protein-fasta", required=True)
    ap.add_argument("--output-metrics-json", required=True)
    ap.add_argument(
        "--copy-through",
        action="store_true",
        help="Copy consensus files unchanged (no miniprot / augMP genome)",
    )
    args = ap.parse_args(argv)
    try:
        return run(args)
    except Exception:
        logger.exception("protein_refine failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
