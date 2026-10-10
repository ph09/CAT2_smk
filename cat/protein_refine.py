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
import shutil
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from cat.convert_miniprot_to_genepred import (
    _exon_frames,
    parse_cigar_to_exons,
    parse_paf_row,
)
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


def load_query_symbols(protein_fasta: str | None, gp_attrs: str | None) -> dict[str, str]:
    """Map protein query IDs (FASTA header token / ENST) to gene symbols."""
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
    if protein_fasta and Path(protein_fasta).exists():
        n_gn = 0
        with open(protein_fasta) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                qid = header.split()[0]
                gn = GN_RE.search(header)
                if gn:
                    out[qid] = gn.group(1)
                    n_gn += 1
                    continue
                if qid in out:
                    continue
                # Human GFF-derived headers are often just ENST ids.
                if qid in out or qid.split(".")[0] in out:
                    continue
        logger.info("  UniProt GN= symbols from protein FASTA: %s", f"{n_gn:,}")
    return out


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


def load_ref_structure(ref_gp: str | None, query_symbols: dict[str, str]) -> dict[str, tuple[int, int]]:
    """Reference symbol → (cds_aa, coding_exons) of its longest coding isoform.

    Completeness is judged against the gene, not the matched peptide, which
    can itself be a short isoform or a UniProt fragment.
    """
    out: dict[str, tuple[int, int]] = {}
    if not ref_gp or not Path(ref_gp).exists():
        return out
    for tx in gene_pred_iterator(ref_gp):
        if tx.cds_size <= 0:
            continue
        sym = query_symbols.get(tx.name) or query_symbols.get(tx.name.split(".")[0])
        if not sym:
            sym = tx.name2 if tx.name2 and tx.name2 not in NA else None
        if not sym:
            continue
        rec = (tx.cds_size // 3, coding_exon_count(tx))
        if rec > out.get(sym, (0, 0)):
            out[sym] = rec
    return out


def parent_structure(sym: str, ref_struct: dict[str, tuple[int, int]]) -> tuple[int, int]:
    return ref_struct.get(sym, (0, 0))


def completeness(h, ref_struct: dict[str, tuple[int, int]]) -> float:
    """Identity × fraction of the reference gene's peptide covered by the hit."""
    p_aa = parent_structure(h.symbol, ref_struct)[0]
    if not p_aa or not h.qlen:
        return h.score
    return h.pid * min(1.0, h.qcov * h.qlen / p_aa)


def is_complete_copy(n_exons: int, aa: int, sym: str, ref_struct: dict[str, tuple[int, int]]) -> bool:
    p_aa, p_ex = parent_structure(sym, ref_struct)
    if not p_aa:
        return False
    if not MIN_COPY_AA_FRAC * p_aa <= aa <= MAX_COPY_AA_FRAC * p_aa:
        return False
    if p_ex >= 2 and n_exons < max(2, math.ceil(MIN_COPY_EXON_FRAC * p_ex)):
        return False
    return True


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
            exons, matched_aa, aligned_aa, _fs = parse_cigar_to_exons(cg, rec["t_start"])
            if not exons:
                continue
            qcov = aligned_aa / qlen
            pid = matched_aa / aligned_aa if aligned_aa else 0.0
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
                "cds_size": tx.cds_size,
                "coding_exons": coding_exon_count(tx),
            }
        )
    return genes, info


def pick_representatives(genes: list[dict]) -> list[int]:
    """Longest-CDS transcript index per gene_id."""
    best: dict[str, int] = {}
    for i, g in enumerate(genes):
        prev = best.get(g["gene_id"])
        if prev is None or g["cds_size"] > genes[prev]["cds_size"]:
            best[g["gene_id"]] = i
        elif g["cds_size"] == genes[prev]["cds_size"] and (g["end"] - g["start"]) > (
            genes[prev]["end"] - genes[prev]["start"]
        ):
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
    if g["cds_size"] <= 0 or g["cds_size"] % 3:
        return False
    if str(g["info"].get("frameshift", "")).lower() in {"true", "1"}:
        return False
    if protein is None:
        return True
    return protein.startswith("M") and protein.endswith("*") and "*" not in protein[:-1]


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
            own_pid = own.pid if own else 0.0
            if rival and rival.pid >= MIN_PRIMARY_PID and rival.pid >= own_pid + CAT_OVERRIDE_MARGIN:
                corrections.append((rival.pid, rival.qcov, gi, rival, own))
            elif own and own.pid >= MIN_NAMED_PID:
                assign[gi] = _assignment(src, own, "confirmed")
            continue
        call, hit = coding_call(g, by_sym, ref_struct, proteins.get(g["transcript_id"]))
        calls[gi] = call
        if call == "real_copy":
            copies.append((completeness(hit, ref_struct), hit.pid, gi, hit))

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

    copies.sort(key=lambda x: (x[0], x[1]), reverse=True)
    for _sc, _pid, gi, hit in copies:
        g = genes[gi]
        span = (int(g["start"]), int(g["end"]))
        key = (hit.symbol, g["chrom"])
        if any(overlap_bp(span, s) > 0 for s in named_spans[key]):
            calls[gi] = "redundant"
            continue
        named_spans[key].append(span)
        label = "named_copy" if present[hit.symbol] > 0 else "named_missing"
        assign[gi] = _assignment(hit.symbol, hit, label)
        present[hit.symbol] += 1

    return assign, calls, present


def rescue_candidates(
    hits: list[Hit],
    present: Counter,
    ref_struct: dict[str, tuple[int, int]],
) -> list[Hit]:
    """Complete, intron-consistent hits to reference genes absent from the genome.

    Up to RESCUE_CANDIDATES_PER_SYMBOL per symbol, most complete first; the
    exon-level placement check in ``place_rescues`` picks at most one.
    """
    by_sym: dict[str, list[Hit]] = defaultdict(list)
    for h in hits:
        if present.get(h.symbol, 0) > 0:
            continue
        p_aa, p_ex = parent_structure(h.symbol, ref_struct)
        if not p_aa or h.pid < MIN_COPY_PID or h.qcov < MIN_COPY_QCOV:
            continue
        covered = h.qcov * h.qlen
        if not MIN_COPY_AA_FRAC * p_aa <= covered <= MAX_COPY_AA_FRAC * p_aa:
            continue
        if p_ex >= 2 and h.n_exons < max(2, math.ceil(MIN_COPY_EXON_FRAC * p_ex)):
            continue
        by_sym[h.symbol].append(h)
    out: list[Hit] = []
    for sym, hs in by_sym.items():
        hs.sort(key=lambda h: (completeness(h, ref_struct), h.pid), reverse=True)
        out.extend(hs[:RESCUE_CANDIDATES_PER_SYMBOL])
    out.sort(key=lambda h: (completeness(h, ref_struct), h.pid), reverse=True)
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


def place_rescues(
    cands: list[Hit],
    exon_lookup: dict[tuple, list[tuple[int, int]]],
    genes: list[dict],
    assign: dict[int, dict],
    ref_struct: dict[str, tuple[int, int]],
) -> tuple[list[tuple[Hit, GenePredTranscript]], set[int]]:
    """Accept one rescue per absent symbol.

    A rescue is blocked when more than RESCUE_MAX_NAMED_OVERLAP of its CDS
    falls on coding exons of a named gene (same strand). Unnamed coding models
    whose exons it overlaps are fragments of it and are superseded.
    """
    by_chrom: dict[str, list[int]] = defaultdict(list)
    for gi, g in enumerate(genes):
        by_chrom[g["chrom"]].append(gi)
    segs_cache: dict[int, list[tuple[int, int]]] = {}
    accepted: list[tuple[Hit, GenePredTranscript]] = []
    taken: dict[tuple[str, str], list[list[tuple[int, int]]]] = defaultdict(list)
    superseded: set[int] = set()
    used: set[str] = set()
    for h in cands:
        if h.symbol in used:
            continue
        exons = exon_lookup.get((h.chrom, h.start, h.end, h.query))
        if not exons:
            continue
        tx = GenePredTranscript(
            gp_line_from_exons(f"MPR-{_safe(h.symbol)}-{h.start}", h.chrom, h.strand, exons, h.symbol).split("\t")
        )
        if not is_complete_copy(coding_exon_count(tx), tx.cds_size // 3, h.symbol, ref_struct):
            continue
        segs = _coding_segments(tx)
        if any(_segment_overlap(segs, t) > 0 for t in taken[(h.chrom, h.strand)]):
            continue
        named_ov = 0
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
            if gi in assign or g["source_name"]:
                named_ov += ov
            else:
                unnamed.append(gi)
        if named_ov > RESCUE_MAX_NAMED_OVERLAP * tx.cds_size:
            continue
        accepted.append((h, tx))
        superseded.update(unnamed)
        taken[(h.chrom, h.strand)].append(segs)
        used.add(h.symbol)
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
) -> tuple[int, int]:
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
        if new_span + 50 >= old_span and new_span - 50 <= old_span:
            continue
        line = gp_line_from_exons(
            g["transcript_id"], chosen.chrom, chosen.strand, exons, a["symbol"]
        )
        new_tx = GenePredTranscript(line.split("\t"))
        genes[gi]["tx"] = new_tx
        genes[gi]["start"] = new_tx.start
        genes[gi]["end"] = new_tx.stop
        genes[gi]["strand"] = new_tx.strand
        a["hit"] = chosen
        if new_span < old_span:
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
        out_info_rows.append(rec)
        attrs = {k: rec[k] for k in rec}
        gene_dict[g["chrom"]].setdefault(g["gene_id"], []).append((g["tx"], attrs))

    for rec, hit in orphans:
        out_info_rows.append(rec)
        tx = rec["_tx"]
        attrs = {k: v for k, v in rec.items() if k != "_tx"}
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
    write_consensus_fastas(gene_dict, args.output_fasta, args.output_protein_fasta, args.fasta)

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
    query_symbols = load_query_symbols(args.protein_fasta, args.ref_gp_attrs)
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
    work = [pc_genes[i] for i in pick_representatives(pc_genes)]
    chrom_hits = hits_by_chrom(hits)
    index_by_chrom = {c: LocusIndex(hs) for c, hs in chrom_hits.items()}
    ref_struct = load_ref_structure(args.ref_gp, query_symbols)
    proteins = load_proteins(
        args.consensus_protein_fasta,
        {g["transcript_id"] for g in work if not g["source_name"]},
    )

    assign, calls, present = assign_names(work, index_by_chrom, ref_symbols, ref_struct, proteins)
    cands = rescue_candidates(hits, present, ref_struct)
    labels = Counter(a["label"] for a in assign.values())
    call_counts = Counter(calls.values())
    logger.info("  names: %s", ", ".join(f"{k}={v}" for k, v in sorted(labels.items())))
    logger.info(
        "  unnamed-model calls: %s",
        ", ".join(f"{k}={v}" for k, v in call_counts.most_common()),
    )
    logger.info(
        "  rescue candidates: %s hits for %s absent symbols",
        len(cands), len({h.symbol for h in cands}),
    )

    needed = set()
    for wi, a in assign.items():
        chosen = choose_protein_span(work[wi], a["hit"], index_by_chrom.get(work[wi]["chrom"]), ref_struct)
        needed.add((chosen.chrom, chosen.start, chosen.end, chosen.query))
    for h in cands:
        needed.add((h.chrom, h.start, h.end, h.query))
    exon_lookup = fetch_exons(args.miniprot_paf, query_symbols, needed)
    logger.info("  exon records fetched: %s / %s", len(exon_lookup), len(needed))

    n_shrink, n_expand = apply_structure(work, assign, index_by_chrom, exon_lookup, ref_struct)
    logger.info("  shrunk %s chimeric; expanded %s truncated", n_shrink, n_expand)

    rescues, superseded = place_rescues(cands, exon_lookup, work, assign, ref_struct)
    logger.info(
        "  rescued %s absent genes; superseded %s unnamed fragment models",
        len(rescues), len(superseded),
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
                keep_tx.discard(g["transcript_id"])
                n_drop_iso += 1

    for wi in superseded:
        for gi in by_gid[work[wi]["gene_id"]]:
            keep_tx.discard(all_genes[gi]["transcript_id"])

    orphan_rows = []
    for h, tx in rescues:
        tx_id = tx.name
        gene_id = tx_id
        template = info.iloc[0].to_dict() if len(info) else {}
        rec = info_row_for_orphan(template, gene_id, tx_id, h, args.genome)
        rec["_tx"] = tx
        rec["protein_refine"] = "rescued"
        rec["protein_refine_call"] = "rescued"
        orphan_rows.append((rec, h))

    stats = {
        "genome": args.genome,
        "names": dict(labels),
        "unnamed_model_calls": dict(call_counts),
        "n_rescued": len(orphan_rows),
        "n_superseded": len(superseded),
        "n_shrink": n_shrink,
        "n_expand": n_expand,
        "n_drop_isoform": n_drop_iso,
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
    write_outputs(all_genes, info, assign_by_gid, calls_by_gid, orphan_rows, keep_tx, args, stats)
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
