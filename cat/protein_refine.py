#!/usr/bin/env python3
"""Correct consensus gene structures and names from miniprot protein hits.

Projection (transMap / augTM) can fuse neighbouring paralogs, drop long first
introns, and carry a wrong ``source_gene_common_name``. Miniprot already mapped
the protein DB; this pass uses those hits to:

1. Reassign exclusive identities (one primary copy per reference gene symbol,
   plus near-identical extras). CAT names lose if another peptide is clearly
   better. Remaining models get unique-best names only when the margin is real.
2. Replace chimeric / truncated CDS with the owned full-protein miniprot locus.
3. Rescue unused reference genes whose protein hits sit in a clean syntenic gap.

Intended to run after ``generate_consensus`` and before ``annotate_novel_genes``.
"""
from __future__ import annotations

import argparse
import json
import logging
import re
import shutil
from collections import defaultdict
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
MIN_EXTRA_PID = 0.80
MIN_ORTHOLOG_PID = 0.68
MIN_ORTHOLOG_QCOV = 0.40
UNIQUE_MARGIN = 0.012
UNIQUE_UNUSED_MARGIN = 0.008
CAT_OVERRIDE_MARGIN = 0.08
FULL_PROTEIN_QCOV = 0.80
MIN_PROTEIN_SPAN = 2_000
MIN_RECOVER_SPAN = 2_000
ASSIGN_MAX_GAP = 20_000
SYNTENY_GAP_MIN = 8_000
SYNTENY_GAP_MAX = 500_000


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


def interval_gap(a: tuple[int, int], b: tuple[int, int]) -> int:
    if a[1] < b[0]:
        return b[0] - a[1]
    if b[1] < a[0]:
        return a[0] - b[1]
    return 0


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


class Hit:
    __slots__ = ("chrom", "start", "end", "strand", "query", "symbol", "pid", "qcov")

    def __init__(self, chrom, start, end, strand, query, symbol, pid, qcov):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.query = query
        self.symbol = symbol
        self.pid = float(pid)
        self.qcov = float(qcov)

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
                rec["q_name"], sym, pid, qcov,
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


def best_gene_for_hit(h: Hit, gene_spans: list[tuple[int, int] | None], genes: list[dict]) -> int | None:
    hb = h.span
    best_gi, best_ov, nearby = None, 0, []
    for gi, gsp in enumerate(gene_spans):
        if gsp is None:
            continue
        if genes[gi]["chrom"] != h.chrom:
            continue
        ov = overlap_bp(gsp, hb)
        if ov > best_ov:
            best_ov, best_gi = ov, gi
        elif ov == 0:
            strand = genes[gi].get("strand") or ""
            if strand and h.strand and strand != h.strand:
                continue
            gap = interval_gap(gsp, hb)
            if gap <= ASSIGN_MAX_GAP:
                nearby.append((gap, gi))
    if best_ov >= 100:
        return best_gi
    if nearby:
        nearby.sort()
        return nearby[0][1]
    return None


def assign_identities(
    genes: list[dict],
    hits: list[Hit],
    ref_symbols: set[str],
) -> tuple[dict[int, dict], list[Hit]]:
    gene_spans: list[tuple[int, int] | None] = []
    for g in genes:
        gene_spans.append((int(g["start"]), int(g["end"])))

    used_gene: set[int] = set()
    used_query: set[str] = set()
    assign: dict[int, dict] = {}

    named = []
    for gi, g in enumerate(genes):
        src = g.get("source_name") or ""
        if src not in ref_symbols:
            continue
        gsp = gene_spans[gi]
        best = None
        rival = None
        for h in hits:
            if h.chrom != g["chrom"]:
                continue
            if overlap_bp(gsp, h.span) < 100:
                continue
            if h.symbol == src:
                if best is None or (h.pid, h.qcov) > (best.pid, best.qcov):
                    best = h
            else:
                if rival is None or h.pid > rival.pid:
                    rival = h
        if best is None or best.pid < MIN_NAMED_PID:
            continue
        if rival is not None and rival.pid >= best.pid + CAT_OVERRIDE_MARGIN:
            continue
        named.append((best.pid, g["cds_size"], gi, src, best))
    named.sort(reverse=True)
    for _pid, _cds, gi, query, h in named:
        if gi in used_gene or query in used_query:
            continue
        assign[gi] = {"symbol": query, "pid": h.pid, "qcov": h.qcov, "hit": h, "extra": False}
        used_gene.add(gi)
        used_query.add(query)

    scores = []
    for h in hits:
        if h.symbol not in ref_symbols:
            continue
        gi = best_gene_for_hit(h, gene_spans, genes)
        if gi is None:
            continue
        src = genes[gi].get("source_name") or ""
        boost = 0.05 if src == h.symbol else 0.0
        scores.append((h.pid * h.qcov + boost, genes[gi]["cds_size"], h.pid, h.qcov, gi, h))
    scores.sort(key=lambda x: (x[0], x[1], x[2]), reverse=True)
    for _sc, _cds, pid, qcov, gi, h in scores:
        if gi in used_gene:
            continue
        query = h.symbol
        if query not in used_query:
            if pid < MIN_PRIMARY_PID or qcov < MIN_PRIMARY_QCOV:
                continue
            assign[gi] = {"symbol": query, "pid": pid, "qcov": qcov, "hit": h, "extra": False}
            used_gene.add(gi)
            used_query.add(query)
        else:
            if pid < MIN_EXTRA_PID or qcov < MIN_PRIMARY_QCOV:
                continue
            assign[gi] = {"symbol": query, "pid": pid, "qcov": qcov, "hit": h, "extra": True}
            used_gene.add(gi)

    for gi, g in enumerate(genes):
        if gi in used_gene:
            continue
        gsp = gene_spans[gi]
        by_q: dict[str, Hit] = {}
        for h in hits:
            if h.chrom != g["chrom"] or h.symbol not in ref_symbols:
                continue
            if overlap_bp(gsp, h.span) < 80:
                continue
            prev = by_q.get(h.symbol)
            if prev is None or h.pid * h.qcov > prev.pid * prev.qcov:
                by_q[h.symbol] = h
        if not by_q:
            continue
        ranked = sorted(by_q.values(), key=lambda h: h.pid * h.qcov, reverse=True)
        best = ranked[0]
        if best.pid < MIN_ORTHOLOG_PID or best.qcov < MIN_ORTHOLOG_QCOV:
            continue
        margin = 1.0
        if len(ranked) > 1:
            margin = best.pid * best.qcov - ranked[1].pid * ranked[1].qcov
        extra = best.symbol in used_query
        if extra:
            if margin < UNIQUE_MARGIN or best.pid < 0.72:
                continue
        else:
            second_used = len(ranked) > 1 and ranked[1].symbol in used_query
            if margin < UNIQUE_UNUSED_MARGIN and not second_used:
                continue
            if margin < 0.004:
                continue
        assign[gi] = {
            "symbol": best.symbol,
            "pid": best.pid,
            "qcov": best.qcov,
            "hit": best,
            "extra": extra,
        }
        used_gene.add(gi)
        used_query.add(best.symbol)

    orphans = synteny_orphans(genes, assign, hits, ref_symbols, used_query)
    return assign, orphans


def synteny_orphans(
    genes: list[dict],
    assign: dict[int, dict],
    hits: list[Hit],
    ref_symbols: set[str],
    used_query: set[str],
) -> list[Hit]:
    """Unused symbols whose best hit sits in a clean gap between assigned neighbours."""
    primary: dict[str, tuple[str, int, int]] = {}
    for gi, a in assign.items():
        if a.get("extra"):
            continue
        g = genes[gi]
        primary[a["symbol"]] = (g["chrom"], int(g["start"]), int(g["end"]))

    occupied: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for gi, a in assign.items():
        g = genes[gi]
        occupied[g["chrom"]].append((int(g["start"]), int(g["end"])))

    # Neighbour walk uses target genomic order of *assigned primaries*, not the
    # human chromosome (which rearranges in NWM). A symbol is rescued only if
    # some pair of already-named genes on the same contig leave a 8–500 kb gap
    # containing a high-coverage hit to that unused symbol and nothing else.
    prim_by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for sym, (chrom, s, e) in primary.items():
        prim_by_chrom[chrom].append((s, e, sym))
    for chrom in prim_by_chrom:
        prim_by_chrom[chrom].sort()

    orphans: list[Hit] = []
    for h in sorted(hits, key=lambda x: (-x.pid, -x.qcov, -x.length)):
        if h.symbol not in ref_symbols or h.symbol in used_query:
            continue
        if h.pid < MIN_NAMED_PID or h.qcov < FULL_PROTEIN_QCOV:
            continue
        if h.length < MIN_RECOVER_SPAN:
            continue
        if any(interval_gap(h.span, oc) == 0 for oc in occupied.get(h.chrom, [])):
            continue
        loc = prim_by_chrom.get(h.chrom) or []
        in_clean_gap = False
        for i in range(len(loc) - 1):
            _s1, e1, _left = loc[i]
            s2, _e2, _right = loc[i + 1]
            gap_lo, gap_hi = e1, s2
            gap = gap_hi - gap_lo
            if gap < SYNTENY_GAP_MIN or gap > SYNTENY_GAP_MAX:
                continue
            mid = (h.start + h.end) // 2
            if gap_lo < mid < gap_hi:
                in_clean_gap = True
                break
        if not in_clean_gap:
            continue
        orphans.append(h)
        used_query.add(h.symbol)
        occupied[h.chrom].append(h.span)
    return orphans


def choose_protein_span(gene: dict, hit: Hit, hits: list[Hit]) -> Hit:
    """Owned full-protein hit: rank by (qcov, span, pid), not identity first."""
    gs, ge = int(gene["start"]), int(gene["end"])
    chrom = gene["chrom"]
    query = hit.symbol
    cands = [
        h
        for h in hits
        if h.chrom == chrom
        and h.symbol == query
        and h.pid >= MIN_ORTHOLOG_PID
        and (
            overlap_bp((gs, ge), h.span) >= 50
            or interval_gap((gs, ge), h.span) <= 40_000
        )
    ]
    if not cands:
        cands = [hit]
    full = [h for h in cands if h.qcov >= FULL_PROTEIN_QCOV and h.length >= MIN_PROTEIN_SPAN]
    pool = full or cands
    return max(pool, key=lambda h: (h.qcov, h.length, h.pid))


def overlaps_other_gene(
    span: tuple[int, int],
    chrom: str,
    self_gi: int,
    genes: list[dict],
    assign: dict[int, dict],
) -> bool:
    for gj, a in assign.items():
        if gj == self_gi or a.get("extra") and genes[gj]["gene_id"] == genes[self_gi]["gene_id"]:
            continue
        g = genes[gj]
        if g["chrom"] != chrom:
            continue
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
    hits: list[Hit],
    exon_lookup: dict[tuple, list[tuple[int, int]]],
) -> tuple[int, int]:
    n_shrink = n_expand = 0
    for gi, a in sorted(assign.items(), key=lambda kv: (genes[kv[0]]["chrom"], genes[kv[0]]["start"])):
        g = genes[gi]
        chosen = choose_protein_span(g, a["hit"], hits)
        if chosen.qcov < FULL_PROTEIN_QCOV or chosen.length < MIN_PROTEIN_SPAN:
            continue
        if overlaps_other_gene(chosen.span, g["chrom"], gi, genes, assign):
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
    assign: dict[int, dict],
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
            rec["source_gene_common_name"] = a["symbol"]
            rec["source_gene"] = a["symbol"]
            rec["transcript_class"] = "ortholog"
            rec["gene_biotype"] = "protein_coding"
            rec["transcript_biotype"] = "protein_coding"
            rec["protein_refine"] = "extra" if a["extra"] else "primary"
        else:
            rec.setdefault("protein_refine", "unchanged")
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

    assign: dict[int, dict] = {}
    orphans: list[Hit] = []
    for chrom, idxs in _group_by_chrom(work):
        local_genes = [work[i] for i in idxs]
        local_hits = chrom_hits.get(chrom, [])
        local_assign, local_orphans = assign_identities(local_genes, local_hits, ref_symbols)
        for gi, a in local_assign.items():
            assign[idxs[gi]] = a
        orphans.extend(local_orphans)
    logger.info(
        "  assigned %s genes (%s extras); %s synteny rescues",
        len(assign),
        sum(1 for a in assign.values() if a["extra"]),
        len(orphans),
    )

    needed = set()
    for wi, a in assign.items():
        chosen = choose_protein_span(work[wi], a["hit"], chrom_hits.get(work[wi]["chrom"], []))
        needed.add((chosen.chrom, chosen.start, chosen.end, chosen.query))
    for h in orphans:
        needed.add((h.chrom, h.start, h.end, h.query))
    exon_lookup = fetch_exons(args.miniprot_paf, query_symbols, needed)
    logger.info("  exon records fetched: %s / %s", len(exon_lookup), len(needed))

    # Map assign indices from `work` back onto `work` genes (already those objects).
    n_shrink, n_expand = apply_structure(work, assign, hits, exon_lookup)
    logger.info("  shrunk %s chimeric; expanded %s truncated", n_shrink, n_expand)

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

    orphan_rows = []
    for h in orphans:
        key = (h.chrom, h.start, h.end, h.query)
        exons = exon_lookup.get(key)
        if not exons:
            continue
        tx_id = f"MPR-{_safe(h.symbol)}-{h.start}"
        gene_id = tx_id
        line = gp_line_from_exons(tx_id, h.chrom, h.strand, exons, h.symbol)
        tx = GenePredTranscript(line.split("\t"))
        template = info.iloc[0].to_dict() if len(info) else {}
        rec = info_row_for_orphan(template, gene_id, tx_id, h, args.genome)
        rec["_tx"] = tx
        rec["protein_refine"] = "rescued"
        orphan_rows.append((rec, h))

    stats = {
        "genome": args.genome,
        "n_assigned": len(assign),
        "n_extra": sum(1 for a in assign.values() if a["extra"]),
        "n_rescued": len(orphan_rows),
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
    write_outputs(all_genes, info, assign_by_gid, orphan_rows, keep_tx, args, stats)
    logger.info("  wrote %s", args.output_gp)
    return 0


def _group_by_chrom(genes: list[dict]) -> list[tuple[str, list[int]]]:
    by = defaultdict(list)
    for i, g in enumerate(genes):
        by[g["chrom"]].append(i)
    return [(c, by[c]) for c in sorted(by)]


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
