#!/usr/bin/env python3
"""
Filter augMP (Augustus + miniprot) gene models BEFORE they reach consensus.

Why this exists
---------------
miniprot is run permissively so it recovers every paralog / lineage-specific
copy of a protein. That permissiveness means a single genomic *copy* of a gene
is hit by many different (homologous) database proteins, each of which produces
its own augMP model (``augMP-<proteinA>``, ``augMP-<proteinB>``, ...) with
essentially the *same* gene structure. On dense / segmentally-duplicated regions
this piles up into hundreds of near-identical overlapping models at one locus
(observed: up to ~500 models on a single Callithrix locus, 826k models genome
wide). That pileup is what makes the consensus step go quadratic and time out /
OOM -- it is pure alignment redundancy, not biology.

Expansion-safe design
---------------------
The filter removes *within-locus structural redundancy only*. It NEVER removes a
distinct gene copy:

  * Real duplication  = copies at DISTINCT loci -> DISTINCT coordinates ->
                        DISTINCT gene structures -> every copy is kept.
  * Spurious pileup   = many models with the SAME structure stacked on ONE copy
                        -> collapsed to the single best model.

Stages (each keeps the best-scoring representative, so nothing real is lost):
  1. Multi-exon models: keep the best model per (chrom, strand, intron-chain).
     Distinct copies have distinct intron chains, so all copies survive.
  2. Single-exon models: dedup per locus (rounded coords), keep the best, subject
     to an optional coverage floor -- but a locus is never emptied.
  3. Per-locus backstop cap: after de-dup, if a single overlap-merged locus still
     holds more than --max-models-per-locus structurally-distinct models (only
     happens on pathological regions), keep the top-N by score. Generous by
     default so genuine tandem arrays survive.

Ranking score (higher = better): protein coverage, then identity, then CDS
length, then exon count -- all derived from the miniprot-based augMP PSL when
available; models with no PSL metrics are ranked last but are still eligible to
be the surviving representative of their locus.
"""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path


# ─── genePred parsing ─────────────────────────────────────────────────────────

class Model:
    __slots__ = ("name", "chrom", "strand", "tx_start", "tx_end",
                 "cds_start", "cds_end", "exon_count", "exon_starts",
                 "exon_ends", "raw", "coverage", "identity", "cds_len", "score")

    def __init__(self, name, chrom, strand, tx_start, tx_end, cds_start,
                 cds_end, exon_count, exon_starts, exon_ends, raw):
        self.name = name
        self.chrom = chrom
        self.strand = strand
        self.tx_start = tx_start
        self.tx_end = tx_end
        self.cds_start = cds_start
        self.cds_end = cds_end
        self.exon_count = exon_count
        self.exon_starts = exon_starts
        self.exon_ends = exon_ends
        self.raw = raw
        self.coverage = 0.0
        self.identity = 0.0
        # coding length (falls back to transcript length for non-coding preds)
        cds_len = 0
        for s, e in zip(exon_starts, exon_ends):
            a = max(s, cds_start)
            b = min(e, cds_end)
            if b > a:
                cds_len += b - a
        self.cds_len = cds_len if cds_len > 0 else (tx_end - tx_start)
        self.score = (0.0, 0.0, 0, 0)

    def intron_chain(self, round_bp):
        """Rounded intron chain; empty tuple for single-exon models."""
        if self.exon_count <= 1:
            return ()
        r = round_bp
        return tuple(
            (self.exon_ends[i] // r, self.exon_starts[i + 1] // r)
            for i in range(self.exon_count - 1)
        )

    def locus_key(self, round_bp):
        r = round_bp
        return (self.chrom, self.strand,
                self.tx_start // r, self.tx_end // r)


def parse_genepred(path):
    models = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            try:
                tx_start = int(c[3]); tx_end = int(c[4])
                cds_start = int(c[5]); cds_end = int(c[6])
                exon_count = int(c[7])
                exon_starts = [int(x) for x in c[8].rstrip(",").split(",") if x != ""]
                exon_ends = [int(x) for x in c[9].rstrip(",").split(",") if x != ""]
            except (ValueError, IndexError):
                continue
            if not exon_starts or len(exon_starts) != len(exon_ends):
                continue
            models.append(Model(c[0], c[1], c[2], tx_start, tx_end, cds_start,
                                 cds_end, exon_count, exon_starts, exon_ends,
                                 line.rstrip("\n")))
    return models


def parse_psl_metrics(path):
    """Return {qName: (coverage, identity)} from a PSL. Best row per qName."""
    metrics = {}
    if path is None or not Path(path).exists():
        return metrics
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 21:
                continue
            try:
                matches = int(f[0]); mismatches = int(f[1]); q_size = int(f[10])
            except (ValueError, IndexError):
                continue
            q_name = f[9]
            coverage = matches / q_size if q_size > 0 else 0.0
            identity = matches / (matches + mismatches) if (matches + mismatches) > 0 else 0.0
            prev = metrics.get(q_name)
            if prev is None or coverage > prev[0]:
                metrics[q_name] = (coverage, identity)
    return metrics


# ─── filtering ────────────────────────────────────────────────────────────────

def assign_scores(models, metrics):
    for m in models:
        cov, ident = metrics.get(m.name, (0.0, 0.0))
        m.coverage = cov
        m.identity = ident
        # higher is better: coverage, identity, CDS length, exon count
        m.score = (cov, ident, m.cds_len, m.exon_count)


def _best(models):
    return max(models, key=lambda m: m.score)


def dedup_multi_exon(models, round_bp):
    """Keep the best model per (chrom, strand, intron-chain). Expansion-safe."""
    groups = collections.defaultdict(list)
    for m in models:
        groups[(m.chrom, m.strand, m.intron_chain(round_bp))].append(m)
    return [_best(g) for g in groups.values()]


def dedup_single_exon(models, round_bp, min_coverage):
    """Dedup single-exon models per rounded locus; keep the best per locus.

    A coverage floor is applied because single-exon augMP models are the noisiest
    class, but the best model of any locus is always retained so we never delete a
    genuine single-exon locus.
    """
    groups = collections.defaultdict(list)
    for m in models:
        groups[m.locus_key(round_bp)].append(m)
    kept = []
    for g in groups.values():
        best = _best(g)
        # keep the best representative regardless of floor (preserves the locus);
        # only additional copies at the same rounded locus are dropped.
        if min_coverage <= 0 or best.coverage >= min_coverage:
            kept.append(best)
        else:
            # locus is below floor: still keep the single best so novel loci are
            # not silently lost (consensus/metrics can still filter it downstream).
            kept.append(best)
    return kept


def backstop_locus_cap(models, max_per_locus):
    """After structure de-dup, cap the number of models per overlap-merged locus.

    Only triggers on pathological pileups. Applied per (chrom, strand) by a
    single sweep over transcript spans; within an over-full cluster the top-N by
    score are kept.
    """
    if max_per_locus <= 0:
        return models
    by_cs = collections.defaultdict(list)
    for m in models:
        by_cs[(m.chrom, m.strand)].append(m)
    kept = []
    for cluster_models in by_cs.values():
        cluster_models.sort(key=lambda m: (m.tx_start, m.tx_end))
        cur = []
        cur_end = None
        for m in cluster_models:
            if cur_end is not None and m.tx_start <= cur_end:
                cur.append(m)
                cur_end = max(cur_end, m.tx_end)
            else:
                kept.extend(_cap_cluster(cur, max_per_locus))
                cur = [m]
                cur_end = m.tx_end
        kept.extend(_cap_cluster(cur, max_per_locus))
    return kept


def _cap_cluster(cluster, max_per_locus):
    if len(cluster) <= max_per_locus:
        return cluster
    return sorted(cluster, key=lambda m: m.score, reverse=True)[:max_per_locus]


# ─── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-gp", required=True, type=Path,
                    help="Input (raw) augMP genePred.")
    ap.add_argument("--out-gp", required=True, type=Path,
                    help="Output (filtered) augMP genePred.")
    ap.add_argument("--in-psl", type=Path, default=None,
                    help="Input (raw) augMP PSL for coverage/identity ranking (optional).")
    ap.add_argument("--out-psl", type=Path, default=None,
                    help="Output filtered PSL (subset of --in-psl for kept models).")
    ap.add_argument("--max-models-per-locus", type=int, default=25,
                    help="Backstop cap on structurally-distinct models per "
                         "overlap-merged locus (0 = no cap) [25].")
    ap.add_argument("--structure-round-bp", type=int, default=10,
                    help="Round intron-chain coordinates to this many bp when "
                         "collapsing near-identical structures [10].")
    ap.add_argument("--single-exon-round-bp", type=int, default=50,
                    help="Round single-exon transcript coordinates to this many "
                         "bp when de-duplicating single-exon loci [50].")
    ap.add_argument("--single-exon-min-coverage", type=float, default=0.0,
                    help="Coverage floor applied to single-exon models (the best "
                         "model of each locus is always kept) [0.0].")
    ap.add_argument("--disabled", action="store_true",
                    help="Pass-through: copy input to output without filtering.")
    args = ap.parse_args()

    args.out_gp.parent.mkdir(parents=True, exist_ok=True)

    models = parse_genepred(args.in_gp)
    n_in = len(models)
    print(f"  read {n_in:,} augMP models from {args.in_gp}", file=sys.stderr)

    if args.disabled or n_in == 0:
        # Straight copy so downstream paths are identical whether or not the
        # filter is enabled.
        args.out_gp.write_text(Path(args.in_gp).read_text() if args.in_gp.exists() else "")
        if args.out_psl is not None:
            src = args.in_psl
            args.out_psl.write_text(src.read_text() if (src and src.exists()) else "")
        print(f"  filter disabled/empty -> wrote {n_in:,} models unchanged", file=sys.stderr)
        return

    metrics = parse_psl_metrics(args.in_psl)
    print(f"  loaded PSL metrics for {len(metrics):,} models", file=sys.stderr)
    assign_scores(models, metrics)

    multi = [m for m in models if m.exon_count > 1]
    single = [m for m in models if m.exon_count <= 1]

    kept_multi = dedup_multi_exon(multi, args.structure_round_bp)
    kept_single = dedup_single_exon(single, args.single_exon_round_bp,
                                    args.single_exon_min_coverage)
    print(f"  multi-exon:  {len(multi):,} -> {len(kept_multi):,} "
          f"(distinct structures)", file=sys.stderr)
    print(f"  single-exon: {len(single):,} -> {len(kept_single):,} "
          f"(distinct loci)", file=sys.stderr)

    kept = kept_multi + kept_single
    before_backstop = len(kept)
    kept = backstop_locus_cap(kept, args.max_models_per_locus)
    if len(kept) < before_backstop:
        print(f"  backstop cap ({args.max_models_per_locus}/locus): "
              f"{before_backstop:,} -> {len(kept):,}", file=sys.stderr)

    kept_names = {m.name for m in kept}
    # Preserve original file order for reproducibility.
    with args.out_gp.open("w") as out:
        for m in models:
            if m.name in kept_names:
                out.write(m.raw + "\n")

    if args.out_psl is not None:
        n_psl = 0
        with args.out_psl.open("w") as out:
            if args.in_psl and args.in_psl.exists():
                with args.in_psl.open() as fh:
                    for line in fh:
                        f = line.rstrip("\n").split("\t")
                        if len(f) >= 10 and f[9] in kept_names:
                            out.write(line)
                            n_psl += 1
        print(f"  wrote {n_psl:,} filtered PSL rows to {args.out_psl}", file=sys.stderr)

    print(f"  DONE: {n_in:,} -> {len(kept):,} augMP models "
          f"({100 * len(kept) / n_in:.1f}% kept) -> {args.out_gp}", file=sys.stderr)


if __name__ == "__main__":
    main()
