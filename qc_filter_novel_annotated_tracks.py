#!/usr/bin/env python3
"""QC + filter novel-annotated consensus GFF3/GP, then (re)build bigGenePred tracks.

Drops clear chimeric / artifactual transcripts while keeping real large genes
(DMD, CNTNAP2, RBFOX1, …):

  - genomic span > 5 Mb
  - max intron > 2 Mb
  - max intron > 1.5 Mb AND exon_bp < 3 kb
  - exon count > 500
  - malformed exon coordinates

Writes per-genome:
  cleaned/*.gp, cleaned/*.gff3
  cleaned/QC_ARTIFACTS.tsv, cleaned/QC_SUMMARY.md
and rebuilds browser_biggenepred/*.bb from the cleaned GP.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

# Reuse autosql / conversion helpers from the track builder
sys.path.insert(0, str(Path(__file__).resolve().parent))
from generate_novel_annotated_biggenepred import (  # noqa: E402
    BIGGENEPRED_AS,
    DEFAULT_SIZE_ROOTS,
    find_chrom_sizes,
    genome_from_stem,
    gp_to_biggenepred_bed,
    load_gp_info,
    trackdb_stanza,
    which,
)

MAX_SPAN = 5_000_000
MAX_INTRON = 2_000_000
SPARSE_INTRON = 1_500_000
SPARSE_EXON_BP = 3_000
MAX_EXONS = 500


def parse_gp_row(line: str) -> dict | None:
    if not line.strip() or line.startswith("#"):
        return None
    f = line.rstrip("\n").split("\t")
    if len(f) < 10:
        return None
    try:
        tx0, tx1 = int(f[3]), int(f[4])
        n_ex = int(f[7])
        starts = [int(x) for x in f[8].rstrip(",").split(",") if x != ""]
        ends = [int(x) for x in f[9].rstrip(",").split(",") if x != ""]
    except ValueError:
        return {
            "fields": f,
            "name": f[0],
            "name2": f[11] if len(f) > 11 else f[0],
            "chrom": f[1] if len(f) > 1 else "?",
            "span": 0,
            "max_intron": 0,
            "n_ex": 0,
            "exon_bp": 0,
            "reasons": ["parse_error"],
            "line": line,
        }
    reasons: list[str] = []
    span = tx1 - tx0
    if len(starts) != len(ends) or n_ex != len(starts):
        reasons.append("exon_count_mismatch")
    exon_bp = 0
    max_intron = 0
    if not reasons:
        for i, (s, e) in enumerate(zip(starts, ends)):
            if e <= s:
                reasons.append("zero_or_neg_exon")
                break
            if i and s < ends[i - 1]:
                reasons.append("overlapping_exons")
                break
            exon_bp += e - s
        if not reasons:
            for i in range(len(starts) - 1):
                max_intron = max(max_intron, starts[i + 1] - ends[i])
    if span > MAX_SPAN:
        reasons.append(f"span>{MAX_SPAN}")
    if max_intron > MAX_INTRON:
        reasons.append(f"intron>{MAX_INTRON}")
    if max_intron > SPARSE_INTRON and exon_bp < SPARSE_EXON_BP:
        reasons.append(f"intron>{SPARSE_INTRON}&exon_bp<{SPARSE_EXON_BP}")
    if n_ex > MAX_EXONS:
        reasons.append(f"n_ex>{MAX_EXONS}")
    return {
        "fields": f,
        "name": f[0],
        "name2": f[11] if len(f) > 11 else f[0],
        "chrom": f[1],
        "span": span,
        "max_intron": max_intron,
        "n_ex": n_ex,
        "exon_bp": exon_bp,
        "reasons": reasons,
        "line": line if line.endswith("\n") else line + "\n",
    }


def filter_gp(
    gp_path: Path,
) -> tuple[list[str], list[dict], Counter, set[str]]:
    """Return keep_lines, dropped rows, reason counts, fully-dropped gene IDs."""
    keep_lines: list[str] = []
    dropped: list[dict] = []
    reason_counts: Counter = Counter()
    gene_txs: dict[str, set[str]] = defaultdict(set)
    dropped_tx: set[str] = set()
    with gp_path.open() as fh:
        for line in fh:
            row = parse_gp_row(line)
            if row is None:
                continue
            gene_txs[row["name2"]].add(row["name"])
            if row["reasons"]:
                dropped.append(row)
                dropped_tx.add(row["name"])
                for r in row["reasons"]:
                    reason_counts[r] += 1
            else:
                keep_lines.append(row["line"])
    drop_genes = {
        g for g, txs in gene_txs.items() if txs and txs <= dropped_tx
    }
    return keep_lines, dropped, reason_counts, drop_genes


def filter_gff3(
    gff3_path: Path, drop_tx: set[str], drop_genes: set[str], out_path: Path
) -> tuple[int, int]:
    """Drop GFF3 features for rejected transcripts / fully-dropped genes.

    Fast path: if nothing to drop, copy the file. Otherwise one compiled
    regex over attribute IDs (drop sets are typically <500).
    """
    import re

    if not drop_tx and not drop_genes:
        if out_path.exists() or out_path.is_symlink():
            out_path.unlink()
        out_path.symlink_to(gff3_path.resolve())
        return 0, 0

    ids = sorted(drop_tx | drop_genes, key=len, reverse=True)
    alt = "|".join(re.escape(i) for i in ids)
    # Match ID/Parent/transcript_id/gene_id values (optional gene:/transcript: prefix)
    pat = re.compile(
        rf"(?:^|\t|;)(?:ID|Parent|transcript_id|gene_id)="
        rf"(?:gene:|transcript:)?(?:{alt})(?:;|$|\t|\n)"
    )

    n_in = n_out = 0
    with gff3_path.open() as inf, out_path.open("w") as outf:
        for line in inf:
            if line.startswith("#") or not line.strip():
                outf.write(line)
                continue
            n_in += 1
            if pat.search(line):
                continue
            outf.write(line)
            n_out += 1
    return n_in, n_out


def build_bb_from_gp(
    gp: Path,
    genome: str,
    out_bb: Path,
    sizes: Path,
    tools: dict[str, str],
    gp_info: Path | None = None,
) -> None:
    out_bb.parent.mkdir(parents=True, exist_ok=True)
    info_by_tx = load_gp_info(gp_info) if gp_info and gp_info.is_file() else {}
    with tempfile.TemporaryDirectory(prefix="bgp_") as tmp:
        tmp_path = Path(tmp)
        as_path = tmp_path / "bigGenePred.as"
        as_path.write_text(BIGGENEPRED_AS)
        bed = tmp_path / "out.bed"
        gp_to_biggenepred_bed(gp, info_by_tx, bed)
        subprocess.run([tools["bedSort"], str(bed), str(bed)], check=True)
        subprocess.run(
            [
                tools["bedToBigBed"],
                "-type=bed12+8",
                "-tab",
                f"-as={as_path}",
                "-extraIndex=name2,geneName,geneName2",
                str(bed),
                str(sizes),
                str(out_bb),
            ],
            check=True,
        )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--gff3-dir",
        type=Path,
        default=Path(
            "/private/groups/cgl/pnhebbar/cat2/panprimate_out_combined/consensus_gene_set"
        ),
    )
    ap.add_argument(
        "--cleaned-dir",
        type=Path,
        default=None,
        help="Default: <gff3-dir>/cleaned_for_browser",
    )
    ap.add_argument(
        "--track-dir",
        type=Path,
        default=None,
        help="Default: <gff3-dir>/browser_biggenepred",
    )
    ap.add_argument("--size-root", type=Path, action="append", default=None)
    ap.add_argument("--skip-tracks", action="store_true")
    ap.add_argument(
        "--write-gff3",
        action="store_true",
        help="Also write cleaned GFF3 (slow on NFS; default: cleaned GP + tracks only)",
    )
    ap.add_argument("--genome", action="append", default=None)
    args = ap.parse_args()

    gff3_dir = args.gff3_dir
    cleaned_dir = args.cleaned_dir or (gff3_dir / "cleaned_for_browser")
    track_dir = args.track_dir or (gff3_dir / "browser_biggenepred")
    cleaned_dir.mkdir(parents=True, exist_ok=True)
    size_roots = args.size_root or DEFAULT_SIZE_ROOTS

    gps = sorted(gff3_dir.glob("*_consensus_novel_annotated.gp"))
    if args.genome:
        want = set(args.genome)
        gps = [p for p in gps if genome_from_stem(p.name[: -len(".gp")]) in want]
    if not gps:
        raise SystemExit("No GP files found")

    tools = None
    if not args.skip_tracks:
        tools = {
            "bedToBigBed": which("bedToBigBed"),
            "bedSort": which("bedSort"),
        }

    artifact_rows: list[dict] = []
    summary_rows: list[dict] = []
    total_in = total_drop = 0

    for gp in gps:
        genome = genome_from_stem(gp.name[: -len(".gp")])
        stem = gp.name[: -len(".gp")]
        gff3 = gff3_dir / f"{stem}.gff3"
        print(f"[qc] {genome}", file=sys.stderr)

        keep_lines, dropped, reason_counts, drop_genes = filter_gp(gp)
        n_in = len(keep_lines) + len(dropped)
        total_in += n_in
        total_drop += len(dropped)

        out_gp = cleaned_dir / gp.name
        with out_gp.open("w") as fh:
            fh.writelines(keep_lines)

        drop_tx = {r["name"] for r in dropped}

        # Always record drop IDs for optional later GFF3 filtering
        drop_id_path = cleaned_dir / f"{stem}.drop_ids.txt"
        drop_id_path.write_text("\n".join(sorted(drop_tx | drop_genes)) + ("\n" if drop_tx or drop_genes else ""))

        n_gff_in = n_gff_out = 0
        if args.write_gff3 and gff3.is_file():
            out_gff3 = cleaned_dir / gff3.name
            n_gff_in, n_gff_out = filter_gff3(gff3, drop_tx, drop_genes, out_gff3)
            if n_gff_in == 0 and n_gff_out == 0 and not drop_tx:
                print("  gff3: copied (no drops)", file=sys.stderr)
            else:
                print(f"  gff3: wrote cleaned ({n_gff_out:,}/{n_gff_in:,} lines)", file=sys.stderr)
        elif gff3.is_file():
            print(
                f"  gff3: skipped rewrite (see {drop_id_path.name}; pass --write-gff3)",
                file=sys.stderr,
            )
        else:
            print(f"  warning: missing {gff3.name}", file=sys.stderr)

        for r in dropped:
            artifact_rows.append(
                {
                    "genome": genome,
                    "transcript_id": r["name"],
                    "gene_id": r["name2"],
                    "chrom": r["chrom"],
                    "span": r["span"],
                    "max_intron": r["max_intron"],
                    "n_exons": r["n_ex"],
                    "exon_bp": r["exon_bp"],
                    "reasons": ";".join(r["reasons"]),
                }
            )

        summary_rows.append(
            {
                "genome": genome,
                "tx_in": n_in,
                "tx_kept": len(keep_lines),
                "tx_dropped": len(dropped),
                "genes_fully_dropped": len(drop_genes),
                "gff3_lines_in": n_gff_in,
                "gff3_lines_out": n_gff_out,
                "top_reasons": ",".join(f"{k}:{v}" for k, v in reason_counts.most_common(5)),
            }
        )
        print(
            f"  kept {len(keep_lines):,}/{n_in:,}  dropped {len(dropped)}  "
            f"genes_gone={len(drop_genes)}",
            file=sys.stderr,
        )

        if tools is not None:
            sizes = find_chrom_sizes(genome, size_roots)
            bb = track_dir / f"{stem}.bb"
            gp_info = gff3_dir / f"{stem}.gp_info"
            print(f"  [bb] {bb.name}", file=sys.stderr)
            build_bb_from_gp(out_gp, genome, bb, sizes, tools, gp_info=gp_info)
            (track_dir / f"{stem}.trackDb.txt").write_text(
                trackdb_stanza(genome, bb.name)
            )

    # write QC tables
    art_path = cleaned_dir / "QC_ARTIFACTS.tsv"
    with art_path.open("w", newline="") as fh:
        fields = [
            "genome",
            "transcript_id",
            "gene_id",
            "chrom",
            "span",
            "max_intron",
            "n_exons",
            "exon_bp",
            "reasons",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in sorted(artifact_rows, key=lambda r: (-r["span"], r["genome"])):
            w.writerow(row)

    sum_path = cleaned_dir / "QC_SUMMARY.tsv"
    with sum_path.open("w", newline="") as fh:
        fields = list(summary_rows[0])
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)

    md = [
        "# Novel-annotated GFF3/GP artifact QC\n\n",
        "Filters (transcript-level; real large genes like DMD/CNTNAP2/RBFOX1 kept):\n\n",
        f"- genomic span > {MAX_SPAN:,} bp\n",
        f"- max intron > {MAX_INTRON:,} bp\n",
        f"- max intron > {SPARSE_INTRON:,} bp **and** exon_bp < {SPARSE_EXON_BP:,}\n",
        f"- exon count > {MAX_EXONS}\n",
        "- malformed exon coordinates\n\n",
        f"**Dropped {total_drop:,} / {total_in:,} transcripts "
        f"({100 * total_drop / total_in:.3f}%).**\n\n",
        "Cleaned GFF3/GP: this directory. Tracks rebuilt in "
        f"`{track_dir.name}/` from cleaned GP.\n\n",
        "| genome | tx in | kept | dropped | genes fully dropped |\n",
        "|---|---:|---:|---:|---:|\n",
    ]
    for r in summary_rows:
        md.append(
            f"| {r['genome']} | {r['tx_in']} | {r['tx_kept']} | "
            f"{r['tx_dropped']} | {r['genes_fully_dropped']} |\n"
        )
    md.append("\nSee `QC_ARTIFACTS.tsv` for every dropped transcript.\n")
    (cleaned_dir / "QC_SUMMARY.md").write_text("".join(md))

    if tools is not None and track_dir.is_dir():
        parts = sorted(track_dir.glob("*.trackDb.txt"))
        # exclude aggregate if re-run
        parts = [p for p in parts if p.name != "trackDb.all.txt"]
        (track_dir / "trackDb.all.txt").write_text(
            "".join(p.read_text() for p in parts)
        )

    print(
        f"Done. Dropped {total_drop:,}/{total_in:,}. "
        f"QC → {cleaned_dir}  tracks → {track_dir}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
