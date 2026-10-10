#!/usr/bin/env python3
"""Build UCSC bigGenePred (.bb) browser tracks from novel-annotated consensus GFF3/GP.

For each ``*_consensus_novel_annotated.gff3`` (paired ``.gp`` + ``.gp_info``), writes:
  - ``{stem}.bb``          bigGenePred bigBed, itemRgb colored by gene class
  - ``{stem}.bgp.bed``     sorted bed12+8 (optional keep)
  - ``{stem}.trackDb.txt`` trackDb stanza (itemRgb on)

Display name (BED ``name``) prefers a human-readable symbol:
  gene/transcript symbol → ``GENE-like`` for novel paralogs → Ensembl ID → CAT gene ID.

Colors (itemRgb):
  protein-coding  blue    76,85,212
  non-coding RNA  green   85,212,76
  pseudogene      pink    255,105,180
  predicted       purple  135,76,212

Chrom sizes are resolved from ``*/genome_files/{genome}.chrom.sizes`` under
the provided search roots (default: panprimate_out{,_pr,_apes,_combined}).
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

BIGGENEPRED_AS = '''table bigGenePred
"bigGenePred gene models"
(
string chrom;       "Reference sequence chromosome or scaffold"
uint   chromStart;  "Start position in chromosome"
uint   chromEnd;    "End position in chromosome"
string name;        "Name. Typically geneName2 for genes or ENST for transcripts"
uint score;         "Score of gene prediction (0-1000)"
char[1] strand;     "+ or - for strand"
uint thickStart;    "Coding region start"
uint thickEnd;      "Coding region end"
uint reserved;      "RGB color string"
int blockCount;     "Number of blocks"
int[blockCount] blockSizes; "Comma separated list of block sizes"
int[blockCount] chromStarts; "Start positions relative to chromStart"
string name2;       "Alternative/human readable name (geneName for genes)"
string cdsStartStat; "Status of CDS start annotation (none, unknown, incomplete, or complete)"
string cdsEndStat;   "Status of CDS end annotation (none, unknown, incomplete, or complete)"
int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame for exon"
string type;        "Transcript type"
string geneName;    "Primary identifier for gene"
string geneName2;   "Alternative/human readable gene name"
string geneType;    "Gene type"
)
'''

DEFAULT_SIZE_ROOTS = [
    Path("/private/groups/cgl/pnhebbar/cat2/panprimate_out"),
    Path("/private/groups/cgl/pnhebbar/cat2/panprimate_out_pr"),
    Path("/private/groups/cgl/pnhebbar/cat2/panprimate_out_apes"),
    Path("/private/groups/cgl/pnhebbar/cat2/panprimate_out_combined"),
]

# Browser RGB (R,G,B). Distinct, colorblind-friendly-ish on white.
COLOR_CODING = "76,85,212"       # blue
COLOR_NCRNA = "85,212,76"        # green
COLOR_PSEUDO = "255,105,180"     # pink
COLOR_PREDICTED = "135,76,212"   # purple

CODING_GENE_BIOTYPES = {
    "protein_coding",
    "mRNA",
    "IG_C_gene",
    "IG_D_gene",
    "IG_J_gene",
    "IG_V_gene",
    "TR_C_gene",
    "TR_D_gene",
    "TR_J_gene",
    "TR_V_gene",
}
PREDICTED_GENE_BIOTYPES = {
    "unknown_likely_coding",
    "TEC",
    "fragment",
}
PREDICTED_TRANSCRIPT_CLASSES = {
    "putative_novel",
}
NA_VALUES = {"", "N/A", "NA", "None", "none", ".", "nan", "NaN"}

# miniprot / Augustus protein IDs like 109500.t1
_PROTEIN_TX_RE = re.compile(r"^\d+\.t\d+$")
_ENSEMBL_RE = re.compile(r"^ENS[A-Z]*[GTP]\d+", re.IGNORECASE)
_CAT_GENE_RE = re.compile(r"(?:^|_)G\d{5,}$")
_CAT_TX_RE = re.compile(r"(?:^|_)T\d{5,}$")

_STANDALONES = Path(__file__).resolve().parent / "standalones"


def which(name: str) -> str:
    if _STANDALONES.is_dir():
        os.environ["PATH"] = f"{_STANDALONES}:{os.environ.get('PATH', '')}"
    path = shutil.which(name)
    if not path:
        raise SystemExit(
            f"Required tool not on PATH: {name}\n"
            "Activate the cat conda env or add kent binaries to PATH."
        )
    return path


def genome_from_stem(stem: str) -> str:
    """T2T_Homo_sapiens.pri_consensus_novel_annotated → T2T_Homo_sapiens.pri"""
    for suffix in (
        "_consensus_novel_annotated",
        "_consensus",
        "_novel_annotated",
    ):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def species_short(genome: str) -> str:
    """PR00030~Macaca_nemestrina.pri → Macaca_nemestrina (UCSC shortLabel ≤17)."""
    s = genome.split("~", 1)[-1]
    s = s.replace(".pri", "").replace("T2T_", "")
    return s[:17]


def find_chrom_sizes(genome: str, roots: list[Path]) -> Path:
    for root in roots:
        cand = root / "genome_files" / f"{genome}.chrom.sizes"
        if cand.is_file():
            return cand
    raise FileNotFoundError(
        f"No chrom.sizes for {genome} under {[str(r) for r in roots]}"
    )


def ensure_gp(gff3: Path, gp: Path, gff3_to_gp: str) -> Path:
    if gp.is_file() and gp.stat().st_mtime >= gff3.stat().st_mtime:
        return gp
    print(f"  [gp] converting {gff3.name} → {gp.name}", file=sys.stderr)
    subprocess.run(
        [
            gff3_to_gp,
            "-rnaNameAttr=transcript_id",
            "-geneNameAttr=gene_id",
            "-honorStartStopCodons",
            str(gff3),
            str(gp),
        ],
        check=True,
    )
    return gp


def trackdb_stanza(genome: str, bb_name: str) -> str:
    short = species_short(genome)
    return (
        f"track {genome}_novel_annotated\n"
        f"shortLabel {short} CAT\n"
        f"longLabel {genome} CAT consensus "
        f"(blue=coding, green=ncRNA, pink=pseudogene, purple=predicted)\n"
        f"type bigGenePred\n"
        f"bigDataUrl {bb_name}\n"
        f"visibility pack\n"
        f"itemRgb on\n"
        f"searchIndex name,name2,geneName,geneName2\n"
        f"labelFields name,name2,geneName\n"
        f"defaultLabelFields name\n"
    )


def _blank(val: object) -> bool:
    return val is None or str(val).strip() in NA_VALUES


def is_human_readable(name: str | None) -> bool:
    """True if *name* looks like a gene/transcript symbol, not an opaque ID."""
    if name is None or _blank(name):
        return False
    s = str(name).strip()
    if _ENSEMBL_RE.match(s):
        return False
    if _PROTEIN_TX_RE.match(s):
        return False
    if _CAT_GENE_RE.search(s) or _CAT_TX_RE.search(s):
        return False
    return True


def _paralog_symbol(description: str | None) -> str | None:
    if description is None or _blank(description):
        return None
    m = re.match(r"^paralog of (.+)$", str(description).strip(), re.IGNORECASE)
    if not m:
        return None
    gene = m.group(1).strip()
    return gene if gene and is_human_readable(gene) else None


def display_names(info: dict, tx_id: str, gene_id: str) -> tuple[str, str]:
    """Return (bed_name, gene_symbol) for browser display.

    bed_name is the main label (isoform symbol when available).
    gene_symbol is the gene-level readable name (name2 / geneName2).
    """
    tx_sym = str(info.get("source_transcript_name") or "").strip()
    gene_sym = str(info.get("source_gene_common_name") or "").strip()
    src_gene = str(info.get("source_gene") or "").strip()
    tx_class = str(info.get("transcript_class") or "").strip()
    novel_class = str(info.get("novel_class") or "").strip()
    gene_from_para = _paralog_symbol(info.get("novel_gene_description"))
    is_novel = (
        tx_class in PREDICTED_TRANSCRIPT_CLASSES
        or novel_class in {"paralog", "lineage_specific"}
    )

    if is_human_readable(gene_sym):
        gene_label = gene_sym
    elif gene_from_para:
        gene_label = f"{gene_from_para}-like"
    elif not _blank(gene_sym):
        gene_label = gene_sym  # ENSG… still better than CAT IDs
    elif not _blank(src_gene):
        gene_label = src_gene.split(".")[0]
    elif is_novel:
        m = re.search(r"(G\d+)$", gene_id)
        gene_label = f"novel-{m.group(1)}" if m else gene_id
    else:
        gene_label = gene_id

    if is_human_readable(tx_sym):
        bed_name = tx_sym
    elif gene_from_para:
        bed_name = f"{gene_from_para}-like"
    elif is_human_readable(gene_sym):
        bed_name = gene_sym
    elif not _blank(tx_sym) and not _PROTEIN_TX_RE.match(tx_sym):
        bed_name = tx_sym
    else:
        bed_name = gene_label

    return _sanitize_name(bed_name), _sanitize_name(gene_label)


def _sanitize_name(name: str) -> str:
    """BED name: no tabs/spaces; keep symbols like HLA-A, BRCA2-206."""
    s = re.sub(r"[\t\n\r]+", " ", name).strip()
    s = s.replace(" ", "_")
    s = s.replace(",", "_")
    return s[:255] if s else "unknown"


def color_class(info: dict) -> tuple[str, str]:
    """Return (rgb, class_label) for a transcript.

    Predicted (purple) wins over coding so novel miniprot/augustus genes
    are visually distinct from lifted orthologs.
    """
    gene_bt = str(info.get("gene_biotype") or "").strip()
    tx_bt = str(info.get("transcript_biotype") or "").strip()
    tx_class = str(info.get("transcript_class") or "").strip()

    is_pseudo = "pseudogene" in gene_bt.lower() or "pseudogene" in tx_bt.lower()
    is_predicted = (
        tx_class in PREDICTED_TRANSCRIPT_CLASSES
        or gene_bt in PREDICTED_GENE_BIOTYPES
        or tx_bt in PREDICTED_GENE_BIOTYPES
    )

    # Classified pseudogenes stay pink even if they were called putative_novel.
    if is_pseudo:
        return COLOR_PSEUDO, "pseudogene"
    if is_predicted:
        return COLOR_PREDICTED, "predicted"
    if gene_bt in CODING_GENE_BIOTYPES:
        return COLOR_CODING, "protein_coding"
    return COLOR_NCRNA, "ncRNA"


def load_gp_info(path: Path) -> dict[str, dict]:
    """transcript_id → selected gp_info fields."""
    keep = {
        "transcript_id",
        "gene_id",
        "gene_biotype",
        "transcript_biotype",
        "transcript_class",
        "source_gene_common_name",
        "source_transcript_name",
        "novel_gene_description",
        "novel_class",
        "source_gene",
    }
    out: dict[str, dict] = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tx = row.get("transcript_id")
            if not tx:
                continue
            out[tx] = {k: row.get(k, "") for k in keep}
    return out


def _int_list(field: str) -> list[int]:
    return [int(x) for x in field.rstrip(",").split(",") if x != ""]


def gp_to_biggenepred_bed(
    gp: Path,
    info_by_tx: dict[str, dict],
    bed: Path,
) -> Counter:
    """Write a bed12+8 bigGenePred file. Returns counts by color class."""
    counts: Counter = Counter()
    n_missing_info = 0
    with gp.open() as inf, bed.open("w") as outf:
        for line in inf:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            tx_id = f[0]
            chrom = f[1]
            strand = f[2]
            chrom_start = int(f[3])
            chrom_end = int(f[4])
            thick_start = int(f[5])
            thick_end = int(f[6])
            block_count = int(f[7])
            exon_starts = _int_list(f[8])
            exon_ends = _int_list(f[9])
            score = int(f[10]) if len(f) > 10 and f[10] not in NA_VALUES else 0
            gene_id = f[11] if len(f) > 11 else tx_id
            cds_start_stat = f[12] if len(f) > 12 else "none"
            cds_end_stat = f[13] if len(f) > 13 else "none"
            if len(f) > 14 and f[14].strip():
                exon_frames = _int_list(f[14])
            else:
                exon_frames = [-1] * block_count

            if len(exon_starts) != block_count or len(exon_ends) != block_count:
                continue
            block_sizes = [e - s for s, e in zip(exon_starts, exon_ends)]
            chrom_starts = [s - chrom_start for s in exon_starts]
            if len(exon_frames) != block_count:
                exon_frames = (exon_frames + [-1] * block_count)[:block_count]

            info = info_by_tx.get(tx_id, {})
            if not info:
                n_missing_info += 1
                info = {
                    "gene_id": gene_id,
                    "gene_biotype": "",
                    "transcript_biotype": "",
                    "transcript_class": "",
                }
            rgb, klass = color_class(info)
            bed_name, gene_label = display_names(info, tx_id, gene_id)
            tx_type = str(info.get("transcript_biotype") or klass)
            gene_type = str(info.get("gene_biotype") or klass)
            gene_name = str(info.get("gene_id") or gene_id)
            score = max(0, min(1000, score))

            row = [
                chrom,
                str(chrom_start),
                str(chrom_end),
                bed_name,
                str(score),
                strand,
                str(thick_start),
                str(thick_end),
                rgb,
                str(block_count),
                ",".join(str(x) for x in block_sizes) + ",",
                ",".join(str(x) for x in chrom_starts) + ",",
                tx_id,  # name2: unique CAT transcript ID (searchable)
                cds_start_stat,
                cds_end_stat,
                ",".join(str(x) for x in exon_frames) + ",",
                tx_type,
                gene_name,
                gene_label,
                gene_type,
            ]
            outf.write("\t".join(row) + "\n")
            counts[klass] += 1
    if n_missing_info:
        print(
            f"  warning: {n_missing_info} transcripts had no gp_info row",
            file=sys.stderr,
        )
    return counts


def resolve_gp(gff3: Path, prefer_cleaned: bool) -> Path:
    stem = gff3.name[: -len(".gff3")] if gff3.name.endswith(".gff3") else gff3.stem
    cleaned = gff3.parent / "cleaned_for_browser" / f"{stem}.gp"
    if prefer_cleaned and cleaned.is_file():
        return cleaned
    gp = gff3.with_suffix(".gp")
    if not gp.exists():
        gp = gff3.parent / f"{stem}.gp"
    return gp


def convert_one(
    gff3: Path,
    out_dir: Path,
    size_roots: list[Path],
    keep_bed: bool,
    tools: dict[str, str],
    prefer_cleaned: bool = True,
) -> Path:
    stem = gff3.name[: -len(".gff3")] if gff3.name.endswith(".gff3") else gff3.stem
    genome = genome_from_stem(stem)
    gp = resolve_gp(gff3, prefer_cleaned)
    gp = ensure_gp(gff3, gp, tools["gff3ToGenePred"])
    gp_info = gff3.parent / f"{stem}.gp_info"
    if not gp_info.is_file():
        raise FileNotFoundError(f"Missing gp_info for {gff3.name}: {gp_info}")
    sizes = find_chrom_sizes(genome, size_roots)

    out_dir.mkdir(parents=True, exist_ok=True)
    bb = out_dir / f"{stem}.bb"
    bed_out = out_dir / f"{stem}.bgp.bed"
    trackdb = out_dir / f"{stem}.trackDb.txt"

    print(f"  [info] {genome}: loading {gp_info.name}", file=sys.stderr)
    info_by_tx = load_gp_info(gp_info)
    src = "cleaned GP" if "cleaned_for_browser" in str(gp) else "GP"
    print(f"  [bed]  {genome}: coloring {src} ({gp.name})", file=sys.stderr)

    with tempfile.TemporaryDirectory(prefix="bgp_") as tmp:
        tmp_path = Path(tmp)
        as_path = tmp_path / "bigGenePred.as"
        as_path.write_text(BIGGENEPRED_AS)
        bed = tmp_path / "out.bed"
        counts = gp_to_biggenepred_bed(gp, info_by_tx, bed)
        subprocess.run([tools["bedSort"], str(bed), str(bed)], check=True)
        print(f"  [bb]   {genome}: bedToBigBed → {bb.name}", file=sys.stderr)
        subprocess.run(
            [
                tools["bedToBigBed"],
                "-type=bed12+8",
                "-tab",
                f"-as={as_path}",
                "-extraIndex=name2,geneName,geneName2",
                str(bed),
                str(sizes),
                str(bb),
            ],
            check=True,
        )
        if keep_bed:
            shutil.copy2(bed, bed_out)

    trackdb.write_text(trackdb_stanza(genome, bb.name))
    n = sum(counts.values())
    parts = ", ".join(f"{k}={v:,}" for k, v in counts.most_common())
    print(f"  [ok]   {n:,} items  {parts}", file=sys.stderr)
    return bb


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--gff3-dir",
        type=Path,
        default=Path(
            "/private/groups/cgl/pnhebbar/cat2/panprimate_out_combined/consensus_gene_set"
        ),
        help="Directory containing *_consensus_novel_annotated.gff3",
    )
    ap.add_argument(
        "--pattern",
        default="*_consensus_novel_annotated.gff3",
        help="Glob under --gff3-dir (default: *_consensus_novel_annotated.gff3)",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: <gff3-dir>/browser_biggenepred)",
    )
    ap.add_argument(
        "--size-root",
        type=Path,
        action="append",
        default=None,
        help="Root containing genome_files/*.chrom.sizes (repeatable)",
    )
    ap.add_argument("--keep-bed", action="store_true", help="Also write .bgp.bed")
    ap.add_argument(
        "--genome",
        action="append",
        default=None,
        help="Restrict to genome id(s) (repeatable)",
    )
    ap.add_argument(
        "--no-cleaned",
        action="store_true",
        help="Use the raw GP next to the GFF3 instead of cleaned_for_browser/*.gp",
    )
    args = ap.parse_args()

    tools = {
        "bedToBigBed": which("bedToBigBed"),
        "bedSort": which("bedSort"),
        "gff3ToGenePred": which("gff3ToGenePred"),
    }
    size_roots = args.size_root or DEFAULT_SIZE_ROOTS
    gffs = sorted(args.gff3_dir.glob(args.pattern))
    if args.genome:
        want = set(args.genome)
        gffs = [g for g in gffs if genome_from_stem(g.name[: -len(".gff3")]) in want]
    if not gffs:
        raise SystemExit(f"No GFF3s matching {args.pattern} in {args.gff3_dir}")

    out_dir = args.out_dir or (args.gff3_dir / "browser_biggenepred")
    print(f"Converting {len(gffs)} novel-annotated GFF3(s) → {out_dir}", file=sys.stderr)
    ok = 0
    for gff3 in gffs:
        try:
            bb = convert_one(
                gff3,
                out_dir,
                size_roots,
                args.keep_bed,
                tools,
                prefer_cleaned=not args.no_cleaned,
            )
            print(f"[ok] {bb}", file=sys.stderr)
            ok += 1
        except Exception as e:
            print(f"[FAIL] {gff3.name}: {e}", file=sys.stderr)
            raise
    # aggregate trackDb
    parts = sorted(out_dir.glob("*.trackDb.txt"))
    parts = [p for p in parts if p.name != "trackDb.all.txt"]
    if parts:
        (out_dir / "trackDb.all.txt").write_text("".join(p.read_text() for p in parts))
    print(f"Done: {ok}/{len(gffs)} bigGenePred tracks", file=sys.stderr)


if __name__ == "__main__":
    main()
