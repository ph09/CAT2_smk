#!/usr/bin/env python3
"""NCBI RefSeq protein.faa -> protein-DB FASTA with UniProt-style GN= symbols.

protein_refine maps miniprot queries to gene symbols through ``GN=`` in the
header, which RefSeq headers lack. The symbol comes from the ``gene=``
attribute of the CDS rows in the annotation GFF. One protein per gene (the
longest) is kept, matching the one-canonical-per-gene UniProt proteomes; LOC
symbols are kept, since uncharacterised loci are where lineage-specific genes
sit.

Header written:
    >{taxid}|ref|{XP_acc}|{XP_acc}_{CODE} {product} OS={organism} OX={taxid} GN={gene}
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

ORG_RE = re.compile(r"\[([^\]]+)\]\s*$")


def cds_gene_symbols(gff: Path) -> dict[str, str]:
    """protein_id -> gene symbol from CDS rows."""
    out: dict[str, str] = {}
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            pid, gene = attrs.get("protein_id"), attrs.get("gene")
            if pid and gene:
                out.setdefault(pid, gene)
    return out


def read_fasta(path: Path):
    name, seq = None, []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--protein-faa", type=Path, required=True)
    ap.add_argument("--gff", type=Path, required=True)
    ap.add_argument("--taxid", type=int, required=True)
    ap.add_argument("--code", required=True, help="5-letter species code, e.g. EULRU")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    sym = cds_gene_symbols(args.gff)
    best: dict[str, tuple[str, str, str, str]] = {}
    n_in = n_nosym = 0
    for header, seq in read_fasta(args.protein_faa):
        n_in += 1
        acc, _, desc = header.partition(" ")
        gene = sym.get(acc)
        if not gene:
            n_nosym += 1
            continue
        m = ORG_RE.search(desc)
        org = m.group(1) if m else ""
        product = ORG_RE.sub("", desc).strip()
        if gene not in best or len(seq) > len(best[gene][3]):
            best[gene] = (acc, product, org, seq)

    with args.out.open("w") as fh:
        for gene in sorted(best):
            acc, product, org, seq = best[gene]
            fh.write(
                f">{args.taxid}|ref|{acc}|{acc}_{args.code} {product} "
                f"OS={org} OX={args.taxid} GN={gene}\n"
            )
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
    print(
        f"{args.protein_faa}: {n_in:,} proteins, {n_nosym:,} without a CDS gene symbol, "
        f"{len(best):,} genes written"
    )


if __name__ == "__main__":
    main()
