"""Classify miniprot-only models as protein_coding vs unknown_likely_coding.

A protein-to-genome alignment can be a real gene. It can also be a retrocopy,
a domain hit, or a fragment of a gene already annotated in this genome. This
module keeps those models in the annotation and splits them on evidence.


  1. transcribed      RNA/IsoSeq support, >=3 exons, >=100 aa
  2. missing gene     paralog of a gene NOT already projected here,
                      >=3 exons, >=100 aa
  3. local CNV        paralog of a gene that IS projected, copy on the SAME
                      chromosome as a projected copy, >=4 exons, >=150 aa
  4. interchrom CNV   paralog of a projected gene on a *different* chromosome,
                      but near-full-length (>=70% of parent exons and peptide,
                      and the local-CNV floors). Real non-local segdups
                      (gorilla JMJD7-PLA2G4B at the chr1 inversion vs ancestral
                      HSA15; chimpanzee extras on other chroms).
  5. strong de novo   no reference homolog (lineage_specific / unknown),
                      >=5 exons, >=200 aa
  6. otherwise        unknown_likely_coding

"""
from __future__ import annotations

from typing import Optional

PROTEIN_CODING = "protein_coding"
UNKNOWN_LIKELY_CODING = "unknown_likely_coding"

# Thresholds. Named so configs / tests can import them.
RNA_MIN_EXONS = 3
RNA_MIN_AA = 100
MISSING_MIN_EXONS = 3
MISSING_MIN_AA = 100
LOCAL_CNV_MIN_EXONS = 4
LOCAL_CNV_MIN_AA = 150
INTERCHROM_CNV_MIN_EXON_FRAC = 0.70
INTERCHROM_CNV_MIN_AA_FRAC = 0.70
DENOVO_MIN_EXONS = 5
DENOVO_MIN_AA = 200


def has_rna_support(attrs: dict) -> bool:
    """True if IsoSeq / RNA flags on a gp_info-like dict indicate transcription."""
    pacbio = attrs.get("pacbio_isoform_supported")
    if pacbio is True or str(pacbio).lower() in {"1", "true"}:
        return True
    for key in ("exon_rna_support", "intron_rna_support"):
        val = attrs.get(key, "")
        if val is True:
            return True
        if isinstance(val, str) and "1" in val.split(","):
            return True
    return False


def parent_symbol_from_description(desc: str) -> str:
    """'paralog of PCCA' -> 'PCCA'; otherwise ''."""
    d = (desc or "").strip()
    if d.lower().startswith("paralog of "):
        return d[11:].split()[0].strip().upper()
    return ""


def copy_is_complete(
    n_exons: int,
    cds_aa: int,
    parent_n_exons: Optional[int] = None,
    parent_cds_aa: Optional[int] = None,
) -> bool:
    """Near-full-length relative to the projected parent.

    Requires the local-CNV floors plus >=70% of parent exons and peptide.
    Returns False when parent structure is unknown so we do not promote
    interchromosomal fragments by default.
    """
    n_exons = int(n_exons or 0)
    cds_aa = int(cds_aa or 0)
    if n_exons < LOCAL_CNV_MIN_EXONS or cds_aa < LOCAL_CNV_MIN_AA:
        return False
    p_ex = int(parent_n_exons or 0)
    p_aa = int(parent_cds_aa or 0)
    if p_ex <= 0 and p_aa <= 0:
        return False
    if p_ex > 0 and n_exons < int(INTERCHROM_CNV_MIN_EXON_FRAC * p_ex):
        return False
    if p_aa > 0 and cds_aa < int(INTERCHROM_CNV_MIN_AA_FRAC * p_aa):
        return False
    return True


def call_miniprot_novel(
    n_exons: int,
    cds_aa: int,
    has_rna: bool = False,
    novel_class: Optional[str] = None,
    parent_already_projected: Optional[bool] = None,
    same_chrom_as_parent: Optional[bool] = None,
    parent_n_exons: Optional[int] = None,
    parent_cds_aa: Optional[int] = None,
    as_coding: bool = False,
) -> tuple[str, str]:
    """Return ``(biotype, reason)`` for one miniprot-only gene.

    ``novel_class`` is ``'paralog'``, ``'lineage_specific'``, or None when
    DIAMOND has not been run yet. ``parent_already_projected`` /
    ``same_chrom_as_parent`` / parent sizes are None in that case.

    ``as_coding`` forces protein_coding (high_recall). Models are still
    labeled with a reason so they can be audited.
    """
    n_exons = int(n_exons or 0)
    cds_aa = int(cds_aa or 0)
    ncl = (novel_class or "").strip().lower()
    if ncl in {"n/a", "na", ""}:
        ncl = ""

    if as_coding:
        return PROTEIN_CODING, "high_recall"

    if has_rna and n_exons >= RNA_MIN_EXONS and cds_aa >= RNA_MIN_AA:
        return PROTEIN_CODING, "transcribed"

    if ncl == "paralog" and parent_already_projected is False:
        if n_exons >= MISSING_MIN_EXONS and cds_aa >= MISSING_MIN_AA:
            return PROTEIN_CODING, "missing_gene"

    if ncl == "paralog" and parent_already_projected is True:
        if (
            same_chrom_as_parent
            and n_exons >= LOCAL_CNV_MIN_EXONS
            and cds_aa >= LOCAL_CNV_MIN_AA
        ):
            return PROTEIN_CODING, "local_cnv"
        if same_chrom_as_parent is False and copy_is_complete(
            n_exons, cds_aa, parent_n_exons, parent_cds_aa
        ):
            return PROTEIN_CODING, "interchrom_cnv"
        return UNKNOWN_LIKELY_CODING, (
            "dispersed_paralog" if same_chrom_as_parent is False else "redundant_paralog"
        )

    # lineage_specific, empty class, or consensus-time (class unknown):
    # only a long multi-exon model is called coding without a named parent.
    if ncl in ("", "lineage_specific") or ncl is None:
        if n_exons >= DENOVO_MIN_EXONS and cds_aa >= DENOVO_MIN_AA:
            return PROTEIN_CODING, "strong_de_novo"

    if n_exons <= 1:
        return UNKNOWN_LIKELY_CODING, "single_exon"
    if n_exons == 2:
        return UNKNOWN_LIKELY_CODING, "two_exon"
    return UNKNOWN_LIKELY_CODING, "weak_structure"


if __name__ == "__main__":
    # Local tandem-like extra (KZNF).
    bt, reason = call_miniprot_novel(
        7, 400, novel_class="paralog",
        parent_already_projected=True, same_chrom_as_parent=True,
        parent_n_exons=8, parent_cds_aa=420,
    )
    assert (bt, reason) == (PROTEIN_CODING, "local_cnv"), (bt, reason)

    # Gorilla ancestral HSA15 PLA2G4B: complete, other chrom than chr1 expansion.
    bt, reason = call_miniprot_novel(
        22, 729, novel_class="paralog",
        parent_already_projected=True, same_chrom_as_parent=False,
        parent_n_exons=21, parent_cds_aa=782,
    )
    assert (bt, reason) == (PROTEIN_CODING, "interchrom_cnv"), (bt, reason)

    # Truncated ancestral JMJD7 fragment: stay unknown.
    bt, reason = call_miniprot_novel(
        4, 113, novel_class="paralog",
        parent_already_projected=True, same_chrom_as_parent=False,
        parent_n_exons=6, parent_cds_aa=135,
    )
    assert (bt, reason) == (UNKNOWN_LIKELY_CODING, "dispersed_paralog"), (bt, reason)

    # PCCA domain hits: incomplete, stay unknown.
    bt, reason = call_miniprot_novel(
        12, 278, novel_class="paralog",
        parent_already_projected=True, same_chrom_as_parent=False,
        parent_n_exons=25, parent_cds_aa=729,
    )
    assert (bt, reason) == (UNKNOWN_LIKELY_CODING, "dispersed_paralog"), (bt, reason)

    # Missing parent stats: do not promote interchromosomal copies.
    bt, reason = call_miniprot_novel(
        22, 729, novel_class="paralog",
        parent_already_projected=True, same_chrom_as_parent=False,
    )
    assert (bt, reason) == (UNKNOWN_LIKELY_CODING, "dispersed_paralog"), (bt, reason)

    print("ok")
