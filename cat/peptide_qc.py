"""Species-agnostic check: is this translated CDS clearly not a protein?

Alignment frameshifts are normal across species; ``GenePredTranscript``
already skips those indels when translating. Stops that remain are therefore
coordinate/frame errors or processed pseudogenes, not ordinary divergence.
One or two internal stops are allowed (selenocysteine, a leftover codon).
A missing start alone is not enough to drop a model.
"""
from __future__ import annotations

import logging

from tools.bio import get_sequence_dict, translate_sequence

logger = logging.getLogger(__name__)

MIN_INTERNAL_STOPS = 3
# Same floor as protein-only novel CDS in consensus. A StringTie model with
# no annotated CDS is still kept as a transcript; this is only whether it
# also counts as a protein.
MIN_NOVEL_ORF_AA = 100


def recode_selenocysteine(aa: str, cds: str | None = None, *, seleno: bool = False) -> str:
    """Internal in-frame TGA → U, only for genes whose parent protein uses Sec.

    ``seleno`` is true when the UniProt / protein-DB sequence for this symbol
    contains ``U``. A premature TGA in JAK1 or any other non-selenoprotein
    is left as ``*``. TAG/TAA and the terminal stop always stay ``*``.
    """
    if not aa or not seleno or not cds:
        return aa or ""
    cds = str(cds).upper().replace("U", "T")
    out = list(aa)
    last = len(out) - 1
    for i, res in enumerate(out):
        if res != "*":
            continue
        if i == last:
            continue
        if cds[i * 3 : (i + 1) * 3] == "TGA":
            out[i] = "U"
    return "".join(out)


def protein_from_transcript(tx_obj, seq_dict, *, seleno: bool = False) -> str:
    """Translate CDS; recode internal TGA as Sec only when ``seleno``."""
    aa = tx_obj.get_protein_sequence(seq_dict) or ""
    if not aa:
        return ""
    get_cds = getattr(tx_obj, "get_cds", None)
    if get_cds is None or not seleno:
        return aa
    try:
        cds = get_cds(seq_dict, ignore_frameshift=True)
    except TypeError:
        cds = get_cds(seq_dict)
    return recode_selenocysteine(aa, cds, seleno=True)


def clearly_broken_peptide(aa: str | None) -> bool:
    """True when a translated CDS is not a plausible protein.

    Frame-aware translation already skips alignment indels, so three or more
    leftover internal stops mean the coordinates or frame are wrong — not
    ordinary cross-species gaps. One or two stops are kept (selenocysteine,
    a residual codon). A missing start by itself is not a drop. An empty
    string is not judged (no CDS to score).
    """
    if not aa:
        return False
    seq = aa[:-1] if aa.endswith("*") else aa
    return bool(seq) and seq.count("*") >= MIN_INTERNAL_STOPS


def coding_attrs(attrs: dict) -> bool:
    bt = str(attrs.get("gene_biotype") or attrs.get("transcript_biotype") or "")
    if bt == "protein_coding":
        return True
    return str(attrs.get("transcript_class") or "").startswith("putative_novel")


def best_orf_peptide(mrna: str, min_aa: int = MIN_NOVEL_ORF_AA) -> str:
    """Longest ATG-to-stop ORF across three frames of a spliced mRNA.

    StringTie models have no CDS. Translating the whole transcript from nt 0
    is almost always 5' UTR in the wrong frame, which looks like hundreds of
    stops. A real novel isoform still has an ORF in the spliced RNA — that is
    the CDS we test. Returns ``''`` if nothing reaches ``min_aa``; the
    transcript is kept, it just is not emitted as a protein.
    """
    if not mrna:
        return ""
    seq = str(mrna).upper()
    best = ""
    best_n = 0
    for frame in range(3):
        aa = translate_sequence(seq[frame:])
        i = 0
        while True:
            m = aa.find("M", i)
            if m < 0:
                break
            stop = aa.find("*", m)
            if stop < 0:
                orf = aa[m:]
                i = m + 1
            else:
                orf = aa[m : stop + 1]
                i = stop + 1
            n = len(orf) - (1 if orf.endswith("*") else 0)
            if n >= min_aa and n > best_n:
                best, best_n = orf, n
    return best


def load_seleno_symbols(protein_fasta: str | None) -> set[str]:
    """Symbols whose protein-DB peptides contain U (UniProt Sec)."""
    out: set[str] = set()
    if not protein_fasta:
        return out
    from pathlib import Path
    import re
    gn_re = re.compile(r"(?:^|\s)GN=([^\s;]+)")
    if not Path(protein_fasta).exists():
        return out
    cur = None
    buf: list[str] = []

    def _flush() -> None:
        if cur and "U" in "".join(buf).upper():
            out.add(cur)

    with open(protein_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                _flush()
                header = line[1:].strip()
                gn = gn_re.search(header)
                cur = gn.group(1) if gn else None
                buf = []
            else:
                buf.append(line.strip())
        _flush()
    return out


def _symbol_of(attrs: dict) -> str:
    for key in ("source_gene_common_name", "source_gene", "gene_name"):
        val = str(attrs.get(key) or "").strip()
        if val and val not in {"N/A", "NA", "None", "nan", "."}:
            return val
    return ""


def peptide_for_transcript(tx_obj, attrs: dict, seq_dict, seleno_symbols: set[str] | None = None) -> str:
    if getattr(tx_obj, "cds_size", 0) > 0:
        return protein_from_transcript(
            tx_obj, seq_dict, seleno=_symbol_of(attrs) in (seleno_symbols or ())
        )
    if str(attrs.get("transcript_class") or "").startswith("putative_novel"):
        return best_orf_peptide(tx_obj.get_mrna(seq_dict))
    return ""


def drop_clearly_broken_coding(final_consensus, tx_dict, fasta: str, seleno_symbols: set[str] | None = None):
    """Remove protein-coding transcripts whose annotated CDS is clearly junk.

    StringTie / putative_novel models with no CDS are not dropped — they are
    the source of novel isoforms. Their protein record is decided later by
    ``best_orf_peptide``. Returns ``(kept, n_dropped)``.
    """
    if not fasta or not final_consensus:
        return list(final_consensus), 0
    seq_dict = get_sequence_dict(fasta)
    kept = []
    n_drop = 0
    for tx, attrs in final_consensus:
        tx_obj = tx_dict.get(tx)
        if tx_obj is None or not coding_attrs(attrs) or getattr(tx_obj, "cds_size", 0) <= 0:
            kept.append((tx, attrs))
            continue
        try:
            aa = peptide_for_transcript(tx_obj, attrs, seq_dict, seleno_symbols)
        except (ValueError, KeyError, AssertionError, RuntimeError):
            kept.append((tx, attrs))
            continue
        if clearly_broken_peptide(aa):
            n_drop += 1
            continue
        kept.append((tx, attrs))
    if n_drop:
        logger.info("  dropped %s clearly broken coding transcripts", n_drop)
    return kept, n_drop
