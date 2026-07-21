#!/usr/bin/env python3
"""
Build a multi-species protein database for CAT2's protein-based gene finding
(miniprot / augMP).

Why this exists
---------------
CAT2 finds genes two complementary ways:

  * transMap projects the *reference* annotation across the HAL alignment. It can
    only recover orthologs of genes that exist in the reference. It structurally
    cannot find a gene that has no reference ortholog (lineage-specific genes, or
    genes lost in the reference lineage but retained elsewhere).
  * augMP aligns a protein database to every target genome with miniprot. This is
    the mode that can discover genes *absent from the reference* -- but only if the
    protein database actually contains proteins from those lineages.

So if you annotate (say) a panprimate alignment using only human proteins, augMP
can never see a non-human, lineage-specific gene. This script assembles a broader
protein set from several species so miniprot has evidence for those genes.

What it does
------------
Given a list of species names and/or NCBI taxon IDs, it:
  1. resolves each to a UniProt *reference* proteome,
  2. downloads that proteome's protein FASTA (cached on disk),
  3. optionally folds in a local ``--base-fasta`` (e.g. proteins derived from your
     reference GFF3),
  4. drops sequences below ``--min-len`` and (by default) collapses byte-identical
     sequences so miniprot does not re-align many near-duplicate copies,
  5. writes a single combined FASTA plus a per-species summary TSV.

Only the Python standard library is used (no requests/biopython) so it runs in any
environment. It talks to the public UniProt REST API, so the machine running it
needs outbound internet access.

Species names are normalised leniently, so HAL genome names paste in directly:
``PR00246~Eulemur_fulvus.pri`` -> ``Eulemur fulvus`` and ``GRCh38_Homo_sapiens``
-> ``Homo sapiens`` (best-effort; give a clean ``Genus species`` or a ``--taxa``
ID if a name does not resolve).

Example
-------
    python scripts/build_protein_db.py \\
        --out panprimate_data/protein_db.fa \\
        --species "Homo sapiens,Macaca mulatta,Callithrix jacchus,Mus musculus" \\
        --base-fasta panprimate_data/GRCh38.prot.fa
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Iterable, Iterator, Optional, Tuple

UNIPROT_REST = "https://rest.uniprot.org"

# UniProt proteomeType values, best first. We prefer a single curated reference
# proteome per species; "Other"/"Redundant" are last-resort fallbacks.
_PROTEOME_TYPE_RANK = {
    "Reference and representative proteome": 0,
    "Reference proteome": 1,
    "Representative proteome": 2,
    "Other proteome": 3,
    "Redundant proteome": 4,
}

# Genome-assembly / haplotype suffixes we strip off HAL-style names before trying
# to read a species name out of them (e.g. "...pri", "...alt", "...hap1").
_ASSEMBLY_SUFFIX_RE = re.compile(
    r"\.(pri|alt|mat|pat|hap\d+|hap|dip|v\d+|\d+)$", re.IGNORECASE
)


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------

def _http_get(url: str, *, timeout: int, retries: int, accept_gzip: bool = True) -> bytes:
    """GET *url*, retrying transient failures with exponential backoff.

    Transparently gunzips gzip-encoded bodies. Raises the last error after
    ``retries`` attempts.
    """
    last_err: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "cat2-build-protein-db"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
                enc = (resp.headers.get("Content-Encoding") or "").lower()
            if accept_gzip and (enc == "gzip" or data[:2] == b"\x1f\x8b"):
                data = gzip.decompress(data)
            return data
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
            last_err = exc
            # 4xx (except 429) are not worth retrying.
            code = getattr(exc, "code", None)
            if code is not None and 400 <= code < 500 and code != 429:
                break
            if attempt < retries:
                delay = min(60.0, 3.0 * (2 ** (attempt - 1)))
                print(f"  [retry {attempt}/{retries}] {url}: {exc}; sleeping {delay:.0f}s",
                      file=sys.stderr)
                time.sleep(delay)
    raise RuntimeError(f"GET failed after {retries} attempt(s): {url}: {last_err}")


# ---------------------------------------------------------------------------
# Name normalisation
# ---------------------------------------------------------------------------

def normalize_species_name(raw: str) -> str:
    """Best-effort clean-up of a species name so common HAL genome names resolve.

    Handles ``PREFIX~Genus_species.pri`` and ``ASM_Genus_species`` shapes. Returns
    a ``Genus species`` string. Non-destructive for already-clean names.
    """
    name = raw.strip()
    if "~" in name:                      # drop assembly/sample prefix, e.g. PR00246~
        name = name.split("~", 1)[1]
    name = _ASSEMBLY_SUFFIX_RE.sub("", name)
    name = name.replace("_", " ").strip()
    # Collapse an assembly token glued to the front, e.g. "GRCh38 Homo sapiens" or
    # "T2T Homo sapiens" -> keep the last two capitalised binomial-looking tokens.
    tokens = name.split()
    if len(tokens) > 2:
        # Find a plausible "Genus species" tail: a Capitalised word followed by a
        # lower-case word. Fall back to the last two tokens.
        for i in range(len(tokens) - 1):
            if tokens[i][:1].isupper() and tokens[i + 1][:1].islower():
                name = " ".join(tokens[i:i + 2])
                break
        else:
            name = " ".join(tokens[-2:])
    return name


# ---------------------------------------------------------------------------
# UniProt proteome resolution
# ---------------------------------------------------------------------------

def _proteome_search(query: str, *, timeout: int, retries: int, size: int = 25) -> list:
    """Run a UniProt proteomes search and return the parsed ``results`` list."""
    url = (
        f"{UNIPROT_REST}/proteomes/search?query="
        f"{urllib.parse.quote(query)}&format=json&size={size}"
    )
    data = _http_get(url, timeout=timeout, retries=retries)
    return json.loads(data.decode("utf-8")).get("results", [])


def _pick_reference(results: list) -> Optional[dict]:
    """Choose the best proteome from search *results* (reference > representative
    > other; ties broken by protein count). Returns the raw record or None."""
    scored = []
    for r in results:
        ptype = r.get("proteomeType", "")
        rank = _PROTEOME_TYPE_RANK.get(ptype, 9)
        pcount = r.get("proteinCount") or 0
        scored.append((rank, -int(pcount), r))
    if not scored:
        return None
    scored.sort(key=lambda t: (t[0], t[1]))
    return scored[0][2]


def resolve_proteome(*, species: Optional[str], taxon: Optional[int],
                     timeout: int, retries: int) -> Optional[dict]:
    """Resolve one species name or taxon id to a chosen proteome record.

    Returns a dict with keys: upid, taxon_id, scientific_name, proteome_type,
    protein_count -- or None if nothing usable was found.
    """
    if taxon is not None:
        query = f"organism_id:{int(taxon)}"
    else:
        query = f'organism_name:"{species}"'
    results = _proteome_search(query, timeout=timeout, retries=retries)
    rec = _pick_reference(results)
    if rec is None:
        return None
    return _rec_to_row(rec, fallback_name=species or str(taxon))


def _busco_score(rec: dict) -> int:
    """BUSCO completeness score (0-100) for a proteome record, or -1 if absent."""
    try:
        return int(rec.get("proteomeCompletenessReport", {})
                   .get("buscoReport", {}).get("score"))
    except (TypeError, ValueError):
        return -1


def _rec_to_row(rec: dict, *, fallback_name: str = "") -> dict:
    """Flatten a raw UniProt proteome record to our summary row shape."""
    tax = rec.get("taxonomy", {}) or {}
    return {
        "upid": rec.get("id"),
        "taxon_id": tax.get("taxonId"),
        "scientific_name": tax.get("scientificName", fallback_name),
        "proteome_type": rec.get("proteomeType", "?"),
        "protein_count": rec.get("proteinCount") or 0,
        "busco": _busco_score(rec),
    }


def resolve_taxon_id(name: str, *, timeout: int, retries: int) -> Optional[Tuple[int, str, str]]:
    """Resolve a clade/organism *name* to ``(taxon_id, scientific_name, rank)``.

    Prefers an exact (case-insensitive) scientific-name match; otherwise takes the
    first result. Returns None if nothing matches.
    """
    url = (f"{UNIPROT_REST}/taxonomy/search?query="
           f"{urllib.parse.quote(name)}&format=json&size=25")
    results = json.loads(_http_get(url, timeout=timeout, retries=retries)
                         .decode("utf-8")).get("results", [])
    if not results:
        return None
    lname = name.strip().lower()
    exact = [r for r in results
             if str(r.get("scientificName", "")).strip().lower() == lname]
    rec = (exact or results)[0]
    tid = rec.get("taxonId")
    if tid is None:
        return None
    return int(tid), rec.get("scientificName", name), rec.get("rank", "")


def proteomes_under_taxon(taxon_id: int, *, max_n: int, include_other: bool,
                          timeout: int, retries: int) -> list:
    """Return the best proteome per organism under a clade *taxon_id* (subtree).

    A ``taxonomy_id:<id>`` proteomes query matches the node and all descendants, so
    for a genus/family this yields every member species that has a proteome. We
    keep one proteome per organism (best proteomeType, then BUSCO, then protein
    count), keep only reference/representative proteomes (unless ``include_other``),
    and return at most ``max_n`` proteomes ordered by quality. Each element is a
    ``_rec_to_row`` dict.
    """
    url = (f"{UNIPROT_REST}/proteomes/search?query="
           f"{urllib.parse.quote(f'(taxonomy_id:{int(taxon_id)})')}"
           f"&format=json&size=500")
    results = json.loads(_http_get(url, timeout=timeout, retries=retries)
                         .decode("utf-8")).get("results", [])
    max_rank = 4 if include_other else 2  # 0-2 = reference/representative
    best_by_org: dict = {}
    for rec in results:
        ptype = rec.get("proteomeType", "")
        rank = _PROTEOME_TYPE_RANK.get(ptype, 9)
        if rank > max_rank:
            continue
        org = (rec.get("taxonomy", {}) or {}).get("taxonId")
        if org is None:
            continue
        key = (rank, -_busco_score(rec), -int(rec.get("proteinCount") or 0))
        prev = best_by_org.get(org)
        if prev is None or key < prev[0]:
            best_by_org[org] = (key, rec)
    ordered = sorted(best_by_org.values(), key=lambda t: t[0])
    return [_rec_to_row(rec) for _key, rec in ordered[:max_n]]


def download_proteome_fasta(upid: str, cache: Path, *, timeout: int, retries: int) -> Path:
    """Download a proteome's protein FASTA to *cache* (cached across runs)."""
    if cache.exists() and cache.stat().st_size > 0:
        print(f"  cached: {cache}", file=sys.stderr)
        return cache
    query = urllib.parse.quote(f"(proteome:{upid})")
    url = f"{UNIPROT_REST}/uniprotkb/stream?query={query}&format=fasta&compressed=true"
    data = _http_get(url, timeout=timeout, retries=retries)
    cache.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache.with_suffix(cache.suffix + ".part")
    tmp.write_bytes(data)
    tmp.replace(cache)
    return cache


# ---------------------------------------------------------------------------
# FASTA handling
# ---------------------------------------------------------------------------

def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield ``(header_without_>, sequence)`` from a (optionally gzipped) FASTA."""
    opener = gzip.open if path.suffix == ".gz" else open
    header = None
    seq_parts: list[str] = []
    with opener(path, "rt") as fh:  # type: ignore[arg-type]
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts)


def _tag_header(header: str, tag: Optional[str]) -> str:
    """Prefix a provenance tag onto a FASTA header for traceability."""
    if not tag:
        return header
    return f"{tag}|{header}"


def write_combined(sources: Iterable[Tuple[Optional[str], Path]], out: Path, *,
                   min_len: int, dedup: bool, line_width: int = 60) -> dict:
    """Merge FASTA *sources* into *out*.

    ``sources`` is an iterable of ``(tag, path)``. Returns per-tag stats
    ``{tag: {"read": n, "kept": n}}``.
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    stats: dict = {}
    n_total_kept = 0
    with open(out, "w") as fo:
        for tag, path in sources:
            key = tag or str(path)
            s = stats.setdefault(key, {"read": 0, "kept": 0})
            for header, seq in iter_fasta(path):
                s["read"] += 1
                seq = seq.replace("*", "").strip()
                if len(seq) < min_len:
                    continue
                if dedup:
                    h = hashlib.sha1(seq.encode("ascii", "ignore")).hexdigest()
                    if h in seen:
                        continue
                    seen.add(h)
                fo.write(">" + _tag_header(header, tag) + "\n")
                for i in range(0, len(seq), line_width):
                    fo.write(seq[i:i + line_width] + "\n")
                s["kept"] += 1
                n_total_kept += 1
    stats["_total_kept"] = n_total_kept
    return stats


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _split_list(values: Optional[list]) -> list[str]:
    """Flatten repeated and/or comma-separated CLI values into a clean list."""
    out: list[str] = []
    for v in values or []:
        out.extend(part.strip() for part in str(v).split(",") if part.strip())
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Build a multi-species protein DB for miniprot/augMP from "
                    "UniProt reference proteomes.")
    ap.add_argument("--out", required=True, help="Output combined protein FASTA.")
    ap.add_argument("--species", action="append", default=[],
                    help="Species names (repeatable and/or comma-separated). "
                         "HAL-style names are accepted and normalised.")
    ap.add_argument("--taxa", action="append", default=[],
                    help="NCBI taxon IDs (repeatable and/or comma-separated). "
                         "Unambiguous; use when a name does not resolve.")
    ap.add_argument("--clades", action="append", default=[],
                    help="Genus/family (or any higher rank) names or taxon IDs "
                         "(repeatable and/or comma-separated). For each clade, one "
                         "proteome per member species is pulled (best BUSCO first, "
                         "capped by --max-per-clade). Use this to get coverage for "
                         "poorly annotated species via their better-annotated "
                         "relatives, e.g. --clades Cercopithecidae,Lemuridae.")
    ap.add_argument("--max-per-clade", type=int, default=25,
                    help="Max proteomes to take per clade, best BUSCO/quality first "
                         "[25]. Guards against pulling dozens of redundant proteomes.")
    ap.add_argument("--clade-include-other", action="store_true",
                    help="When expanding clades, also accept 'Other'/'Redundant' "
                         "proteomes (default: reference/representative only).")
    ap.add_argument("--base-fasta", default=None,
                    help="Optional local protein FASTA to include (e.g. proteins "
                         "derived from your reference GFF3).")
    ap.add_argument("--source", default="uniprot", choices=("uniprot",),
                    help="Protein source (only 'uniprot' is implemented).")
    ap.add_argument("--cache-dir", default=None,
                    help="Directory for cached per-proteome downloads "
                         "(default: <out dir>/_cache).")
    ap.add_argument("--min-len", type=int, default=20,
                    help="Drop proteins shorter than this many residues [20].")
    ap.add_argument("--no-dedup", action="store_true",
                    help="Keep byte-identical duplicate sequences (default: drop).")
    ap.add_argument("--timeout", type=int, default=180, help="HTTP timeout (s) [180].")
    ap.add_argument("--retries", type=int, default=4, help="HTTP retries [4].")
    ap.add_argument("--strict", action="store_true",
                    help="Fail if any species fails to resolve/download "
                         "(default: warn and continue).")
    ap.add_argument("--summary", default=None,
                    help="Path for the per-species summary TSV "
                         "(default: <out>.summary.tsv).")
    args = ap.parse_args()

    species = [normalize_species_name(s) for s in _split_list(args.species)]
    taxa: list[int] = []
    for t in _split_list(args.taxa):
        try:
            taxa.append(int(t))
        except ValueError:
            print(f"WARNING: ignoring non-integer taxon '{t}'", file=sys.stderr)

    clades = _split_list(args.clades)

    if not species and not taxa and not clades and not args.base_fasta:
        ap.error("nothing to do: give --species, --taxa, --clades, and/or --base-fasta")

    out = Path(args.out)
    cache_dir = Path(args.cache_dir) if args.cache_dir else out.parent / "_cache"
    summary_path = Path(args.summary) if args.summary else Path(str(out) + ".summary.tsv")

    print(f"build_protein_db: {len(species)} species, {len(taxa)} taxa, "
          f"{len(clades)} clade(s), base_fasta={'yes' if args.base_fasta else 'no'}",
          file=sys.stderr)

    # Resolve every requested proteome to a flat record, deduplicating by UPID so a
    # species named directly and also reached via its clade is only downloaded once.
    sources: list[Tuple[Optional[str], Path]] = []
    resolved_rows: list[dict] = []
    failures: list[str] = []
    recs_by_upid: dict = {}

    def _add_rec(rec: Optional[dict], label: str) -> None:
        if rec is None or not rec.get("upid"):
            failures.append(label)
            print(f"  no proteome found for {label}", file=sys.stderr)
            return
        recs_by_upid.setdefault(rec["upid"], rec)

    targets: list[Tuple[Optional[str], Optional[int]]] = (
        [(s, None) for s in species] + [(None, t) for t in taxa]
    )
    for sp, tx in targets:
        label = sp if sp is not None else f"taxon:{tx}"
        print(f"resolving {label} ...", file=sys.stderr)
        try:
            rec = resolve_proteome(species=sp, taxon=tx, timeout=args.timeout,
                                   retries=args.retries)
        except Exception as exc:  # noqa: BLE001 - want to keep going per species
            rec = None
            print(f"  ERROR resolving {label}: {exc}", file=sys.stderr)
        _add_rec(rec, label)

    # Clades: expand each to one proteome per member species (best BUSCO first).
    for clade in clades:
        print(f"expanding clade {clade} ...", file=sys.stderr)
        try:
            if clade.isdigit():
                tid, cname = int(clade), f"taxon:{clade}"
            else:
                resolved = resolve_taxon_id(clade, timeout=args.timeout, retries=args.retries)
                if resolved is None:
                    _add_rec(None, f"clade:{clade}")
                    continue
                tid, cname, rank = resolved
                print(f"  -> {cname} (taxon {tid}, rank {rank})", file=sys.stderr)
            recs = proteomes_under_taxon(
                tid, max_n=args.max_per_clade,
                include_other=args.clade_include_other,
                timeout=args.timeout, retries=args.retries)
        except Exception as exc:  # noqa: BLE001
            recs = []
            print(f"  ERROR expanding clade {clade}: {exc}", file=sys.stderr)
        if not recs:
            _add_rec(None, f"clade:{clade}")
            continue
        print(f"  {len(recs)} proteome(s) under {clade}", file=sys.stderr)
        for rec in recs:
            _add_rec(rec, f"clade:{clade}:{rec.get('scientific_name')}")

    # Download every unique proteome.
    for upid, rec in recs_by_upid.items():
        print(f"  -> {upid} ({rec['scientific_name']}, {rec['proteome_type']}, "
              f"{rec['protein_count']} proteins, BUSCO {rec.get('busco', -1)})",
              file=sys.stderr)
        cache = cache_dir / f"{upid}.fasta"
        try:
            download_proteome_fasta(upid, cache, timeout=args.timeout, retries=args.retries)
        except Exception as exc:  # noqa: BLE001
            failures.append(f"{rec['scientific_name']} ({upid})")
            print(f"  ERROR downloading {upid}: {exc}", file=sys.stderr)
            continue
        tag = str(rec.get("taxon_id") or upid)
        sources.append((tag, cache))
        resolved_rows.append({**rec, "cache": str(cache)})

    if failures and args.strict:
        raise SystemExit(f"ERROR: failed to obtain proteomes for: {', '.join(failures)}")

    if args.base_fasta:
        bf = Path(args.base_fasta)
        if not bf.exists():
            raise SystemExit(f"ERROR: --base-fasta not found: {bf}")
        sources.append(("ref", bf))

    if not sources:
        raise SystemExit("ERROR: no protein sources were obtained; nothing written.")

    print(f"merging {len(sources)} source(s) -> {out}", file=sys.stderr)
    stats = write_combined(sources, out, min_len=args.min_len, dedup=not args.no_dedup)

    # Write summary TSV.
    with open(summary_path, "w") as sf:
        sf.write("source\ttaxon_id\tupid\tproteome_type\tbusco\tproteins_downloaded\t"
                 "proteins_read\tproteins_kept\n")
        for row in resolved_rows:
            tag = str(row.get("taxon_id") or row["upid"])
            s = stats.get(tag, {"read": 0, "kept": 0})
            sf.write(f"{row['scientific_name']}\t{row.get('taxon_id','')}\t{row['upid']}\t"
                     f"{row['proteome_type']}\t{row.get('busco','')}\t{row['protein_count']}\t"
                     f"{s['read']}\t{s['kept']}\n")
        if args.base_fasta:
            s = stats.get("ref", {"read": 0, "kept": 0})
            sf.write(f"(base_fasta) {args.base_fasta}\t\t\tlocal\t\t\t"
                     f"{s['read']}\t{s['kept']}\n")

    total = stats.get("_total_kept", 0)
    print(f"DONE: wrote {total:,} proteins to {out}", file=sys.stderr)
    print(f"      summary: {summary_path}", file=sys.stderr)
    if failures:
        print(f"WARNING: {len(failures)} species/taxa did not resolve: "
              f"{', '.join(failures)}", file=sys.stderr)
    if total == 0:
        raise SystemExit("ERROR: 0 proteins written (all filtered?).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
