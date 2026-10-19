"""Reference GFF3 normalisation to the GENCODE-style layout CAT2 reads.

CAT2 reads the reference through gff3ToGenePred -attrsOut, keyed by feature ID,
and needs on every transcript: gene_id, gene_name, gene_biotype|gene_type,
transcript_id, transcript_name, transcript_biotype|transcript_type, with
ID == transcript_id and Parent == gene_id. GENCODE GFF3 has all of that.

  ensembl  Native Ensembl GFF3 (any species, ftp.ensembl.org): strips the
           gene:/transcript:/CDS: ID prefixes, copies gene name/biotype onto
           transcripts, biotype -> gene_biotype/transcript_biotype, collapses
           ncRNA_gene/pseudogene -> gene and every transcript type -> transcript.
  refseq   NCBI RefSeq / Gnomon GFF3 (GCF_*_genomic.gff): genes keyed by NCBI
           GeneID, transcripts by accession, biotypes mapped to Ensembl names,
           readthrough / MANE tags carried over, optional seqid renaming.

Duplicate IDs at a second locus (PAR copies, alt scaffolds) get a _locN suffix:
it survives CAT2's version (.N / -N) and txTM-copy (_N) stripping.

usage:
  python -m cat.reference_gff3 detect  IN.gff3[.gz]
  python -m cat.reference_gff3 ensembl IN.gff3[.gz] OUT.gff3
  python -m cat.reference_gff3 refseq  IN.gff[.gz]  OUT.gff3 [--assembly-report R] [--chrom-map M] ...
"""
import argparse
import collections
import gzip
import re
import sys
from urllib.parse import unquote

SUBFEATURES = {"exon", "CDS", "five_prime_UTR", "three_prime_UTR", "start_codon", "stop_codon", "UTR"}
REQUIRED_TX = ("gene_id", "gene_name", "transcript_id", "transcript_name")


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_attrs(s):
    d = collections.OrderedDict()
    for kv in s.strip().rstrip(";").split(";"):
        if not kv:
            continue
        k, _, v = kv.strip().partition("=")
        d[k] = v
    return d


def format_attrs(d):
    return ";".join(f"{k}={v}" for k, v in d.items() if v is not None and v != "")


def esc(v):
    """GFF3-escape a value built from unescaped text."""
    return str(v).replace("%", "%25").replace(";", "%3B").replace("=", "%3D").replace("&", "%26").replace(",", "%2C")


def iter_features(path):
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) == 9:
                yield f


# ---------------------------------------------------------------------------
# format detection


def detect_format(path, max_lines=200000):
    """'gencode' (usable as is), 'ensembl' (native Ensembl), or 'refseq' (NCBI)."""
    n = 0
    for f in iter_features(path):
        n += 1
        if n > max_lines:
            break
        a = parse_attrs(f[8])
        if f[2] in SUBFEATURES or "Parent" not in a:
            continue
        if a["ID"].startswith("rna-") or a["Parent"].startswith("gene-") or ("GeneID:" in a.get("Dbxref", "") and "gene_id" not in a):
            return "refseq"
        if a["ID"].startswith("transcript:") or a["Parent"].startswith("gene:"):
            return "ensembl"
        has_bt = ("gene_biotype" in a or "gene_type" in a) and ("transcript_biotype" in a or "transcript_type" in a)
        if all(k in a for k in REQUIRED_TX) and has_bt:
            return "gencode"
        return "ensembl" if "biotype" in a else "unknown"
    return "unknown"


# ---------------------------------------------------------------------------
# shared writer


class Model:
    __slots__ = ("seqid", "source", "start", "end", "strand", "attrs", "children")

    def __init__(self, seqid, source, start, end, strand, attrs):
        self.seqid, self.source, self.start, self.end, self.strand = seqid, source, int(start), int(end), strand
        self.attrs = attrs
        self.children = []          # gene: transcripts; transcript: (type, start, end, phase, attrs)


def write_gff3(genes, out):
    """genes: list of gene Models whose children are transcript Models."""
    with open(out, "w") as fh:
        fh.write("##gff-version 3\n")
        for g in sorted(genes, key=lambda g: (g.seqid, g.start, g.end)):
            if not g.children:
                continue
            g.start = min(g.start, min(t.start for t in g.children))
            g.end = max(g.end, max(t.end for t in g.children))
            fh.write("\t".join([g.seqid, g.source, "gene", str(g.start), str(g.end), ".", g.strand, ".",
                                format_attrs(g.attrs)]) + "\n")
            for t in sorted(g.children, key=lambda t: (t.start, t.end)):
                fh.write("\t".join([t.seqid, t.source, "transcript", str(t.start), str(t.end), ".", t.strand, ".",
                                    format_attrs(t.attrs)]) + "\n")
                for typ, s, e, ph, a in sorted(t.children, key=lambda c: (c[1], c[0] != "exon")):
                    fh.write("\t".join([t.seqid, t.source, typ, str(s), str(e), ".", t.strand, ph,
                                        format_attrs(a)]) + "\n")
            fh.write("###\n")


# ---------------------------------------------------------------------------
# Ensembl


def _strip_prefix(x):
    return x.split(":", 1)[1] if ":" in x and x.split(":", 1)[0] in ("gene", "transcript", "CDS", "exon") else x


def normalize_ensembl(inp, out, log=sys.stderr):
    genes, tx = collections.OrderedDict(), {}
    subs = []
    n = collections.Counter()
    for f in iter_features(inp):
        a = parse_attrs(f[8])
        if f[2] in SUBFEATURES:
            subs.append((f, a))
        elif "Parent" not in a:
            if "gene_id" in a or a.get("ID", "").startswith("gene:"):
                gid = a.get("gene_id") or _strip_prefix(a["ID"])
                bt = a.get("gene_biotype") or a.get("gene_type") or a.get("biotype") or "unknown"
                ga = collections.OrderedDict(ID=gid, gene_id=gid, gene_name=a.get("gene_name") or a.get("Name") or gid,
                                             gene_biotype=bt)
                for k, v in a.items():
                    if k not in ("ID", "Name", "biotype", "gene_id", "gene_name", "gene_biotype", "gene_type"):
                        ga[k] = v
                genes[gid] = Model(f[0], f[1], f[3], f[4], f[6], ga)
            else:
                n[f"dropped top-level {f[2]}"] += 1
        else:
            pid = _strip_prefix(a["Parent"].split(",")[0])
            tid = a.get("transcript_id") or _strip_prefix(a["ID"])
            tx[tid] = (pid, f, a)
    for tid, (pid, f, a) in tx.items():
        g = genes.get(pid)
        if g is None:
            n["transcript without gene"] += 1
            continue
        bt = a.get("transcript_biotype") or a.get("transcript_type") or a.get("biotype") or g.attrs["gene_biotype"]
        ta = collections.OrderedDict(ID=tid, Parent=g.attrs["gene_id"], gene_id=g.attrs["gene_id"],
                                     transcript_id=tid, gene_name=g.attrs["gene_name"],
                                     gene_biotype=g.attrs["gene_biotype"], transcript_biotype=bt,
                                     transcript_name=a.get("transcript_name") or a.get("Name") or tid)
        if "description" in g.attrs:
            ta["description"] = g.attrs["description"]
        for k, v in a.items():
            if k not in ta and k not in ("Name", "biotype", "transcript_type", "gene_type"):
                ta[k] = v
        m = Model(f[0], f[1], f[3], f[4], f[6], ta)
        g.children.append(m)
        tx[tid] = m
    for f, a in subs:
        for p in a["Parent"].split(","):
            m = tx.get(_strip_prefix(p))
            if not isinstance(m, Model):
                n[f"{f[2]} without transcript"] += 1
                continue
            sa = collections.OrderedDict(Parent=m.attrs["transcript_id"])
            if f[2] == "CDS":
                sa["ID"] = "CDS:" + m.attrs["transcript_id"]
            for k in ("exon_id", "protein_id", "rank"):
                if k in a:
                    sa[k] = a[k]
            m.children.append((f[2], int(f[3]), int(f[4]), f[7], sa))
    write_gff3(list(genes.values()), out)
    nt = sum(len(g.children) for g in genes.values())
    print(f"ensembl -> CAT2 GFF3: {sum(1 for g in genes.values() if g.children):,} genes, {nt:,} transcripts"
          + "".join(f"; {k} {v:,}" for k, v in sorted(n.items())), file=log)
    return n


# ---------------------------------------------------------------------------
# RefSeq

# NCBI gene_biotype -> Ensembl gene biotype (IG/TR resolved from the symbol)
REFSEQ_GENE_BIOTYPE = {
    "protein_coding": "protein_coding", "lncRNA": "lncRNA", "antisense_RNA": "lncRNA",
    "pseudogene": "pseudogene", "transcribed_pseudogene": "transcribed_pseudogene",
    "miRNA": "miRNA", "snoRNA": "snoRNA", "snRNA": "snRNA", "scaRNA": "scaRNA", "rRNA": "rRNA",
    "tRNA": "tRNA", "misc_RNA": "misc_RNA", "ncRNA": "ncRNA", "vault_RNA": "vault_RNA", "Y_RNA": "Y_RNA",
    "RNase_P_RNA": "misc_RNA", "RNase_MRP_RNA": "misc_RNA", "SRP_RNA": "misc_RNA",
    "telomerase_RNA": "misc_RNA", "guide_RNA": "misc_RNA", "ncRNA_pseudogene": "pseudogene",
    "other": "other",
}
SEGMENT = {"V_segment": "V", "D_segment": "D", "J_segment": "J", "C_region": "C",
           "V_segment_pseudogene": "V", "D_segment_pseudogene": "D", "J_segment_pseudogene": "J",
           "C_region_pseudogene": "C"}
TX_TYPE_BIOTYPE = {"lnc_RNA": "lncRNA", "antisense_RNA": "lncRNA", "primary_transcript": "miRNA",
                   "snoRNA": "snoRNA", "snRNA": "snRNA", "scaRNA": "scaRNA", "rRNA": "rRNA", "tRNA": "tRNA",
                   "vault_RNA": "vault_RNA", "Y_RNA": "Y_RNA", "RNase_P_RNA": "misc_RNA",
                   "RNase_MRP_RNA": "misc_RNA", "SRP_RNA": "misc_RNA", "telomerase_RNA": "misc_RNA"}
REFSEQ_GENE_TYPES = {"gene", "pseudogene"}
ACC_RE = re.compile(r"^[A-Z]{2}_\d+(\.\d+)?$")
KEEP_GENE_ATTRS = ("description", "gene_synonym", "Dbxref", "partial", "exception")
KEEP_TX_ATTRS = ("product", "Dbxref", "partial", "exception", "model_evidence", "inference")


def refseq_gene_biotype(raw, symbol, pseudo):
    if raw in SEGMENT:
        locus = "TR" if symbol.upper().startswith("TR") else "IG"
        return f"{locus}_{SEGMENT[raw]}_{'pseudogene' if 'pseudogene' in raw else 'gene'}"
    bt = REFSEQ_GENE_BIOTYPE.get(raw, raw or "other")
    if pseudo and "pseudogene" not in bt:
        bt = "pseudogene"
    return bt


def read_seq_map(assembly_report=None, chrom_map=None, style="ucsc",
                 drop_roles=("alt-scaffold", "fix-patch", "novel-patch")):
    """RefSeq accession -> output seqid; None means keep seqids as they are."""
    if chrom_map:
        m = {}
        for line in _open(chrom_map):
            p = line.split()
            if len(p) >= 2 and not line.startswith("#"):
                m[p[0]] = p[1]
        return m
    if not assembly_report:
        return None
    col = {"ucsc": 9, "genbank": 4, "refseq": 6, "sequence": 0}[style]
    m = {}
    for line in _open(assembly_report):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n\r").split("\t")
        if len(p) < 10 or p[1] in drop_roles:
            continue
        name = p[col] if p[col] not in ("na", "") else p[6]
        for key in (p[6], p[4], p[0]):
            if key and key != "na":
                m[key] = name
    return m


def _uniq(base, seen):
    seen[base] += 1
    return base if seen[base] == 1 else f"{base}_loc{seen[base]}"


def convert_refseq(inp, out, seq_map=None, curated_only=False, log=sys.stderr):
    feats = collections.OrderedDict()     # ID -> [type, rows, attrs, parent]
    anon = 0
    n = collections.Counter()
    for f in iter_features(inp):
        a = parse_attrs(f[8])
        fid = a.get("ID")
        if fid is None:
            anon += 1
            fid = f"__anon{anon}"
        rec = feats.get(fid)
        if rec is None:
            feats[fid] = [f[2], [f], a, a.get("Parent", "").split(",")[0] or None]
        else:
            rec[1].append(f)
    kids = collections.defaultdict(list)
    for fid, rec in feats.items():
        if rec[3]:
            kids[rec[3]].append(fid)

    gene_seen, tx_seen = collections.Counter(), collections.Counter()
    genes = []
    for fid, (typ, rows, a, parent) in feats.items():
        if parent is not None or typ not in REFSEQ_GENE_TYPES:
            if parent is None:
                n[f"dropped top-level {typ}"] += 1
            continue
        g0 = rows[0]
        seqid = g0[0]
        if seq_map is not None:
            if seqid not in seq_map:
                n["genes on dropped/unmapped sequences"] += 1
                continue
            seqid = seq_map[seqid]
        if "trans-splicing" in a.get("exception", "") or len(rows) > 1:
            n["trans-spliced / multi-row genes dropped"] += 1
            continue
        xref = dict(x.split(":", 1) for x in unquote(a.get("Dbxref", "")).split(",") if ":" in x)
        symbol = unquote(a.get("gene") or a.get("Name") or xref.get("GeneID") or fid)
        gbt = refseq_gene_biotype(a.get("gene_biotype", ""), symbol, a.get("pseudo") == "true")
        gid = _uniq(xref.get("GeneID") or re.sub(r"^gene-", "", fid), gene_seen)
        desc = unquote(a.get("description", ""))
        readthrough = "readthrough" in desc.lower()
        ga = collections.OrderedDict(ID=gid, gene_id=gid, gene_name=esc(symbol), gene_biotype=gbt)
        if readthrough:
            ga["tag"] = "readthrough_gene"
        for k in KEEP_GENE_ATTRS:
            if k in a:
                ga[k] = a[k]
        ga["ncbi_gene_id"] = xref.get("GeneID", "")
        g = Model(seqid, g0[1], g0[3], g0[4], g0[6], ga)

        direct = [c for c in kids[fid] if feats[c][0] in ("exon", "CDS")]
        txs = [c for c in kids[fid] if feats[c][0] not in SUBFEATURES]
        units = [(c, feats[c]) for c in txs]
        if not txs and (direct or gbt != "protein_coding"):
            # gene-direct exon/CDS, or a single-line pseudogene / ncRNA gene: one transcript
            units = [(None, ("transcript", [g0], a, fid))]
        for cid, (ttyp, trows, ta0, _) in units:
            if len(trows) > 1:
                n["multi-row transcripts dropped"] += 1
                continue
            sub = [feats[x] for x in (kids[cid] if cid else direct)]
            if cid and any(feats[x][0] not in SUBFEATURES for x in kids[cid]):
                n["nested transcripts (mature miRNA etc.) dropped"] += sum(feats[x][0] not in SUBFEATURES for x in kids[cid])
            exons = [(int(r[3]), int(r[4])) for s in sub if s[0] == "exon" for r in s[1]]
            cds = [(int(r[3]), int(r[4]), r[7]) for s in sub if s[0] == "CDS" for r in s[1]]
            if not exons:
                exons = [(s, e) for s, e, _ in cds] or ([(int(g0[3]), int(g0[4]))] if cid is None else [])
            if not exons or any(r[0] != g0[0] or r[6] != g0[6] for s in sub for r in s[1]):
                n["transcripts without exons or off-locus dropped"] += 1
                continue
            raw_tid = ta0.get("transcript_id") or (re.sub(r"^rna-", "", cid) if cid else "")
            if curated_only and (trows[0][1] == "Gnomon" or raw_tid[:3] in ("XM_", "XR_", "XP_")):
                n["Gnomon models dropped (--curated-only)"] += 1
                continue
            if not ACC_RE.match(raw_tid):
                raw_tid = f"{gid}_t"            # symbolic NCBI names (trnaa-agc) collide under version stripping
            tid = _uniq(raw_tid, tx_seen)
            if cds:
                tbt = "protein_coding" if gbt == "protein_coding" else gbt
            elif gbt == "protein_coding":
                tbt = "processed_transcript"
            elif "pseudogene" in gbt or gbt.startswith(("IG_", "TR_")):
                tbt = gbt
            else:
                tbt = TX_TYPE_BIOTYPE.get(ttyp, gbt)
            tags = [t.strip().replace(" ", "_") for t in unquote(ta0.get("tag", "")).split(",") if t.strip()]
            if readthrough:
                tags.append("readthrough_transcript")
            tname = unquote(ta0.get("Name", "")) if cid else ""
            ta = collections.OrderedDict(ID=tid, Parent=gid, gene_id=gid, transcript_id=tid, gene_name=esc(symbol),
                                         gene_biotype=gbt, transcript_biotype=tbt,
                                         transcript_name=esc(f"{symbol}-{tname}" if tname else tid))
            if tags:
                ta["tag"] = ",".join(esc(t) for t in tags)
            if desc:
                ta["description"] = a["description"]
            for k in KEEP_TX_ATTRS:
                if k in ta0:
                    ta[k] = ta0[k]
            pid = next((s[2].get("protein_id") for s in sub if s[0] == "CDS" and s[2].get("protein_id")), None)
            if pid:
                ta["protein_id"] = pid
            r0 = trows[0]
            t = Model(seqid, r0[1], min(s for s, _ in exons), max(e for _, e in exons), g0[6], ta)
            t.children = [("exon", s, e, ".", collections.OrderedDict(Parent=tid)) for s, e in sorted(set(exons))]
            t.children += [("CDS", s, e, ph, collections.OrderedDict(ID=f"CDS:{tid}", Parent=tid))
                           for s, e, ph in sorted(set(cds))]
            g.children.append(t)
            n[f"transcripts {tbt}"] += 1
        if g.children:
            genes.append(g)
            n[f"genes {gbt}"] += 1
        else:
            n["genes without usable transcripts dropped"] += 1
    write_gff3(genes, out)
    print(f"refseq -> CAT2 GFF3: {len(genes):,} genes, {sum(len(g.children) for g in genes):,} transcripts", file=log)
    for k, v in sorted(n.items()):
        print(f"  {k:55s} {v:>9,}", file=log)
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    d = sp.add_parser("detect")
    d.add_argument("gff3")
    e = sp.add_parser("ensembl")
    e.add_argument("gff3")
    e.add_argument("out")
    r = sp.add_parser("refseq")
    r.add_argument("gff3")
    r.add_argument("out")
    r.add_argument("--assembly-report", help="NCBI *_assembly_report.txt, renames seqids (see --seq-names)")
    r.add_argument("--seq-names", default="ucsc", choices=("ucsc", "genbank", "refseq", "sequence"),
                   help="assembly-report column to use as seqid (default ucsc: chr1...; na falls back to refseq)")
    r.add_argument("--keep-alt", action="store_true", help="keep alt-scaffold / patch sequences (dropped by default)")
    r.add_argument("--chrom-map", help="2-column file: RefSeq seqid -> genome FASTA seqid (overrides the report)")
    r.add_argument("--curated-only", action="store_true", help="drop Gnomon (XM_/XR_) models, keep NM_/NR_")
    args = ap.parse_args(argv)
    if args.cmd == "detect":
        print(detect_format(args.gff3))
    elif args.cmd == "ensembl":
        normalize_ensembl(args.gff3, args.out)
    else:
        seq_map = read_seq_map(args.assembly_report, args.chrom_map, args.seq_names,
                               () if args.keep_alt else ("alt-scaffold", "fix-patch", "novel-patch"))
        convert_refseq(args.gff3, args.out, seq_map, args.curated_only)


if __name__ == "__main__":
    main()
