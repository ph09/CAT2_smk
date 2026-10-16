"""Close-relative evidence copies and protein-coding precision rules.

Runs after annotate_novel_genes on its gp / gp_info / GFF3 and rewrites them.

Evidence loci
    miniprot alignments of close-relative proteins (``--close-taxids``) to this
    genome that are full length (query coverage >= 90%), >= 90% identical and
    free of frameshifts and stops. Overlapping same-strand alignments are one
    locus. A locus is a retrocopy when it is intronless while its proteins or
    gene symbols align with introns elsewhere in the genome.

Evidence copies (a locus with no protein_coding gene covering >= 50% of it)
    The intact non-coding model already at the locus is promoted, otherwise the
    intact miniprot_consensus model is inserted and the non-coding models it
    overlaps are removed. A candidate must cover >= 80% of the locus, the locus
    >= 80% of its CDS, and share <= 20% of its CDS with any protein_coding gene.
    Skipped: transposon-like or uncharacterized evidence, retrocopies, and loci
    the projected reference calls non-coding (pseudogene, lncRNA, TR/IG, TEC;
    unitary pseudogenes excepted).

Precision rules (genes not projected from a reference protein-coding gene)
    R1 no intact coding transcript             -> unknown_likely_coding
    R2 projected from a reference non-coding gene -> that gene's biotype
    R3 IsoSeq novel isoform sharing CDS with a reference-sourced gene
                                               -> isoform of that gene
    R4 novel augMP/miniprot gene, CDS > 50% soft-masked, no evidence locus
       matching it both ways >= 80%            -> unknown_likely_coding
    R5 protein_refine copy named by symbol, < 70% of that symbol's longest
       projected protein                       -> unknown_likely_coding

Evidence extensions (protein_coding genes projected from a reference gene)
    A gene whose longest protein is < 90% of its reference gene's longest
    protein gets, as an extra isoform, an intact miniprot model built from a
    full-length close-relative alignment of the same gene symbol that is >= 10%
    longer and at most 5% longer than the reference protein (close-relative
    predictions often carry spurious extensions). The model must cover >= 50%
    of the gene's CDS and share <= 20% of its CDS with any other protein_coding
    gene.

Demoted genes stay in the annotation; the reason is written to
miniprot_novel_call. Evidence copies are labelled 'evidence_copy', added
isoforms 'evidence_extension'.
"""
from __future__ import annotations

import argparse
import bisect
import collections
import csv
import logging
import re
from pathlib import Path

from cat.consensus_parallel import _merged_cds_blocks, _shared_bases
from cat.convert_miniprot_to_genepred import paf_identity, parse_cigar_to_exons, parse_paf_row
from tools.bio import get_sequence_dict
from tools.transcripts import GenePredTranscript

logger = logging.getLogger(__name__)

PROTEIN_CODING = "protein_coding"
UNKNOWN_LIKELY_CODING = "unknown_likely_coding"
PROJ_MODES = ("transMap", "transMap_pairwise", "augTM", "augTM_pairwise")

EV_MIN_QCOV = 0.9
EV_MIN_PID = 0.9
COPY_MIN_AA = 100
COPY_LOCUS_COV = 0.8
COPY_MODEL_COV = 0.8
COPY_MAX_PC_SHARE = 0.2
RECOVERED_COV = 0.5
MASK_MAX = 0.5
NAMED_COPY_MIN_FRAC = 0.7
EXT_MIN_GAIN = 1.1
EXT_MIN_OLD_COV = 0.5
EXT_MATCH_COV = 0.9
EXT_TRUNCATED = 0.9
EXT_MAX_OVER_REF = 1.05

TE_RE = re.compile(r"retrovir|transpos|reverse transcriptase|LINE-1|\bL1\b|ORF1p|ORF2p|\bgag\b|\bpol\b|"
                   r"\benv\b|envelope|integrase|piggybac|tigger|jerky|pogo|endogenous", re.I)
UNCHAR_RE = re.compile(r"(PREDICTED: )?uncharacterized", re.I)


# ── evidence loci ─────────────────────────────────────────────────────────────

def read_protein_db_headers(path, taxa):
    """protein id -> (SYMBOL, description) for proteins of ``taxa`` with a GN= symbol."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if not (line.startswith(">") and line[1:].startswith(taxa)):
                continue
            m = re.search(r"GN=(\S+)", line)
            if not m:
                continue
            parts = line[1:].rstrip("\n").split(None, 1)
            desc = parts[1].split(" OS=")[0].strip() if len(parts) > 1 else ""
            out[parts[0]] = (m.group(1).upper(), desc)
    return out


def read_evidence_alignments(paf, prot_info, taxa):
    """Full-length, high-identity, intact alignments: (chrom, strand, exons, symbol, desc, protein)."""
    alns = []
    with open(paf) as fh:
        for line in fh:
            if not line.startswith(taxa):
                continue
            rec = parse_paf_row(line)
            if rec is None or rec["q_name"] not in prot_info:
                continue
            t = rec["tags"]
            if int(t.get("fs", 0)) or int(t.get("st", 0)):
                continue
            if (rec["q_end"] - rec["q_start"]) / rec["q_len"] < EV_MIN_QCOV or paf_identity(rec) < EV_MIN_PID:
                continue
            exons, *_ = parse_cigar_to_exons(t["cg"], rec["t_start"], rec["t_end"], rec["strand"])
            if exons:
                sym, desc = prot_info[rec["q_name"]]
                alns.append((rec["t_name"], rec["strand"], exons, sym, desc, rec["q_name"]))
    return alns


def build_evidence_loci(alns):
    """Merge overlapping same-strand alignments into loci and flag retrocopies."""
    by = collections.defaultdict(list)
    for a in alns:
        by[(a[0], a[1])].append(a)
    groups = []
    for (chrom, strand), al in by.items():
        al.sort(key=lambda a: a[2][0][0])
        cur, end = [], -1
        for a in al:
            if cur and a[2][0][0] < end:
                cur.append(a)
                end = max(end, a[2][-1][1])
            else:
                if cur:
                    groups.append((chrom, strand, cur))
                cur, end = [a], a[2][-1][1]
        if cur:
            groups.append((chrom, strand, cur))
    prot_introns = collections.defaultdict(int)
    sym_introns = collections.defaultdict(int)
    for a in alns:
        prot_introns[a[5]] = max(prot_introns[a[5]], len(a[2]) - 1)
        sym_introns[a[3]] = max(sym_introns[a[3]], len(a[2]) - 1)
    loci = []
    for chrom, strand, al in groups:
        symbol = collections.Counter(a[3] for a in al).most_common(1)[0][0]
        desc = collections.Counter(a[4] for a in al if a[3] == symbol).most_common(1)[0][0]
        introns = max(len(a[2]) - 1 for a in al)
        merged = []
        for s, e in sorted(b for a in al for b in a[2]):
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        retro = introns == 0 and (any(prot_introns[a[5]] >= 1 for a in al)
                                  or any(sym_introns[a[3]] >= 1 for a in al))
        loci.append({"chrom": chrom, "strand": strand, "blocks": merged,
                     "len": sum(e - s for s, e in merged), "symbol": symbol, "desc": desc,
                     "introns": introns, "retro": retro})
    return loci


# ── helpers ───────────────────────────────────────────────────────────────────

def cds_blocks(tx):
    return [(max(e.start, tx.thick_start), min(e.stop, tx.thick_stop)) for e in tx.exon_intervals
            if min(e.stop, tx.thick_stop) > max(e.start, tx.thick_start)]


def intact_peptide(tx, seqs):
    p = tx.get_protein_sequence(seqs).upper().rstrip("*")
    return bool(p) and "*" not in p


def masked_fraction(blocks, chrom_seq):
    tot = sum(b - a for a, b in blocks)
    if not tot:
        return 1.0
    low = sum(sum(1 for ch in chrom_seq[a:b] if ch.islower()) for a, b in blocks)
    return low / tot


class Index:
    """Interval index of (start, stop, payload) per (chrom, strand)."""

    def __init__(self):
        self.d = collections.defaultdict(list)

    def add(self, chrom, strand, start, stop, payload):
        self.d[(chrom, strand)].append((start, stop, payload))

    def finish(self):
        self.starts, self.maxlen = {}, {}
        for k, v in self.d.items():
            v.sort(key=lambda x: x[0])
            self.starts[k] = [x[0] for x in v]
            self.maxlen[k] = max(x[1] - x[0] for x in v)
        return self

    def overlapping(self, chrom, strand, start, stop):
        k = (chrom, strand)
        if k not in self.d:
            return
        v, st = self.d[k], self.starts[k]
        i = bisect.bisect_left(st, start - self.maxlen[k])
        while i < len(v) and v[i][0] < stop:
            if v[i][1] > start:
                yield v[i][2]
            i += 1


# ── the annotation ────────────────────────────────────────────────────────────

class Annotation:
    def __init__(self, gp_path, gp_info_path, ref_biotypes):
        with open(gp_info_path) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            self.fields = list(rd.fieldnames)
            self.info = collections.OrderedDict((x["transcript_id"], x) for x in rd)
        self.gp_rows = collections.OrderedDict()
        with open(gp_path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                self.gp_rows[f[0]] = f
        self.ref_bt = ref_biotypes
        self.ref_pc = {k for k, v in ref_biotypes.items() if v == PROTEIN_CODING}
        self.genes = collections.defaultdict(list)
        for t, x in self.info.items():
            if t in self.gp_rows:
                self.genes[x["gene_id"]].append(t)
        self._tx = {}

    def tx(self, t):
        if t not in self._tx:
            self._tx[t] = GenePredTranscript(self.gp_rows[t])
        return self._tx[t]

    def src(self, t):
        return self.info[t]["source_gene"].split(".")[0]

    def gene_x(self, g):
        return self.info[self.genes[g][0]]

    def is_pc(self, g):
        return self.gene_x(g)["gene_biotype"] == PROTEIN_CODING

    def ref_sourced(self, g):
        return any(self.src(t) in self.ref_pc for t in self.genes[g])

    def coding(self, g):
        return [self.tx(t) for t in self.genes[g] if int(self.gp_rows[t][6]) > int(self.gp_rows[t][5])]

    def ref_noncoding(self, t):
        x = self.info[t]
        b = x.get("gene_biotype", "")
        return x.get("alignment_mode") in PROJ_MODES and b != PROTEIN_CODING and "unitary" not in b and (
            "pseudogene" in b or b in ("lncRNA", "TEC") or b.startswith(("TR_", "IG_")))


def read_ref_protein_aa(ref_gp, gp_attrs):
    """reference gene (versionless) -> longest protein length (aa)."""
    tx2g = {}
    with open(gp_attrs) as fh:
        for line in fh:
            k, a, v = line.rstrip("\n").split("\t")[:3]
            if k.startswith("ENST") and a == "gene_id":
                tx2g[k] = v.split(".")[0]
    out = collections.defaultdict(int)
    with open(ref_gp) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            g = tx2g.get(f[0])
            cs, ce = int(f[5]), int(f[6])
            if not g or ce <= cs:
                continue
            n = sum(max(0, min(int(e), ce) - max(int(s), cs))
                    for s, e in zip(f[8].rstrip(",").split(","), f[9].rstrip(",").split(",")))
            out[g] = max(out[g], n // 3 - 1)
    return out


def read_ref_biotypes(gp_attrs):
    bt = {}
    with open(gp_attrs) as fh:
        for line in fh:
            k, a, v = line.rstrip("\n").split("\t")[:3]
            if k.startswith("ENSG") and a in ("gene_biotype", "gene_type"):
                bt[k.split(".")[0]] = v
    return bt


# ── evidence copies ───────────────────────────────────────────────────────────

def evidence_copies(ann, loci, miniprot_rows, seqs, stats):
    """Return (promote_genes, inserts, drop_genes). inserts: list of (gp_row, locus)."""
    pc_idx, cand_idx, nc_idx = Index(), Index(), Index()
    for g, ts in ann.genes.items():
        txs = ann.coding(g)
        if ann.is_pc(g) and txs:
            m = _merged_cds_blocks(txs)
            pc_idx.add(txs[0].chromosome, txs[0].strand, m[0][0], m[-1][1], m)
    for t in ann.info:
        if t not in ann.gp_rows:
            continue
        tx = ann.tx(t)
        if ann.ref_noncoding(t):
            nc_idx.add(tx.chromosome, tx.strand, tx.start, tx.stop, [(e.start, e.stop) for e in tx.exon_intervals])
        elif ann.info[t]["gene_biotype"] != PROTEIN_CODING and tx.cds_size > 0:
            cand_idx.add(tx.chromosome, tx.strand, tx.thick_start, tx.thick_stop, ("final", t, tx))
    for row in miniprot_rows:
        tx = GenePredTranscript(row)
        if tx.cds_size > 0:
            cand_idx.add(tx.chromosome, tx.strand, tx.thick_start, tx.thick_stop, ("miniprot", row, tx))
    for ix in (pc_idx, cand_idx, nc_idx):
        ix.finish()

    chosen = {}
    for r in loci:
        c, s = r["chrom"], r["strand"]
        lo, hi = r["blocks"][0][0], r["blocks"][-1][1]
        if any(_shared_bases(r["blocks"], m) >= RECOVERED_COV * r["len"] for m in pc_idx.overlapping(c, s, lo, hi)):
            continue
        if TE_RE.search(r["desc"]):
            stats["copy: transposon-like evidence"] += 1
            continue
        if UNCHAR_RE.match(r["desc"]):
            stats["copy: uncharacterized evidence"] += 1
            continue
        if r["retro"]:
            stats["copy: retrocopy"] += 1
            continue
        why = "copy: no candidate model"
        best = None
        for src, ref, tx in sorted(cand_idx.overlapping(c, s, lo, hi), key=lambda x: x[0] != "final"):
            b = cds_blocks(tx)
            L = sum(y - x for x, y in b)
            sh = _shared_bases(r["blocks"], b)
            if sh < COPY_LOCUS_COV * r["len"] or sh < COPY_MODEL_COV * L:
                why = "copy: candidate does not match evidence span"
                continue
            if tx.cds_size < 3 * COPY_MIN_AA or not intact_peptide(tx, seqs):
                why = "copy: candidate peptide broken/short"
                continue
            if any(_shared_bases(b, ex) > 0 for ex in nc_idx.overlapping(c, s, tx.thick_start, tx.thick_stop)):
                why = "copy: reference calls the locus non-coding"
                continue
            if any(_shared_bases(b, m) > COPY_MAX_PC_SHARE * L for m in pc_idx.overlapping(c, s, tx.thick_start, tx.thick_stop)):
                why = "copy: candidate shares CDS with a protein_coding gene"
                continue
            best = (src, ref, tx)
            break
        if best is None:
            stats[why] += 1
            continue
        key = best[1] if best[0] == "final" else best[1][0]
        chosen.setdefault(key, (best, r))

    promote = {ann.info[k]["gene_id"] for k, ((src, _, _), _) in chosen.items() if src == "final"}
    inserts = [(ref, r) for k, ((src, ref, tx), r) in chosen.items() if src == "miniprot"]
    drop = set()
    for row, r in inserts:
        tx = GenePredTranscript(row)
        b = cds_blocks(tx)
        for src, t2, tx2 in cand_idx.overlapping(tx.chromosome, tx.strand, tx.thick_start, tx.thick_stop):
            if src == "final" and ann.info[t2]["gene_id"] not in promote and _shared_bases(b, cds_blocks(tx2)) > 0:
                drop.add(ann.info[t2]["gene_id"])
    stats["copy: genes promoted in place"] = len(promote)
    stats["copy: miniprot models inserted"] = len(inserts)
    stats["copy: non-coding genes replaced"] = len(drop)
    return promote, inserts, drop


# ── precision rules ───────────────────────────────────────────────────────────

def precision_rules(ann, loci, seqs, stats, skip_genes=()):
    """Return (set_biotype: gene -> (biotype, reason), regene: gene -> reference-sourced gene)."""
    known_idx = Index()
    longest_aa = collections.defaultdict(int)
    for g in ann.genes:
        if ann.is_pc(g) and ann.ref_sourced(g):
            txs = ann.coding(g)
            if txs:
                m = _merged_cds_blocks(txs)
                known_idx.add(txs[0].chromosome, txs[0].strand, m[0][0], m[-1][1], (m, g))
                name = ann.gene_x(g)["source_gene_common_name"].upper()
                longest_aa[name] = max(longest_aa[name], max(t.cds_size for t in txs) // 3)
    known_idx.finish()
    ev_idx = Index()
    for r in loci:
        if not UNCHAR_RE.match(r["desc"]):
            ev_idx.add(r["chrom"], r["strand"], r["blocks"][0][0], r["blocks"][-1][1], r)
    ev_idx.finish()

    def evidence_match(txs):
        for tx in txs:
            b = _merged_cds_blocks([tx])
            L = sum(y - x for x, y in b)
            for r in ev_idx.overlapping(tx.chromosome, tx.strand, b[0][0], b[-1][1]):
                sh = _shared_bases(r["blocks"], b)
                if sh >= COPY_LOCUS_COV * r["len"] and sh >= COPY_MODEL_COV * L:
                    return True
        return False

    set_bt, regene = {}, {}
    for g in ann.genes:
        if g in skip_genes or not ann.is_pc(g) or ann.ref_sourced(g):
            continue
        txs = ann.coding(g)
        if not txs:
            continue
        x = ann.gene_x(g)
        s = x["source_gene"].split(".")[0]
        if s in ann.ref_bt and ann.ref_bt[s] != PROTEIN_CODING:
            set_bt[g] = (ann.ref_bt[s], "noncoding_source")
            stats["R2 projected from a reference non-coding gene"] += 1
            continue
        if not any(intact_peptide(t, seqs) for t in txs):
            set_bt[g] = (UNKNOWN_LIKELY_CODING, "internal_stop")
            stats["R1 no intact coding transcript"] += 1
            continue
        m = _merged_cds_blocks(txs)
        if x["transcript_class"] == "putative_novel_isoform":
            best = max(((_shared_bases(m, km), kg) for km, kg in
                        known_idx.overlapping(txs[0].chromosome, txs[0].strand, m[0][0], m[-1][1])),
                       default=(0, None))
            if best[0] > 0:
                regene[g] = best[1]
                stats["R3 IsoSeq isoform merged into its gene"] += 1
                continue
        if x["transcript_class"] == "putative_novel" and x["alignment_mode"] in ("augMP", "miniprot") \
                and masked_fraction(m, seqs[txs[0].chromosome]) > MASK_MAX:
            if evidence_match(txs):
                stats["R4 masked, kept by close-relative evidence"] += 1
            else:
                set_bt[g] = (UNKNOWN_LIKELY_CODING, "masked_no_close_relative")
                stats["R4 masked, no close-relative evidence"] += 1
                continue
        call = x.get("protein_refine_call", "")
        if call in ("real_copy", "named_copy", "rescued") and not s.startswith("ENS"):
            name = x["source_gene_common_name"].upper()
            aa = max(t.cds_size for t in txs) // 3
            if longest_aa.get(name) and aa < NAMED_COPY_MIN_FRAC * longest_aa[name]:
                set_bt[g] = (UNKNOWN_LIKELY_CODING, "short_named_copy")
                stats["R5 short named copy"] += 1
    return set_bt, regene


# ── evidence extensions ───────────────────────────────────────────────────────

def evidence_extensions(ann, alns, miniprot_rows, seqs, stats, ref_aa, skip_genes=()):
    """Return [(gp_row, gene_id)]: full-length miniprot models to add as isoforms."""
    mp_idx = Index()
    for row in miniprot_rows:
        tx = GenePredTranscript(row)
        if tx.cds_size > 0:
            mp_idx.add(tx.chromosome, tx.strand, tx.thick_start, tx.thick_stop, (row, tx))
    mp_idx.finish()
    pc_idx = Index()
    pc = {}
    for g in ann.genes:
        if g in skip_genes or not ann.is_pc(g):
            continue
        txs = ann.coding(g)
        if not txs:
            continue
        m = _merged_cds_blocks(txs)
        pc[g] = (txs, m)
        pc_idx.add(txs[0].chromosome, txs[0].strand, m[0][0], m[-1][1], (g, m))
    pc_idx.finish()
    by_sym = collections.defaultdict(list)
    for a in alns:
        by_sym[(a[0], a[1], a[3])].append(a)

    out = []
    for g, (txs, m) in pc.items():
        name = ann.gene_x(g)["source_gene_common_name"].upper()
        if not name or name in ("N/A", "NONE", "NAN"):
            continue
        ref_max = max((ref_aa.get(ann.src(t), 0) for t in ann.genes[g]), default=0)
        old_aa = max(t.cds_size for t in txs) // 3 - 1
        if not ref_max or old_aa >= EXT_TRUNCATED * ref_max:
            continue
        c, s_ = txs[0].chromosome, txs[0].strand
        old_len = sum(b - a for a, b in m)
        best = None
        for a in by_sym.get((c, s_, name), []):
            ex = a[2]
            if ex[-1][1] <= m[0][0] or ex[0][0] >= m[-1][1]:
                continue
            ev_len = sum(y - x for x, y in ex)
            if ev_len // 3 < EXT_MIN_GAIN * old_aa or ev_len // 3 > EXT_MAX_OVER_REF * ref_max \
                    or _shared_bases(ex, m) < EXT_MIN_OLD_COV * old_len:
                continue
            for row, tx in mp_idx.overlapping(c, s_, ex[0][0], ex[-1][1]):
                b = cds_blocks(tx)
                L = sum(y - x for x, y in b)
                sh = _shared_bases(ex, b)
                if sh < EXT_MATCH_COV * ev_len or sh < EXT_MATCH_COV * L:
                    continue
                new_aa = L // 3 - 1
                if new_aa < EXT_MIN_GAIN * old_aa or new_aa > EXT_MAX_OVER_REF * ref_max \
                        or _shared_bases(b, m) < EXT_MIN_OLD_COV * old_len:
                    continue
                if any(og != g and _shared_bases(b, om) > COPY_MAX_PC_SHARE * L
                       for og, om in pc_idx.overlapping(c, s_, tx.thick_start, tx.thick_stop)):
                    stats["ext: would fuse with another protein_coding gene"] += 1
                    continue
                if not intact_peptide(tx, seqs):
                    continue
                if best is None or L > best[1]:
                    best = (row, L)
        if best:
            out.append((best[0], g))
            stats[f"ext: isoform added ({ann.gene_x(g)['alignment_mode']} gene)"] += 1
    stats["ext: genes extended"] = len(out)
    return out


# ── writing ───────────────────────────────────────────────────────────────────

def _attrs(s):
    d = collections.OrderedDict()
    for item in s.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def _render(d):
    return ";".join(f"{k}={v}" for k, v in d.items())


def _prefix(ann):
    m = re.match(r"^(.*_)G\d+$", next(iter(ann.genes)))
    return m.group(1) if m else "EV_"


def _next_index(ids, tag):
    """1 + the largest number already used after `tag` in ids (so reruns never collide)."""
    nums = [int(m.group(1)) for i in ids for m in [re.search(tag + r"(\d+)$", i)] if m]
    return max(nums, default=0) + 1


def new_ids(ann, n):
    """n fresh (gene_id, transcript_id) pairs in the annotation's naming scheme."""
    prefix = _prefix(ann)
    start = max(_next_index(ann.genes, "GEVC"), _next_index(ann.info, "TEVC"))
    return [(f"{prefix}GEVC{i:06d}", f"{prefix}TEVC{i:06d}") for i in range(start, start + n)]


def gff3_records(row, gene_id, tx_id, attrs, with_gene=True):
    tx = GenePredTranscript(row)
    chrom, strand = tx.chromosome, tx.strand
    base = dict(attrs, gene_id=gene_id, transcript_id=tx_id)
    lines = [[chrom, "CAT", "gene", tx.start + 1, tx.stop, ".", strand, ".",
              _render(collections.OrderedDict(gene_id=gene_id, gene_biotype=attrs["gene_biotype"],
                                              source_gene_common_name=attrs["source_gene_common_name"],
                                              gene_sources="miniprot", ID=gene_id, Name=gene_id))],
             [chrom, "CAT", "transcript", tx.start + 1, tx.stop, ".", strand, ".",
              _render(collections.OrderedDict(base, Parent=gene_id, ID=tx_id, Name=gene_id))]]
    if not with_gene:
        lines = lines[1:]
    exons = [(e.start, e.stop) for e in tx.exon_intervals]
    order = exons if strand == "+" else exons[::-1]
    for i, (a, b) in enumerate(order):
        lines.append([chrom, "CAT", "exon", a + 1, b, ".", strand, ".",
                      _render(collections.OrderedDict(base, Parent=tx_id, ID=f"exon:{tx_id}:{i}", Name=gene_id))])
    done = 0
    for i, (a, b) in enumerate(order):
        a2, b2 = max(a, tx.thick_start), min(b, tx.thick_stop)
        if b2 <= a2:
            continue
        phase = (3 - done % 3) % 3
        lines.append([chrom, "CAT", "CDS", a2 + 1, b2, ".", strand, str(phase),
                      _render(collections.OrderedDict(base, Parent=tx_id, ID=f"CDS:{tx_id}:{i}", Name=gene_id))])
        done += b2 - a2
    return ["\t".join(map(str, l)) for l in lines]


def write_outputs(ann, args, set_bt, regene, promote, inserts, drop, extensions=()):
    tx_set = {}       # transcript -> dict of gp_info/GFF3 attribute overrides
    for g in promote:
        for t in ann.genes[g]:
            tx_set[t] = {"gene_biotype": PROTEIN_CODING, "transcript_biotype": PROTEIN_CODING,
                         "miniprot_novel_call": "evidence_copy"}
    for g, (bt, why) in set_bt.items():
        for t in ann.genes[g]:
            tx_set[t] = {"gene_biotype": bt, "transcript_biotype": bt, "miniprot_novel_call": why}
    for g, kg in regene.items():
        for t in ann.genes[g]:
            tx_set[t] = {"gene_id": kg}
    drop_tx = {t for g in drop for t in ann.genes[g]}
    ids = new_ids(ann, len(inserts))
    new_rows = []
    for (row, r), (gid, tid) in zip(inserts, ids):
        x = {k: "N/A" for k in ann.fields}
        x.update(gene_id=gid, transcript_id=tid, gene_biotype=PROTEIN_CODING, transcript_biotype=PROTEIN_CODING,
                 alignment_mode="miniprot", transcript_class="putative_novel", source_transcript=row[0],
                 alignment_id=row[0], source_gene_common_name=r["symbol"], miniprot_novel_call="evidence_copy")
        new_rows.append((row, gid, tid, x))
    ext_rows = []
    tprefix = _prefix(ann)
    for i, (row, g) in enumerate(extensions, _next_index(ann.info, "TEVX")):
        x = dict(ann.gene_x(g))
        x.update(transcript_id=f"{tprefix}TEVX{i:06d}", transcript_biotype=PROTEIN_CODING,
                 alignment_mode="miniprot", alignment_id=row[0], source_transcript=row[0],
                 miniprot_novel_call="evidence_extension")
        ext_rows.append((row, g, x["transcript_id"], x))

    with open(args.output_gp, "w") as fh:
        for t, f in ann.gp_rows.items():
            if t not in drop_tx:
                fh.write("\t".join(f) + "\n")
        for row, gid, tid, _ in new_rows + ext_rows:
            fh.write("\t".join([tid] + row[1:]) + "\n")
    with open(args.output_gp_info, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=ann.fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for t, x in ann.info.items():
            if t in drop_tx:
                continue
            w.writerow(dict(x, **tx_set.get(t, {})))
        for _, _, _, x in new_rows + ext_rows:
            w.writerow(x)

    # gene spans after R3 merges
    span = {}
    for g, kg in regene.items():
        for t in ann.genes[g] + ann.genes[kg]:
            tx = ann.tx(t)
            a, b = span.get(kg, (tx.start, tx.stop))
            span[kg] = (min(a, tx.start), max(b, tx.stop))
    for row, g, _, _ in ext_rows:
        tx = GenePredTranscript(row)
        a, b = span.get(g, (min(ann.tx(t).start for t in ann.genes[g]), max(ann.tx(t).stop for t in ann.genes[g])))
        span[g] = (min(a, tx.start), max(b, tx.stop))
    merged_away = set(regene)
    gene_bt = {g: v[0] for g, v in set_bt.items()}
    gene_bt.update({g: PROTEIN_CODING for g in promote})
    with open(args.gff3) as inf, open(args.output_gff3, "w") as out:
        for line in inf:
            if line.startswith("#") or not line.strip():
                out.write(line)
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                out.write(line)
                continue
            a = _attrs(p[8])
            gid = a.get("gene_id") or a.get("ID", "")
            if gid in drop:
                continue
            if p[2] == "gene":
                if gid in merged_away:
                    continue
                if gid in span:
                    p[3], p[4] = str(span[gid][0] + 1), str(span[gid][1])
                if gid in gene_bt:
                    a["gene_biotype"] = gene_bt[gid]
            else:
                o = tx_set.get(a.get("transcript_id", ""))
                if o:
                    for k, v in o.items():
                        if k in a or k == "miniprot_novel_call":
                            a[k] = v
                    if "gene_id" in o and p[2] in ("transcript", "mRNA"):
                        a["Parent"] = o["gene_id"]
            p[8] = _render(a)
            out.write("\t".join(p) + "\n")
        for row, gid, tid, x in new_rows:
            attrs = collections.OrderedDict((k, x[k]) for k in (
                "source_gene_common_name", "gene_biotype", "transcript_biotype", "alignment_mode",
                "alignment_id", "source_transcript", "transcript_class", "miniprot_novel_call"))
            for l in gff3_records(row, gid, tid, attrs):
                out.write(l + "\n")
        for row, g, tid, x in ext_rows:
            attrs = collections.OrderedDict((k, x.get(k, "N/A")) for k in (
                "source_gene", "source_gene_common_name", "gene_biotype", "transcript_biotype", "alignment_mode",
                "alignment_id", "source_transcript", "transcript_class", "miniprot_novel_call"))
            for l in gff3_records(row, g, tid, attrs, with_gene=False):
                out.write(l + "\n")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gp", required=True)
    ap.add_argument("--gp-info", required=True)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--output-gp", required=True)
    ap.add_argument("--output-gp-info", required=True)
    ap.add_argument("--output-gff3", required=True)
    ap.add_argument("--genome-fasta", required=True, help="soft-masked genome FASTA")
    ap.add_argument("--miniprot-paf", required=True)
    ap.add_argument("--miniprot-gp", required=True, help="miniprot_consensus.gp (candidate models)")
    ap.add_argument("--protein-db", required=True, help="protein FASTA miniprot was run with")
    ap.add_argument("--close-taxids", required=True, help="comma-separated taxids counted as close relatives")
    ap.add_argument("--ref-gp-attrs", required=True)
    ap.add_argument("--ref-gp", required=True, help="reference genePred (reference protein lengths)")
    ap.add_argument("--report", default=None)
    args = ap.parse_args()

    taxa = tuple(t.strip() + "|" for t in args.close_taxids.split(",") if t.strip())
    prot_info = read_protein_db_headers(args.protein_db, taxa)
    logger.info(f"close-relative proteins with a gene symbol: {len(prot_info):,}")
    alns = read_evidence_alignments(args.miniprot_paf, prot_info, taxa)
    loci = build_evidence_loci(alns)
    logger.info(f"evidence alignments {len(alns):,} -> loci {len(loci):,} "
                f"({sum(r['retro'] for r in loci):,} retrocopies)")
    ann = Annotation(args.gp, args.gp_info, read_ref_biotypes(args.ref_gp_attrs))
    seqs = get_sequence_dict(args.genome_fasta, upper=False)
    miniprot_rows = [l.rstrip("\n").split("\t") for l in open(args.miniprot_gp)]
    stats = collections.Counter()
    n_pc_before = sum(ann.is_pc(g) for g in ann.genes)
    promote, inserts, drop = evidence_copies(ann, loci, miniprot_rows, seqs, stats)
    set_bt, regene = precision_rules(ann, loci, seqs, stats, skip_genes=promote | drop)
    extensions = evidence_extensions(ann, alns, miniprot_rows, seqs, stats,
                                     read_ref_protein_aa(args.ref_gp, args.ref_gp_attrs), skip_genes=promote | drop | set(set_bt) | set(regene))
    write_outputs(ann, args, set_bt, regene, promote, inserts, drop, extensions)
    n_pc_after = n_pc_before + len(promote) + len(inserts) - len(set_bt) - len(regene) \
        - sum(ann.is_pc(g) for g in drop)
    lines = [f"protein_coding genes {n_pc_before:,} -> {n_pc_after:,}"] + \
        [f"  {v:7,}  {k}" for k, v in sorted(stats.items())]
    for l in lines:
        logger.info(l)
    if args.report:
        Path(args.report).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
