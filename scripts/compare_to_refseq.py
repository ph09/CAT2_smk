#!/usr/bin/env python3
"""
Compare a CAT2 consensus annotation to a RefSeq GFF on the SAME assembly.

The two annotations use different sequence names (GenBank CM*/JA* vs RefSeq
NC_*/NW_*), but the underlying assembly is identical, so sequences are matched
1:1 by exact length (falling back to file order for any ambiguous lengths).

Comparison is done at the protein-coding gene level using CDS-exon overlap on the
same strand:

  * recall     = RefSeq PC genes whose CDS overlaps any of our PC CDS
  * precision  = our PC genes whose CDS overlaps any RefSeq PC CDS
  * novel validation = for our protein-only-novel genes (split into paralog vs
    lineage_specific), how many are independently annotated (CDS overlap) in
    RefSeq -- i.e. corroborated real genes -- vs unique to our set.

Usage:
  compare_to_refseq.py --name <genome> --gp <consensus.gp> \
      --gp-info <consensus_novel_annotated.gp_info> --fai <genome.fa.fai> \
      --refseq-gff <refseq.gff>
"""
import argparse
import bisect
import collections
import sys


def parse_fai(fai):
    """our_seqname -> length"""
    out = {}
    with open(fai) as fh:
        for ln in fh:
            p = ln.split('\t')
            out[p[0]] = int(p[1])
    return out




def build_seq_map(fai_lengths, refseq_lengths):
    """our_seqname -> refseq_seqname, matched by exact length (unique lengths)."""
    ref_by_len = collections.defaultdict(list)
    for name, ln in refseq_lengths.items():
        ref_by_len[ln].append(name)
    our_by_len = collections.defaultdict(list)
    for name, ln in fai_lengths.items():
        our_by_len[ln].append(name)

    seq_map = {}
    ambiguous = 0
    unmatched = 0
    for name, ln in fai_lengths.items():
        cand = ref_by_len.get(ln, [])
        if len(cand) == 1 and len(our_by_len[ln]) == 1:
            seq_map[name] = cand[0]
        elif len(cand) == 0:
            unmatched += 1
        else:
            ambiguous += 1
    return seq_map, ambiguous, unmatched


def _attr(s, key):
    for item in s.split(';'):
        if item.startswith(key + '='):
            return item[len(key) + 1:]
    return None


def parse_refseq(gff, pc_only=True):
    """
    Returns:
      gene_cds[gene_id] -> list of (seq, strand, start, end)  [0-based half-open]
      pc_gene_ids -> set
    """
    gene_biotype = {}
    gene_span = {}
    gene_cds = collections.defaultdict(list)
    ref_lengths = {}
    with open(gff) as fh:
        for ln in fh:
            if ln.startswith('#'):
                # ##sequence-region lines are interleaved with features (one per
                # sequence), so scan the whole file for them.
                if ln.startswith('##sequence-region'):
                    q = ln.split()
                    ref_lengths[q[1]] = int(q[3])
                continue
            p = ln.rstrip('\n').split('\t')
            if len(p) < 9:
                continue
            ftype = p[2]
            if ftype in ('gene', 'pseudogene'):
                gid = _attr(p[8], 'gene') or _attr(p[8], 'ID')
                bt = _attr(p[8], 'gene_biotype') or ('pseudogene' if ftype == 'pseudogene' else None)
                if gid:
                    gene_biotype[gid] = bt
                    gene_span[gid] = (p[0], int(p[3]) - 1, int(p[4]))
            elif ftype == 'CDS':
                gid = _attr(p[8], 'gene')
                if not gid:
                    continue
                gene_cds[gid].append((p[0], p[6], int(p[3]) - 1, int(p[4])))
    pc_gene_ids = {g for g, bt in gene_biotype.items() if bt == 'protein_coding'}
    return gene_cds, pc_gene_ids, ref_lengths, gene_biotype, gene_span


def parse_our_gp_info(gp_info):
    """transcript_id -> (gene_id, gene_biotype, is_novel, novel_class)"""
    tx = {}
    with open(gp_info) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        idx = {c: i for i, c in enumerate(header)}
        for ln in fh:
            p = ln.rstrip('\n').split('\t')
            if len(p) < len(header):
                continue
            gid = p[idx['gene_id']]
            bt = p[idx['gene_biotype']]
            is_novel = p[idx.get('protein_only_novel', -1)] == 'True' if 'protein_only_novel' in idx else False
            nc = 'N/A'
            if 'novel_class' in idx:
                nc = p[idx['novel_class']]
            if nc == 'N/A' and is_novel:
                desc = p[idx['novel_gene_description']] if 'novel_gene_description' in idx else 'N/A'
                nc = 'paralog' if desc.startswith('paralog') else 'lineage_specific'
            tx[p[idx['transcript_id']]] = (gid, bt, is_novel, nc)
    return tx


def parse_our_gp(gp, tx_info, seq_map):
    """
    Aggregate our transcripts into genes with CDS-exon intervals mapped into
    RefSeq sequence space.
    Returns gene[gene_id] = {'bt','novel','nc','cds':[(refseq_seq,strand,s,e)]}
    """
    gene = {}
    with open(gp) as fh:
        for ln in fh:
            p = ln.rstrip('\n').split('\t')
            txid = p[0]
            info = tx_info.get(txid)
            if info is None:
                continue
            gid, bt, is_novel, nc = info
            chrom = seq_map.get(p[1])
            if chrom is None:
                continue  # sequence not mapped to RefSeq
            strand = p[2]
            cds_s = int(p[5]); cds_e = int(p[6])
            if cds_e <= cds_s:
                continue
            exs = [int(x) for x in p[8].strip(',').split(',') if x]
            exe = [int(x) for x in p[9].strip(',').split(',') if x]
            g = gene.get(gid)
            if g is None:
                g = {'bt': bt, 'novel': is_novel, 'nc': nc, 'cds': []}
                gene[gid] = g
            for s, e in zip(exs, exe):
                cs = max(s, cds_s); ce = min(e, cds_e)
                if ce > cs:
                    g['cds'].append((chrom, strand, cs, ce))
    return gene


def build_merged(cds_iter):
    """(seq,strand) -> (starts[], stops[]) merged sorted intervals."""
    raw = collections.defaultdict(list)
    for seq, strand, s, e in cds_iter:
        raw[(seq, strand)].append((s, e))
    merged = {}
    for key, ivs in raw.items():
        ivs.sort()
        starts, stops = [], []
        cs, ce = ivs[0]
        for s, e in ivs[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                starts.append(cs); stops.append(ce); cs, ce = s, e
        starts.append(cs); stops.append(ce)
        merged[key] = (starts, stops)
    return merged


def overlaps(merged, seq, strand, s, e):
    f = merged.get((seq, strand))
    if not f:
        return False
    starts, stops = f
    j = bisect.bisect_left(starts, e) - 1
    return j >= 0 and stops[j] > s


def build_merged_spans(span_iter):
    """seq -> (starts[], stops[]) merged, strand-agnostic (for gene-locus tests)."""
    raw = collections.defaultdict(list)
    for seq, s, e in span_iter:
        raw[seq].append((s, e))
    merged = {}
    for key, ivs in raw.items():
        ivs.sort()
        starts, stops = [], []
        cs, ce = ivs[0]
        for s, e in ivs[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                starts.append(cs); stops.append(ce); cs, ce = s, e
        starts.append(cs); stops.append(ce)
        merged[key] = (starts, stops)
    return merged


def overlaps_seq(merged, seq, s, e):
    f = merged.get(seq)
    if not f:
        return False
    starts, stops = f
    j = bisect.bisect_left(starts, e) - 1
    return j >= 0 and stops[j] > s


def _biotype_group(bt):
    if bt == 'protein_coding':
        return 'protein_coding'
    if bt and 'pseudogene' in bt:
        return 'pseudogene'
    if bt == 'lncRNA':
        return 'lncRNA'
    return 'other_ncRNA'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--name', required=True)
    ap.add_argument('--gp', required=True)
    ap.add_argument('--gp-info', required=True)
    ap.add_argument('--fai', required=True)
    ap.add_argument('--refseq-gff', required=True)
    args = ap.parse_args()

    fai_len = parse_fai(args.fai)
    ref_gene_cds, ref_pc, ref_len, ref_biotype, ref_span = parse_refseq(args.refseq_gff)
    seq_map, amb, unm = build_seq_map(fai_len, ref_len)

    # Restrict to sequences present in both assemblies (exclude MT / RefSeq-only
    # scaffolds, which are unrecoverable by construction).
    mapped = set(seq_map.values())
    ref_pc = {g for g in ref_pc
              if ref_gene_cds.get(g) and ref_gene_cds[g][0][0] in mapped}

    print(f"===== {args.name} =====")
    print(f"seq mapping: {len(seq_map)} of {len(fai_len)} our sequences matched to RefSeq "
          f"by length (ambiguous={amb}, unmatched={unm})")
    print(f"RefSeq: {len(ref_pc):,} protein_coding genes (on shared sequences; MT/unmapped excluded)")

    tx_info = parse_our_gp_info(args.gp_info)
    our_gene = parse_our_gp(args.gp, tx_info, seq_map)
    our_pc = {g: d for g, d in our_gene.items() if d['bt'] == 'protein_coding'}
    print(f"Ours:   {len(our_pc):,} protein_coding genes (with CDS on mapped sequences)")

    # Merged CDS sets
    ref_merged = build_merged(iv for g in ref_pc for iv in ref_gene_cds.get(g, []))
    our_merged = build_merged(iv for d in our_pc.values() for iv in d['cds'])

    # Recall: RefSeq PC gene recovered if any CDS overlaps our PC CDS
    recovered = 0
    for g in ref_pc:
        for (seq, strand, s, e) in ref_gene_cds.get(g, []):
            if overlaps(our_merged, seq, strand, s, e):
                recovered += 1
                break
    recall = 100.0 * recovered / max(1, len(ref_pc))

    # Precision + novel validation: our PC gene overlaps a RefSeq PC CDS?
    our_matched = 0
    novel_tot = collections.Counter()
    novel_val = collections.Counter()
    for gid, d in our_pc.items():
        hit = any(overlaps(ref_merged, seq, strand, s, e) for (seq, strand, s, e) in d['cds'])
        if hit:
            our_matched += 1
        if d['novel']:
            key = d['nc'] if d['nc'] in ('paralog', 'lineage_specific') else 'unclassified'
            novel_tot[key] += 1
            novel_tot['ALL'] += 1
            if hit:
                novel_val[key] += 1
                novel_val['ALL'] += 1
    precision = 100.0 * our_matched / max(1, len(our_pc))

    print(f"RECALL   : {recovered:,}/{len(ref_pc):,} RefSeq PC genes recovered ({recall:.1f}%)")
    print(f"PRECISION: {our_matched:,}/{len(our_pc):,} of our PC genes overlap a RefSeq PC gene ({precision:.1f}%)")

    # Characterise what RefSeq calls at the loci of our protein-only-novel genes.
    # Priority: PC CDS (same strand) > pseudogene > lncRNA > other ncRNA > nothing.
    span_merged = {}
    for grp in ('pseudogene', 'lncRNA', 'other_ncRNA'):
        span_merged[grp] = build_merged_spans(
            (ref_span[g][0], ref_span[g][1], ref_span[g][2])
            for g, bt in ref_biotype.items()
            if g in ref_span and _biotype_group(bt) == grp
        )

    def classify_novel(d):
        # d['cds'] = list of (seq,strand,s,e)
        if any(overlaps(ref_merged, seq, strand, s, e) for (seq, strand, s, e) in d['cds']):
            return 'refseq_protein_coding'
        # strand-agnostic span test against other biotypes
        by_seq = collections.defaultdict(lambda: [1 << 62, -1])
        for (seq, strand, s, e) in d['cds']:
            b = by_seq[seq]
            b[0] = min(b[0], s); b[1] = max(b[1], e)
        for grp in ('pseudogene', 'lncRNA', 'other_ncRNA'):
            for seq, (s, e) in by_seq.items():
                if overlaps_seq(span_merged[grp], seq, s, e):
                    return 'refseq_' + grp
        return 'refseq_intergenic'

    cats = ('refseq_protein_coding', 'refseq_pseudogene', 'refseq_lncRNA',
            'refseq_other_ncRNA', 'refseq_intergenic')
    breakdown = {k: collections.Counter() for k in ('ALL', 'paralog', 'lineage_specific')}
    for gid, d in our_pc.items():
        if not d['novel']:
            continue
        cat = classify_novel(d)
        breakdown['ALL'][cat] += 1
        key = d['nc'] if d['nc'] in ('paralog', 'lineage_specific') else None
        if key:
            breakdown[key][cat] += 1

    print("NOVEL genes — what RefSeq annotates at their locus:")
    for grpkey in ('ALL', 'paralog', 'lineage_specific'):
        tot = sum(breakdown[grpkey].values())
        if not tot:
            continue
        print(f"  {grpkey} (n={tot:,}):")
        for cat in cats:
            n = breakdown[grpkey][cat]
            if n:
                print(f"     {cat:24s}: {n:,} ({100.0*n/tot:.1f}%)")
    print()


if __name__ == '__main__':
    main()
