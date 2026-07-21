"""Postprocessing for the consensus gene set.
"""
import logging
import os
from collections import defaultdict

import pandas as pd

logger = logging.getLogger(__name__)


# ---- ref loading -----------------------------------------------------------

def _find_ref_gp_attrs(ref_gp):
    """Look for a ref_gp_attrs file next to the ref_gp."""
    candidates = [
        ref_gp + '_attrs',
        ref_gp.replace('.gp', '.gp_attrs'),
        ref_gp.replace('.gp', '_attrs'),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return None


def _load_ref_pc_info(ref_gp, ref_gp_attrs=None):
    """Return gene_name -> {'max_exons', 'span'} restricted to pc."""
    if ref_gp_attrs is None:
        ref_gp_attrs = _find_ref_gp_attrs(ref_gp)
    if not ref_gp_attrs or not os.path.exists(ref_gp_attrs):
        logger.warning(f"ref_gp_attrs not found near {ref_gp}; postprocess skipped")
        return None
    tx_to_biotype = {}
    tx_to_name = {}
    with open(ref_gp_attrs) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 3:
                continue
            if f[1] == 'gene_biotype':
                tx_to_biotype[f[0]] = f[2]
            elif f[1] == 'gene_name':
                tx_to_name[f[0]] = f[2]
    out = defaultdict(lambda: {'max_exons': 0, 'span_s': None, 'span_e': None})
    with open(ref_gp) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 12:
                continue
            if tx_to_biotype.get(f[0]) != 'protein_coding':
                continue
            name = tx_to_name.get(f[0], f[11])
            s = int(f[3]); e = int(f[4])
            d = out[name]
            d['max_exons'] = max(d['max_exons'], int(f[7]))
            d['span_s'] = s if d['span_s'] is None else min(d['span_s'], s)
            d['span_e'] = e if d['span_e'] is None else max(d['span_e'], e)
    return {n: {'max_exons': d['max_exons'],
                'span': (d['span_e'] - d['span_s']) if d['span_s'] is not None else 0}
            for n, d in out.items()}


# ---- helpers to inspect a gene record -------------------------------------

def _gene_name(entries):
    """Get gene_name from the first transcript's attrs."""
    for tx_obj, attrs in entries:
        name = attrs.get('source_gene_common_name') or attrs.get('source_gene')
        if name and name != 'N/A':
            return name
    return getattr(entries[0][0], 'name2', '') if entries else ''


def _tx_features(attrs):
    """Extract (mode, n_exons, n_introns_total, n_introns_supported)."""
    mode = attrs.get('alignment_mode', '')
    esup = attrs.get('exon_annotation_support', '') or ''
    isup = attrs.get('intron_annotation_support', '') or ''
    n_exons = (esup.count(',') + 1) if esup else 0
    n_introns_total = (isup.count(',') + 1) if isup else 0
    n_introns_supported = sum(1 for x in isup.split(',') if x == '1')
    return mode, n_exons, n_introns_total, n_introns_supported


# ---- Step 1: split runaway pc genes ---------------------------------------

def split_runaway_pc_genes(consensus_gene_dict,
                            ref_pc_info,
                            runaway_ratio=3.0,
                            runaway_min_bp=100_000,
                            cluster_gap_mult=5.0,
                            cluster_gap_floor=100_000,
                            id_template='{base}_dup{idx}'):
    """For each pc gene record whose transcripts span far more than the
    reference, cluster the transcripts and split into per-cluster records.

    Returns: number of source genes split.
    """
    n_split = 0
    for chrom in list(consensus_gene_dict.keys()):
        chrom_dict = consensus_gene_dict[chrom]
        for gene_id in list(chrom_dict.keys()):
            entries = chrom_dict[gene_id]
            if not entries:
                continue
            # Only act on pc records (use first transcript's biotype as proxy)
            first_attrs = entries[0][1]
            if first_attrs.get('gene_biotype') and first_attrs['gene_biotype'] != 'protein_coding':
                continue
            name = _gene_name(entries)
            ref = ref_pc_info.get(name)
            if not ref or ref['span'] <= 0:
                continue
            ref_span = ref['span']
            # Compute gene-level span on target
            starts = [t.start for t, _ in entries]
            stops  = [t.stop  for t, _ in entries]
            gene_span = max(stops) - min(starts)
            if gene_span < runaway_min_bp or gene_span / ref_span < runaway_ratio:
                continue
            # Cluster transcripts by genomic gap
            sorted_entries = sorted(entries, key=lambda e: e[0].start)
            gap_threshold = max(int(cluster_gap_mult * ref_span), cluster_gap_floor)
            clusters = []
            cur = [sorted_entries[0]]
            cur_e = sorted_entries[0][0].stop
            for e in sorted_entries[1:]:
                if e[0].start - cur_e > gap_threshold:
                    clusters.append(cur)
                    cur = [e]
                    cur_e = e[0].stop
                else:
                    cur.append(e)
                    cur_e = max(cur_e, e[0].stop)
            clusters.append(cur)
            if len(clusters) < 2:
                continue
            # Sort clusters by size desc; cluster 0 keeps the original id
            clusters.sort(key=lambda cl: -len(cl))
            chrom_dict[gene_id] = clusters[0]
            for idx, cl in enumerate(clusters[1:], start=1):
                new_gene_id = id_template.format(base=gene_id, idx=idx)
                chrom_dict[new_gene_id] = cl
                for tx_obj, attrs in cl:
                    attrs['gene_id'] = new_gene_id
            n_split += 1
            logger.debug(f"split runaway gene {name} ({gene_id}) into "
                         f"{len(clusters)} clusters")
    return n_split


# ---- Step 2: drop / reclassify ---------------------------------------------

STRONG_MODES = {'transMap', 'transMap_pairwise',
                'augTM', 'augTM_pairwise',
                'augPB', 'strg'}


def _collect_modes(entries):
    return set(_tx_features(attrs)[0] for _, attrs in entries)


def _collect_introns(entries):
    total = supported = 0
    max_exons = 0
    for _, attrs in entries:
        _, ne, nt, ns = _tx_features(attrs)
        total += nt
        supported += ns
        max_exons = max(max_exons, ne)
    return total, supported, max_exons


def reclassify_and_drop_pc_records(consensus_gene_dict,
                                    ref_pc_info,
                                    min_introns_for_low_support=3,
                                    augpb_chimera_exon_ratio=1.5,
                                    low_support_fraction=0.3,
                                    protect_strong_modes=True):
    """Reclassify processed pseudogenes and drop redundant copies.

    Mutates consensus_gene_dict and the attrs dicts in place.

    ``protect_strong_modes`` (default True) prevents the low-intron-support rule
    from dropping a locus that is backed by a strong mode (transMap/augTM/augPB/
    strg): a diverged paralog can legitimately have low *reference* intron support
    yet be a real gene copy, so we only drop weak-mode (e.g. augMP/txTM-only)
    duplicates. ``low_support_fraction`` is the supported/total intron ratio below
    which a low-support duplicate is eligible for dropping.

    Returns (n_reclassified, n_dropped, report_rows).
    """
    # Index gene records by gene_name to find alt loci.
    by_name = defaultdict(list)  # name -> [(chrom, gene_id)]
    for chrom, chrom_dict in consensus_gene_dict.items():
        for gene_id, entries in chrom_dict.items():
            if not entries:
                continue
            first_attrs = entries[0][1]
            if first_attrs.get('gene_biotype') and \
               first_attrs['gene_biotype'] != 'protein_coding':
                continue
            name = _gene_name(entries)
            by_name[name].append((chrom, gene_id))

    n_drop = 0
    n_reclass = 0
    report_rows = []
    drop_keys = []

    for name, loci in by_name.items():
        ref = ref_pc_info.get(name)
        if not ref:
            continue
        ref_exons = ref['max_exons']

        # Per-locus stats
        loci_info = []
        for chrom, gid in loci:
            entries = consensus_gene_dict[chrom][gid]
            modes = _collect_modes(entries)
            total, supported, max_exons = _collect_introns(entries)
            loci_info.append({
                'chrom': chrom, 'gid': gid, 'entries': entries,
                'modes': modes, 'max_tx_exons': max_exons,
                'introns_total': total, 'introns_supported': supported,
            })
        all_modes = set().union(*(li['modes'] for li in loci_info))
        # Strong alt locus: any locus has at least one STRONG mode
        any_strong = any(li['modes'] & STRONG_MODES for li in loci_info)

        for li in loci_info:
            entries = li['entries']
            modes = li['modes']
            max_exons = li['max_tx_exons']
            total = li['introns_total']
            supported = li['introns_supported']
            chrom = li['chrom']
            gid = li['gid']
            has_alt = len(loci_info) > 1
            # Strong alt: an *other* locus has a STRONG mode
            other_strong = any(
                (other['modes'] & STRONG_MODES)
                for other in loci_info if other is not li
            )

            # Rule A: reclassify processed pseudogenes (single-exon + ref multi-exon + alt exists)
            if max_exons == 1 and ref_exons >= 2 and has_alt:
                for _, attrs in entries:
                    attrs['gene_biotype'] = 'processed_pseudogene'
                    attrs['transcript_biotype'] = 'transcribed_pseudogene'
                n_reclass += 1
                report_rows.append((gid, name, 'reclassify', 'processed_pseudogene',
                                    f"{chrom}:?", f"max_tx_exons=1, ref_exons={ref_exons}"))
                continue

            # Rules B: drop weak duplicates when alt has strong support
            if not (has_alt and other_strong):
                continue

            drop_reason = None
            if modes == {'txTM'}:
                drop_reason = 'lone_txTM_with_strong_alt'
            elif modes == {'augPB'} and ref_exons > 0 and \
                 max_exons >= ref_exons * augpb_chimera_exon_ratio:
                drop_reason = (f'augpb_chimera_with_alt (tx_exons={max_exons},'
                               f' ref_exons={ref_exons})')
            elif total >= min_introns_for_low_support and \
                 supported / total < low_support_fraction and \
                 not (protect_strong_modes and (modes & STRONG_MODES)):
                drop_reason = (f'low_intron_support_with_alt '
                               f'({supported}/{total})')

            if drop_reason:
                drop_keys.append((chrom, gid))
                report_rows.append((gid, name, 'drop', drop_reason,
                                    f"{chrom}:?",
                                    f"modes={','.join(sorted(modes))}"))

    for chrom, gid in drop_keys:
        if gid in consensus_gene_dict[chrom]:
            del consensus_gene_dict[chrom][gid]
            n_drop += 1

    return n_reclass, n_drop, report_rows


# ---- file rewrite ---------------------------------------------------------

def rewrite_consensus_gp_and_info(consensus_gp, consensus_gp_info,
                                  consensus_gene_dict, genome):
    """Rewrite the consensus .gp and .gp_info files from the (post-processed)
    consensus_gene_dict. We preserve the existing transcript_id values stored
    in attrs['transcript_id']; we just re-emit lines in iteration order.
    """
    rows = []
    with open(consensus_gp, 'w') as out_gp:
        for chrom in consensus_gene_dict:
            for gene_id, entries in consensus_gene_dict[chrom].items():
                for tx_obj, attrs in entries:
                    name = attrs.get('transcript_id', tx_obj.name)
                    score = int(round(attrs.get('score', 0))) if isinstance(
                        attrs.get('score'), (int, float)) else 0
                    out_gp.write('\t'.join(
                        tx_obj.get_gene_pred(name=name, name2=gene_id,
                                              score=score)) + '\n')
                    gp_attrs = {k: v for k, v in attrs.items()
                                if not k.startswith('_')}
                    # Make sure gene_id matches the (possibly updated) key
                    gp_attrs['gene_id'] = gene_id
                    rows.append(gp_attrs)
    df = pd.DataFrame(rows).set_index(['gene_id', 'transcript_id'])
    if 'alternative_source_transcripts' not in df.columns:
        df['alternative_source_transcripts'] = 'N/A'
    with open(consensus_gp_info, 'w') as fh:
        df.to_csv(fh, sep='\t', na_rep='N/A')


# ---- public driver --------------------------------------------------------

def apply_postprocess(consensus_gene_dict, ref_gp,
                       consensus_gp=None, consensus_gp_info=None,
                       genome=None, report_path=None,
                       split_runaway=True, clean_pseudo=True,
                       min_introns_for_low_support=3,
                       augpb_chimera_exon_ratio=1.5,
                       low_support_fraction=0.3,
                       protect_strong_modes=True):
    """Run full postprocess on consensus_gene_dict and (optionally) rewrite
    the .gp / .gp_info files. Returns a stats dict.
    """
    stats = {'split': 0, 'reclassify': 0, 'drop': 0, 'report': []}
    ref_info = _load_ref_pc_info(ref_gp)
    if ref_info is None:
        logger.warning("Skipping consensus postprocess (no ref pc info)")
        return stats

    if split_runaway:
        n_split = split_runaway_pc_genes(consensus_gene_dict, ref_info)
        stats['split'] = n_split
        if n_split:
            logger.info(f"  consensus postprocess: split {n_split:,} "
                        f"runaway pc genes")

    if clean_pseudo:
        n_rec, n_drop, rows = reclassify_and_drop_pc_records(
            consensus_gene_dict, ref_info,
            min_introns_for_low_support=min_introns_for_low_support,
            augpb_chimera_exon_ratio=augpb_chimera_exon_ratio,
            low_support_fraction=low_support_fraction,
            protect_strong_modes=protect_strong_modes,
        )
        stats['reclassify'] = n_rec
        stats['drop']       = n_drop
        stats['report']     = rows
        if n_rec or n_drop:
            logger.info(f"  consensus postprocess: reclassified {n_rec:,} "
                        f"as processed_pseudogene; dropped {n_drop:,} "
                        f"redundant pc records")

    if report_path and stats['report']:
        try:
            with open(report_path, 'w') as fh:
                fh.write("gene_id\tgene_name\taction\tdetail\tlocation\tcontext\n")
                for row in sorted(stats['report'], key=lambda r: (r[2], r[1])):
                    fh.write('\t'.join(str(x) for x in row) + '\n')
            logger.info(f"  consensus postprocess: wrote report to {report_path}")
        except Exception as exc:
            logger.warning(f"  failed to write report: {exc}")

    if consensus_gp and (stats['split'] or stats['reclassify'] or stats['drop']):
        rewrite_consensus_gp_and_info(consensus_gp, consensus_gp_info,
                                       consensus_gene_dict, genome)
        logger.info(f"  consensus postprocess: rewrote {consensus_gp}")

    return stats
