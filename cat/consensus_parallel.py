import tools.bio
import argparse
import collections
import json
import logging
import os
import pickle
import re
import subprocess
import sys
import time
import warnings
import pandas as pd
import multiprocessing as mp
from functools import partial
from pathlib import Path
import tools.fileOps
import tools.intervals
import tools.mathOps
import tools.misc
import tools.nameConversions
import tools.procOps
import tools.sqlInterface
import tools.transcripts
from tools.defaultOrderedDict import DefaultOrderedDict
from sqlalchemy import inspect

# Suppress warnings for better performance
warnings.filterwarnings('ignore')
warnings.filterwarnings('ignore', category=FutureWarning, module='pandas')

# De novo mode prefixes recognized by the consensus pipeline
DENOVO_PREFIXES = ('augPB-', 'strg-')

# Higher value wins when resolving overlapping genes at the same assembly locus.
# augMP is lowest among map-based modes (rescue only when txTM/transMap absent).
OVERLAP_RESOLUTION_SOURCE_PRIORITY = {
    'transMap': 3,
    'transMap_pairwise': 3,
    'txTM': 2,
    'augTM': 1,
    'augTMR': 1,
    'augTM_pairwise': 1,
    'augTMR_pairwise': 1,
    'augMP': 0,
}

# Reference gene biotypes eligible for non-PC rescue (transMap / txTM only).
RESCUE_REF_NONCODING_BIOTYPES = frozenset({
    'pseudogene',
    'transcribed_pseudogene',
    'lncRNA',
    'unprocessed_pseudogene',
})


def _is_ref_noncoding_gene(gene_id, gene_biotype_map):
    return _lookup_ref_gene_biotype(gene_id, gene_biotype_map) in RESCUE_REF_NONCODING_BIOTYPES


def _is_augmp_attrs(attrs):
    if attrs.get('alignment_mode') == 'augMP':
        return True
    for key in ('alignment_id', 'source_transcript', 'transcript_id'):
        val = attrs.get(key)
        if val and tools.nameConversions.aln_id_is_augustus_mp(str(val)):
            return True
    return False


def _effective_gene_biotypes(tx_list):
    """Biotypes used for overlap conflict resolution (augMP with CDS is protein_coding)."""
    biotypes = set()
    for tx_obj, _tx_id, attrs in tx_list:
        bt = attrs.get('gene_biotype') or attrs.get('transcript_biotype') or 'unknown'
        if _is_augmp_attrs(attrs) and getattr(tx_obj, 'cds_size', 0) > 0:
            bt = 'protein_coding'
        biotypes.add(bt)
    return biotypes


def _interval_is_protein_coding(tx_list):
    return 'protein_coding' in _effective_gene_biotypes(tx_list)


def _different_biotypes_allow_coexistence(tx_list_a, tx_list_b):
    """True only when overlap of unlike kinds (e.g. PC + lncRNA) is intentional."""
    if _interval_is_protein_coding(tx_list_a) and _interval_is_protein_coding(tx_list_b):
        return False
    a = _effective_gene_biotypes(tx_list_a)
    b = _effective_gene_biotypes(tx_list_b)
    return len(a & b) == 0


def _is_denovo(aln_id):
    """Check if an alignment ID belongs to a de novo mode (augPB or strg)."""
    return isinstance(aln_id, str) and aln_id.startswith(DENOVO_PREFIXES)


def norm_ensg(gene_id):
    if gene_id is None or (isinstance(gene_id, float) and pd.isna(gene_id)):
        return None
    return re.sub(r'\.\d+$', '', str(gene_id))


def _lookup_ref_gene_biotype(gene_id, gene_biotype_map):
    if gene_id is None or (isinstance(gene_id, float) and pd.isna(gene_id)):
        return None
    gid = str(gene_id)
    biotype = gene_biotype_map.get(gid)
    if biotype is None:
        biotype = gene_biotype_map.get(norm_ensg(gid))
    return biotype


def _is_ref_protein_coding_gene(gene_id, gene_biotype_map):
    return _lookup_ref_gene_biotype(gene_id, gene_biotype_map) == 'protein_coding'


def apply_reference_gene_biotype_policy(final_consensus, gene_biotype_map, metrics=None):
    """Align consensus biotypes with reference; fix spurious PC and mis-called pseudogenes."""
    if not gene_biotype_map:
        return final_consensus
    out = []
    demoted = 0
    restored_pc = 0
    for aln_id, attrs in final_consensus:
        attrs = dict(attrs)
        if _is_denovo(aln_id):
            out.append((aln_id, attrs))
            continue
        sg = attrs.get('source_gene')
        ref_biotype = _lookup_ref_gene_biotype(sg, gene_biotype_map)
        if not ref_biotype:
            if _is_augmp_attrs(attrs):
                attrs['gene_biotype'] = 'protein_coding'
                attrs['transcript_biotype'] = 'protein_coding'
            out.append((aln_id, attrs))
            continue
        if ref_biotype != 'protein_coding' and attrs.get('gene_biotype') == 'protein_coding':
            attrs['gene_biotype'] = ref_biotype
            if attrs.get('transcript_biotype') == 'protein_coding':
                attrs['transcript_biotype'] = ref_biotype
            demoted += 1
        elif (
            ref_biotype == 'protein_coding'
            and attrs.get('transcript_class') == 'processed_pseudogene'
        ):
            attrs['transcript_class'] = 'ortholog'
            attrs['gene_biotype'] = 'protein_coding'
            attrs['transcript_biotype'] = 'protein_coding'
            restored_pc += 1
        out.append((aln_id, attrs))
    if demoted:
        logger.info(
            f"  Demoted {demoted} transcripts from protein_coding using reference gene biotypes"
        )
        if metrics is not None:
            metrics['Demoted spurious protein_coding by ref biotype'] = (
                metrics.get('Demoted spurious protein_coding by ref biotype', 0) + demoted
            )
    if restored_pc:
        logger.info(
            f"  Restored {restored_pc} ref protein_coding genes misclassified as processed_pseudogene"
        )
        if metrics is not None:
            metrics['Restored ref PC from processed_pseudogene'] = (
                metrics.get('Restored ref PC from processed_pseudogene', 0) + restored_pc
            )
    return out


def _ref_gene_span_coords(norm_gid, ref_gene_coords):
    if norm_gid in ref_gene_coords:
        return ref_gene_coords[norm_gid]
    for gid, coords in ref_gene_coords.items():
        if norm_ensg(gid) == norm_gid:
            return coords
    return None


def norm_ref_transcript_id(tx_id):
    """Reference / txTM transcript accession for matching (versionless)."""
    return re.sub(r'\.\d+$', '', str(tx_id))


def norm_match_transcript_id(tx_id):
    """Accession key for rt2t comparison and isoform rescue (keep version suffix)."""
    raw = str(tx_id)
    raw = re.sub(r'^augMP-rna-', '', raw)
    raw = re.sub(r'^rna-', '', raw)
    raw = re.sub(r'^txTM-', '', raw)
    raw = re.sub(r'_\d+$', '', raw)
    return raw


def build_alignment_coverage_map(mrna_metrics_df, tx_dict=None, alignment_source_map=None):
    """AlignmentId -> AlnCoverage_mRNA (percent) from metrics DB load."""
    cov = {}

    def _add(aid, val):
        try:
            v = float(val)
        except (TypeError, ValueError):
            return
        if pd.isna(v):
            return
        aid = str(aid)
        cov[aid] = v
        cov[norm_ref_transcript_id(aid)] = v

    if mrna_metrics_df is not None and len(mrna_metrics_df) > 0:
        if 'AlnCoverage_mRNA' in mrna_metrics_df.columns:
            for _, row in mrna_metrics_df.iterrows():
                _add(row['AlignmentId'], row['AlnCoverage_mRNA'])
        if 'AlnCoverage' in mrna_metrics_df.columns:
            for _, row in mrna_metrics_df.iterrows():
                _add(row['AlignmentId'], row['AlnCoverage'])
        elif 'classifier' in mrna_metrics_df.columns and 'value' in mrna_metrics_df.columns:
            sub = mrna_metrics_df[mrna_metrics_df['classifier'] == 'AlnCoverage']
            for _, row in sub.iterrows():
                _add(row['AlignmentId'], row['value'])

    # Map txTM CNV copy IDs (e.g. XM_…_1) to base metrics rows (XM_…).
    if tx_dict and alignment_source_map:
        for aln_id in tx_dict:
            mode = alignment_source_map.get(aln_id, '')
            if mode != 'txTM':
                continue
            base = norm_ref_transcript_id(normalize_alignment_id(aln_id, mode))
            if base in cov and aln_id not in cov:
                cov[aln_id] = cov[base]
    return cov


def _lookup_tx_coverage(coverage_map, aln_id, tx_obj, mode=None):
    if not coverage_map:
        return None
    keys = [aln_id, getattr(tx_obj, 'name', None)]
    if mode:
        keys.append(normalize_alignment_id(aln_id, mode))
    for key in keys:
        if not key:
            continue
        val = coverage_map.get(str(key))
        if val is not None:
            return val
        val = coverage_map.get(norm_ref_transcript_id(str(key)))
        if val is not None:
            return val
        if mode == 'txTM':
            val = coverage_map.get(
                norm_ref_transcript_id(normalize_alignment_id(str(key), mode))
            )
            if val is not None:
                return val
    return None


def _index_rescue_candidates(tx_dict, alignment_source_map):
    """Map versionless transcript name -> [(aln_id, tx_obj, mode)]."""
    by_tx = collections.defaultdict(list)
    rescue_modes = {'transMap', 'transMap_pairwise', 'txTM'}
    for aln_id, tx_obj in tx_dict.items():
        mode = alignment_source_map.get(aln_id, '')
        if mode not in rescue_modes:
            continue
        base_id = normalize_alignment_id(aln_id, mode) if mode == 'txTM' else tx_obj.name
        base = norm_match_transcript_id(base_id)
        by_tx[base].append((aln_id, tx_obj, mode))
    return by_tx


def _pick_rescue_candidate(candidates, present_aln, coverage_map, min_txTM_cov, mode_priority):
    """Best transMap/txTM candidate not already in consensus."""
    best = None
    for aln_id, tx_obj, mode in candidates:
        if aln_id in present_aln:
            continue
        if mode == 'txTM' and min_txTM_cov > 0:
            cov = _lookup_tx_coverage(coverage_map, aln_id, tx_obj, mode=mode)
            if cov is None or cov < min_txTM_cov:
                continue
        tx_len = max(1, int(tx_obj.stop) - int(tx_obj.start))
        pri = mode_priority[mode]
        key = (pri, -tx_len)
        if best is None or key < best[0]:
            best = (key, aln_id, tx_obj, mode)
    return best


_REF_NAME_LOOKUP_CACHE = {}


def _ref_name_lookups(ref_df):
    """Build (and cache) O(1) name lookups from ref_df so per-transcript rescue
    doesn't re-scan the whole reference DataFrame with regex on every call.

    Returns (gene_name_by_ensg, tx_name_exact, tx_name_base) where the keys exactly
    reproduce the previous per-call matching:
      * gene key   = version-stripped GeneId  (matched against norm_ensg(gene_id))
      * tx exact   = full TranscriptId
      * tx base    = version-stripped TranscriptId
    Cached per ref_df object (one ref_df per consensus process).
    """
    if ref_df is None:
        return {}, {}, {}
    key = id(ref_df)
    cached = _REF_NAME_LOOKUP_CACHE.get(key)
    if cached is not None:
        return cached
    gene_name_by_ensg = {}
    tx_name_exact = {}
    tx_name_base = {}
    if 'GeneId' in ref_df.columns and 'GeneName' in ref_df.columns:
        gids_base = ref_df['GeneId'].astype(str).str.replace(r'\.\d+$', '', regex=True)
        for gb, gn in zip(gids_base, ref_df['GeneName']):
            if gb not in gene_name_by_ensg:
                gene_name_by_ensg[gb] = gn
    if 'TranscriptId' in ref_df.columns and 'TranscriptName' in ref_df.columns:
        tids = ref_df['TranscriptId'].astype(str)
        tids_base = tids.str.replace(r'\.[0-9]+$', '', regex=True)
        for te, tb, tn in zip(tids, tids_base, ref_df['TranscriptName']):
            if te not in tx_name_exact:
                tx_name_exact[te] = tn
            if tb not in tx_name_base:
                tx_name_base[tb] = tn
    result = (gene_name_by_ensg, tx_name_exact, tx_name_base)
    _REF_NAME_LOOKUP_CACHE[key] = result
    return result


def _build_rescue_support_row(
    aln_id,
    tx_obj,
    gene_id,
    ref_df,
    coverage_map=None,
    mode=None,
    gene_biotype='protein_coding',
    transcript_biotype=None,
):
    """Minimal support row for create_transcript_attributes on rescued transMap models."""
    if transcript_biotype is None:
        transcript_biotype = gene_biotype
    gene_name = None
    transcript_name = None
    if ref_df is not None and gene_id:
        gene_name_by_ensg, tx_name_exact, tx_name_base = _ref_name_lookups(ref_df)
        gene_name = gene_name_by_ensg.get(norm_ensg(gene_id))
        tx_key = str(tx_obj.name)
        transcript_name = tx_name_exact.get(tx_key)
        if transcript_name is None:
            tx_base = re.sub(r'\.[0-9]+$', '', tx_key)
            transcript_name = tx_name_base.get(tx_base)
    return pd.Series({
        'AlignmentId': aln_id,
        'GeneId': gene_id,
        'TranscriptId': tx_obj.name,
        'TranscriptBiotype': transcript_biotype,
        'GeneBiotype': gene_biotype,
        'GeneName': gene_name,
        'TranscriptName': transcript_name,
        'AlnGoodness_mRNA': 1.0,
        'AlnCoverage_mRNA': float(
            _lookup_tx_coverage(coverage_map, aln_id, tx_obj, mode=mode) or 100.0
        ) if coverage_map else 100.0,
        'AlnIdentity_mRNA': 100.0,
        'TranscriptScore': 0,
        'Frameshift': False,
        'ValidStart': True,
        'ValidStop': True,
        'ProperOrf': True,
        'IntronRnaSupportPercent': 0,
        'ExonRnaSupportPercent': 0,
        'IntronAnnotSupportPercent': 0,
        'ExonAnnotSupportPercent': 0,
        'ExonAnnotSupport': [],
        'IntronAnnotSupport': [],
        'ExonRnaSupport': [],
        'IntronRnaSupport': [],
    })


def rescue_missing_reference_pc_genes(
    final_consensus,
    tx_dict,
    alignment_source_map,
    ref_gene_coords,
    gene_biotype_map,
    ref_df,
    args,
    mrna_metrics_df=None,
    metrics=None,
):
    """
    Re-add reference protein_coding genes that have a passing transMap/txTM model in inputs
    but were dropped during overlap resolution or never selected.
    """
    ref_pc = {norm_ensg(g) for g, b in gene_biotype_map.items() if b == 'protein_coding'}
    present_pc = set()
    for _, attrs in final_consensus:
        if attrs.get('gene_biotype') != 'protein_coding':
            continue
        sg = attrs.get('source_gene')
        if sg and sg != 'N/A':
            present_pc.add(norm_ensg(sg))

    missing = sorted(ref_pc - present_pc)
    if not missing:
        return final_consensus

    min_ratio = float(getattr(args, 'min_pc_len_ratio_vs_reference', 0.0) or 0.0)
    mode_priority = {'transMap': 0, 'transMap_pairwise': 1, 'txTM': 2}
    coverage_map = build_alignment_coverage_map(
        mrna_metrics_df, tx_dict=tx_dict, alignment_source_map=alignment_source_map
    )
    rescued = []

    # Index transMap/txTM transcripts by ENSG (avoid O(missing * len(tx_dict)) scan).
    logger.info(
        f"  Rescue pass: {len(missing)} ref PC genes missing after filtering; "
        f"indexing transMap/txTM candidates..."
    )
    by_ensg = collections.defaultdict(list)
    for aln_id, tx_obj in tx_dict.items():
        mode = alignment_source_map.get(aln_id, '')
        if mode not in mode_priority:
            continue
        g = norm_ensg(tx_obj.name2 or '')
        if g:
            by_ensg[g].append((aln_id, tx_obj, mode))
    logger.info(f"  Rescue index: {sum(len(v) for v in by_ensg.values())} candidate transcripts")

    gene_id_by_norm = {}
    if ref_df is not None and 'GeneId' in ref_df.columns:
        for gid in ref_df['GeneId'].dropna().unique():
            gene_id_by_norm[norm_ensg(gid)] = gid

    for norm_gid in missing:
        ref_coords = _ref_gene_span_coords(norm_gid, ref_gene_coords)
        if ref_coords:
            ref_chrom, ref_start, ref_end = ref_coords
            ref_len = max(1, int(ref_end) - int(ref_start))
        else:
            ref_len = None

        candidates = by_ensg.get(norm_gid, [])
        has_transmap = any(m in ('transMap', 'transMap_pairwise') for _, _, m in candidates)
        orphan_min_ratio = float(
            getattr(args, 'min_pc_len_ratio_txTM_only_rescue', min_ratio) or min_ratio
        )
        ratio_threshold = min_ratio if has_transmap else orphan_min_ratio

        best = None
        for aln_id, tx_obj, mode in candidates:
            tx_len = max(1, int(tx_obj.stop) - int(tx_obj.start))
            if ref_len is not None and ratio_threshold > 0:
                ratio = tx_len / ref_len
                if ratio < ratio_threshold:
                    continue
            pri = mode_priority[mode]
            ratio_for_key = (tx_len / ref_len) if ref_len else 1.0
            key = (pri, -ratio_for_key, -tx_len)
            if best is None or key < best[0]:
                gene_id = gene_id_by_norm.get(norm_gid, tx_obj.name2)
                best = (key, aln_id, tx_obj, mode, gene_id)

        if best is None:
            continue
        _, aln_id, tx_obj, mode, gene_id = best
        row = _build_rescue_support_row(aln_id, tx_obj, gene_id, ref_df, coverage_map, mode=mode)
        norm_tx_id = tools.nameConversions.strip_alignment_numbers(aln_id)
        attrs = create_transcript_attributes(
            row, mode, norm_tx_id, tx_dict=tx_dict, args=args, ref_df=ref_df
        )
        if attrs.get('transcript_class') is None:
            attrs['transcript_class'] = 'ortholog'
        attrs['gene_biotype'] = 'protein_coding'
        attrs['transcript_biotype'] = 'protein_coding'
        attrs['source_gene'] = gene_id
        attrs['alignment_mode'] = mode
        rescued.append((aln_id, attrs))

    if not rescued:
        return final_consensus

    logger.info(
        f"  Rescued {len(rescued)} reference protein_coding genes from transMap/txTM inputs "
        f"({len(missing)} were missing after filtering)"
    )
    if metrics is not None:
        metrics['Rescued missing ref PC genes'] = len(rescued)
    return list(final_consensus) + rescued


def rescue_missing_reference_noncoding_genes(
    final_consensus,
    tx_dict,
    alignment_source_map,
    ref_gene_coords,
    gene_biotype_map,
    ref_df,
    args,
    mrna_metrics_df=None,
    metrics=None,
):
    """
    Re-add reference pseudogene / lncRNA genes with a passing transMap/txTM model that were
    dropped during selection or overlap resolution. Does not use augMP.
    """
    if not getattr(args, 'rescue_reference_noncoding_genes', True):
        return final_consensus
    if not gene_biotype_map:
        return final_consensus

    ref_nc = {
        norm_ensg(g)
        for g, b in gene_biotype_map.items()
        if b in RESCUE_REF_NONCODING_BIOTYPES
    }
    present_nc = set()
    for _, attrs in final_consensus:
        sg = attrs.get('source_gene')
        if not sg or sg == 'N/A':
            continue
        if _is_ref_noncoding_gene(sg, gene_biotype_map):
            present_nc.add(norm_ensg(sg))

    missing = sorted(ref_nc - present_nc)
    if not missing:
        return final_consensus

    min_ratio = float(getattr(args, 'min_nc_len_ratio_vs_reference', 0.0) or 0.0)
    orphan_min_ratio = float(
        getattr(args, 'min_nc_len_ratio_txTM_only_rescue', min_ratio) or min_ratio
    )
    min_txTM_cov = float(getattr(args, 'rescue_min_txTM_coverage_noncoding', 50) or 0)
    mode_priority = {'transMap': 0, 'transMap_pairwise': 1, 'txTM': 2}
    coverage_map = build_alignment_coverage_map(
        mrna_metrics_df, tx_dict=tx_dict, alignment_source_map=alignment_source_map
    )
    rescued = []

    logger.info(
        f"  Non-PC rescue: {len(missing)} ref non-coding genes missing; "
        f"indexing transMap/txTM candidates..."
    )
    by_ensg = collections.defaultdict(list)
    for aln_id, tx_obj in tx_dict.items():
        mode = alignment_source_map.get(aln_id, '')
        if mode not in mode_priority:
            continue
        g = norm_ensg(tx_obj.name2 or '')
        if g:
            by_ensg[g].append((aln_id, tx_obj, mode))

    gene_id_by_norm = {}
    biotype_by_norm = {}
    if ref_df is not None and 'GeneId' in ref_df.columns:
        for gid in ref_df['GeneId'].dropna().unique():
            ng = norm_ensg(gid)
            gene_id_by_norm[ng] = gid
    for g, b in gene_biotype_map.items():
        ng = norm_ensg(g)
        if b in RESCUE_REF_NONCODING_BIOTYPES:
            biotype_by_norm[ng] = b

    for norm_gid in missing:
        ref_biotype = biotype_by_norm.get(norm_gid, 'pseudogene')
        ref_coords = _ref_gene_span_coords(norm_gid, ref_gene_coords)
        if ref_coords:
            ref_chrom, ref_start, ref_end = ref_coords
            ref_len = max(1, int(ref_end) - int(ref_start))
        else:
            ref_len = None

        candidates = by_ensg.get(norm_gid, [])
        has_transmap = any(m in ('transMap', 'transMap_pairwise') for _, _, m in candidates)
        ratio_threshold = min_ratio if has_transmap else orphan_min_ratio

        best = None
        for aln_id, tx_obj, mode in candidates:
            if mode == 'txTM' and min_txTM_cov > 0:
                cov = _lookup_tx_coverage(coverage_map, aln_id, tx_obj, mode=mode)
                if cov is None or cov < min_txTM_cov:
                    continue
            tx_len = max(1, int(tx_obj.stop) - int(tx_obj.start))
            if ref_len is not None and ratio_threshold > 0:
                if tx_len / ref_len < ratio_threshold:
                    continue
            pri = mode_priority[mode]
            ratio_for_key = (tx_len / ref_len) if ref_len else 1.0
            key = (pri, -ratio_for_key, -tx_len)
            if best is None or key < best[0]:
                gene_id = gene_id_by_norm.get(norm_gid, tx_obj.name2)
                best = (key, aln_id, tx_obj, mode, gene_id)

        if best is None:
            continue
        _, aln_id, tx_obj, mode, gene_id = best
        row = _build_rescue_support_row(
            aln_id,
            tx_obj,
            gene_id,
            ref_df,
            coverage_map,
            mode=mode,
            gene_biotype=ref_biotype,
            transcript_biotype=ref_biotype,
        )
        norm_tx_id = tools.nameConversions.strip_alignment_numbers(aln_id)
        attrs = create_transcript_attributes(
            row, mode, norm_tx_id, tx_dict=tx_dict, args=args, ref_df=ref_df
        )
        if attrs.get('transcript_class') is None:
            attrs['transcript_class'] = 'ortholog'
        attrs['gene_biotype'] = ref_biotype
        attrs['transcript_biotype'] = ref_biotype
        attrs['source_gene'] = gene_id
        attrs['alignment_mode'] = mode
        rescued.append((aln_id, attrs))

    if not rescued:
        return final_consensus

    logger.info(
        f"  Rescued {len(rescued)} reference non-coding genes from transMap/txTM "
        f"({len(missing)} were missing after filtering)"
    )
    if metrics is not None:
        metrics['Rescued missing ref non-coding genes'] = len(rescued)
    return list(final_consensus) + rescued


def rescue_missing_reference_transcripts(
    final_consensus,
    tx_dict,
    alignment_source_map,
    ref_df,
    gene_biotype_map,
    args,
    mrna_metrics_df=None,
    metrics=None,
):
    """
    Re-add reference PC transcript isoforms present in transMap/txTM inputs but dropped
    during per-locus selection (not a global dedup of paralogs).

    Uses a separate, typically lower, txTM coverage floor than the pre-score filter so
    isoforms in the 70–79% band can be recovered without reopening global noise.
    """
    if ref_df is None or len(ref_df) == 0:
        return final_consensus
    if not getattr(args, 'rescue_reference_isoforms', True):
        return final_consensus

    min_txTM_cov = float(getattr(args, 'rescue_min_txTM_coverage', 80) or 0)
    mode_priority = {'transMap': 0, 'transMap_pairwise': 1, 'txTM': 2}
    coverage_map = build_alignment_coverage_map(
        mrna_metrics_df, tx_dict=tx_dict, alignment_source_map=alignment_source_map
    )

    present_aln = {aln for aln, _ in final_consensus}
    present_tx = set()
    for _, attrs in final_consensus:
        st = attrs.get('source_transcript') or attrs.get('normalized_transcript_id')
        if st:
            present_tx.add(norm_match_transcript_id(st))

    ref_pc = ref_df[ref_df['GeneBiotype'] == 'protein_coding'].copy()
    if 'TranscriptBiotype' in ref_pc.columns:
        ref_pc = ref_pc[
            ref_pc['TranscriptBiotype'].fillna('protein_coding') == 'protein_coding'
        ]
    ref_pc['tx_base'] = ref_pc['TranscriptId'].astype(str).map(norm_match_transcript_id)

    missing_tx = ref_pc[~ref_pc['tx_base'].isin(present_tx)]
    if len(missing_tx) == 0:
        return final_consensus

    by_tx = _index_rescue_candidates(tx_dict, alignment_source_map)
    gene_id_by_tx = dict(
        zip(
            ref_pc['tx_base'],
            ref_pc['GeneId'].astype(str),
        )
    )

    rescued = []
    skipped_no_pick = 0
    skipped_no_candidate = 0

    for tx_base in missing_tx['tx_base'].unique():
        candidates = by_tx.get(tx_base, [])
        if not candidates:
            skipped_no_candidate += 1
            continue

        pick = _pick_rescue_candidate(
            candidates, present_aln, coverage_map, min_txTM_cov, mode_priority
        )
        if pick is None:
            skipped_no_pick += 1
            continue

        _, aln_id, tx_obj, mode = pick
        gene_id = gene_id_by_tx.get(tx_base, tx_obj.name2)
        row = _build_rescue_support_row(aln_id, tx_obj, gene_id, ref_df, coverage_map, mode=mode)
        norm_tx_id = tools.nameConversions.strip_alignment_numbers(aln_id)
        attrs = create_transcript_attributes(
            row, mode, norm_tx_id, tx_dict=tx_dict, args=args, ref_df=ref_df
        )
        attrs['transcript_class'] = attrs.get('transcript_class') or 'ortholog'
        attrs['gene_biotype'] = 'protein_coding'
        attrs['transcript_biotype'] = 'protein_coding'
        attrs['source_gene'] = str(gene_id)
        attrs['source_transcript'] = norm_match_transcript_id(tx_obj.name)
        attrs['normalized_transcript_id'] = attrs['source_transcript']
        attrs['alignment_mode'] = mode
        rescued.append((aln_id, attrs))
        present_aln.add(aln_id)
        present_tx.add(tx_base)

    if not rescued:
        return final_consensus

    logger.info(
        f"  Rescued {len(rescued)} reference PC transcript isoforms "
        f"({len(missing_tx)} were missing from consensus; "
        f"{skipped_no_candidate} had no transMap/txTM candidate, "
        f"{skipped_no_pick} had candidates but none selected "
        f"(coverage <{min_txTM_cov} for txTM-only, or already present))"
    )
    if metrics is not None:
        metrics['Rescued missing ref PC transcripts'] = len(rescued)
    return list(final_consensus) + rescued


def rescue_alternative_source_isoforms(
    final_consensus,
    tx_dict,
    alignment_source_map,
    ref_df,
    args,
    mrna_metrics_df=None,
    metrics=None,
):
    """
    Emit separate transcript records for reference isoforms listed only in
    alternative_source_transcripts on a kept gene (common when augMP wins a locus).
    """
    if not getattr(args, 'rescue_reference_isoforms', True):
        return final_consensus
    if not getattr(args, 'rescue_alternative_isoforms', True):
        return final_consensus

    min_txTM_cov = float(getattr(args, 'rescue_min_txTM_coverage', 70) or 0)
    mode_priority = {'transMap': 0, 'transMap_pairwise': 1, 'txTM': 2}
    coverage_map = build_alignment_coverage_map(
        mrna_metrics_df, tx_dict=tx_dict, alignment_source_map=alignment_source_map
    )
    by_tx = _index_rescue_candidates(tx_dict, alignment_source_map)

    present_aln = {aln for aln, _ in final_consensus}
    present_tx = set()
    for _, attrs in final_consensus:
        st = attrs.get('source_transcript') or attrs.get('normalized_transcript_id')
        if st:
            present_tx.add(norm_match_transcript_id(st))

    gene_id_by_tx = {}
    if ref_df is not None and 'TranscriptId' in ref_df.columns and 'GeneId' in ref_df.columns:
        ref_pc = ref_df[ref_df['GeneBiotype'] == 'protein_coding']
        gene_id_by_tx = dict(
            zip(
                ref_pc['TranscriptId'].astype(str).map(norm_match_transcript_id),
                ref_pc['GeneId'].astype(str),
            )
        )

    rescued = []
    skipped_no_candidate = 0
    skipped_no_pick = 0

    for _, host_attrs in final_consensus:
        alt_field = host_attrs.get('alternative_source_transcripts', '')
        if not alt_field or alt_field == 'N/A':
            continue
        host_gene = host_attrs.get('source_gene')
        host_gene_name = host_attrs.get('source_gene_common_name') or host_attrs.get('gene_name')

        for raw in str(alt_field).split(','):
            raw = raw.strip()
            if not raw:
                continue
            tx_base = norm_match_transcript_id(raw)
            if tx_base in present_tx:
                continue

            candidates = by_tx.get(tx_base, [])
            if not candidates:
                skipped_no_candidate += 1
                continue

            pick = _pick_rescue_candidate(
                candidates, present_aln, coverage_map, min_txTM_cov, mode_priority
            )
            if pick is None:
                skipped_no_pick += 1
                continue

            _, aln_id, tx_obj, mode = pick
            gene_id = gene_id_by_tx.get(tx_base, host_gene or tx_obj.name2)
            row = _build_rescue_support_row(aln_id, tx_obj, gene_id, ref_df, coverage_map, mode=mode)
            norm_tx_id = tools.nameConversions.strip_alignment_numbers(aln_id)
            attrs = create_transcript_attributes(
                row, mode, norm_tx_id, tx_dict=tx_dict, args=args, ref_df=ref_df
            )
            attrs['transcript_class'] = attrs.get('transcript_class') or 'ortholog'
            attrs['gene_biotype'] = 'protein_coding'
            attrs['transcript_biotype'] = 'protein_coding'
            attrs['source_gene'] = str(gene_id)
            attrs['source_gene_common_name'] = host_gene_name
            attrs['source_transcript'] = norm_match_transcript_id(tx_obj.name)
            attrs['normalized_transcript_id'] = attrs['source_transcript']
            attrs['alignment_mode'] = mode
            rescued.append((aln_id, attrs))
            present_aln.add(aln_id)
            present_tx.add(tx_base)

    if not rescued:
        return final_consensus

    logger.info(
        f"  Promoted {len(rescued)} isoforms from alternative_source_transcripts "
        f"({skipped_no_candidate} had no transMap/txTM candidate, "
        f"{skipped_no_pick} had candidates but none selected)"
    )
    if metrics is not None:
        metrics['Promoted alternative ref isoforms'] = len(rescued)
    return list(final_consensus) + rescued


warnings.filterwarnings('ignore', category=DeprecationWarning, module='pandas')

logger = logging.getLogger(__name__)
ID_TEMPLATE = '{genome:.10}_{tag_type}{unique_id:07d}'


def main():
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    args = parser.parse_args()
    
    logger.info("="*80)
    logger.info(f"Starting consensus generation for {args.genome}")
    logger.info("="*80)
    
    # Run the main consensus logic
    metrics = generate_consensus(args)

    # Write the final metrics file
    logger.info(f"Writing metrics to {args.metrics_json}")
    with open(args.metrics_json, 'w') as outf:
        json.dump(metrics, outf, indent=4)
    logger.info(f"Successfully generated consensus gene set and metrics for {args.genome}.")
    logger.info("="*80)


def add_arguments(parser):
    """Adds all command-line arguments to the argparse parser."""
    # --- Input Files ---
    parser.add_argument("--gp-list", nargs='+', required=True, help="Space-separated list of all genePred files to consider.")
    parser.add_argument("--db-path", required=True, help="Path to the genome's primary database.")
    parser.add_argument("--ref-db-path", required=True, help="Path to the reference genome's database.")
    parser.add_argument("--fasta", required=True, help="Path to the target genome FASTA file.")
    parser.add_argument("--genome", required=True, help="Name of the target genome.")
    
    # --- Output Files ---
    parser.add_argument("--consensus-gp", required=True, help="Output path for consensus genePred.")
    parser.add_argument("--consensus-gp-info", required=True, help="Output path for consensus gp_info TSV.")
    parser.add_argument("--consensus-gff3", required=True, help="Output path for consensus GFF3.")
    parser.add_argument("--consensus-fasta", required=True, help="Output path for consensus transcript FASTA.")
    parser.add_argument("--protein-fasta", required=True, help="Output path for consensus protein FASTA.")
    parser.add_argument("--metrics-json", required=True, help="Output path for consensus metrics JSON.")
    
    # --- Control and Filtering Parameters ---
    parser.add_argument("--intron-rnaseq-support", type=float, default=0.0, help="Percent of introns that must be supported by RNA-seq.")
    parser.add_argument("--exon-rnaseq-support", type=float, default=0.0, help="Percent of exons supported by RNA-seq.")
    parser.add_argument("--intron-annot-support", type=float, default=0.0, help="Percent of introns supported by reference annotation.")
    parser.add_argument("--exon-annot-support", type=float, default=0.0, help="Percent of exons supported by reference annotation.")
    parser.add_argument("--original-intron-support", type=float, default=0.0, help="Percent of original introns that must be preserved.")
    parser.add_argument("--in-species-rna-support-only", action="store_true", help="Use in-species RNA-seq support only, ignoring cross-species evidence.")
    parser.add_argument("--filter-overlapping-genes", action="store_true", help="Filter out overlapping CDS intervals from different genes.")
    parser.add_argument("--overlapping-ignore-bases", type=int, default=0, help="Number of bases to ignore when clustering for overlap.")
    parser.add_argument("--cnv-score-similarity", type=float, default=0.80, help="Keep multi-locus transcripts if scores are within this fraction of max (default: 0.80 = 80%)")
    parser.add_argument(
        "--min-pc-len-ratio-vs-reference",
        type=float,
        default=0.4,
        help=(
            "Remove transcripts when their span is below this fraction of the reference source_gene span. "
            "Applies to all gene biotypes. Set to 0 to disable. Default: 0.4."
        ),
    )
    parser.add_argument(
        "--filter-spurious-pc-overlaps-not-in-reference",
        action="store_true",
        help=(
            "If set, remove overlapping protein-coding gene loci when the overlap is not supported by the "
            "reference gene coordinates (including different reference chromosomes). This is a strict cleanup "
            "pass to reduce spurious loci."
        ),
    )
    
    # --- Parallelization Parameters ---
    parser.add_argument("--num-workers", type=int, default=None, help="Number of parallel workers")
    
    # --- De Novo Parameters ---
    parser.add_argument("--denovo-tx-modes", nargs='*', default=[], help="List of de novo modes to consider (e.g., augPB).")
    parser.add_argument("--denovo-num-introns", type=int, default=1, help="A de novo isoform must have at least this many introns.")
    parser.add_argument("--denovo-splice-support", type=float, default=1.0, help="Percent of de novo splices that must be RNA-seq supported.")
    parser.add_argument("--denovo-exon-support", type=float, default=1.0, help="Percent of de novo exons that must be RNA-seq supported.")
    parser.add_argument("--denovo-ignore-novel-genes", action="store_true", help="If set, only incorporate de novo transcripts as novel isoforms, not novel genes.")
    parser.add_argument("--denovo-only-novel-genes", action="store_true", help="If set, only incorporate de novo transcripts if they are novel genes.")
    parser.add_argument("--denovo-allow-unsupported", action="store_true", help="Allow de novo transcripts with novel splices that lack RNA-seq support.")
    parser.add_argument("--denovo-allow-bad-annot-or-tm", action="store_true", help="Allow de novo models flagged as 'badAnnotOrTm'.")
    parser.add_argument("--denovo-allow-novel-ends", action="store_true", help="Allow de novo models with novel 5' or 3' ends.")
    parser.add_argument("--denovo-novel-end-distance", type=int, default=0)

    # --- PacBio Parameters ---
    parser.add_argument("--require-pacbio-support", action="store_true", help="If set, remove any consensus transcript not validated by Iso-Seq data.")
    parser.add_argument("--hints-db-has-rnaseq", action="store_true", help="Flag if the hints DB contains RNA-seq, for tagging purposes.")


def normalize_alignment_id(aln_id, mode):
    """
    Normalize alignment IDs according to mode:
    - transMap: keep as-is
    - txTM: strip underscore and numbers after (e.g., ENST00000123_1 -> ENST00000123)
    - augTM/augTMR: strip prefix (e.g., augTM-ENST00000123 -> ENST00000123)
    """
    if mode == 'transMap':
        return aln_id
    elif mode == 'txTM':
        # Strip cross-mode _cp suffix first, then txTM CNV copy _N (e.g. XM_….1_14_cp9 → XM_….1).
        import re
        base = re.sub(r'_cp\d*$', '', aln_id)
        base = re.sub(r'_\d+$', '', base)
        return base
    elif mode in ['augTM', 'augTMR', 'augMP']:
        # Strip the prefix
        if aln_id.startswith('augTM-'):
            return aln_id[6:]  # len('augTM-') = 6
        elif aln_id.startswith('augTMR-'):
            return aln_id[7:]  # len('augTMR-') = 7
        elif aln_id.startswith('augMP-'):
            return aln_id[6:]  # len('augMP-') = 6
        return aln_id
    elif mode in ('augPB', 'strg'):
        # Keep denovo IDs as-is
        return aln_id
    return aln_id


def normalize_gene_id(gene_id, mode):
    """
    Normalize gene IDs according to mode (same as transcript ID normalization).
    This is important for txTM CNV copies which append _N to both transcript and gene IDs.
    
    - transMap: keep as-is
    - txTM: strip _N suffix (where N is 1-20) OR _N_M pattern
    - augTM/augTMR/augMP: keep as-is
    - augPB: keep as-is
    
    Examples for txTM:
    - ENSG00000026103.25_1 → ENSG00000026103.25 (CNV copy, strip _1)
    - hg002_chrY_paternal_691_1 → hg002_chrY_paternal_691 (CNV copy, strip _1)
    - hg002_chrY_paternal_691 → hg002_chrY_paternal_691 (original, keep as-is)
    - hg002_chrY_paternal_157_10 → hg002_chrY_paternal_157 (CNV copy 10, strip _10)
    """
    if mode == 'txTM':
        import re
        # First, check for double _N_M pattern (e.g., _691_1 or _157_10)
        match = re.search(r'_\d+_(\d+)$', gene_id)
        if match:
            # Has double _N_M pattern, strip the last _M
            base = re.sub(r'_\d+$', '', gene_id)
            return base
        # Second, check for simple _N suffix (txTM CNV copies)
        # This handles ENSG00000026103.25_1 → ENSG00000026103.25
        # Match _<digits> where there's a version number (dot + digits) before it
        # This ensures we only strip CNV suffixes, not gene names ending in numbers
        match = re.search(r'\.\d+_(\d+)$', gene_id)
        if match:
            # Has version.number_CNV pattern (e.g., ENSG...25_21), strip the _CNV part
            base = gene_id[:match.end() - len(match.group(1)) - 1]
            return base
    # For all other modes or no matching pattern, keep gene ID as-is
    return gene_id


def identify_mode(aln_id, gp_file=None):
    """Identify the alignment mode from the alignment ID or file path"""
    if aln_id.startswith('augPB-'):
        return 'augPB'
    elif aln_id.startswith('strg-'):
        return 'strg'
    elif aln_id.startswith('augTMR-'):
        # Check if it's pairwise by looking at the file path
        if gp_file and 'augTMR_pairwise' in gp_file:
            return 'augTMR_pairwise'
        return 'augTMR'
    elif aln_id.startswith('augTM-'):
        # Check if it's pairwise by looking at the file path
        if gp_file and 'augTM_pairwise' in gp_file:
            return 'augTM_pairwise'
        return 'augTM'
    elif gp_file and 'txTM' in gp_file:
        return 'txTM'
    elif gp_file and 'transMap_pairwise' in gp_file:
        return 'transMap_pairwise'
    elif gp_file and 'transMap' in gp_file:
        return 'transMap'
    # Default fallback
    return 'transMap'


def backfill_transmap_multimapped_metrics(mrna_metrics_df, cds_metrics_df, valid_aln_ids, alignment_source_map):
    """
    Backfill metrics for transMap multi-mapped transcripts with _N suffix.
    TransMap stores metrics by base transcript ID (e.g., NM_004043.3) but
    alignment IDs in tx_dict have a _N suffix for each mapping location
    (e.g., NM_004043.3_1). Duplicate base metrics for all _N variants.
    """
    if len(mrna_metrics_df) == 0 and len(cds_metrics_df) == 0:
        return mrna_metrics_df, cds_metrics_df

    existing_mrna_ids = set(mrna_metrics_df['AlignmentId'].values) if len(mrna_metrics_df) > 0 else set()
    existing_cds_ids = set(cds_metrics_df['AlignmentId'].values) if len(cds_metrics_df) > 0 else set()

    missing_mrna = []
    missing_cds = []

    for aln_id in valid_aln_ids:
        mode = alignment_source_map.get(aln_id, 'unknown')
        if mode in ('transMap', 'transMap_pairwise'):
            m = re.match(r'^(.+)_(\d+)$', aln_id)
            if m:
                base_id = m.group(1)
                if aln_id not in existing_mrna_ids and base_id in existing_mrna_ids:
                    missing_mrna.append((aln_id, base_id))
                if aln_id not in existing_cds_ids and base_id in existing_cds_ids:
                    missing_cds.append((aln_id, base_id))

    if missing_mrna or missing_cds:
        n_unique = len(set([x[0] for x in missing_mrna + missing_cds]))
        logger.info(f"  Backfilling metrics for {n_unique} transMap multi-mapped transcripts...")

    if missing_mrna and len(mrna_metrics_df) > 0:
        new_rows = []
        for derived_id, base_id in missing_mrna:
            base_rows = mrna_metrics_df[mrna_metrics_df['AlignmentId'] == base_id]
            for _, row in base_rows.iterrows():
                new_row = row.copy()
                new_row['AlignmentId'] = derived_id
                new_rows.append(new_row)
        if new_rows:
            mrna_metrics_df = pd.concat([mrna_metrics_df, pd.DataFrame(new_rows)], ignore_index=True)
            logger.info(f"    ✓ Added {len(new_rows)} mRNA metric records for transMap multi-mapped transcripts")

    if missing_cds and len(cds_metrics_df) > 0:
        new_rows = []
        for derived_id, base_id in missing_cds:
            base_rows = cds_metrics_df[cds_metrics_df['AlignmentId'] == base_id]
            for _, row in base_rows.iterrows():
                new_row = row.copy()
                new_row['AlignmentId'] = derived_id
                new_rows.append(new_row)
        if new_rows:
            cds_metrics_df = pd.concat([cds_metrics_df, pd.DataFrame(new_rows)], ignore_index=True)
            logger.info(f"    ✓ Added {len(new_rows)} CDS metric records for transMap multi-mapped transcripts")

    return mrna_metrics_df, cds_metrics_df


def backfill_cnv_metrics(mrna_metrics_df, cds_metrics_df, valid_aln_ids, alignment_source_map, db_path):
    """
    Backfill metrics for txTM CNV copies that don't have metrics in the database.
    TxTM computes metrics for the base transcript but not each _N copy.
    We duplicate the base metrics for all CNV copies.
    """
    if len(mrna_metrics_df) == 0 and len(cds_metrics_df) == 0:
        return mrna_metrics_df, cds_metrics_df
    
    # Find txTM transcripts with _N suffix that are missing metrics
    existing_mrna_ids = set(mrna_metrics_df['AlignmentId'].values) if len(mrna_metrics_df) > 0 else set()
    existing_cds_ids = set(cds_metrics_df['AlignmentId'].values) if len(cds_metrics_df) > 0 else set()
    
    missing_mrna = []
    missing_cds = []
    
    for aln_id in valid_aln_ids:
        mode = alignment_source_map.get(aln_id, 'unknown')
        if mode == 'txTM' and '_' in aln_id.split('.')[-1]:
            # This is a txTM CNV copy (_N suffix)
            base_id = normalize_alignment_id(aln_id, mode)  # Strips _N
            
            # Check if this copy is missing metrics but base has them
            if aln_id not in existing_mrna_ids and base_id in existing_mrna_ids:
                missing_mrna.append((aln_id, base_id))
            
            if aln_id not in existing_cds_ids and base_id in existing_cds_ids:
                missing_cds.append((aln_id, base_id))
    
    if missing_mrna or missing_cds:
        logger.info(f"  Backfilling metrics for {len(set([x[0] for x in missing_mrna + missing_cds]))} txTM CNV copies...")
    
    # Duplicate metrics for missing CNV copies
    if missing_mrna and len(mrna_metrics_df) > 0:
        new_rows = []
        for cnv_id, base_id in missing_mrna:
            base_rows = mrna_metrics_df[mrna_metrics_df['AlignmentId'] == base_id]
            for _, row in base_rows.iterrows():
                new_row = row.copy()
                new_row['AlignmentId'] = cnv_id
                new_rows.append(new_row)
        if new_rows:
            mrna_metrics_df = pd.concat([mrna_metrics_df, pd.DataFrame(new_rows)], ignore_index=True)
            logger.info(f"    ✓ Added {len(new_rows)} mRNA metric records for CNV copies")
    
    if missing_cds and len(cds_metrics_df) > 0:
        new_rows = []
        for cnv_id, base_id in missing_cds:
            base_rows = cds_metrics_df[cds_metrics_df['AlignmentId'] == base_id]
            for _, row in base_rows.iterrows():
                new_row = row.copy()
                new_row['AlignmentId'] = cnv_id
                new_rows.append(new_row)
        if new_rows:
            cds_metrics_df = pd.concat([cds_metrics_df, pd.DataFrame(new_rows)], ignore_index=True)
            logger.info(f"    ✓ Added {len(new_rows)} CDS metric records for CNV copies")
    
    return mrna_metrics_df, cds_metrics_df


def generate_consensus(args):
    """
    Main consensus finding logic with chromosome-based processing.
    """
    start_time = time.time()
    
    logger.info("\n" + "="*80)
    logger.info("STEP 1: Loading Input Data")
    logger.info("="*80)
    logger.info(f"Input genePred files: {len(args.gp_list)}")
    for gp_file in args.gp_list:
        logger.info(f"  - {gp_file}")
    
    logger.info("\nBuilding alignment source map...")
    
    # Create a mapping from alignment IDs to their source files
    alignment_source_map = {}
    tx_by_mode = collections.defaultdict(dict)
    
    # Load genePreds and track their sources
    for gp_idx, gp_file in enumerate(args.gp_list, 1):
        logger.info(f"  Loading {gp_idx}/{len(args.gp_list)}: {gp_file}...")
        tx_count = 0
        
        # Determine mode from filename (check more specific patterns first)
        if 'transMap_pairwise' in gp_file:
            mode = 'transMap_pairwise'
        elif 'transMap' in gp_file:
            mode = 'transMap'
        elif 'txTM' in gp_file:
            mode = 'txTM'
        elif 'augTMR_pairwise' in gp_file:
            mode = 'augTMR_pairwise'
        elif 'augTM_pairwise' in gp_file:
            mode = 'augTM_pairwise'
        elif 'augTMR' in gp_file:
            mode = 'augTMR'
        elif 'augTM' in gp_file:
            mode = 'augTM'
        elif 'augPB' in gp_file:
            mode = 'augPB'
        elif '_strg' in gp_file:
            mode = 'strg'
        else:
            mode = 'unknown'
        
        # Load transcripts - handle duplicates within the same mode (e.g., paralogs, CNV copies)
        for t in tools.transcripts.gene_pred_iterator(gp_file):
            tx_count += 1
            # Use pre-determined mode or infer from ID if unknown
            if mode == 'unknown':
                inferred_mode = identify_mode(t.name, gp_file)
                actual_mode = inferred_mode
            else:
                actual_mode = mode
            
            # Normalize txTM IDs BEFORE duplicate detection to catch CNV copies at different loci
            # TxTM uses _N suffix (e.g., gene_1, gene_2) for CNV copies that should be preserved
            if actual_mode == 'txTM':
                normalized_name = normalize_alignment_id(t.name, actual_mode)
                if normalized_name != t.name:
                    logger.debug(f"    Normalized txTM ID: {t.name} → {normalized_name}")
                    t.name = normalized_name
                # Also normalize the gene ID (name2 field) to strip txTM's _N suffix
                if t.name2:
                    normalized_gene = normalize_gene_id(t.name2, actual_mode)
                    if normalized_gene != t.name2:
                        logger.debug(f"    Normalized txTM gene ID: {t.name2} → {normalized_gene}")
                        t.name2 = normalized_gene
            
            # Handle duplicates within the same mode (isoforms vs paralogs at different loci)
            original_name = t.name
            original_gene_name = t.name2  # Save original gene ID
            if original_name in tx_by_mode[actual_mode]:
                existing = tx_by_mode[actual_mode][original_name]
                # Check if transcripts are identical (same hash AND same location)
                same_location = (existing.chromosome == t.chromosome and 
                               existing.start == t.start and 
                               existing.stop == t.stop)
                if hash(existing) == hash(t) and same_location:
                    # Truly identical transcript at same location - skip
                    continue
                
                # Check if transcripts overlap (isoforms) or are at different loci (paralogs)
                overlaps = (existing.chromosome == t.chromosome and
                           t.start <= existing.stop and t.stop >= existing.start)
                
                # Different transcript - add with unique suffix to transcript ID
                suffix_num = 2
                new_name = f"{original_name}_cp{suffix_num}"
                while new_name in tx_by_mode[actual_mode]:
                    suffix_num += 1
                    new_name = f"{original_name}_cp{suffix_num}"
                t.name = new_name
            
            # CROSS-MODE duplicate handling: If this transcript ID exists in alignment_source_map 
            # (meaning another mode already has it), add a suffix to keep both versions
            # This ensures ALL transcripts from ALL modes are included in consensus generation
            if t.name in alignment_source_map and alignment_source_map[t.name] != actual_mode:
                # Transcript ID exists in another mode - add suffix to distinguish
                original_name = t.name
                suffix_num = 2
                new_name = f"{original_name}_cp{suffix_num}"
                # Check across ALL modes to find an unused name
                while any(new_name in mode_txs for mode_txs in tx_by_mode.values()) or new_name in alignment_source_map:
                    suffix_num += 1
                    new_name = f"{original_name}_cp{suffix_num}"
                logger.debug(f"    Cross-mode duplicate: {t.name} exists in {alignment_source_map[t.name]}, renaming {actual_mode} version to {new_name}")
                t.name = new_name
            
            # Store the transcript
            tx_by_mode[actual_mode][t.name] = t
            
            # Track in alignment_source_map (now each transcript has a unique ID)
            alignment_source_map[t.name] = actual_mode
        
        logger.info(f"    ✓ Loaded {tx_count} transcripts (mode: {mode})")
    
    # Report cross-mode duplicate handling
    cross_mode_duplicates = sum(1 for tx_id in alignment_source_map.keys() if re.search(r'_cp\d+$', tx_id))
    if cross_mode_duplicates > 0:
        logger.info(f"\n  ℹ️  Found {cross_mode_duplicates} cross-mode duplicate transcript IDs")
        logger.info(f"     (renamed with _cp suffix to include all versions in consensus)")
    
    # Set the global mapping for use by alignment_type function
    tools.nameConversions.set_alignment_source_map(alignment_source_map)
    
    # Load all genePreds
    logger.info("\nLoading all genePred transcripts into memory...")
    tx_dict = tools.transcripts.load_gps(args.gp_list)
    
    logger.info(f"\n✓ Loaded {len(tx_dict)} total transcripts")
    
    # Load .gp_attrs files for biotype information (especially for txTM)
    logger.info("\nLoading genePred attribute files (.gp_attrs) for biotype information...")
    gp_attrs_biotypes = {}  # transcript_id -> biotype
    for gp_file in args.gp_list:
        attrs_file = gp_file + '_attrs'
        if os.path.exists(attrs_file):
            logger.info(f"  Loading {attrs_file}...")
            with open(attrs_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        transcript_id, attr_name, attr_value = parts[0], parts[1], parts[2]
                        if attr_name == 'gene_biotype':
                            gp_attrs_biotypes[transcript_id] = attr_value
    
    logger.info(f"✓ Loaded biotype info for {len(gp_attrs_biotypes)} transcripts from .gp_attrs files")
    logger.info(f"  Modes found: {list(tx_by_mode.keys())}")
    for mode, txs in tx_by_mode.items():
        logger.info(f"    {mode}: {len(txs)} transcripts")
    
    # Update alignment_source_map for transcripts that were renamed during loading
    # (e.g., cross-mode conflicts where load_gps added "_cp", "_cp2", "_cp3" suffix)
    # Strategy: load_gps returns transcripts with _cp suffix already applied,
    # so tx_by_mode contains the _cp name directly for the mode that loaded it
    renamed_count = 0
    debug_234_238 = []
    for tx_id in tx_dict.keys():
        # Check for _cp, _cp2, _cp3, etc. suffix
        if re.search(r'_cp\d*$', tx_id) and tx_id not in alignment_source_map:
            # Check which mode has this _cp transcript
            found_mode = None
            for mode, mode_txs in tx_by_mode.items():
                if tx_id in mode_txs:
                    found_mode = mode
                    break
            
            if found_mode:
                alignment_source_map[tx_id] = found_mode
                renamed_count += 1
                # Debug genes 234/238
                if '234' in tx_id or '238' in tx_id:
                    tx = tx_dict[tx_id]
                    debug_234_238.append(f"{tx_id} → {found_mode} at {tx.start}-{tx.stop}")
    if renamed_count > 0:
        logger.info(f"  Updated alignment_source_map for {renamed_count} renamed transcripts (cross-mode conflicts)")
    
    # Track unique genes from each source
    logger.info("\n  Gene counts by source:")
    genes_by_source = {}
    for mode, txs in tx_by_mode.items():
        mode_genes = set()
        for tx_id, tx_obj in txs.items():
            # Extract gene from transcript ID
            # For transMap/txTM, the transcript ID itself can be used to get gene from ref_df later
            # For now, just track transcript IDs as proxy
            mode_genes.add(tx_id)
        genes_by_source[mode] = mode_genes
        logger.info(f"    {mode}: {len(mode_genes)} transcript IDs")
    
    # Load reference annotation information
    logger.info("\nLoading reference annotation...")
    ref_df = tools.sqlInterface.load_annotation(args.ref_db_path)
    logger.info(f"✓ Loaded {len(ref_df)} reference annotations")
    ref_biotype_counts = collections.Counter(ref_df.TranscriptBiotype)
    coding_count = ref_biotype_counts['protein_coding']
    non_coding_count = sum(y for x, y in ref_biotype_counts.items() if x != 'protein_coding')
    
    # Create readthrough gene map from ExtraTags (check both readthrough_gene and readthrough_transcript)
    readthrough_gene_set = set()
    if 'ExtraTags' in ref_df.columns:
        # Some annotations have readthrough_gene tag, some have readthrough_transcript tag
        readthrough_mask = (ref_df['ExtraTags'].str.contains('readthrough_gene', na=False) | 
                           ref_df['ExtraTags'].str.contains('readthrough_transcript', na=False))
        readthrough_gene_set = set(ref_df[readthrough_mask]['GeneId'].unique())
        logger.info(f"  Found {len(readthrough_gene_set)} readthrough genes in reference")
    
    # Build reference gene coordinate map for checking overlaps
    # This will be used to check if overlapping genes in the target also overlap in the reference
    logger.info("  Building reference gene coordinate map...")
    ref_gene_coords = {}
    genes_with_overlaps_in_ref = set()  # Track which genes overlap with others in reference
    
    if 'GeneId' in ref_df.columns and 'ExtraTags' in ref_df.columns:
        # Extract coordinates from ExtraTags using vectorized string operations
        # ExtraTags format: "...;Seqid=chr1;...;Start=12345;...;End=67890;..."
        ref_df_with_coords = ref_df[ref_df['ExtraTags'].notna()].copy()
        
        # Extract Seqid, Start, End using regex
        ref_df_with_coords['Seqid'] = ref_df_with_coords['ExtraTags'].str.extract(r'Seqid=([^;]+)', expand=False)
        ref_df_with_coords['Start'] = ref_df_with_coords['ExtraTags'].str.extract(r'Start=(\d+)', expand=False)
        ref_df_with_coords['End'] = ref_df_with_coords['ExtraTags'].str.extract(r'End=(\d+)', expand=False)
        
        # Filter rows where all three fields are present
        ref_df_with_coords = ref_df_with_coords[
            ref_df_with_coords['Seqid'].notna() & 
            ref_df_with_coords['Start'].notna() & 
            ref_df_with_coords['End'].notna()
        ]
        
        if len(ref_df_with_coords) > 0:
            # Convert Start and End to integers
            ref_df_with_coords['Start'] = ref_df_with_coords['Start'].astype(int)
            ref_df_with_coords['End'] = ref_df_with_coords['End'].astype(int)
            
            # Group by GeneId and get min Start, max End
            gene_coords_grouped = ref_df_with_coords.groupby('GeneId').agg({
                'Seqid': 'first',  # Chromosome should be same for all transcripts of a gene
                'Start': 'min',
                'End': 'max'
            })
            
            # Convert to dictionary (vectorized - avoid iterrows)
            ref_gene_coords = {
                gene_id: (row.Seqid, row.Start, row.End)
                for gene_id, row in zip(gene_coords_grouped.index, gene_coords_grouped.itertuples(index=False))
            }
            
            # Find which genes overlap with others in the reference
            # This helps identify gene families that naturally overlap
            logger.info("  Identifying genes that overlap in reference...")
            ref_genes_by_chrom = collections.defaultdict(list)
            for gene_id, (chrom, start, end) in ref_gene_coords.items():
                ref_genes_by_chrom[chrom].append((start, end, gene_id))
            
            for chrom, gene_list in ref_genes_by_chrom.items():
                gene_list.sort()  # Sort by start position
                for i in range(len(gene_list)):
                    start_i, end_i, gene_i = gene_list[i]
                    for j in range(i + 1, len(gene_list)):
                        start_j, end_j, gene_j = gene_list[j]
                        if start_j >= end_i:  # No more overlaps possible
                            break
                        # Check if they overlap
                        if not (end_i <= start_j or end_j <= start_i):
                            genes_with_overlaps_in_ref.add(gene_i)
                            genes_with_overlaps_in_ref.add(gene_j)
            
            logger.info(f"  Found {len(genes_with_overlaps_in_ref)} genes that overlap with others in reference")
            
        logger.info(f"  Built reference gene coordinates for {len(ref_gene_coords)} genes")
    
    gene_biotype_map = tools.sqlInterface.get_gene_biotype_map(args.ref_db_path)
    transcript_biotype_map = tools.sqlInterface.get_transcript_biotype_map(args.ref_db_path)
    logger.info(f"  Reference has {len(gene_biotype_map)} genes, {len(transcript_biotype_map)} transcripts")
    logger.info(f"  Biotype breakdown: {dict(ref_biotype_counts)}")
    
    # Map transcript IDs to gene IDs for tracking (vectorized for speed)
    ref_df['TranscriptId_base'] = ref_df['TranscriptId'].astype(str).str.replace(r"\.[0-9]+$", "", regex=True)
    tx_to_gene_map = dict(zip(ref_df['TranscriptId_base'], ref_df['GeneId']))
    
    # Now track actual gene IDs from each source
    logger.info("\n  Mapping to actual gene IDs:")
    gene_ids_by_source = {}
    pc_gene_ids_by_source = {}
    for mode, tx_ids in genes_by_source.items():
        mode_gene_ids = set()
        mode_pc_gene_ids = set()
        for tx_id in tx_ids:
            # Normalize the tx_id
            normalized_id = normalize_alignment_id(tx_id, mode)
            # Strip version to match tx_to_gene_map (which is versionless)
            normalized_id_base = re.sub(r'\.[0-9]+$', '', normalized_id)
            gene_id = tx_to_gene_map.get(normalized_id_base)
            if gene_id:
                mode_gene_ids.add(gene_id)
                # Check if protein coding
                if gene_biotype_map.get(gene_id) == 'protein_coding':
                    mode_pc_gene_ids.add(gene_id)
        gene_ids_by_source[mode] = mode_gene_ids
        pc_gene_ids_by_source[mode] = mode_pc_gene_ids
        logger.info(f"    {mode}: {len(mode_gene_ids)} genes ({len(mode_pc_gene_ids)} protein-coding)")
    
    # Load transMap evaluation data
    logger.info("\nLoading transMap evaluation data...")
    tm_eval_df = load_transmap_evals(args.db_path)
    logger.info(f"✓ Loaded {len(tm_eval_df)} transMap evaluations")
    
    # Determine which modes are available
    tx_modes = list(tx_by_mode.keys())
    tx_modes_with_metrics = [x for x in tx_modes if x not in ('augPB', 'strg')]
    
    logger.info("\n" + "="*80)
    logger.info("STEP 2: Loading Alignment Metrics")
    logger.info("="*80)
    # Load metrics (filter to actual transcripts for speed)
    valid_aln_ids = set(tx_dict.keys())
    # Build extended set including transMap base IDs for multi-mapped transcripts.
    # e.g. NM_004043.3_1 in tx_dict but metrics DB has NM_004043.3 — include both.
    transmap_base_ids = set()
    for aln_id in valid_aln_ids:
        if alignment_source_map.get(aln_id, '') in ('transMap', 'transMap_pairwise'):
            m = re.match(r'^(.+)_(\d+)$', aln_id)
            if m:
                transmap_base_ids.add(m.group(1))
    valid_aln_ids_extended = valid_aln_ids | transmap_base_ids
    if transmap_base_ids:
        logger.info(f"  Extended metrics lookup with {len(transmap_base_ids)} transMap base IDs for multi-mapped transcripts")

    if tx_modes_with_metrics:
        logger.info(f"Loading metrics for modes: {tx_modes_with_metrics}")

        mrna_dfs = []
        for tx_mode in tx_modes_with_metrics:
            df = load_metrics_from_db(args.db_path, tx_mode, 'mRNA')
            # For transMap, include base IDs so multi-mapped variants can be backfilled
            filter_ids = valid_aln_ids_extended if tx_mode in ('transMap', 'transMap_pairwise') else valid_aln_ids
            df = df[df['AlignmentId'].isin(filter_ids)]
            mrna_dfs.append(df)
        mrna_metrics_df = pd.concat(mrna_dfs) if mrna_dfs else pd.DataFrame()
        logger.info(f"✓ Loaded {len(mrna_metrics_df)} mRNA metrics (filtered to actual transcripts)")

        cds_dfs = []
        for tx_mode in tx_modes_with_metrics:
            df = load_metrics_from_db(args.db_path, tx_mode, 'CDS')
            filter_ids = valid_aln_ids_extended if tx_mode in ('transMap', 'transMap_pairwise') else valid_aln_ids
            df = df[df['AlignmentId'].isin(filter_ids)]
            cds_dfs.append(df)
        cds_metrics_df = pd.concat(cds_dfs) if cds_dfs else pd.DataFrame()
        logger.info(f"✓ Loaded {len(cds_metrics_df)} CDS metrics (filtered to actual transcripts)")

        # Backfill metrics for txTM CNV copies (_N suffix) that don't have metrics
        mrna_metrics_df, cds_metrics_df = backfill_cnv_metrics(
            mrna_metrics_df, cds_metrics_df, valid_aln_ids, alignment_source_map, args.db_path
        )
        # Backfill metrics for transMap multi-mapped transcripts (_N suffix variants)
        mrna_metrics_df, cds_metrics_df = backfill_transmap_multimapped_metrics(
            mrna_metrics_df, cds_metrics_df, valid_aln_ids, alignment_source_map
        )
        
        eval_dfs = []
        for tx_mode in tx_modes_with_metrics:
            df = load_evaluations_from_db(args.db_path, tx_mode)
            df = df[df.index.isin(valid_aln_ids)]
            eval_dfs.append(df)
        eval_df = pd.concat(eval_dfs).reset_index() if eval_dfs else pd.DataFrame()
        logger.info(f"✓ Loaded {len(eval_df)} evaluation entries (filtered to actual transcripts)")
    else:
        logger.warning("No modes with metrics found")
        mrna_metrics_df = pd.DataFrame()
        cds_metrics_df = pd.DataFrame()
        eval_df = pd.DataFrame()
    
    # Create support dataframe
    logger.info("\nCreating support dataframe...")
    logger.info(f"  args.denovo_tx_modes = {args.denovo_tx_modes} (type: {type(args.denovo_tx_modes)})")
    support_df = create_support_dataframe(tx_dict, args.db_path, ref_df, alignment_source_map, denovo_tx_modes=args.denovo_tx_modes)
    
    # Process by chromosome
    logger.info("\n" + "="*80)
    logger.info("STEP 3: Grouping Transcripts by Chromosome")
    logger.info("="*80)
    tx_by_chrom = collections.defaultdict(list)
    for aln_id, tx_obj in tx_dict.items():
        tx_by_chrom[tx_obj.chromosome].append(aln_id)
    
    logger.info(f"✓ Found {len(tx_by_chrom)} chromosomes")
    for chrom in sorted(tx_by_chrom.keys()):
        logger.info(f"  {chrom}: {len(tx_by_chrom[chrom])} transcripts")
    
    logger.info("\n" + "="*80)
    logger.info("STEP 4: Processing Chromosomes")
    logger.info("="*80)

    # Multiprocessing mode (only mode supported)
    max_workers = args.num_workers if args.num_workers else min(mp.cpu_count(), len(tx_by_chrom), 23)
    logger.info(f"Using multiprocessing with {max_workers} parallel workers for {len(tx_by_chrom)} chromosomes")

    all_consensus_transcripts = []
    metrics = initialize_metrics()

    # Prepare chromosome tasks
    chrom_tasks = []
    for chrom_num, chrom in enumerate(sorted(tx_by_chrom.keys()), 1):
        chrom_tx_ids = tx_by_chrom[chrom]
        chrom_tx_set = set(chrom_tx_ids)

        # Filter dataframes to this chromosome
        chrom_support_df = support_df[support_df['AlignmentId'].isin(chrom_tx_set)].copy()
        chrom_mrna_df = mrna_metrics_df[mrna_metrics_df['AlignmentId'].isin(chrom_tx_set)].copy() if len(mrna_metrics_df) > 0 else pd.DataFrame()
        chrom_cds_df = cds_metrics_df[cds_metrics_df['AlignmentId'].isin(chrom_tx_set)].copy() if len(cds_metrics_df) > 0 else pd.DataFrame()
        chrom_eval_df = eval_df[eval_df['AlignmentId'].isin(chrom_tx_set)].copy() if len(eval_df) > 0 else pd.DataFrame()

        chrom_tasks.append((
            chrom, chrom_num, len(tx_by_chrom), chrom_tx_ids, tx_dict, chrom_support_df,
            chrom_mrna_df, chrom_cds_df, chrom_eval_df, tm_eval_df,
            ref_df, alignment_source_map, args, readthrough_gene_set, ref_gene_coords,
            genes_with_overlaps_in_ref, gp_attrs_biotypes, gene_biotype_map
        ))

    # Process chromosomes in parallel
    if max_workers > 1:
        logger.info("Processing chromosomes in parallel...")
        with mp.Pool(processes=max_workers) as pool:
            results = pool.starmap(process_chromosome_wrapper, chrom_tasks)
    else:
        logger.info("Processing chromosomes sequentially...")
        results = [process_chromosome_wrapper(*task) for task in chrom_tasks]

    # Aggregate results
    for chrom_consensus, chrom_metrics in results:
        all_consensus_transcripts.extend(chrom_consensus)
        merge_metrics(metrics, chrom_metrics)

    logger.info(f"\n✓ Total consensus transcripts across all chromosomes: {len(all_consensus_transcripts)}")
    elapsed = time.time() - start_time
    logger.info(f"  Chromosome processing took {elapsed:.1f} seconds")
    
    # Convert to consensus dict format
    logger.info("\n" + "="*80)
    logger.info("STEP 5: Final Filtering and Cleanup")
    logger.info("="*80)
    logger.info("Converting to consensus dictionary format...")
    consensus_dict = {}
    for aln_id, attrs in all_consensus_transcripts:
        consensus_dict[aln_id] = attrs
    logger.info(f"✓ Consensus dict has {len(consensus_dict)} entries")
    
    # Deduplication is now done per-chromosome (see process_chromosome Step 9)
    # This prevents genes on different chromosomes from being treated as duplicates
    logger.info("\nSkipping global deduplication (already done per-chromosome)...")
    deduplicated_consensus = consensus_dict
    if metrics['Duplicate transcripts']:
        logger.info(f"  Total duplicates removed across all chromosomes: {dict(metrics['Duplicate transcripts'])}")
    
    args.ref_pc_ensg = {norm_ensg(g) for g, b in gene_biotype_map.items() if b == 'protein_coding'}
    final_consensus = finalize_consensus_after_source_gene_resolution(
        consensus_dict,
        tx_dict,
        metrics,
        args,
        readthrough_gene_set,
        ref_gene_coords,
        genes_with_overlaps_in_ref,
        run_resolve_overlapping_different_genes=True,
    )

    final_consensus = apply_reference_gene_biotype_policy(
        final_consensus, gene_biotype_map, metrics=metrics
    )

    final_consensus = rescue_missing_reference_pc_genes(
        final_consensus,
        tx_dict,
        alignment_source_map,
        ref_gene_coords,
        gene_biotype_map,
        ref_df,
        args,
        mrna_metrics_df=mrna_metrics_df,
        metrics=metrics,
    )
    final_consensus = rescue_missing_reference_noncoding_genes(
        final_consensus,
        tx_dict,
        alignment_source_map,
        ref_gene_coords,
        gene_biotype_map,
        ref_df,
        args,
        mrna_metrics_df=mrna_metrics_df,
        metrics=metrics,
    )
    final_consensus = rescue_missing_reference_transcripts(
        final_consensus,
        tx_dict,
        alignment_source_map,
        ref_df,
        gene_biotype_map,
        args,
        mrna_metrics_df=mrna_metrics_df,
        metrics=metrics,
    )
    final_consensus = rescue_alternative_source_isoforms(
        final_consensus,
        tx_dict,
        alignment_source_map,
        ref_df,
        args,
        mrna_metrics_df=mrna_metrics_df,
        metrics=metrics,
    )

    # Count genes and protein coding genes
    unique_genes = set()
    unique_pc_genes = set()
    unique_pc_genes_by_source = collections.defaultdict(set)
    
    for aln_id, attrs in final_consensus:
        gene = attrs.get('source_gene')
        mode = alignment_source_map.get(aln_id, 'unknown')
        if gene and gene != 'N/A':
            unique_genes.add(gene)
            if attrs.get('gene_biotype') == 'protein_coding':
                unique_pc_genes.add(gene)
                unique_pc_genes_by_source[mode].add(gene)
    
    logger.info(f"  Total unique genes: {len(unique_genes)}")
    logger.info(f"  Protein-coding genes: {len(unique_pc_genes)}")
    logger.info(f"  Protein-coding genes by source: {dict((k, len(v)) for k, v in unique_pc_genes_by_source.items())}")
    
    # Compare input to output for each source
    logger.info("\n  Protein-coding gene comparison (input vs output):")
    for mode in pc_gene_ids_by_source.keys():
        input_genes = pc_gene_ids_by_source[mode]
        output_genes = unique_pc_genes_by_source.get(mode, set())
        missing_genes = input_genes - output_genes
        logger.info(f"    {mode}: {len(input_genes)} input → {len(output_genes)} output")
        if missing_genes:
            logger.info(f"      Missing {len(missing_genes)} genes: {sorted(list(missing_genes))[:10]}{'...' if len(missing_genes) > 10 else ''}")
    
    # Calculate final metrics
    logger.info("\n" + "="*80)
    logger.info("STEP 6: Calculating Metrics")
    logger.info("="*80)
    calculate_completeness(final_consensus, metrics, gene_biotype_map, transcript_biotype_map)
    
    # Log summary metrics
    if 'Multi-locus mappings' in metrics and metrics['Multi-locus mappings']:
        logger.info("\nMulti-locus mapping summary:")
        for mode, count in metrics['Multi-locus mappings'].items():
            kept = metrics['Multi-locus kept'].get(mode, 0)
            logger.info(f"  {mode}: {count} transcripts map to multiple loci, kept {kept} total copies")
    
    if 'AugPB Classes' in metrics and metrics['AugPB Classes']:
        logger.info("\nAugPB classification summary:")
        for cls, count in metrics['AugPB Classes'].items():
            logger.info(f"  {cls}: {count}")
    
    # Write outputs
    logger.info("\n" + "="*80)
    logger.info("STEP 7: Writing Output Files")
    logger.info("="*80)
    
    logger.info(f"Writing consensus genePred to {args.consensus_gp}")
    consensus_gene_dict = write_consensus_gps(args.consensus_gp, args.consensus_gp_info,
                                              final_consensus, tx_dict, args.genome)
    logger.info(f"✓ Wrote {len(final_consensus)} transcripts to genePred")
    
    logger.info(f"\nWriting GFF3 to {args.consensus_gff3}")
    write_consensus_gff3(consensus_gene_dict, args.consensus_gff3)
    logger.info("✓ Wrote GFF3 file")
    
    logger.info(f"\nWriting FASTA files...")
    logger.info(f"  Transcript FASTA: {args.consensus_fasta}")
    logger.info(f"  Protein FASTA: {args.protein_fasta}")
    write_consensus_fastas(consensus_gene_dict, args.consensus_fasta, args.protein_fasta, args.fasta)
    logger.info("✓ Wrote FASTA files")
    
    total_time = time.time() - start_time
    logger.info(f"\n✓ Total consensus generation time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")
    
    return metrics


def load_alt_gene_assignments(db_path, denovo_tx_modes):
    """Load alternative gene assignments (AssignedGeneId) for augPB transcripts"""
    session = tools.sqlInterface.start_session(db_path)
    r = []
    for tx_mode in denovo_tx_modes:
        try:
            table = tools.sqlInterface.tables['alt_names'][tx_mode]
            
            # Check if the table exists in the database
            inspector = inspect(session.bind)
            if inspector.has_table(table.__tablename__):
                alt_data = tools.sqlInterface.load_alternatives(table, session)
                if not alt_data.empty:
                    r.append(alt_data)
        except Exception as e:
            logger.warning(f"Error loading alternative names for {tx_mode}: {e}")
    
    session.close()
    
    if not r:
        # If no valid data, create an empty dataframe with required columns
        empty_df = pd.DataFrame(columns=['TranscriptId', 'AssignedGeneId', 'AlternativeGeneIds', 'ResolutionMethod'])
        # rename TranscriptId to AlignmentId
        empty_df.columns = [x if x != 'TranscriptId' else 'AlignmentId' for x in empty_df.columns]
        return empty_df
    
    df = pd.concat(r, ignore_index=True)
    df.columns = [x if x != 'TranscriptId' else 'AlignmentId' for x in df.columns]
    return df


def create_support_dataframe(tx_dict, db_path, ref_df, alignment_source_map, denovo_tx_modes=None):
    """Create support dataframe with RNA-seq and annotation support"""
    logger.info("  Loading alignment evaluation data...")
    
    # Load evaluation data
    tm_eval = tools.sqlInterface.load_alignment_evaluation(db_path)
    tm_filter_eval = tools.sqlInterface.load_filter_evaluation(db_path)
    
    # Build a mapping of original IDs to renamed IDs for cross-mode conflicts
    # (e.g., when load_gps renamed transcripts by adding "_cp")
    id_rename_map = {}
    for tx_id in tx_dict.keys():
        if tx_id.endswith('_cp'):
            original_id = tx_id[:-3]
            # Check if the original ID exists in tx_dict and is from a different source
            if original_id in tx_dict:
                # We have a conflict: original_id exists in tx_dict (likely from a different mode)
                # and tx_id is the renamed version. The database likely has data under original_id
                # that should be associated with the renamed version.
                # We need to determine which mode each belongs to
                original_mode = alignment_source_map.get(original_id)
                renamed_mode = alignment_source_map.get(tx_id)
                
                # If they're from different modes, the database entry for original_id from renamed_mode
                # should be updated to use tx_id
                if original_mode != renamed_mode:
                    # The database has entries under original_id for renamed_mode's data
                    # We need to map: (original_id, renamed_mode) -> tx_id
                    id_rename_map[(original_id, renamed_mode)] = tx_id
    
    # Update AlignmentIds in database DataFrames to match tx_dict's naming
    if id_rename_map:
        logger.info(f"  Updating AlignmentIds for {len(id_rename_map)} cross-mode conflicts...")
        
        # For tm_eval and tm_filter_eval, we need to map based on the mode
        # But we don't have a mode column in these DataFrames, so we need to infer it
        # from the AlignmentId itself using alignment_source_map
        
        def update_alignment_ids(df):
            """Update AlignmentIds in a DataFrame based on the rename map"""
            if 'AlignmentId' not in df.columns:
                return df
            
            updated_count = 0
            new_ids = []
            for aln_id in df['AlignmentId']:
                # Determine the mode for this alignment
                # The database was populated before renaming, so we need to check if this ID
                # should be renamed based on which mode it belongs to
                # 
                # Problem: We can't easily determine the mode from the database alone.
                # Solution: Check if aln_id is NOT in tx_dict but aln_id+"_cp" IS in tx_dict
                if aln_id not in tx_dict and f"{aln_id}_cp" in tx_dict:
                    new_ids.append(f"{aln_id}_cp")
                    updated_count += 1
                else:
                    new_ids.append(aln_id)
            
            df['AlignmentId'] = new_ids
            if updated_count > 0:
                logger.info(f"    Updated {updated_count} AlignmentIds to match tx_dict naming")
            return df
        
        tm_eval = update_alignment_ids(tm_eval)
        tm_filter_eval = update_alignment_ids(tm_filter_eval)
    
    # Filter to only transcripts we have before merging (faster)
    valid_aln_ids = set(tx_dict.keys())
    tm_eval = tm_eval[tm_eval['AlignmentId'].isin(valid_aln_ids)]
    tm_filter_eval = tm_filter_eval[tm_filter_eval['AlignmentId'].isin(valid_aln_ids)]
    
    tm_eval_df = pd.merge(tm_eval, tm_filter_eval, on=['TranscriptId', 'AlignmentId'])
    logger.info(f"    Loaded {len(tm_eval_df)} evaluation records (filtered to actual transcripts)")
    
    # Create base support_df
    support_df = tm_eval_df[['GeneId', 'TranscriptId', 'AlignmentId']].copy()
    
    # CRITICAL: Update GeneId from tx_dict to get _cp suffixed gene IDs for paralogs
    # This must happen BEFORE normalization, as tx_dict has the corrected gene IDs
    def get_gene_id_from_tx_dict(aln_id):
        tx_obj = tx_dict.get(aln_id)
        return tx_obj.name2 if tx_obj and tx_obj.name2 else None
    
    support_df['GeneId_from_tx'] = support_df['AlignmentId'].apply(get_gene_id_from_tx_dict)
    # Update GeneId where we have a value from tx_dict (this includes _cp suffixed paralogs)
    has_tx_gene = support_df['GeneId_from_tx'].notna()
    if has_tx_gene.any():
        support_df.loc[has_tx_gene, 'GeneId'] = support_df.loc[has_tx_gene, 'GeneId_from_tx']
        logger.info(f"    Updated GeneId from tx_dict for {has_tx_gene.sum()} transcripts (includes _cp paralogs)")
    support_df.drop(columns=['GeneId_from_tx'], inplace=True)
    
    # Normalize GeneId for txTM transcripts to strip _N suffix (but NOT _cpN suffix!)
    # This ensures consistency with the normalized gene IDs in tx_dict
    # Note: _cp suffixes are preserved because they don't match the _\d+$ pattern
    txTM_mask = support_df['AlignmentId'].apply(lambda x: alignment_source_map.get(x, '') == 'txTM')
    if txTM_mask.any():
        support_df.loc[txTM_mask, 'GeneId'] = support_df.loc[txTM_mask, 'GeneId'].apply(
            lambda x: normalize_gene_id(x, 'txTM') if pd.notna(x) else x
        )
        logger.info(f"    Normalized GeneId for {txTM_mask.sum()} txTM transcripts in support_df")
    
    # Add missing transcripts (e.g., augPB, txTM without metrics)
    logger.info("  Adding missing transcripts to support dataframe...")
    existing_aln_ids = set(support_df['AlignmentId'].values)
    missing_aln_ids = set(tx_dict.keys()) - existing_aln_ids
    
    if missing_aln_ids:
        logger.info(f"    Found {len(missing_aln_ids)} transcripts not in evaluation data (likely augPB/txTM)")
        
        # Build a lookup dict for faster reference gene finding (avoid repeated DataFrame queries)
        ref_gene_lookup = dict(zip(ref_df['TranscriptId'], ref_df['GeneId']))
        
        missing_transcripts = []
        txTM_examples = []
        for tx_id in missing_aln_ids:
            mode = alignment_source_map.get(tx_id, 'unknown')
            
            # Get gene ID directly from the genePred object (preserves _N suffixes for txTM)
            tx_obj = tx_dict.get(tx_id)
            if tx_obj and tx_obj.name2:
                gene_id = tx_obj.name2
            else:
                gene_id = f'UNKNOWN_GENE_{tx_id}'
            
            # Normalize based on mode for TranscriptId lookup
            if mode == 'txTM':
                # Strip txTM's _N suffix first, then strip_alignment_numbers handles the rest
                normalized_id = normalize_alignment_id(tx_id, mode)
                base_tx_id = tools.nameConversions.strip_alignment_numbers(normalized_id)
                if len(txTM_examples) < 3:
                    txTM_examples.append(f"{tx_id} → gene_id={gene_id}, base_tx={base_tx_id}")
            else:
                base_tx_id = tools.nameConversions.strip_alignment_numbers(tx_id)
            
            missing_transcripts.append({
                'GeneId': gene_id,
                'TranscriptId': base_tx_id,
                'AlignmentId': tx_id
            })
        
        missing_df = pd.DataFrame(missing_transcripts)
        support_df = pd.concat([support_df, missing_df], ignore_index=True)
        logger.info(f"    ✓ Added {len(missing_transcripts)} missing transcripts")
    
    # Add support vectors (optimized - batch process)
    # Pre-compute dimensions for all transcripts
    logger.info("  Creating support vectors...")
    aln_ids = support_df['AlignmentId'].values
    
    # Batch create support vectors
    intron_supports = []
    exon_supports = []
    intron_annot_supports = []
    exon_annot_supports = []
    cds_annot_supports = []
    
    for aln_id in aln_ids:
        tx = tx_dict[aln_id]
        num_exons = len(tx.exon_frames)
        num_introns = num_exons - 1
        is_denovo = _is_denovo(aln_id)
        
        # For denovo modes: full RNA support (PacBio/StringTie-validated), for others: no support
        support_val = 1 if is_denovo else 0
        
        intron_supports.append([support_val] * num_introns)
        exon_supports.append([support_val] * num_exons)
        intron_annot_supports.append([0] * num_introns)
        exon_annot_supports.append([0] * num_exons)
        cds_annot_supports.append([0] * num_exons)
    
    support_df['AllSpeciesIntronRnaSupport'] = intron_supports
    support_df['AllSpeciesExonRnaSupport'] = exon_supports
    support_df['IntronRnaSupport'] = intron_supports
    support_df['ExonRnaSupport'] = exon_supports
    support_df['IntronAnnotSupport'] = intron_annot_supports
    support_df['CdsAnnotSupport'] = cds_annot_supports
    support_df['ExonAnnotSupport'] = exon_annot_supports
    
    # Calculate percentages (vectorized for speed)
    logger.info("  Calculating support percentages...")
    
    def calc_percent_batch(support_vectors):
        """Vectorized percentage calculation"""
        return [100.0 * sum(1 for x in vec if x > 0) / len(vec) if len(vec) > 0 else 0.0 
                for vec in support_vectors]
    
    support_df['IntronAnnotSupportPercent'] = calc_percent_batch(intron_annot_supports)
    support_df['ExonAnnotSupportPercent'] = calc_percent_batch(exon_annot_supports)
    support_df['CdsAnnotSupportPercent'] = calc_percent_batch(cds_annot_supports)
    support_df['ExonRnaSupportPercent'] = calc_percent_batch(exon_supports)
    support_df['IntronRnaSupportPercent'] = calc_percent_batch(intron_supports)
    support_df['AllSpeciesExonRnaSupportPercent'] = calc_percent_batch(exon_supports)
    support_df['AllSpeciesIntronRnaSupportPercent'] = calc_percent_batch(intron_supports)
    
    # Mark denovo transcripts (augPB or strg)
    support_df['IsAugPB'] = [_is_denovo(aln_id) for aln_id in aln_ids]
    
    # Load alternative gene assignments for augPB transcripts (AssignedGeneId, AlternativeGeneIds)
    if denovo_tx_modes:
        logger.info(f"  Loading alternative gene assignments for augPB transcripts (denovo_tx_modes={denovo_tx_modes})...")
        alt_assignments_df = load_alt_gene_assignments(db_path, denovo_tx_modes)
        if len(alt_assignments_df) > 0:
            logger.info(f"    Loaded {len(alt_assignments_df)} alternative gene assignments")
            # Check how many have non-empty AssignedGeneId
            non_empty = alt_assignments_df[alt_assignments_df['AssignedGeneId'].notna() & (alt_assignments_df['AssignedGeneId'] != '')].shape[0]
            logger.info(f"    {non_empty} have non-empty AssignedGeneId")
            # Merge with support_df
            support_df = pd.merge(support_df, alt_assignments_df, on='AlignmentId', how='left')
            logger.info(f"    ✓ Merged alternative gene assignments into support dataframe")
            # Verify the merge worked
            merged_with_assigned = support_df[(support_df['AssignedGeneId'].notna()) & (support_df['AssignedGeneId'] != '')].shape[0]
            logger.info(f"    After merge: {merged_with_assigned} rows have AssignedGeneId")
            
            # Look up source gene biotype from reference annotation
            logger.info(f"    Looking up source gene biotypes from reference annotation...")
            # Create a mapping of GeneId -> GeneBiotype from ref_df
            gene_biotype_map = dict(zip(ref_df['GeneId'], ref_df['GeneBiotype']))
            # Map AssignedGeneId to SourceGeneBiotype
            support_df['SourceGeneBiotype'] = support_df['AssignedGeneId'].map(gene_biotype_map)
            source_biotypes_found = support_df['SourceGeneBiotype'].notna().sum()
            logger.info(f"    Found biotypes for {source_biotypes_found} source genes")
        else:
            logger.info("    No alternative gene assignments found")
            # Add empty columns
            support_df['AssignedGeneId'] = None
            support_df['AlternativeGeneIds'] = None
            support_df['ResolutionMethod'] = None
            support_df['SourceGeneBiotype'] = None
    else:
        logger.info("  No denovo_tx_modes specified - skipping alternative gene assignments")
        # No denovo modes, add empty columns
        support_df['AssignedGeneId'] = None
        support_df['AlternativeGeneIds'] = None
        support_df['ResolutionMethod'] = None
        support_df['SourceGeneBiotype'] = None
    
    logger.info(f"Created support dataframe with {len(support_df)} entries")
    return support_df


def merge_metrics(target_metrics, source_metrics):
    """Merge metrics from one chromosome into the global metrics"""
    for key in source_metrics:
        if isinstance(source_metrics[key], collections.Counter):
            target_metrics[key].update(source_metrics[key])
        elif isinstance(source_metrics[key], collections.defaultdict):
            for subkey in source_metrics[key]:
                if isinstance(source_metrics[key][subkey], list):
                    target_metrics[key][subkey].extend(source_metrics[key][subkey])
                elif isinstance(source_metrics[key][subkey], (dict, collections.defaultdict)):
                    for subsubkey in source_metrics[key][subkey]:
                        if isinstance(source_metrics[key][subkey][subsubkey], list):
                            target_metrics[key][subkey][subsubkey].extend(source_metrics[key][subkey][subsubkey])
                        else:
                            target_metrics[key][subkey][subsubkey] += source_metrics[key][subkey][subsubkey]
                else:
                    target_metrics[key][subkey] += source_metrics[key][subkey]
        elif isinstance(source_metrics[key], dict) and not isinstance(source_metrics[key], collections.Counter):
            # Plain dict (e.g. 'denovo' with nested mode dicts)
            for subkey in source_metrics[key]:
                if subkey not in target_metrics.get(key, {}):
                    continue
                if isinstance(source_metrics[key][subkey], dict):
                    for subsubkey in source_metrics[key][subkey]:
                        target_metrics[key][subkey][subsubkey] += source_metrics[key][subkey][subsubkey]
                elif isinstance(source_metrics[key][subkey], (int, float)):
                    target_metrics[key][subkey] += source_metrics[key][subkey]
        elif isinstance(source_metrics[key], int):
            target_metrics[key] += source_metrics[key]
        elif isinstance(source_metrics[key], list):
            target_metrics[key].extend(source_metrics[key])


def process_chromosome_wrapper(chrom, chrom_num, total_chroms, chrom_tx_ids, tx_dict, support_df, 
                              mrna_metrics_df, cds_metrics_df, eval_df, tm_eval_df, ref_df, 
                              alignment_source_map, args, readthrough_gene_set, ref_gene_coords,
                              genes_with_overlaps_in_ref, gp_attrs_transcript_biotypes,
                              gp_attrs_gene_biotypes, ref_gene_biotype_map):
    """
    Wrapper for process_chromosome that handles logging and metrics initialization.
    Designed for parallel execution.
    """
    # Set up logging for this process
    logger = logging.getLogger(__name__)
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Processing chromosome {chrom_num}/{total_chroms}: {chrom}")
    logger.info(f"{'='*60}")
    logger.info(f"  Transcripts on {chrom}: {len(chrom_tx_ids)}")
    
    # Initialize metrics for this chromosome
    chrom_metrics = initialize_metrics()
    
    # Process the chromosome
    chrom_consensus = process_chromosome(
        chrom, chrom_tx_ids, tx_dict, support_df,
        mrna_metrics_df, cds_metrics_df, eval_df, tm_eval_df,
        ref_df, alignment_source_map, args, chrom_metrics, readthrough_gene_set,
        ref_gene_coords, genes_with_overlaps_in_ref, gp_attrs_transcript_biotypes,
        gp_attrs_gene_biotypes, ref_gene_biotype_map
    )
    
    logger.info(f"  ✓ Selected {len(chrom_consensus)} consensus transcripts for {chrom}")
    
    return chrom_consensus, chrom_metrics


def process_chromosome(chrom, chrom_tx_ids, tx_dict, support_df, mrna_metrics_df, cds_metrics_df,
                       eval_df, tm_eval_df, ref_df, alignment_source_map, args, metrics, readthrough_gene_set=None,
                       ref_gene_coords=None, genes_with_overlaps_in_ref=None, gp_attrs_transcript_biotypes=None,
                       gp_attrs_gene_biotypes=None, ref_gene_biotype_map=None):
    """
    Process a single chromosome to select consensus transcripts.
    
    Key logic:
    1. Group transcripts by normalized ID and mode
    2. Track multi-locus mappings per transcript-mode combination
    3. Keep CNV copies if alignment scores are similar (within threshold)
    4. Select best transcript per gene
    """
    chrom_start_time = time.time()
    
    # PRE-FILL biotypes from .gp_attrs BEFORE merging with reference
    # This ensures txTM/transMap/augTM transcripts without database entries still get biotypes
    if gp_attrs_transcript_biotypes or gp_attrs_gene_biotypes:
        logger.info(f"  Step 0: Pre-filling biotypes from .gp_attrs...")
        support_df['TranscriptBiotype'] = None
        support_df['GeneBiotype'] = None
        
        prefilled_count = 0

        for idx in support_df.index:
            aln_id = support_df.loc[idx, 'AlignmentId']

            # Set transcript biotype and gene biotype separately
            if aln_id in gp_attrs_transcript_biotypes:
                support_df.loc[idx, 'TranscriptBiotype'] = gp_attrs_transcript_biotypes[aln_id]
                prefilled_count += 1
            if aln_id in gp_attrs_gene_biotypes:
                support_df.loc[idx, 'GeneBiotype'] = gp_attrs_gene_biotypes[aln_id]
            # Also try without _cp suffix
            else:
                original_id = re.sub(r'_cp\d*$', '', aln_id)
                if original_id != aln_id:
                    if original_id in gp_attrs_transcript_biotypes:
                        support_df.loc[idx, 'TranscriptBiotype'] = gp_attrs_transcript_biotypes[original_id]
                        prefilled_count += 1
                    if original_id in gp_attrs_gene_biotypes:
                        support_df.loc[idx, 'GeneBiotype'] = gp_attrs_gene_biotypes[original_id]

        logger.info(f"    ✓ Pre-filled biotypes for {prefilled_count}/{len(support_df)} transcripts from .gp_attrs")

    # Merge with reference to get biotypes
    logger.info(f"  Step 1: Merging with reference annotation...")
    logger.info(f"    {len(support_df)} transcripts in support_df for this chromosome")
    
    # Add versionless TranscriptId column (avoid full copy)
    # Also strip _cp and _cpN suffixes for paralog copies so they match the reference
    # First strip _cp\d* (e.g., _cp, _cp2, _cp3), then strip .\d+ (version number)
    support_df.loc[:, 'TranscriptId_base'] = (support_df['TranscriptId'].astype(str)
                                               .str.replace(r"_cp\d*$", "", regex=True)  # Strip _cp, _cp2, _cp3, etc.
                                               .str.replace(r"\.[0-9]+$", "", regex=True))  # Strip .1, .2, etc.
    
    # Only copy ref_df columns we need (reduces memory)
    ref_cols_needed = ['TranscriptId', 'GeneId', 'TranscriptBiotype', 'GeneBiotype']
    # Add optional columns if they exist
    if 'GeneName' in ref_df.columns:
        ref_cols_needed.append('GeneName')
    if 'TranscriptName' in ref_df.columns:
        ref_cols_needed.append('TranscriptName')
    ref_df_subset = ref_df[ref_cols_needed].copy()
    ref_df_subset['TranscriptId_base'] = ref_df_subset['TranscriptId'].astype(str).str.replace(r"\.[0-9]+$", "", regex=True)
    
    support_ref_df = pd.merge(support_df, ref_df_subset, on='TranscriptId_base', how='left', suffixes=('', '_ref'))
    
    # Clean up columns - use pre-filled biotypes where available, otherwise use reference
    if 'GeneId_ref' in support_ref_df.columns:
        support_ref_df['GeneId'] = support_ref_df['GeneId'].fillna(support_ref_df['GeneId_ref'])
        support_ref_df = support_ref_df.drop('GeneId_ref', axis=1)
    if 'GeneName_ref' in support_ref_df.columns:
        support_ref_df['GeneName'] = support_ref_df.get('GeneName', pd.Series([None]*len(support_ref_df))).fillna(support_ref_df['GeneName_ref'])
        support_ref_df = support_ref_df.drop('GeneName_ref', axis=1)
    
    # For biotypes: keep pre-filled values (from .gp_attrs), fill missing with reference values
    if 'TranscriptBiotype_ref' in support_ref_df.columns:
        # Pre-filled biotypes take precedence
        support_ref_df['TranscriptBiotype'] = support_ref_df['TranscriptBiotype'].fillna(support_ref_df['TranscriptBiotype_ref'])
        support_ref_df = support_ref_df.drop('TranscriptBiotype_ref', axis=1)
    if 'GeneBiotype_ref' in support_ref_df.columns:
        support_ref_df['GeneBiotype'] = support_ref_df['GeneBiotype'].fillna(support_ref_df['GeneBiotype_ref'])
        support_ref_df = support_ref_df.drop('GeneBiotype_ref', axis=1)

    if ref_gene_biotype_map:
        ref_biotype_applied = 0
        for idx in support_ref_df.index:
            ref_biotype = _lookup_ref_gene_biotype(support_ref_df.at[idx, 'GeneId'], ref_gene_biotype_map)
            if ref_biotype:
                support_ref_df.at[idx, 'GeneBiotype'] = ref_biotype
                ref_biotype_applied += 1
        logger.info(
            f"    Applied reference gene biotypes to {ref_biotype_applied} rows (authoritative over CDS/.gp_attrs)"
        )
    
    if 'TranscriptId_base' in support_ref_df.columns:
        support_ref_df = support_ref_df.drop('TranscriptId_base', axis=1)
    
    # Merge with tm_eval_df
    support_ref_tm_df = pd.merge(support_ref_df, tm_eval_df, on=['GeneId', 'TranscriptId'], how='left')
    
    # Count transcripts that matched reference
    logger.info(f"  Step 1b: Checking reference matching...")
    matched_ref = support_ref_tm_df[support_ref_tm_df['TranscriptBiotype'].notna()].copy()
    unmatched_ref = support_ref_tm_df[support_ref_tm_df['TranscriptBiotype'].isna()].copy()
    logger.info(f"    Matched reference: {len(matched_ref)} transcripts")
    logger.info(f"    Unmatched reference: {len(unmatched_ref)} transcripts")

    # Show unmatched breakdown by mode
    if len(unmatched_ref) > 0:
        unmatched_by_mode = {}
        for mode in unmatched_ref['AlignmentId'].apply(lambda x: alignment_source_map.get(x, 'unknown')).unique():
            mode_df = unmatched_ref[unmatched_ref['AlignmentId'].apply(lambda x: alignment_source_map.get(x, 'unknown') == mode)]
            unmatched_by_mode[mode] = len(mode_df)
        logger.info(f"    Unmatched by source: {unmatched_by_mode}")
        # Show some examples
        txTM_unmatched = unmatched_ref[unmatched_ref['AlignmentId'].apply(lambda x: alignment_source_map.get(x, 'unknown') == 'txTM')]
        if len(txTM_unmatched) > 0:
            examples = txTM_unmatched['AlignmentId'].head(5).tolist()
            logger.info(f"    Example unmatched txTM IDs: {examples}")
    
    # Fill in biotypes from .gp_attrs for unmatched transcripts
    if (gp_attrs_transcript_biotypes or gp_attrs_gene_biotypes) and len(unmatched_ref) > 0:
        logger.info(f"  Step 1c: Filling biotypes from .gp_attrs for unmatched transcripts...")
        filled_count = 0
        not_filled = []
        for idx in unmatched_ref.index:
            aln_id = support_ref_tm_df.loc[idx, 'AlignmentId']
            # Try direct lookup first - set transcript and gene biotypes separately
            tx_found = False
            gene_found = False
            if aln_id in gp_attrs_transcript_biotypes:
                support_ref_tm_df.loc[idx, 'TranscriptBiotype'] = gp_attrs_transcript_biotypes[aln_id]
                tx_found = True
            if aln_id in gp_attrs_gene_biotypes:
                support_ref_tm_df.loc[idx, 'GeneBiotype'] = gp_attrs_gene_biotypes[aln_id]
                gene_found = True
            if tx_found or gene_found:
                filled_count += 1
            # If not found and has _cp suffix (with or without number), try without suffix
            else:
                # Strip _cp, _cp2, _cp3, etc. to get original ID
                original_id = re.sub(r'_cp\d*$', '', aln_id)
                if original_id != aln_id:
                    if original_id in gp_attrs_transcript_biotypes:
                        support_ref_tm_df.loc[idx, 'TranscriptBiotype'] = gp_attrs_transcript_biotypes[original_id]
                        tx_found = True
                    if original_id in gp_attrs_gene_biotypes:
                        support_ref_tm_df.loc[idx, 'GeneBiotype'] = gp_attrs_gene_biotypes[original_id]
                        gene_found = True
                    if tx_found or gene_found:
                        filled_count += 1
                    else:
                        not_filled.append(aln_id)
                else:
                    not_filled.append(aln_id)
        logger.info(f"    Filled biotypes for {filled_count} transcripts from .gp_attrs")
        if not_filled:
            logger.info(f"    Failed to fill {len(not_filled)} transcripts: {not_filled[:20]}")

    # Split into coding and non-coding
    logger.info(f"  Step 2: Splitting by biotype...")
    
    # Check for missing biotypes - this should be rare if .gp_attrs files are complete
    # Note: denovo (augPB/strg) transcripts are expected to lack biotypes (assigned during classification)
    if 'TranscriptBiotype' in support_ref_tm_df.columns:
        biotype_missing = support_ref_tm_df['TranscriptBiotype'].isna()
        is_denovo_mask = support_ref_tm_df['AlignmentId'].astype(str).apply(_is_denovo)
        
        # Only count non-denovo transcripts with missing biotypes as errors
        non_denovo_missing = biotype_missing & ~is_denovo_mask
        
        if non_denovo_missing.any():
            missing_count = non_denovo_missing.sum()
            missing_pct = 100 * missing_count / len(support_ref_tm_df[~is_denovo_mask])
            
            # Get examples of missing transcripts (non-denovo only)
            missing_examples = support_ref_tm_df[non_denovo_missing][['AlignmentId', 'GeneId']].head(10)
            
            if missing_pct > 1.0:
                logger.error(f"    ERROR: {missing_count} ({missing_pct:.1f}%) non-denovo transcripts missing biotype!")
                logger.error(f"    This indicates incomplete .gp_attrs files or reference annotation.")
                logger.error(f"    Examples: {list(missing_examples['AlignmentId'])}")
                logger.error(f"    These transcripts will be EXCLUDED from protein-coding gene set.")
            else:
                logger.warning(f"    Warning: {missing_count} ({missing_pct:.2f}%) non-denovo transcripts missing biotype")
                logger.warning(f"    Examples: {list(missing_examples['AlignmentId'])}")
                logger.warning(f"    These transcripts will be EXCLUDED from protein-coding gene set.")
        
        # Report denovo count separately (expected to have no biotype initially)
        denovo_count = is_denovo_mask.sum()
        if denovo_count > 0:
            logger.info(f"    Found {denovo_count} denovo transcripts (biotypes assigned during classification)")
    
    # Match consensus.combine_and_filter_dfs: txTM transcripts with CDS but missing or
    # non-coding biotype in the support table still belong in the protein_coding pool.
    if alignment_source_map:
        txTM_alns = {aln for aln, src in alignment_source_map.items() if src == 'txTM'}
        txTM_mask = support_ref_tm_df['AlignmentId'].isin(txTM_alns)
        fixed_biotype = 0
        for idx in support_ref_tm_df[txTM_mask].index:
            aln_id = support_ref_tm_df.at[idx, 'AlignmentId']
            if aln_id in tx_dict:
                tx_obj = tx_dict[aln_id]
                if hasattr(tx_obj, 'cds_size') and tx_obj.cds_size > 0:
                    gid = support_ref_tm_df.at[idx, 'GeneId']
                    if ref_gene_biotype_map and not _is_ref_protein_coding_gene(gid, ref_gene_biotype_map):
                        continue
                    current = support_ref_tm_df.at[idx, 'TranscriptBiotype']
                    if pd.isna(current) or current != 'protein_coding':
                        support_ref_tm_df.at[idx, 'TranscriptBiotype'] = 'protein_coding'
                        fixed_biotype += 1
        if fixed_biotype:
            logger.info(
                f"    Fixed biotype for {fixed_biotype} txTM transcripts with CDS (ref PC genes only)"
            )

    # augMP transcripts often have no reference/attrs biotype annotation. If they have a CDS, treat them
    # as protein_coding so they are not excluded from the coding pool solely due to missing biotype.
    augmp_fixed_biotype = 0
    augmp_mask = support_ref_tm_df['AlignmentId'].apply(tools.nameConversions.aln_id_is_augustus_mp)
    for idx in support_ref_tm_df[augmp_mask].index:
        aln_id = support_ref_tm_df.at[idx, 'AlignmentId']
        tx_obj = tx_dict.get(aln_id)
        if tx_obj is None:
            continue
        if getattr(tx_obj, 'cds_size', 0) > 0:
            gid = support_ref_tm_df.at[idx, 'GeneId']
            if ref_gene_biotype_map and gid and not _is_ref_protein_coding_gene(gid, ref_gene_biotype_map):
                continue
            current = support_ref_tm_df.at[idx, 'TranscriptBiotype']
            if pd.isna(current) or current != 'protein_coding':
                support_ref_tm_df.at[idx, 'TranscriptBiotype'] = 'protein_coding'
                augmp_fixed_biotype += 1
    if augmp_fixed_biotype:
        logger.info(
            f"    Fixed biotype for {augmp_fixed_biotype} augMP transcripts with CDS (ref PC or unknown gene)"
        )

    # transMap / transMap_pairwise: same CDS rule as txTM so full-length map models are not
    # dropped from the coding pool when attrs/reference biotype is missing.
    transmap_fixed_biotype = 0
    for idx in support_ref_tm_df.index:
        aln_id = support_ref_tm_df.at[idx, 'AlignmentId']
        mode = alignment_source_map.get(aln_id, '')
        if mode not in ('transMap', 'transMap_pairwise'):
            continue
        tx_obj = tx_dict.get(aln_id)
        if tx_obj is None:
            continue
        if getattr(tx_obj, 'cds_size', 0) > 0:
            gid = support_ref_tm_df.at[idx, 'GeneId']
            if ref_gene_biotype_map and not _is_ref_protein_coding_gene(gid, ref_gene_biotype_map):
                continue
            current = support_ref_tm_df.at[idx, 'TranscriptBiotype']
            if pd.isna(current) or current != 'protein_coding':
                support_ref_tm_df.at[idx, 'TranscriptBiotype'] = 'protein_coding'
                transmap_fixed_biotype += 1
    if transmap_fixed_biotype:
        logger.info(
            f"    Fixed biotype for {transmap_fixed_biotype} transMap transcripts with CDS (ref PC genes only)"
        )

    # Split into coding and non-coding
    # IMPORTANT: denovo transcripts don't have biotypes yet (assigned during classification)
    # so we keep them in the coding set. Only exclude non-denovo transcripts with missing biotypes.
    has_biotype = support_ref_tm_df['TranscriptBiotype'].notna()
    is_denovo_mask = support_ref_tm_df['AlignmentId'].astype(str).apply(_is_denovo)

    # Coding: ref protein_coding gene + protein_coding transcript, or denovo (even without biotype)
    if ref_gene_biotype_map:
        ref_pc_gene = support_ref_tm_df['GeneId'].apply(
            lambda g: _is_ref_protein_coding_gene(g, ref_gene_biotype_map)
        )
        pc_tx_mask = (
            has_biotype
            & ref_pc_gene
            & (support_ref_tm_df['TranscriptBiotype'] == 'protein_coding')
        )
    else:
        pc_tx_mask = has_biotype & (support_ref_tm_df['TranscriptBiotype'] == 'protein_coding')
    coding_df = support_ref_tm_df[pc_tx_mask | is_denovo_mask].copy()

    # Non-coding: transcripts with non-protein-coding biotypes (excluding denovo)
    non_coding_df = support_ref_tm_df[
        has_biotype & (support_ref_tm_df['TranscriptBiotype'] != 'protein_coding') & ~is_denovo_mask
    ].copy()
    
    # Excluded: non-denovo transcripts without biotypes (data quality issue)
    excluded_df = support_ref_tm_df[~has_biotype & ~is_denovo_mask].copy()

    denovo_in_coding = coding_df['AlignmentId'].astype(str).apply(_is_denovo).sum()
    logger.info(f"    Coding: {len(coding_df)} transcripts (including {denovo_in_coding} denovo)")
    logger.info(f"    Non-coding: {len(non_coding_df)} transcripts")
    if len(excluded_df) > 0:
        logger.warning(f"    Excluded {len(excluded_df)} transcripts with missing biotype")

    # Count unique genes in coding by source
    if len(coding_df) > 0:
        unique_coding_genes = coding_df['GeneId'].nunique()
        logger.info(f"    Unique protein-coding genes after biotype split: {unique_coding_genes}")
        # Break down by mode
        mode_gene_counts = {}
        mode_tx_counts = {}
        for mode in coding_df['AlignmentId'].apply(lambda x: alignment_source_map.get(x, 'unknown')).unique():
            mode_df = coding_df[coding_df['AlignmentId'].apply(lambda x: alignment_source_map.get(x, 'unknown') == mode)]
            mode_gene_counts[mode] = mode_df['GeneId'].nunique()
            mode_tx_counts[mode] = len(mode_df)
        logger.info(f"    By source (genes): {mode_gene_counts}")
        logger.info(f"    By source (transcripts): {mode_tx_counts}")
    
    # Merge with metrics
    logger.info(f"  Step 3: Merging with alignment metrics...")
    if len(mrna_metrics_df) > 0 and len(cds_metrics_df) > 0:
        metrics_df = pd.merge(mrna_metrics_df, cds_metrics_df, on='AlignmentId', how='outer', suffixes=['_mRNA', '_CDS'])
        coding_df = pd.merge(coding_df, metrics_df, on='AlignmentId', how='left')
        logger.info(f"    Merged {len(metrics_df)} metric records")
    
    # Merge with eval_df
    if len(eval_df) > 0:
        logger.info(f"  Step 4: Merging with evaluation data...")
        coding_df = pd.merge(coding_df, eval_df, on='AlignmentId', how='left')
        logger.info(f"    Merged {len(eval_df)} evaluation records")
    
    # Deduplicate by AlignmentId (many-to-many joins can create duplicates)
    logger.info(f"  Step 4b: Deduplicating by AlignmentId...")
    initial_count = len(coding_df)
    # Prioritize protein_coding biotype when deduplicating
    # (Some transcripts match multiple reference entries with different biotypes)
    coding_df['_is_protein_coding'] = (coding_df.get('TranscriptBiotype', 'protein_coding') == 'protein_coding').astype(int)
    coding_df = coding_df.sort_values('_is_protein_coding', ascending=False)
    coding_df = coding_df.drop_duplicates('AlignmentId', keep='first')
    coding_df = coding_df.drop('_is_protein_coding', axis=1)
    if initial_count != len(coding_df):
        logger.info(f"    Removed {initial_count - len(coding_df)} duplicate AlignmentId entries")
    
    # Biotype downgrade (fragment / processed_pseudogene) is now decided per-row
    # in create_transcript_attributes(), based on the cross-species-lenient policy:
    # lift the source biotype as-is unless (a) the alignment is so weak it is clearly
    # a fragment (coverage < 30% AND identity < 30%) or (b) a multi-exon source maps
    # to a near-intronless target (processed pseudogene signature).
    coding_df['is_fragment'] = False
    non_coding_df['is_fragment'] = False

    # Fill missing values
    logger.info(f"  Step 5: Filling missing values...")
    # augMP records now carry REAL miniprot-derived AlnCoverage / AlnIdentity
    # via the ``generate_augMP_psl`` rule. Any augMP record still showing NaN
    # at this point is one whose miniprot PSL row could not be matched (rare:
    # Augustus refined the locus into a model that diverges from the source
    # protein). We deliberately keep those as NaN here — the score / fragment
    # filtering downstream will honestly mark them low-confidence — rather
    # than falsely promoting them to perfect 100 % metrics.
    augmp_mask_now = coding_df['AlignmentId'].apply(tools.nameConversions.aln_id_is_augustus_mp)
    non_augmp_mask = ~augmp_mask_now
    augmp_na_before = 0
    if 'AlnCoverage_mRNA' in coding_df.columns:
        augmp_na_before = coding_df.loc[augmp_mask_now, 'AlnCoverage_mRNA'].isna().sum()

    # Apply txTM coverage filter BEFORE fillna(100) on missing metrics (otherwise strict mode is a no-op).
    txTM_min_cov = float(getattr(args, 'txTM_min_coverage', 0.0) or 0.0)
    txTM_strict = bool(getattr(args, 'txTM_strict_metrics', False))
    txTM_min_cov_orphan = getattr(args, 'txTM_min_coverage_no_transmap', None)
    if txTM_min_cov_orphan is not None:
        txTM_min_cov_orphan = float(txTM_min_cov_orphan)
    txTM_strict_orphan = getattr(args, 'txTM_strict_metrics_no_transmap', None)
    txTM_min_cov_nc = getattr(args, 'txTM_min_coverage_noncoding', None)
    if txTM_min_cov_nc is not None:
        txTM_min_cov_nc = float(txTM_min_cov_nc)
    txTM_min_cov_orphan_nc = getattr(args, 'txTM_min_coverage_no_transmap_noncoding', None)
    if txTM_min_cov_orphan_nc is not None:
        txTM_min_cov_orphan_nc = float(txTM_min_cov_orphan_nc)
    txTM_strict_orphan_nc = getattr(args, 'txTM_strict_metrics_no_transmap_noncoding', None)
    coding_df, non_coding_df = filter_txTM_by_coverage(
        coding_df, non_coding_df, alignment_source_map,
        txTM_min_coverage=txTM_min_cov,
        txTM_strict_metrics=txTM_strict,
        txTM_min_coverage_no_transmap=txTM_min_cov_orphan,
        txTM_strict_metrics_no_transmap=txTM_strict_orphan,
        txTM_min_coverage_noncoding=txTM_min_cov_nc,
        txTM_min_coverage_no_transmap_noncoding=txTM_min_cov_orphan_nc,
        txTM_strict_metrics_no_transmap_noncoding=txTM_strict_orphan_nc,
    )
    logger.info(
        f"    After txTM coverage filter: {len(coding_df)} coding, "
        f"{len(non_coding_df)} non-coding"
    )

    augMP_min_cov_orphan = float(getattr(args, 'augMP_min_coverage_no_anchor', 80) or 0)
    augMP_strict_orphan = getattr(args, 'augMP_strict_metrics_no_anchor', None)
    if augMP_strict_orphan is None:
        augMP_strict_orphan = True
    coding_df, non_coding_df = filter_augMP_by_coverage(
        coding_df, non_coding_df, alignment_source_map,
        augMP_min_coverage_no_anchor=augMP_min_cov_orphan,
        augMP_strict_metrics_no_anchor=augMP_strict_orphan,
    )
    logger.info(
        f"    After augMP coverage filter: {len(coding_df)} coding, "
        f"{len(non_coding_df)} non-coding"
    )

    for col in ['AlnCoverage_mRNA', 'AlnIdentity_mRNA', 'AlnCoverage_CDS', 'AlnIdentity_CDS']:
        if col in coding_df.columns:
            coding_df.loc[non_augmp_mask, col] = coding_df.loc[non_augmp_mask, col].fillna(100.0)
    for col in ['AlnGoodness_mRNA', 'AlnGoodness_CDS']:
        if col in coding_df.columns:
            coding_df.loc[non_augmp_mask, col] = coding_df.loc[non_augmp_mask, col].fillna(1.0)
    if 'Frameshift' in coding_df.columns:
        coding_df['Frameshift'] = coding_df['Frameshift'].fillna(False)
    if augmp_mask_now.any():
        logger.info(
            f"    augMP records: {augmp_mask_now.sum()} total; "
            f"{augmp_na_before} still NaN after miniprot-PSL metrics merge "
            f"(left as NaN; consensus filters will treat them honestly)"
        )
    
    logger.info(f"    After processing: {len(coding_df)} coding, {len(non_coding_df)} non-coding")

    # Score transcripts
    logger.info(f"  Step 6: Scoring transcripts...")

    scored_coding_df = score_transcripts(coding_df, args.in_species_rna_support_only)
    scored_non_coding_df = score_transcripts(non_coding_df, args.in_species_rna_support_only, is_coding=False)
    logger.info(f"    Scored {len(scored_coding_df)} coding + {len(scored_non_coding_df)} non-coding = {len(scored_coding_df) + len(scored_non_coding_df)} total")

    # Track genes after scoring
    if len(scored_coding_df) > 0:
        unique_scored_pc_genes = scored_coding_df['GeneId'].nunique()
        logger.info(f"    Unique protein-coding genes after scoring: {unique_scored_pc_genes}")
        # Break down by mode
        mode_gene_counts_scored = {}
        for mode in scored_coding_df['AlignmentId'].apply(lambda x: alignment_source_map.get(x, 'unknown')).unique():
            mode_df = scored_coding_df[scored_coding_df['AlignmentId'].apply(lambda x: alignment_source_map.get(x, 'unknown') == mode)]
            mode_gene_counts_scored[mode] = mode_df['GeneId'].nunique()
        logger.info(f"    By source after scoring: {mode_gene_counts_scored}")
    
    # Combine
    scored_df = pd.concat([scored_coding_df, scored_non_coding_df], ignore_index=True)
    
    if len(scored_df) == 0:
        logger.warning(f"  No transcripts passed scoring for {chrom}")
        return []
    
    # OPTIMIZATION: Two-pass approach for denovo classification
    # Pass 1: Select consensus from non-denovo sources only (fast)
    # Pass 2: Classify denovo against selected consensus (much fewer comparisons)
    
    # Split scored_df into denovo vs non-denovo
    is_denovo_mask = scored_df['AlignmentId'].apply(_is_denovo)
    denovo_df = scored_df[is_denovo_mask].copy()
    non_denovo_df = scored_df[~is_denovo_mask].copy()
    
    logger.info(f"  Step 7: Tracking multi-locus mappings (non-denovo only)...")
    multi_locus_tracker = track_multi_locus_mappings(non_denovo_df, tx_dict, alignment_source_map)
    logger.info(f"    Found {len(multi_locus_tracker)} unique transcript-mode combinations (non-denovo)")

    txTM_anchor_ov = float(getattr(args, 'txTM_transmap_anchor_overlap', 0.0) or 0.0)
    if txTM_anchor_ov > 0:
        prune_txTM_against_transmap(
            multi_locus_tracker, tx_dict, alignment_source_map,
            min_overlap=txTM_anchor_ov,
        )

    # PASS 1: Select consensus from non-denovo sources
    logger.info(f"  Step 8: Selecting consensus from non-denovo sources (CNV threshold: {args.cnv_score_similarity})...")
    logger.info(f"    Processing {len(non_denovo_df)} non-denovo transcripts...")
    consensus_transcripts = select_consensus_with_cnv(
        non_denovo_df, tx_dict, alignment_source_map, multi_locus_tracker, 
        args.cnv_score_similarity, metrics, consensus_dict={}, args=args, ref_df=ref_df, support_df=None
    )
    logger.info(f"    ✓ Selected {len(consensus_transcripts)} consensus transcripts from non-denovo sources")
    
    # PASS 2: Classify and select denovo against the smaller consensus set
    denovo_modes_active = set(args.denovo_tx_modes) if args and hasattr(args, 'denovo_tx_modes') and args.denovo_tx_modes else set()
    if len(denovo_df) > 0 and denovo_modes_active:
        logger.info(f"  Step 9: Processing {len(denovo_df)} denovo transcripts against consensus (modes: {denovo_modes_active})...")
        
        # Build consensus_dict from selected consensus (much smaller than support_df!)
        consensus_dict = {aln_id: attrs for aln_id, attrs in consensus_transcripts}
        logger.info(f"    Built consensus_dict with {len(consensus_dict)} transcripts")
        
        # Process denovo transcripts
        denovo_consensus = process_augpb_transcripts(
            denovo_df, tx_dict, alignment_source_map, consensus_dict, 
            args, metrics, ref_df
        )
        logger.info(f"    ✓ Selected {len(denovo_consensus)} denovo transcripts")
        
        # Combine with non-denovo consensus
        consensus_transcripts.extend(denovo_consensus)
        logger.info(f"    ✓ Total consensus: {len(consensus_transcripts)} transcripts")
    elif len(denovo_df) > 0:
        logger.info(f"  Step 9: Skipping {len(denovo_df)} denovo transcripts (denovo_tx_modes not configured)")
    else:
        logger.info(f"  Step 9: No denovo transcripts to process")
    
    elapsed = time.time() - chrom_start_time
    logger.info(f"  ✓ Chromosome {chrom} processed in {elapsed:.1f}s")
    logger.info(f"    Returned {len(consensus_transcripts)} consensus transcripts")
    
    # Debug: show which sources contributed
    source_counts = collections.Counter()
    pc_gene_ids_in_consensus = set()
    gene_ids_by_mode_in_consensus = collections.defaultdict(set)
    
    for aln_id, attrs in consensus_transcripts:
        # Get mode from attrs (more reliable than alignment_source_map)
        mode = attrs.get('alignment_mode', alignment_source_map.get(aln_id, 'unknown'))
        source_counts[mode] += 1
        
        # Track gene IDs
        gene_id = attrs.get('source_gene')
        gene_biotype = attrs.get('gene_biotype')
        if gene_id and gene_id != 'N/A':
            gene_ids_by_mode_in_consensus[mode].add(gene_id)
            if gene_biotype == 'protein_coding':
                pc_gene_ids_in_consensus.add(gene_id)
    
    # Track protein-coding genes by source (for accurate reporting)
    pc_gene_ids_by_mode_in_consensus = collections.defaultdict(set)
    for aln_id, attrs in consensus_transcripts:
        mode = attrs.get('alignment_mode', alignment_source_map.get(aln_id, 'unknown'))
        gene_id = attrs.get('source_gene')
        gene_biotype = attrs.get('gene_biotype')
        if gene_id and gene_id != 'N/A' and gene_biotype == 'protein_coding':
            pc_gene_ids_by_mode_in_consensus[mode].add(gene_id)
    
    logger.info(f"    By source (transcripts): {dict(source_counts)}")
    logger.info(f"    Protein-coding genes in consensus for this chrom: {len(pc_gene_ids_in_consensus)}")
    logger.info(f"    By source (protein-coding genes): {dict((k, len(v)) for k, v in pc_gene_ids_by_mode_in_consensus.items())}")
    logger.info(f"    By source (all genes): {dict((k, len(v)) for k, v in gene_ids_by_mode_in_consensus.items())}")
    
    # Perform deduplication per-chromosome (critical for multi-chromosome genomes like HG03456 with 2 Y chromosomes)
    logger.info(f"  Step 9: Deduplicating within chromosome...")
    consensus_dict_per_chrom = {aln_id: attrs for aln_id, attrs in consensus_transcripts}
    deduplicated_dict = deduplicate_consensus(consensus_dict_per_chrom, tx_dict, metrics)
    deduplicated_transcripts = list(deduplicated_dict.items())
    logger.info(f"    After per-chromosome deduplication: {len(deduplicated_transcripts)} transcripts (removed {len(consensus_transcripts) - len(deduplicated_transcripts)})")
    
    return deduplicated_transcripts


def process_augpb_transcripts(augpb_df, tx_dict, alignment_source_map, consensus_dict, 
                              args, metrics, ref_df):
    """
    Process augPB transcripts against selected consensus AND all input transcripts.
    
    CRITICAL: augPB transcripts can only be novel isoforms if they overlap with known genes,
    and can only be novel genes if NO known genes exist at that locus (even if those known
    genes didn't make it into consensus due to missing metrics).
    
    Args:
        augpb_df: DataFrame of augPB transcripts to process
        tx_dict: Dictionary of all transcript objects (includes ALL input transcripts)
        alignment_source_map: Mapping of alignment IDs to sources
        consensus_dict: Dictionary of selected consensus transcripts {aln_id: attrs}
        args: Command-line arguments
        metrics: Metrics dictionary
        ref_df: Reference annotation dataframe
    
    Returns:
        List of (aln_id, attrs) tuples for selected augPB transcripts
    """
    import logging
    import time
    import collections
    
    logger = logging.getLogger(__name__)
    
    # Build spatial index from consensus_dict (for finding overlapping genes)
    logger.info(f"      Building consensus index ({len(consensus_dict)} transcripts)...")
    index_start = time.time()
    pc_gene_index = collections.defaultdict(list)
    
    for consensus_aln_id, consensus_attrs in consensus_dict.items():
        if consensus_aln_id in tx_dict:
            tx = tx_dict[consensus_aln_id]
            # Store: (start, end, aln_id, tx_obj, attrs)
            pc_gene_index[tx.chromosome].append((tx.start, tx.stop, consensus_aln_id, tx, consensus_attrs))
    
    # Sort by start position for efficient overlap checking
    for chrom in pc_gene_index:
        pc_gene_index[chrom].sort(key=lambda x: x[0])
    
    # CRITICAL: Also build index of ALL non-denovo input transcripts
    # This ensures we don't create denovo genes where txTM/transMap/augTM exist,
    # even if those transcripts didn't make it into consensus
    logger.info(f"      Building ALL input transcripts index...")
    all_known_tx_index = collections.defaultdict(list)
    
    for aln_id, tx in tx_dict.items():
        if not _is_denovo(aln_id):
            # This is a known source transcript (txTM/transMap/augTM)
            all_known_tx_index[tx.chromosome].append((tx.start, tx.stop, aln_id, tx))
    
    # Sort by start position
    for chrom in all_known_tx_index:
        all_known_tx_index[chrom].sort(key=lambda x: x[0])
    
    total_known_tx = sum(len(txs) for txs in all_known_tx_index.values())
    
    index_elapsed = time.time() - index_start
    logger.info(f"      ✓ Built indexes in {index_elapsed:.1f}s ({len(consensus_dict)} consensus + {total_known_tx} known input)")
    
    # Process each augPB transcript
    logger.info(f"      Classifying {len(augpb_df)} augPB transcripts...")
    selected_augpb = []
    classified_count = 0
    discarded_count = 0
    discarded_overlap_known = 0
    
    for idx, row in augpb_df.iterrows():
        aln_id = row['AlignmentId']
        mode = alignment_source_map.get(aln_id, 'augPB')
        
        # Get normalized ID
        normalized_id = tools.nameConversions.strip_alignment_numbers(aln_id)
        
        # CRITICAL CHECK: Does this augPB overlap with ANY known transcript in the input?
        tx_obj = tx_dict.get(aln_id)
        overlaps_known_input = False
        if tx_obj:
            chrom_known_txs = all_known_tx_index.get(tx_obj.chromosome, [])
            for start, end, known_aln_id, known_tx in chrom_known_txs:
                if start > tx_obj.stop:
                    break  # Sorted, no more overlaps
                if tx_obj.start <= end and tx_obj.stop >= start and tx_obj.strand == known_tx.strand:
                    # Overlaps with a known transcript on same strand
                    overlaps_known_input = True
                    break
        
        # Create transcript attributes (includes classification)
        attrs = create_transcript_attributes(
            row, mode, normalized_id,
            tx_dict=tx_dict,
            consensus_dict=consensus_dict,
            args=args,
            metrics=metrics,
            pc_gene_index=pc_gene_index,
            ref_df=ref_df
        )
        
        # Skip if classification returned None (doesn't meet criteria)
        if attrs is None or attrs.get('transcript_class') is None:
            discarded_count += 1
            continue
        
        # CRITICAL: If this augPB overlaps with known input but is NOT a novel isoform, discard it
        # This prevents augPB from replacing txTM/transMap/augTM genes
        if overlaps_known_input and attrs.get('transcript_class') != 'putative_novel_isoform':
            discarded_count += 1
            discarded_overlap_known += 1
            continue
        
        classified_count += 1
        
        # Store per-transcript metrics (same as non-denovo in select_consensus_with_cnv)
        cov = row.get('AlnCoverage_mRNA', 0)
        ident = row.get('AlnIdentity_mRNA', 0)
        attrs['_metrics'] = {
            'Coverage': float(cov) if pd.notna(cov) else 0,
            'Identity': float(ident) if pd.notna(ident) else 0,
            'Splice Support': float(row.get('IntronRnaSupportPercent', 0) or 0),
            'Exon Support': float(row.get('ExonRnaSupportPercent', 0) or 0),
            'Original Introns': float(row.get('OriginalIntronsPercent_mRNA', 0) or 0),
            'Splice Annotation Support': float(row.get('IntronAnnotSupportPercent', 0) or 0),
            'Exon Annotation Support': float(row.get('ExonAnnotSupportPercent', 0) or 0),
        }
        
        selected_augpb.append((aln_id, attrs))
    
    logger.info(f"      ✓ Classified: {classified_count}, Discarded: {discarded_count}")
    if discarded_overlap_known > 0:
        logger.info(f"         (including {discarded_overlap_known} that overlap known input genes)")
    
    return selected_augpb


def _tx_exon_intervals(tx_obj):
    """Return sorted (start, end) exon intervals for a genePred."""
    if tx_obj is None:
        return []
    return sorted((iv.start, iv.stop) for iv in tx_obj.exon_intervals)


def _exon_overlap_fraction(tx_a, tx_b):
    """Fraction of tx_a exonic length overlapping tx_b exons (same strand required)."""
    if tx_a is None or tx_b is None or tx_a.strand != tx_b.strand:
        return 0.0
    a_ivs = _tx_exon_intervals(tx_a)
    b_ivs = _tx_exon_intervals(tx_b)
    if not a_ivs:
        return 0.0
    a_len = sum(e - s for s, e in a_ivs)
    if a_len <= 0:
        return 0.0
    overlap = 0
    for s1, e1 in a_ivs:
        for s2, e2 in b_ivs:
            overlap += max(0, min(e1, e2) - max(s1, s2))
    return overlap / a_len


def build_transmap_anchored_transcript_ids(alignment_source_map):
    """Normalized transcript accessions that have a transMap (or pairwise) model."""
    import re
    anchored = set()
    for aln_id, mode in alignment_source_map.items():
        if mode not in ('transMap', 'transMap_pairwise'):
            continue
        for variant in (
            aln_id,
            aln_id.replace('rna-', ''),
            normalize_alignment_id(aln_id, mode),
        ):
            anchored.add(variant)
            anchored.add(re.sub(r'\.\d+$', '', variant))
    return anchored


def _txtm_id_variants(aln_id):
    """Normalization keys for matching a txTM row to a transMap anchor."""
    import re
    base = normalize_alignment_id(aln_id, 'txTM')
    variants = {base, base.replace('rna-', '')}
    for v in list(variants):
        variants.add(re.sub(r'\.\d+$', '', v))
    return variants


def txTM_has_transmap_anchor(aln_id, anchored_ids):
    return bool(_txtm_id_variants(aln_id) & anchored_ids)


def _augmp_id_variants(aln_id):
    """Normalization keys for matching an augMP row to transMap/txTM anchors."""
    import re
    base = normalize_alignment_id(aln_id, 'augMP')
    variants = {base, base.replace('rna-', '')}
    for v in list(variants):
        variants.add(re.sub(r'\.\d+$', '', v))
        variants.add(re.sub(r'_\d+$', '', v))
    return variants


def build_reliable_anchor_ids(alignment_source_map):
    """Transcript accessions with transMap or txTM (modes used to validate orphans)."""
    import re
    anchored = build_transmap_anchored_transcript_ids(alignment_source_map)
    for aln_id, mode in alignment_source_map.items():
        if mode != 'txTM':
            continue
        for variant in _txtm_id_variants(aln_id):
            anchored.add(variant)
            anchored.add(re.sub(r'\.\d+$', '', variant))
    return anchored


def augMP_has_transmap_anchor(aln_id, transmap_anchor_ids):
    """True when transMap exists for this transcript (txTM alone does not exempt augMP)."""
    return bool(_augmp_id_variants(aln_id) & transmap_anchor_ids)


def filter_txTM_by_coverage(
    coding_df,
    non_coding_df,
    alignment_source_map,
    txTM_min_coverage=0.0,
    txTM_strict_metrics=False,
    txTM_min_coverage_no_transmap=None,
    txTM_strict_metrics_no_transmap=None,
    txTM_min_coverage_noncoding=None,
    txTM_min_coverage_no_transmap_noncoding=None,
    txTM_strict_metrics_no_transmap_noncoding=None,
):
    """Drop low-coverage txTM rows before scoring.

    Uses strict thresholds when the same transcript accession has a transMap
    anchor; relaxed thresholds for txTM-only loci (no transMap isoform on this genome).
    Non-coding rows may use separate, typically lower, thresholds.
    """
    if txTM_min_coverage_no_transmap is None:
        txTM_min_coverage_no_transmap = 80.0
    if txTM_strict_metrics_no_transmap is None:
        txTM_strict_metrics_no_transmap = True

    nc_min_cov = txTM_min_coverage if txTM_min_coverage_noncoding is None else txTM_min_coverage_noncoding
    nc_orphan_cov = (
        txTM_min_coverage_no_transmap
        if txTM_min_coverage_no_transmap_noncoding is None
        else txTM_min_coverage_no_transmap_noncoding
    )
    nc_orphan_strict = (
        txTM_strict_metrics_no_transmap
        if txTM_strict_metrics_no_transmap_noncoding is None
        else txTM_strict_metrics_no_transmap_noncoding
    )

    if (
        txTM_min_coverage <= 0
        and not txTM_strict_metrics
        and txTM_min_coverage_no_transmap <= 0
        and not txTM_strict_metrics_no_transmap
        and nc_min_cov <= 0
        and not txTM_strict_metrics
        and nc_orphan_cov <= 0
        and not nc_orphan_strict
    ):
        return coding_df, non_coding_df

    txTM_aln_ids = {aln for aln, src in alignment_source_map.items() if src == 'txTM'}
    if not txTM_aln_ids:
        return coding_df, non_coding_df

    anchored_ids = build_transmap_anchored_transcript_ids(alignment_source_map)

    def _apply(
        df,
        label,
        min_cov,
        strict_metrics,
        min_cov_orphan,
        strict_orphan,
    ):
        if len(df) == 0 or 'AlnCoverage_mRNA' not in df.columns:
            return df
        is_txTM = df['AlignmentId'].isin(txTM_aln_ids)
        if not is_txTM.any():
            return df
        cov = df['AlnCoverage_mRNA']
        has_anchor = df['AlignmentId'].apply(
            lambda aid: txTM_has_transmap_anchor(aid, anchored_ids) if aid in txTM_aln_ids else False
        )
        anchored_txTM = is_txTM & has_anchor
        orphan_txTM = is_txTM & ~has_anchor

        keep = ~is_txTM
        if anchored_txTM.any():
            if strict_metrics:
                keep_anchored = anchored_txTM & cov.notna() & (cov > min_cov)
            else:
                keep_anchored = anchored_txTM & (cov.fillna(100.0) > min_cov)
            keep = keep | keep_anchored
        if orphan_txTM.any():
            if strict_orphan:
                keep_orphan = orphan_txTM & cov.notna() & (cov >= min_cov_orphan)
            else:
                # Allow missing metrics; only drop explicit low-coverage fragments.
                keep_orphan = orphan_txTM & (cov.isna() | (cov >= min_cov_orphan))
            keep = keep | keep_orphan

        n_drop = int(is_txTM.sum() - keep[is_txTM].sum())
        n_drop_anchored = int(anchored_txTM.sum() - keep[anchored_txTM].sum())
        n_drop_orphan = int(orphan_txTM.sum() - keep[orphan_txTM].sum())
        if n_drop > 0:
            logger.info(
                f"    Filtered {n_drop} txTM {label} transcripts "
                f"(anchored: {n_drop_anchored} dropped with "
                f"min_cov={min_cov}, strict={strict_metrics}; "
                f"txTM-only: {n_drop_orphan} dropped with "
                f"min_cov={min_cov_orphan}, "
                f"strict={strict_orphan}; "
                f"{int(keep[orphan_txTM].sum())} txTM-only kept)"
            )
        return df[keep].copy()

    return (
        _apply(
            coding_df,
            'coding',
            txTM_min_coverage,
            txTM_strict_metrics,
            txTM_min_coverage_no_transmap,
            txTM_strict_metrics_no_transmap,
        ),
        _apply(
            non_coding_df,
            'non-coding',
            nc_min_cov,
            txTM_strict_metrics,
            nc_orphan_cov,
            nc_orphan_strict,
        ),
    )


def filter_augMP_by_coverage(
    coding_df,
    non_coding_df,
    alignment_source_map,
    augMP_min_coverage_no_anchor=80.0,
    augMP_strict_metrics_no_anchor=True,
):
    """Drop low-coverage augMP rows when no transMap anchor exists for that transcript."""
    if augMP_min_coverage_no_anchor <= 0 and not augMP_strict_metrics_no_anchor:
        return coding_df, non_coding_df

    augmp_aln_ids = {
        aln for aln, src in alignment_source_map.items()
        if src == 'augMP' or tools.nameConversions.aln_id_is_augustus_mp(aln)
    }
    if not augmp_aln_ids:
        return coding_df, non_coding_df

    transmap_anchor_ids = build_transmap_anchored_transcript_ids(alignment_source_map)

    def _apply(df, label):
        if len(df) == 0 or 'AlnCoverage_mRNA' not in df.columns:
            return df
        is_augmp = df['AlignmentId'].apply(
            lambda aid: aid in augmp_aln_ids or tools.nameConversions.aln_id_is_augustus_mp(aid)
        )
        if not is_augmp.any():
            return df
        cov = df['AlnCoverage_mRNA']
        has_anchor = df['AlignmentId'].apply(
            lambda aid: augMP_has_transmap_anchor(aid, transmap_anchor_ids)
            if (aid in augmp_aln_ids or tools.nameConversions.aln_id_is_augustus_mp(aid))
            else False
        )
        anchored_augmp = is_augmp & has_anchor
        orphan_augmp = is_augmp & ~has_anchor

        keep = ~is_augmp | anchored_augmp
        if orphan_augmp.any():
            if augMP_strict_metrics_no_anchor:
                keep_orphan = orphan_augmp & cov.notna() & (cov >= augMP_min_coverage_no_anchor)
            else:
                keep_orphan = orphan_augmp & (cov.isna() | (cov >= augMP_min_coverage_no_anchor))
            keep = keep | keep_orphan

        n_drop = int(is_augmp.sum() - keep[is_augmp].sum())
        n_drop_orphan = int(orphan_augmp.sum() - keep[orphan_augmp].sum())
        if n_drop > 0:
            logger.info(
                f"    Filtered {n_drop} augMP {label} transcripts "
                f"(augMP-only: {n_drop_orphan} dropped with "
                f"min_cov>={augMP_min_coverage_no_anchor}, "
                f"strict={augMP_strict_metrics_no_anchor}; "
                f"{int(keep[orphan_augmp].sum())} augMP-only kept)"
            )
        return df[keep].copy()

    return _apply(coding_df, 'coding'), _apply(non_coding_df, 'non-coding')


def prune_txTM_against_transmap(multi_locus_tracker, tx_dict, alignment_source_map,
                                min_overlap=0.25):
    """Remove txTM loci that disagree with an existing transMap copy of the same transcript.

    True paralog/CNV copies at separate loci are kept when they lack a transMap anchor
    or when their exon structure overlaps the transMap locus (same gene copy).
    """
    if min_overlap <= 0:
        return 0

    transmap_modes = {'transMap', 'transMap_pairwise'}
    tm_loci = {}
    for key, loci in multi_locus_tracker.items():
        norm_id, mode, chrom = key
        if mode not in transmap_modes:
            continue
        tm_loci[(norm_id, chrom)] = loci

    pruned = 0
    for key in list(multi_locus_tracker.keys()):
        norm_id, mode, chrom = key
        if mode != 'txTM':
            continue
        anchor = tm_loci.get((norm_id, chrom))
        if not anchor:
            continue
        anchor_txs = [tx_dict.get(aln_id) for aln_id, _, _ in anchor]
        anchor_txs = [t for t in anchor_txs if t is not None]
        if not anchor_txs:
            continue
        kept = []
        for aln_id, locus, score in multi_locus_tracker[key]:
            tx_obj = tx_dict.get(aln_id)
            if tx_obj is None:
                continue
            best_ov = max(_exon_overlap_fraction(tx_obj, a) for a in anchor_txs)
            if best_ov >= min_overlap:
                kept.append((aln_id, locus, score))
            else:
                pruned += 1
        multi_locus_tracker[key] = kept

    if pruned > 0:
        logger.info(
            f"      Pruned {pruned} txTM loci with <{min_overlap:.0%} exon overlap vs transMap "
            f"(same transcript id)"
        )
    return pruned


def track_multi_locus_mappings(scored_df, tx_dict, alignment_source_map):
    """
    Track how many loci each transcript-mode combination maps to.
    
    Returns a dict: (normalized_tx_id, mode) -> [(aln_id, locus_info, score), ...]
    """
    multi_locus_map = collections.defaultdict(list)
    
    # Vectorized approach: precompute modes and normalized IDs
    aln_ids = scored_df['AlignmentId'].values
    modes = [alignment_source_map.get(aln_id, 'transMap') for aln_id in aln_ids]
    normalized_ids = [normalize_alignment_id(aln_id, mode) for aln_id, mode in zip(aln_ids, modes)]
    scores = scored_df.get('TranscriptScore', pd.Series([0] * len(scored_df))).values
    
    # Build multi_locus_map
    for aln_id, mode, normalized_id, score in zip(aln_ids, modes, normalized_ids, scores):
        tx_obj = tx_dict.get(aln_id)
        if not tx_obj:
            continue
        
        locus_info = {
            'chromosome': tx_obj.chromosome,
            'start': tx_obj.start,
            'stop': tx_obj.stop,
            'strand': tx_obj.strand
        }
        
        # Include chromosome in key to prevent genes on different chromosomes
        # from being treated as CNV copies (important for individuals with multiple Y chromosomes)
        key = (normalized_id, mode, tx_obj.chromosome)
        multi_locus_map[key].append((aln_id, locus_info, score))
    
    # Log statistics
    multi_locus_count = sum(1 for v in multi_locus_map.values() if len(v) > 1)
    single_locus_count = sum(1 for v in multi_locus_map.values() if len(v) == 1)
    logger.info(f"      Single-locus: {single_locus_count}, Multi-locus: {multi_locus_count}")

    return multi_locus_map


def select_consensus_with_cnv(scored_df, tx_dict, alignment_source_map, multi_locus_tracker, 
                              score_threshold, metrics, consensus_dict=None, args=None, ref_df=None, support_df=None):
    """
    Select consensus transcripts, keeping CNV copies if scores are similar.
    
    Logic:
    1. For each transcript-mode combination, find max score from CORE sources only
    2. Keep all copies with core_score >= threshold * max_core_score
    3. Core sources: txTM, transMap, augTM (augPB doesn't count for CNV threshold)
    4. For single-locus transcripts, just pick the best
    5. CRITICAL: Prioritize known sources (txTM/transMap/augTM) over augPB for same gene
    
    Args:
        support_df: Full dataframe with ALL transcripts (for building pc_gene_index)
    """
    consensus_transcripts = []
    total_kept = 0
    total_discarded = 0
    augpb_classified = 0
    augpb_discarded = 0
    
    import logging
    import time
    logger = logging.getLogger(__name__)
    
    # OPTIMIZATION: Skip spatial index building in Pass 1 (no augPB transcripts present)
    # Spatial index is only needed when classifying augPB transcripts
    pc_gene_index = None
    # Skip index building - not needed for non-augPB transcripts
    
    # OPTIMIZATION: Build AlignmentId index for fast row lookup
    # This replaces slow DataFrame filtering (scored_df[scored_df['AlignmentId'] == aln_id])
    logger.info(f"      Building AlignmentId index for fast lookups...")
    aln_id_to_row = {}
    # Track source priority (lower is better): transMap=0, txTM=1, augTM=2, augMP=2, augPB=9
    source_priority = {'transMap': 0, 'transMap_pairwise': 0, 'txTM': 1,
                      'augTM': 2, 'augTMR': 2, 'augMP': 2, 'augPB': 9, 'strg': 9}
    for idx, row in scored_df.iterrows():
        aln_id = row['AlignmentId']
        mode = alignment_source_map.get(aln_id, 'transMap')
        priority = source_priority.get(mode, 5)
        aln_id_to_row[aln_id] = (row, priority)
    logger.info(f"      ✓ Indexed {len(aln_id_to_row)} transcripts")

    # Group by normalized transcript ID, mode, and chromosome
    for (normalized_id, mode, chrom), loci_list in multi_locus_tracker.items():
        if len(loci_list) == 0:
            continue

        # Calculate core scores (txTM/transMap/augTM only) for CNV threshold
        # Don't let augPB bias the multi-locus selection
        loci_with_core_scores = []
        for aln_id, locus, total_score in loci_list:
            # OPTIMIZATION: Use dict lookup instead of DataFrame filtering
            if aln_id not in aln_id_to_row:
                loci_with_core_scores.append((aln_id, locus, total_score, 0))
                continue
            
            # Calculate score from core sources only
            # A locus can have multiple transcripts (txTM + augPB + augTM), sum core sources
            locus_alignment_ids = [aln_id]
            # Check if there are other sources at this locus (look at original_gene_id)
            # For now, just calculate score from this alignment's mode
            if mode in ['txTM', 'transMap', 'augTM', 'augTMR', 'augMP']:
                core_score = total_score  # This mode is a core source
            else:
                core_score = 0  # augPB doesn't count for CNV threshold
            
            loci_with_core_scores.append((aln_id, locus, total_score, core_score))
        
        # Find max core score (ignore augPB-only loci)
        core_scores = [core_score for _, _, _, core_score in loci_with_core_scores if core_score > 0]
        if core_scores:
            max_core_score = max(core_scores)
            min_score_threshold = max_core_score * score_threshold
        else:
            # All loci are augPB-only, use total score
            max_core_score = max(total_score for _, _, total_score, _ in loci_with_core_scores)
            min_score_threshold = max_core_score * score_threshold
        
        # Filter to keep similar-scoring copies (based on core scores)
        kept_loci = [(aln_id, locus, total_score) for aln_id, locus, total_score, core_score in loci_with_core_scores
                     if (core_score >= min_score_threshold if core_scores else total_score >= min_score_threshold)]
        
        if len(loci_list) > 1:
            metrics['Multi-locus mappings'][mode] += 1
            metrics['Multi-locus kept'][mode] += len(kept_loci)
            
            if len(kept_loci) < len(loci_list):
                total_discarded += (len(loci_list) - len(kept_loci))
                # Log NOTCH2NLB for debugging
                if 'ENST00000593495' in str(normalized_id) or 'ENST00000850899' in str(normalized_id) or '286019' in str(normalized_id):
                    logger.info(f"      Multi-locus {normalized_id} ({mode}): kept {len(kept_loci)}/{len(loci_list)} copies")
                    logger.info(f"        Max core score: {max_core_score if core_scores else max(total_score for _, _, total_score, _ in loci_with_core_scores):.2f}, threshold: {min_score_threshold:.2f}")
                    for aln_id, locus, total_score, core_score in loci_with_core_scores:
                        kept_status = "KEPT" if (aln_id, locus, total_score) in kept_loci else "DISCARDED"
                        logger.info(f"        {aln_id} at {locus}: total={total_score:.2f}, core={core_score:.2f} [{kept_status}]")
        
        # Convert kept loci to consensus format
        for aln_id, locus, score in kept_loci:
            # OPTIMIZATION: Use dict lookup instead of DataFrame filtering
            row_data = aln_id_to_row.get(aln_id)
            if row_data is None:
                continue
            
            # Unpack tuple (row, priority) - we only need the row
            row, priority = row_data
            
            is_augpb = _is_denovo(str(aln_id))
            
            attrs = create_transcript_attributes(row, mode, normalized_id, 
                                                tx_dict=tx_dict, 
                                                consensus_dict=consensus_dict if consensus_dict else {},
                                                args=args,
                                                metrics=metrics,
                                                pc_gene_index=pc_gene_index,
                                                ref_df=ref_df)
            # Skip if classification returned None (doesn't meet criteria)
            if attrs.get('transcript_class') is None:
                if is_augpb:
                    augpb_discarded += 1
                continue
            
            if is_augpb:
                augpb_classified += 1
            
            total_kept += 1
            
            # Store per-transcript plotting metrics in attrs (collected into final metrics later)
            cov = row.get('AlnCoverage_mRNA', 0)
            ident = row.get('AlnIdentity_mRNA', 0)
            attrs['_metrics'] = {
                'Coverage': float(cov) if pd.notna(cov) else 0,
                'Identity': float(ident) if pd.notna(ident) else 0,
                'Splice Support': float(row.get('IntronRnaSupportPercent', 0) or 0),
                'Exon Support': float(row.get('ExonRnaSupportPercent', 0) or 0),
                'Original Introns': float(row.get('OriginalIntronsPercent_mRNA', 0) or 0),
                'Splice Annotation Support': float(row.get('IntronAnnotSupportPercent', 0) or 0),
                'Exon Annotation Support': float(row.get('ExonAnnotSupportPercent', 0) or 0),
            }
            consensus_transcripts.append((aln_id, attrs))
    
    logger.info(f"      Selected {total_kept} transcripts, discarded {total_discarded} low-scoring CNV copies")
    if augpb_classified > 0 or augpb_discarded > 0:
        logger.info(f"      AugPB: {augpb_classified} classified, {augpb_discarded} discarded")

    # Count by mode for debugging
    mode_counts = collections.Counter()
    for aln_id, attrs in consensus_transcripts:
        mode_counts[attrs.get('alignment_mode', 'unknown')] += 1
    logger.info(f"      Selected by mode: {dict(mode_counts)}")
    
    return consensus_transcripts


def _get_ref_exon_count_map(args):
    """Return a dict mapping reference TranscriptId -> exon (block) count.

    Lazily parsed from args.ref_gp on first call and cached on the args object.
    Used to detect processed-pseudogene signatures (multi-exon source -> near-intronless target)
    during biotype assignment. Returns None if no reference genePred is configured.
    """
    cached = getattr(args, "_ref_exon_count_map_cache", None)
    if cached is not None:
        return cached
    ref_gp_path = getattr(args, "ref_gp", None)
    if not ref_gp_path:
        try:
            args._ref_exon_count_map_cache = {}
        except Exception:
            pass
        return {}
    counts = {}
    try:
        with open(ref_gp_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                tx_id = parts[0]
                try:
                    counts[tx_id] = int(parts[7])
                except (ValueError, IndexError):
                    continue
    except OSError:
        counts = {}
    try:
        args._ref_exon_count_map_cache = counts
    except Exception:
        pass
    return counts


def create_transcript_attributes(row, mode, normalized_id, tx_dict=None, consensus_dict=None, args=None, metrics=None, pc_gene_index=None, ref_df=None):
    """Create attribute dictionary for a consensus transcript"""
    import re
    aln_id = str(row['AlignmentId'])
    is_augpb = _is_denovo(aln_id)
    
    # Classify denovo transcripts
    transcript_class = 'ortholog'
    if is_augpb and tx_dict is not None and args is not None:
        denovo_mode = mode if mode in ('augPB', 'strg') else 'augPB'
        transcript_class = classify_augpb_transcript(row, tx_dict, consensus_dict if consensus_dict else {}, args, metrics, pc_gene_index=pc_gene_index, denovo_mode=denovo_mode)
    elif is_augpb:
        transcript_class = 'putative_novel'
    
    # Get the original gene ID from the tx_dict (from genePred name2 field)
    # Normalize to strip _N suffix for txTM CNV copies
    # For augPB, this preserves the unique gene ID (e.g., augPB-1, augPB-2, etc.)
    tx_obj = tx_dict.get(aln_id)
    raw_original_gene_id = tx_obj.name2 if tx_obj else row.get('GeneId', 'N/A')
    # For _cp copies, keep the _cp suffix in original_gene_id to preserve the paralog identity
    # but strip it from source_gene (see below) to point to the reference gene
    original_gene_id = normalize_gene_id(raw_original_gene_id, mode)
    
    # Determine source_gene and source_gene_biotype for augPB transcripts
    # For putative_novel_isoform, use the AssignedGeneId (the known gene it's an isoform of)
    # If no AssignedGeneId, find the overlapping protein-coding gene from consensus_dict
    # For other augPB classes, source_gene is None
    if is_augpb and transcript_class == 'putative_novel_isoform':
        source_gene = None
        source_gene_biotype = None
        # CRITICAL: Always find overlapping gene first (more reliable than AssignedGeneId)
        # This ensures we use the correct gene (MACF1) instead of wrong AssignedGeneId (POTEM)
        # IMPORTANT: Accept ALL biotypes (lncRNA, protein_coding, etc.) to preserve txTM biotypes
        # Try pc_gene_index first (faster), fall back to consensus_dict
        if tx_obj:
            if pc_gene_index is not None:
                # Fast lookup using spatial index
                chrom_txs = pc_gene_index.get(tx_obj.chromosome, [])
                for start, end, consensus_aln_id, consensus_tx, consensus_attrs in chrom_txs:
                    if start > tx_obj.stop:
                        break  # No more overlaps possible
                    # CRITICAL: Only consider same-strand overlaps for novel isoforms
                    if tx_obj.start <= end and tx_obj.stop >= start and consensus_tx.strand == tx_obj.strand:
                        # Get source_gene from consensus attrs (lowercase keys)
                        source_gene = consensus_attrs.get('source_gene', None)
                        if source_gene:
                            source_gene = normalize_gene_id(source_gene, mode)
                            # Also get the biotype to preserve it
                            source_gene_biotype = consensus_attrs.get('gene_biotype', None)
                            break
            elif consensus_dict:
                # Fallback: slow full scan
                for consensus_aln_id, consensus_attrs in consensus_dict.items():
                    if consensus_aln_id in tx_dict:
                        consensus_tx = tx_dict[consensus_aln_id]
                        # Check if same chromosome and overlapping
                        # CRITICAL: Only consider same-strand overlaps for novel isoforms
                        if (consensus_tx.chromosome == tx_obj.chromosome and
                            consensus_tx.strand == tx_obj.strand and
                            tx_obj.start <= consensus_tx.stop and tx_obj.stop >= consensus_tx.start):
                            # Accept ANY biotype (lncRNA, protein_coding, etc.)
                            source_gene = consensus_attrs.get('source_gene')
                            source_gene_biotype = consensus_attrs.get('gene_biotype')
                            if source_gene:
                                break
        
        # If no overlapping gene found, fall back to AssignedGeneId (shouldn't happen for putative_novel_isoform)
        if not source_gene:
            assigned_gene = row.get('AssignedGeneId', None)
            if assigned_gene:
                source_gene = normalize_gene_id(assigned_gene, 'txTM')
                # Try to get biotype from reference if needed
                if not source_gene_biotype:
                    source_gene_biotype = row.get('SourceGeneBiotype', None)
        # If we still don't have a biotype and have a source_gene, try to find it in consensus_dict
        if not source_gene_biotype and source_gene and consensus_dict:
            # Try to get biotype from overlapping gene
            for consensus_aln_id, consensus_attrs in consensus_dict.items():
                if consensus_attrs.get('source_gene') == source_gene:
                    source_gene_biotype = consensus_attrs.get('gene_biotype')
                    break
        # As last resort, try SourceGeneBiotype from row (reference lookup)
        if not source_gene_biotype:
            source_gene_biotype = row.get('SourceGeneBiotype', None)
    elif is_augpb:
        source_gene = None
        source_gene_biotype = None
    else:
        raw_source_gene = row.get('GeneId', 'N/A')
        # CRITICAL: Keep _cp suffix for paralogs - it's part of their unique identity
        # Each paralog at a different locus should have a distinct source_gene
        source_gene = normalize_gene_id(raw_source_gene, mode)
        # Also try normalizing as txTM to strip _N suffixes (for cases where gene IDs
        # have txTM-style suffixes regardless of transcript mode)
        if source_gene != normalize_gene_id(raw_source_gene, 'txTM'):
            source_gene = normalize_gene_id(raw_source_gene, 'txTM')
        source_gene_biotype = row.get('GeneBiotype', None)
    
    # Strip _cp suffix from source_transcript (added for multi-copy transcripts but not in actual source)
    source_transcript_raw = row.get('TranscriptId', 'N/A')
    source_transcript_name_raw = row.get('TranscriptName', row.get('TranscriptId', 'N/A'))
    
    # Ensure we have strings, not NaN or None
    if pd.isna(source_transcript_raw):
        source_transcript_raw = 'N/A'
    else:
        source_transcript_raw = str(source_transcript_raw)
    
    if pd.isna(source_transcript_name_raw):
        source_transcript_name_raw = source_transcript_raw
    else:
        source_transcript_name_raw = str(source_transcript_name_raw)
    
    if not is_augpb and source_transcript_raw != 'N/A':
        # Remove _cp, _cp2, _cp3, etc. suffixes that were added for multi-copy handling
        source_transcript_raw = re.sub(r'_cp\d*$', '', source_transcript_raw)
        source_transcript_name_raw = re.sub(r'_cp\d*$', '', source_transcript_name_raw)
    
    # Strip _cp suffix from source_gene as well (added internally for tracking, not in source)
    if source_gene and isinstance(source_gene, str):
        source_gene = re.sub(r'_cp\d*$', '', source_gene)
    
    # Strip _cp suffix from alignment_id and normalized_transcript_id for clean output
    # (these suffixes are added internally for cross-mode duplicate tracking)
    clean_alignment_id = re.sub(r'_cp\d*$', '', aln_id)
    clean_normalized_id = re.sub(r'_cp\d*$', '', normalized_id)

    # Strip _cp suffix from original_gene_id as well
    if original_gene_id and isinstance(original_gene_id, str):
        original_gene_id = re.sub(r'_cp\d*$', '', original_gene_id)

    # Biotype lift-over policy for cross-species annotation transfer:
    # default to the reference biotype, only downgrade in two narrow cases that
    # are unambiguous signals regardless of evolutionary distance.
    #   1. Processed pseudogene: a multi-exon source transcript mapping to a
    #      near-intronless target. Classic retrotransposed-mRNA signature.
    #   2. Fragment: both alignment coverage AND identity below 30%. Below this
    #      joint threshold we cannot meaningfully claim ortholog status.
    # Frameshifts / improper ORFs alone are NOT used to demote, since real
    # cross-species orthologs frequently have indels in their alignments.
    if not is_augpb and transcript_class == 'ortholog':
        ref_exon_count_map = _get_ref_exon_count_map(args) if args is not None else None

        src_tx_id = source_transcript_raw
        src_exon_count = None
        if ref_exon_count_map and src_tx_id and src_tx_id != 'N/A':
            src_exon_count = ref_exon_count_map.get(src_tx_id)
            if src_exon_count is None:
                base = re.sub(r"\.[0-9]+$", "", src_tx_id)
                src_exon_count = ref_exon_count_map.get(base)
        tgt_exon_count = tx_obj.block_count if tx_obj is not None else None

        tx_biotype = row.get('TranscriptBiotype', None)
        gene_biotype = row.get('GeneBiotype', None)
        src_is_pc = (
            tx_biotype in ('protein_coding', 'mRNA')
            or gene_biotype == 'protein_coding'
        )

        cov = row.get('AlnCoverage_mRNA', None)
        ident = row.get('AlnIdentity_mRNA', None)
        cov_val = float(cov) if pd.notna(cov) else None
        ident_val = float(ident) if pd.notna(ident) else None

        is_processed_pseudo = (
            src_is_pc
            and src_exon_count is not None and src_exon_count >= 3
            and tgt_exon_count is not None and tgt_exon_count <= 1
        )
        ref_pc_ensg = getattr(args, 'ref_pc_ensg', None) if args is not None else None
        if is_processed_pseudo and ref_pc_ensg:
            gid = row.get('GeneId')
            if gid and norm_ensg(gid) in ref_pc_ensg:
                is_processed_pseudo = False
        frag_max_cov = getattr(args, 'fragment_max_coverage', 30.0) if args is not None else 30.0
        frag_max_id = getattr(args, 'fragment_max_identity', 30.0) if args is not None else 30.0
        # A threshold of 0 disables fragment reclassification (nothing scores < 0).
        is_fragment_call = (
            cov_val is not None and ident_val is not None
            and cov_val < frag_max_cov and ident_val < frag_max_id
        )

        if is_processed_pseudo:
            transcript_class = 'processed_pseudogene'
        elif is_fragment_call:
            transcript_class = 'fragment'

    attrs = {
        'source_transcript': source_transcript_raw if not is_augpb else 'N/A',
        'source_transcript_name': source_transcript_name_raw,
        'source_gene': source_gene,  # For augPB novel isoforms, this is the known gene they're an isoform of
        'source_gene_biotype': source_gene_biotype,  # Biotype of the source gene
        'original_gene_id': original_gene_id,  # Original gene ID (cleaned of _cp suffixes)
        'score': int(10 * round(row.get('AlnGoodness_mRNA', 1.0), 3)) if pd.notna(row.get('AlnGoodness_mRNA')) else 100,
        'gene_biotype': (
            # fragment / processed_pseudogene overrides the reference biotype for non-denovo transcripts
            transcript_class if transcript_class in ('fragment', 'processed_pseudogene') else
            # For augPB novel isoforms, inherit the source gene's biotype (preserves lncRNA, protein_coding, etc.)
            (source_gene_biotype if (is_augpb and transcript_class == 'putative_novel_isoform' and source_gene_biotype) else
            # For augPB novel genes, use unknown_likely_coding
            ('unknown_likely_coding' if is_augpb else
             ('unknown' if (pd.isna(row.get('GeneBiotype')) or row.get('GeneBiotype') in ['N/A', 'unknown', None, '']) else row.get('GeneBiotype'))))
        ),
        'transcript_biotype': (
            # fragment / processed_pseudogene overrides the reference biotype for non-denovo transcripts
            transcript_class if transcript_class in ('fragment', 'processed_pseudogene') else
            # For augPB novel isoforms, inherit the source gene's biotype (preserves lncRNA, protein_coding, etc.)
            (source_gene_biotype if (is_augpb and transcript_class == 'putative_novel_isoform' and source_gene_biotype) else
            # For augPB novel genes, use unknown_likely_coding
            ('unknown_likely_coding' if is_augpb else
             ('unknown' if (pd.isna(row.get('TranscriptBiotype')) or row.get('TranscriptBiotype') in ['N/A', 'unknown', None, '']) else row.get('TranscriptBiotype'))))
        ),
        'alignment_id': aln_id,  # Keep internal ID with _cp for source lookup
        'alignment_mode': mode,
        'normalized_transcript_id': normalized_id,  # Keep internal ID with _cp
        'transcript_score': row.get('TranscriptScore', 0),
        'frameshift': str(row.get('Frameshift', 'N/A')) if pd.notna(row.get('Frameshift')) else 'N/A',
        'transcript_class': transcript_class,
        'valid_start': bool(row.get('ValidStart', False)),
        'valid_stop': bool(row.get('ValidStop', False)),
        'proper_orf': bool(row.get('ProperOrf', False)),
        'exon_annotation_support': ','.join(map(str, row.get('ExonAnnotSupport', []))),
        'intron_annotation_support': ','.join(map(str, row.get('IntronAnnotSupport', []))),
        'exon_rna_support': ','.join(map(str, row.get('ExonRnaSupport', []))),
        'intron_rna_support': ','.join(map(str, row.get('IntronRnaSupport', []))),
        'source_gene_common_name': row.get('GeneName', None) if not is_augpb else None
    }
    
    # For augPB putative_novel_isoform, look up gene name from reference if we have source_gene
    if is_augpb and transcript_class == 'putative_novel_isoform' and attrs.get('source_gene') and ref_df is not None:
        source_gene = attrs.get('source_gene')
        # Look up gene name from reference database
        if 'GeneId' in ref_df.columns and 'GeneName' in ref_df.columns:
            ref_gene_match = ref_df[ref_df['GeneId'] == source_gene]
            if len(ref_gene_match) > 0:
                gene_name = ref_gene_match.iloc[0]['GeneName']
                if pd.notna(gene_name):
                    attrs['source_gene_common_name'] = gene_name
    
    # Add novel flags for augPB transcripts
    if is_augpb:
        attrs['novel_5p_cap'] = True
        attrs['novel_poly_a'] = True
        attrs['pacbio_isoform_supported'] = True
    
    return attrs


def _build_pc_gene_index(consensus_dict, tx_dict):
    """
    Build a spatial index of protein-coding transcripts by chromosome.
    Returns: dict[chromosome] = list of (start, end, aln_id, tx_obj, attrs) sorted by start
    """
    index = collections.defaultdict(list)
    for consensus_aln_id, consensus_attrs in consensus_dict.items():
        gene_biotype = consensus_attrs.get('gene_biotype', 'unknown')
        if gene_biotype == 'protein_coding' and consensus_aln_id in tx_dict:
            tx = tx_dict[consensus_aln_id]
            index[tx.chromosome].append((tx.start, tx.stop, consensus_aln_id, tx, consensus_attrs))
    
    # Sort by start position for efficient overlap checking
    for chrom in index:
        index[chrom].sort(key=lambda x: x[0])
    
    return index


def classify_augpb_transcript(row, tx_dict, consensus_dict, args, metrics=None, pc_gene_index=None, all_tx_dict=None, denovo_mode='augPB'):
    """
    Classify de novo (augPB/strg) transcripts into categories:
    - putative_novel: truly novel gene (no overlapping genes from other methods)
    - putative_novel_isoform: novel isoform of known gene (overlaps with gene from other methods, different exon structure)
    - possible_paralog: potential gene family expansion
    - poor_alignment: has some annotation support
    
    Args:
        pc_gene_index: Pre-built spatial index (all transcript biotypes, not just protein-coding)
        all_tx_dict: Dictionary of ALL input transcripts (to check for txTM/transMap/augTM at this locus)
    """
    aln_id = row['AlignmentId']
    tx_obj = tx_dict.get(aln_id)
    if not tx_obj:
        return None
    
    num_introns = len(tx_obj.exon_intervals) - 1

    # Require at least denovo_num_introns spliced introns (default 1 → ≥2 exons) for
    # all de novo modes, including StringTie (strg).
    if num_introns < args.denovo_num_introns:
        if metrics:
            metrics['denovo'][denovo_mode]['Discarded'] += 1
        return None

    # Get support percentages
    intron_rna_support = row.get('IntronRnaSupportPercent', 0)
    exon_rna_support = row.get('ExonRnaSupportPercent', 0)
    intron_annot_support = row.get('IntronAnnotSupportPercent', 0)
    exon_annot_support = row.get('ExonAnnotSupportPercent', 0)
    cds_annot_support = row.get('CdsAnnotSupportPercent', 0)
    
    # Check support thresholds
    if args.in_species_rna_support_only:
        rna_intron = intron_rna_support
        rna_exon = exon_rna_support
    else:
        rna_intron = row.get('AllSpeciesIntronRnaSupportPercent', intron_rna_support)
        rna_exon = row.get('AllSpeciesExonRnaSupportPercent', exon_rna_support)
    
    # Check minimum support thresholds
    if rna_intron < args.denovo_splice_support or rna_exon < args.denovo_exon_support:
        if not args.denovo_allow_unsupported:
            if metrics:
                metrics['denovo'][denovo_mode]['Discarded'] += 1
            return None
    
    # Check if this has alternative gene IDs (gene family/paralog)
    has_alternatives = pd.notna(row.get('AlternativeGeneIds')) and row.get('AlternativeGeneIds') not in ['', 'None']
    if has_alternatives:
        row['AssignedGeneId'] = None
        if metrics:
            metrics['denovo'][denovo_mode]['Possible paralog'] += 1
            metrics['AugPB Classes']['possible_paralog'] += 1
        return 'possible_paralog'
    
    # NEW LOGIC: Check if there are overlapping genes from other methods in the region
    overlapping_genes = []
    assigned_gene = row.get('AssignedGeneId')

    # Bound the number of overlapping models collected per de novo transcript.
    # Normal loci have only a handful of isoforms, but in pathological, highly
    # duplicated regions (e.g. chr1 NBPF / segmental-duplication clusters) a single
    # denovo transcript can span thousands of overlapping models. Collecting them all
    # (and unioning their introns in check_novel_splices) makes this O(N^2) per locus
    # and can hang for hours. The classification outcome only depends on whether a
    # known (non-denovo) overlap exists and on nearby splice structure, both of which
    # are captured well within this generous cap, so truncating dense loci does not
    # change results for real genes.
    MAX_OVERLAPPING_GENES = 500

    # Use spatial index if available, otherwise fall back to full scan (slower)
    if pc_gene_index is not None:
        # Fast lookup using spatial index
        chrom_txs = pc_gene_index.get(tx_obj.chromosome, [])
        # Binary search for transcripts that might overlap
        # Find transcripts that start before our end and end after our start
        for start, end, consensus_aln_id, consensus_tx, consensus_attrs in chrom_txs:
            if start > tx_obj.stop:
                break  # No more overlaps possible (sorted by start)
            # CRITICAL: Only consider overlaps on the SAME STRAND
            # Opposite strand overlaps should be treated as novel genes
            if tx_obj.start <= end and tx_obj.stop >= start and consensus_tx.strand == tx_obj.strand:
                overlapping_genes.append((consensus_aln_id, consensus_tx, consensus_attrs))
                if len(overlapping_genes) >= MAX_OVERLAPPING_GENES:
                    break  # Dense locus: enough context collected to classify
    else:
        # Fallback: full scan (slow but works if index not provided)
        for consensus_aln_id, consensus_attrs in consensus_dict.items():
            if consensus_aln_id in tx_dict:
                consensus_tx = tx_dict[consensus_aln_id]
                # Check if same chromosome and overlapping (include all biotypes)
                # CRITICAL: Only consider overlaps on the SAME STRAND
                if (consensus_tx.chromosome == tx_obj.chromosome and
                    consensus_tx.strand == tx_obj.strand and
                    tx_obj.start <= consensus_tx.stop and tx_obj.stop >= consensus_tx.start):
                    overlapping_genes.append((consensus_aln_id, consensus_tx, consensus_attrs))
                    if len(overlapping_genes) >= MAX_OVERLAPPING_GENES:
                        break  # Dense locus: enough context collected to classify
    
    # If there are overlapping genes from other methods, check if exon structure differs
    if overlapping_genes:
        # Check if it has novel splices compared to overlapping genes
        # Build a temporary consensus_dict from overlapping genes for splice comparison
        # overlapping_genes format: [(consensus_aln_id, consensus_tx, consensus_attrs), ...]
        temp_consensus_dict = {consensus_aln_id: attrs for consensus_aln_id, consensus_tx, attrs in overlapping_genes}
        has_novel_splices = check_novel_splices(tx_obj, temp_consensus_dict, tx_dict)

        if has_novel_splices:
            # Truly novel isoform with NOVEL SPLICE SITES
            # CRITICAL: Must have overlapping gene from NON-augPB source
            # Check if overlapping genes are from known sources (not augPB)
            has_known_overlap = any(
                not _is_denovo(consensus_aln_id)
                for consensus_aln_id, consensus_tx, consensus_attrs in overlapping_genes
            )
            
            if not has_known_overlap:
                if metrics:
                    metrics['denovo'][denovo_mode]['Discarded'] += 1
                return None
            
            # Validate assigned_gene against overlapping genes - use overlapping gene if assigned_gene doesn't match
            if overlapping_genes:
                # Get source_gene from the FIRST NON-augPB overlapping gene
                overlapping_gene_source = None
                for consensus_aln_id, consensus_tx, consensus_attrs in overlapping_genes:
                    if not _is_denovo(consensus_aln_id):
                        overlapping_gene_source = consensus_attrs.get('source_gene') if isinstance(consensus_attrs, dict) else None
                        if overlapping_gene_source:
                            break
                
                if overlapping_gene_source:
                    # Normalize for comparison
                    overlapping_gene_source_norm = normalize_gene_id(overlapping_gene_source, 'txTM')
                    assigned_gene_norm = normalize_gene_id(assigned_gene, 'txTM') if pd.notna(assigned_gene) and assigned_gene not in ['', 'None'] else None
                    
                    # ALWAYS use the overlapping gene's source_gene (more reliable than AssignedGeneId)
                    row['AssignedGeneId'] = overlapping_gene_source
                    
                    if metrics:
                        metrics['denovo'][denovo_mode]['Putative novel isoform'] += 1
                        metrics['AugPB Classes']['putative_novel_isoform'] += 1
                    return 'putative_novel_isoform'
            
            if metrics:
                metrics['denovo'][denovo_mode]['Discarded'] += 1
            return None
        
        if metrics:
            metrics['denovo'][denovo_mode]['Discarded'] += 1
        return None
    
    # No overlapping genes from other methods - this is a novel gene
    # CRITICAL: Clear AssignedGeneId so novel genes get source_gene=None and unknown_likely_coding biotype
    row['AssignedGeneId'] = None
    
    # Check if it has annotation support (poor alignment)
    if intron_annot_support > 0 or exon_annot_support > 0 or cds_annot_support > 0:
        if metrics:
            metrics['denovo'][denovo_mode]['Poor alignment'] += 1
            metrics['AugPB Classes']['poor_alignment'] += 1
        return 'poor_alignment'
    
    if metrics:
        metrics['denovo'][denovo_mode]['Putative novel'] += 1
        metrics['AugPB Classes']['putative_novel'] += 1
    return 'putative_novel'


def check_novel_splices(augpb_tx_obj, consensus_dict, tx_dict):
    """Check if augPB transcript has novel splices compared to consensus"""
    augpb_introns = set(augpb_tx_obj.intron_intervals)
    
    if len(augpb_introns) == 0:
        return False  # Single-exon transcript
    
    # Get all consensus splices
    consensus_splices = set()
    for consensus_aln_id in consensus_dict.keys():
        if consensus_aln_id in tx_dict:
            consensus_tx = tx_dict[consensus_aln_id]
            consensus_splices.update(consensus_tx.intron_intervals)
    
    # If no consensus splices found (empty consensus_dict), cannot determine novelty
    # Return False to avoid false positives
    if len(consensus_splices) == 0:
        return False
    
    # Check if augPB has any novel introns
    novel_introns = augpb_introns - consensus_splices
    
    return len(novel_introns) > 0


def score_transcripts(df, in_species_rna_support_only, is_coding=True):
    """Score transcripts based on alignment quality and support"""
    if len(df) == 0:
        return df
    
    df = df.copy()
    
    if is_coding:
        # For coding transcripts
        aln_id = df.get('AlnIdentity_CDS', df.get('AlnIdentity_mRNA', pd.Series([0]*len(df)))).fillna(0)
        aln_cov = df.get('AlnCoverage_CDS', df.get('AlnCoverage_mRNA', pd.Series([0]*len(df)))).fillna(0)
        orig_intron = df.get('OriginalIntronsPercent_mRNA', pd.Series([0]*len(df))).fillna(0)
    else:
        # For non-coding transcripts
        aln_id = df.get('TransMapIdentity', df.get('AlnIdentity_mRNA', pd.Series([0]*len(df)))).fillna(0)
        aln_cov = df.get('TransMapCoverage', df.get('AlnCoverage_mRNA', pd.Series([0]*len(df)))).fillna(0)
        orig_intron = df.get('TransMapOriginalIntronsPercent', pd.Series([0]*len(df))).fillna(0)
    
    if in_species_rna_support_only:
        rna_support = df.get('ExonRnaSupportPercent', pd.Series([0]*len(df))).fillna(0) + \
                      df.get('IntronRnaSupportPercent', pd.Series([0]*len(df))).fillna(0)
    else:
        rna_support = df.get('AllSpeciesExonRnaSupportPercent', pd.Series([0]*len(df))).fillna(0) + \
                      df.get('AllSpeciesIntronRnaSupportPercent', pd.Series([0]*len(df))).fillna(0)
    
    intron_annot = df.get('IntronAnnotSupportPercent', pd.Series([0]*len(df))).fillna(0)
    exon_annot = df.get('ExonAnnotSupportPercent', pd.Series([0]*len(df))).fillna(0)
    
    df['TranscriptScore'] = aln_id + aln_cov + intron_annot + exon_annot + orig_intron + rna_support
    
    return df


def initialize_metrics():
    """Initialize metrics dictionary"""
    return {
        'Transcript Missing': collections.Counter(),
        'Gene Missing': collections.Counter(),
        'Transcript Modes': collections.Counter(),
        'Duplicate transcripts': collections.Counter(),
        'Discarded by strand resolution': 0,
        'Multi-locus mappings': collections.Counter(),
        'Multi-locus kept': collections.Counter(),
        'Coverage': collections.defaultdict(list),
        'Identity': collections.defaultdict(list),
        'Splice Support': collections.defaultdict(list),
        'Exon Support': collections.defaultdict(list),
        'Original Introns': collections.defaultdict(list),
        'Splice Annotation Support': collections.defaultdict(list),
        'Exon Annotation Support': collections.defaultdict(list),
        'AugPB Classes': collections.Counter(),
        'denovo': {
            'augPB': {
                'Possible paralog': 0,
                'Poor alignment': 0,
                'Putative novel': 0,
                'Putative novel isoform': 0,
                'Discarded': 0
            },
            'strg': {
                'Possible paralog': 0,
                'Poor alignment': 0,
                'Putative novel': 0,
                'Putative novel isoform': 0,
                'Discarded': 0
            }
        }
    }


# Import helper functions from original consensus.py
def load_transmap_evals(db_path):
    tm_eval = tools.sqlInterface.load_alignment_evaluation(db_path)
    tm_filter_eval = tools.sqlInterface.load_filter_evaluation(db_path)
    tm_pw_filter_eval = tools.sqlInterface.load_pairwise_filter_evaluation(db_path)
    
    # Merge regular transMap filter eval with evaluation metrics
    tm_eval_df = pd.merge(tm_eval, tm_filter_eval, on=['TranscriptId', 'AlignmentId'], how='outer')
    
    # Also merge pairwise transMap filter eval if it exists
    if len(tm_pw_filter_eval) > 0:
        # Merge pairwise filter eval with evaluation metrics (may not have metrics for all pairwise alignments)
        tm_pw_merged = pd.merge(tm_eval, tm_pw_filter_eval, on=['TranscriptId', 'AlignmentId'], how='right')
        # Combine both regular and pairwise, keeping all unique combinations
        tm_eval_df = pd.concat([tm_eval_df, tm_pw_merged], ignore_index=True).drop_duplicates(subset=['TranscriptId', 'AlignmentId'], keep='first')
    
    return tm_eval_df.drop('AlignmentId', axis=1)


def load_metrics_from_db(db_path, tx_mode, aln_mode):
    session = tools.sqlInterface.start_session(db_path)
    metrics_table = tools.sqlInterface.tables[aln_mode][tx_mode]['metrics']
    metrics_df = tools.sqlInterface.load_metrics(metrics_table, session)
    # Empty metrics table: txTM/transMap_pairwise/augTM_pairwise/etc.
    # legitimately produce zero rows on small or pairwise-disabled genomes.
    # Returning an empty DataFrame with the expected columns lets the
    # consensus runner concat it harmlessly downstream instead of KeyError'ing
    # on the unstack/numeric-cast block below.
    if len(metrics_df) == 0:
        session.close()
        return pd.DataFrame(columns=[
            'AlignmentId', 'AlnCoverage', 'AlnGoodness', 'AlnIdentity',
            'PercentUnknownBases', 'AdjStart', 'AdjStop', 'OriginalIntrons',
            'OriginalIntronsPercent', 'ProperOrf', 'ValidStart', 'ValidStop',
        ])
    # Some DBs can contain duplicate rows for the same (AlignmentId, classifier).
    # Pandas can't unstack with duplicates, so de-duplicate deterministically.
    if metrics_df.duplicated(subset=['AlignmentId', 'classifier']).any():
        dup_n = int(metrics_df.duplicated(subset=['AlignmentId', 'classifier']).sum())
        logger.warning(
            f"Duplicate metrics rows detected for {tx_mode}/{aln_mode}: "
            f"{dup_n} duplicates on (AlignmentId, classifier). Keeping last occurrence."
        )
        metrics_df = metrics_df.drop_duplicates(subset=['AlignmentId', 'classifier'], keep='last')
    metrics_df = metrics_df.set_index(['AlignmentId', 'classifier']).unstack('classifier')
    metrics_df.columns = [col[1] for col in metrics_df.columns]
    metrics_df = metrics_df.reset_index()
    # augMP PSL metrics (store_psl_metrics.py) only persist AlnCoverage/AlnIdentity.
    # Other modes may also omit classifiers when a table is partially populated.
    numeric_cols = ['AlnCoverage', 'AlnGoodness', 'AlnIdentity', 'PercentUnknownBases']
    optional_cols = numeric_cols + [
        'AdjStart', 'AdjStop', 'OriginalIntrons', 'ProperOrf', 'ValidStart', 'ValidStop',
    ]
    for col in optional_cols:
        if col not in metrics_df.columns:
            metrics_df[col] = pd.NA
    metrics_df[numeric_cols] = metrics_df[numeric_cols].apply(pd.to_numeric, errors='coerce')
    metrics_df['OriginalIntrons'] = metrics_df['OriginalIntrons'].fillna('').astype(str)
    metrics_df['OriginalIntrons'] = [list(map(int, x)) if len(x[0]) > 0 else []
                                     for x in metrics_df['OriginalIntrons'].str.split(',').tolist()]
    metrics_df['OriginalIntronsPercent'] = metrics_df['OriginalIntrons'].apply(
        lambda s: 100 * len([x for x in s if x > 0]) / len(s) if len(s) > 0 else 100
    )
    session.close()
    return metrics_df


def load_evaluations_from_db(db_path, tx_mode):
    def aggfunc(s):
        if s.value_CDS.any():
            c = set(s[s.value_CDS > 0].name)
        else:
            c = set(s[s.value_mRNA > 0].name)
        cols = ['Frameshift', 'CodingInsertion', 'CodingDeletion', 'CodingMult3Indel']
        return pd.Series((
            'CodingDeletion' in c or 'CodingInsertion' in c,
            'CodingInsertion' in c, 
            'CodingDeletion' in c,
            'CodingMult3Deletion' in c or 'CodingMult3Insertion' in c
        ), index=cols)
    
    session = tools.sqlInterface.start_session(db_path)
    cds_table = tools.sqlInterface.tables['CDS'][tx_mode]['evaluation']
    mrna_table = tools.sqlInterface.tables['mRNA'][tx_mode]['evaluation']
    cds_df = tools.sqlInterface.load_evaluation(cds_table, session)
    mrna_df = tools.sqlInterface.load_evaluation(mrna_table, session)
    cds_df = cds_df.set_index('AlignmentId')
    mrna_df = mrna_df.set_index('AlignmentId')
    merged = mrna_df.reset_index().merge(cds_df.reset_index(), how='outer', 
                                         on=['AlignmentId', 'name'], suffixes=['_mRNA', '_CDS'])
    # Pandas compatibility: include_groups was added in newer pandas.
    # Older versions will raise TypeError, so fall back to apply() without it.
    try:
        eval_df = merged.groupby('AlignmentId').apply(aggfunc, include_groups=False)
    except TypeError:
        eval_df = merged.groupby('AlignmentId').apply(aggfunc)
    return eval_df


def deduplicate_consensus(consensus_dict, tx_dict, metrics):
    """Remove duplicate transcripts (same exon structure)"""
    import re
    
    def get_source_priority(aln_id):
        """Return priority for source types (lower is better)"""
        if _is_denovo(aln_id):
            return 100  # Lowest priority - novel predictions
        elif 'transMap' in consensus_dict[aln_id].get('alignment_mode', ''):
            return 1  # Prefer full-length transMap over fragmented txTM
        elif aln_id.startswith('txTM-') or 'txTM' in consensus_dict[aln_id].get('alignment_mode', ''):
            return 2
        elif aln_id.startswith('augTM-'):
            return 3  # Augustus with transMap hints
        elif aln_id.startswith('augTMR-'):
            return 4
        else:
            return 5
    
    def select_best_duplicate(tx_list):
        """Select best transcript from duplicates, preferring known genes over denovo
        BUT keep denovo transcripts if they are novel isoforms"""
        # Separate known genes from denovo
        known_genes = [tx for tx in tx_list if not _is_denovo(tx)]
        augpb_genes = [tx for tx in tx_list if _is_denovo(tx)]
        
        # Check if any denovo transcripts are novel isoforms (should be kept!)
        augpb_novel_isoforms = [
            tx for tx in augpb_genes 
            if consensus_dict[tx].get('transcript_class') == 'putative_novel_isoform'
        ]
        
        # If we have augPB novel isoforms AND known genes, keep BOTH
        # Novel isoforms represent validated PacBio structures that should coexist with known genes
        if augpb_novel_isoforms and known_genes:
            # Return all novel isoforms + best known gene
            # We'll return as list and handle below
            best_known = max(known_genes, key=lambda x: (
                -get_source_priority(x),
                consensus_dict[x].get('score', 0)
            ))
            return [best_known] + augpb_novel_isoforms
        
        # If we have known genes but no novel isoforms, only consider known genes (ignore augPB)
        if known_genes:
            candidates = known_genes
        else:
            candidates = tx_list
        
        # Among candidates, pick by priority then score
        best_tx = max(candidates, key=lambda x: (
            -get_source_priority(x),  # Negative because lower priority number is better
            consensus_dict[x].get('score', 0)
        ))
        
        return best_tx
    
    duplicates = collections.defaultdict(list)
    
    for aln_id in consensus_dict:
        tx = tx_dict[aln_id]
        # Include chromosome in deduplication key to avoid removing genes on different chromosomes
        # (critical for cases like HG03456 with 2 Y chromosomes)
        location_key = (tx.chromosome, tx.start, tx.stop)
        interval_key = (frozenset(tx.exon_intervals), location_key)
        duplicates[interval_key].append(aln_id)

    deduplicated_consensus = {}
    for tx_list in duplicates.values():
        if len(tx_list) > 1:
            metrics['Duplicate transcripts'][len(tx_list)] += 1

            # Pick best using priority + score
            best_result = select_best_duplicate(tx_list)

            # Handle case where multiple transcripts are kept (novel isoforms + known gene)
            if isinstance(best_result, list):
                kept_txs = best_result
                for kept_tx in kept_txs:
                    deduplicated_consensus[kept_tx] = consensus_dict[kept_tx]
                # Store alternatives on the first transcript (strip _cp suffix from multi-copy handling)
                alt_txs = [re.sub(r'_cp\d*$', '', tools.nameConversions.strip_alignment_numbers(x))
                          for x in tx_list if x not in kept_txs]
                if alt_txs:
                    deduplicated_consensus[kept_txs[0]]['alternative_source_transcripts'] = ','.join(alt_txs)
            else:
                # Single best transcript
                best_tx = best_result
                deduplicated_consensus[best_tx] = consensus_dict[best_tx]
                # Strip _cp suffix from alternative transcripts (from multi-copy handling)
                alt_txs = [re.sub(r'_cp\d*$', '', tools.nameConversions.strip_alignment_numbers(x)) 
                          for x in tx_list if x != best_tx]
                if alt_txs:
                    deduplicated_consensus[best_tx]['alternative_source_transcripts'] = ','.join(alt_txs)
        else:
            tx_id = tx_list[0]
            deduplicated_consensus[tx_id] = consensus_dict[tx_id]

    return deduplicated_consensus


def filter_gene_fragments(deduplicated_consensus, tx_dict, ref_gene_coords, metrics, min_length_ratio=0.7):
    """
    Remove genes with invalid source_gene IDs (not in reference).
    This catches cases where source_gene is a transcript ID or other invalid data.
    
    Args:
        deduplicated_consensus: Dictionary of {tx_id: attrs}
        tx_dict: Dictionary of transcript objects
        ref_gene_coords: Dict mapping gene IDs to (chrom, start, end) in reference genome
        metrics: Metrics dictionary
        min_length_ratio: Ignored (kept for compatibility)
    
    Returns:
        Dictionary of {tx_id: attrs} with invalid genes removed
    """
    # Group transcripts by source_gene and chromosome
    genes_by_source = collections.defaultdict(list)
    for tx_id, attrs in deduplicated_consensus.items():
        tx_obj = tx_dict[tx_id]
        source_gene = attrs.get('source_gene', 'N/A')
        alignment_mode = attrs.get('alignment_mode', 'unknown')
        if source_gene and source_gene != 'N/A':
            genes_by_source[(source_gene, tx_obj.chromosome, alignment_mode)].append((tx_obj, tx_id, attrs))
    
    genes_to_remove = set()
    invalid_source_genes = set()
    
    for (source_gene, chrom, mode), tx_list in genes_by_source.items():
        # Skip augMP - user requested exemption
        if mode == 'augMP':
            continue
        
        # Skip denovo modes - they can have novel genes not in reference
        if mode in ('augPB', 'strg'):
            continue
        
        # Skip txTM - can have species-specific genes not in reference
        if mode == 'txTM':
            continue
        
        # Check if source_gene is a valid gene ID in reference
        # Also check if it looks like a transcript ID (starts with ENST)
        is_transcript_id = source_gene.startswith('ENST')
        is_valid_gene = source_gene in ref_gene_coords
        
        if not is_valid_gene or is_transcript_id:
            # Invalid source_gene (likely a transcript ID or bad data) - remove for non-augPB/augMP
            source_name = tx_list[0][2].get('source_gene_common_name', source_gene)
            reason = "transcript ID" if is_transcript_id else "not found in reference"
            logger.info(f"  Removing gene with invalid source_gene: {source_name} ({source_gene}) on {chrom}")
            logger.info(f"    Mode: {mode}, Transcripts: {len(tx_list)} - {reason}")
            
            for _, tx_id, _ in tx_list:
                genes_to_remove.add(tx_id)
                invalid_source_genes.add(source_gene)
    
    metrics['Discarded gene fragments'] = len(genes_to_remove)
    if invalid_source_genes:
        logger.info(f"  Removed {len(invalid_source_genes)} genes with invalid source_gene IDs")
    
    # Return filtered consensus
    return {tx_id: attrs for tx_id, attrs in deduplicated_consensus.items() if tx_id not in genes_to_remove}


def resolve_conflicting_source_genes(
    deduplicated_consensus,
    tx_dict,
    metrics,
    overlap_threshold=0.5,
    readthrough_gene_set=None,
    ref_gene_coords=None,
    genes_with_overlaps_in_ref=None,
    skip_spurious_reference_overlap_removals=False,
    disregard_long_mode_ratio=2.0,
    disregard_long_mode_min_bp=50_000,
    spurious_overlap_min_assembly_overlap_bp=2000,
    spurious_overlap_min_reciprocal=0.02,
    ref_pc_ensg=None,
):
    """
    Resolve cases where different source genes map to overlapping loci on the target assembly.
    Keeps one protein-coding gene per overlapping assembly interval unless the specific gene
    pair also overlaps in the reference genome (nested genes / gene families). Genes that
    overlap another gene in the reference are exempt from multi-way collapse.
    
    Args:
        deduplicated_consensus: Dictionary of {tx_id: attrs}
        tx_dict: Dictionary of transcript objects
        metrics: Metrics dictionary
        overlap_threshold: Minimum reciprocal overlap to consider genes conflicting (default 0.5)
        readthrough_gene_set: Set of gene IDs that are readthrough genes (should coexist with overlaps)
        ref_gene_coords: Dict mapping gene IDs to (chrom, start, end) in reference genome
        genes_with_overlaps_in_ref: Set of gene IDs that overlap with ANY gene in reference
        skip_spurious_reference_overlap_removals: If True, skip N-way / 2-way removals when genes overlap
            on the target but not in reference (normally treated as collapsed paralogs). Default False.
            augPB vs known, reciprocal overlap resolution, and other rules still apply.
    
    Returns:
        Dictionary of {tx_id: attrs} with conflicts resolved
    """
    import intervaltree
    import re
    
    # Group transcripts by chromosome and source_gene
    # Also track genes WITHOUT source_gene separately to check against overlapping families
    genes_by_chrom = collections.defaultdict(lambda: collections.defaultdict(list))
    genes_without_source = collections.defaultdict(lambda: collections.defaultdict(list))  # chrom -> unique_id -> tx_list
    
    for tx_id, attrs in deduplicated_consensus.items():
        tx_obj = tx_dict[tx_id]
        source_gene = attrs.get('source_gene', 'N/A')
        if source_gene and source_gene != 'N/A':
            # Strip _cp suffix for grouping so that genes from different modes
            # with the same base ID are grouped together for conflict resolution
            source_gene_key = re.sub(r'_cp\d+$', '', source_gene)
            # Also strip txTM CNV suffixes (e.g., ENSG00000026103.25_21 → ENSG00000026103.25)
            # Match .<version>_<CNV> pattern to avoid stripping gene names ending in numbers
            source_gene_key = re.sub(r'(\.\d+)_\d+$', r'\1', source_gene_key)
            genes_by_chrom[tx_obj.chromosome][source_gene_key].append((tx_obj, tx_id, attrs))
        else:
            # Track genes without source_gene by their assigned gene ID
            assigned_gene_id = attrs.get('original_gene_id', tx_id)
            genes_without_source[tx_obj.chromosome][assigned_gene_id].append((tx_obj, tx_id, attrs))
    
    # For each chromosome, find overlapping genes from different sources
    genes_to_remove = set()
    conflicts_resolved = 0
    
    # Track competing sources: tx_id -> set of competing modes
    competing_sources = collections.defaultdict(set)
    
    # FIRST: Remove augPB transcripts that overlap with known genes (TxTM/transMap/augTM)
    # UNLESS they are putative_novel_isoform with novel splice sites
    # This ensures de novo transcripts only supplement existing genes, never replace them
    for chrom in genes_by_chrom:
        # Build interval tree of NON-denovo genes (known genes from txTM/transMap/augTM)
        known_gene_tree = intervaltree.IntervalTree()
        augpb_genes_to_check = []  # List of (source_gene, tx_list) for denovo genes
        
        for source_gene, tx_list in genes_by_chrom[chrom].items():
            # Check if ANY transcript in this gene is denovo
            has_denovo = any(_is_denovo(tx_id) for _, tx_id, _ in tx_list)
            has_non_denovo = any(not _is_denovo(tx_id) for _, tx_id, _ in tx_list)
            
            if has_denovo and not has_non_denovo:
                # Pure denovo gene - check it later
                augpb_genes_to_check.append((source_gene, tx_list))
            elif not has_denovo:
                # Pure non-denovo gene - add to tree
                starts = [tx_obj.start for tx_obj, _, _ in tx_list]
                ends = [tx_obj.stop for tx_obj, _, _ in tx_list]
                gene_start = min(starts)
                gene_end = max(ends)
                strand = tx_list[0][0].strand
                known_gene_tree.addi(gene_start, gene_end, (source_gene, strand))
            # Mixed denovo/non-denovo genes are kept (novel isoforms successfully merged)
        
        # Check each denovo gene against known genes
        for source_gene, tx_list in augpb_genes_to_check:
            # Get the span of this denovo gene
            starts = [tx_obj.start for tx_obj, _, _ in tx_list]
            ends = [tx_obj.stop for tx_obj, _, _ in tx_list]
            gene_start = min(starts)
            gene_end = max(ends)
            strand = tx_list[0][0].strand
            
            # Check if it overlaps with any known gene on the same strand
            overlapping_known_genes = known_gene_tree.overlap(gene_start, gene_end)
            overlapping_same_strand = [
                (source_gene_known, gene_strand) for (source_gene_known, gene_strand) 
                in (item.data for item in overlapping_known_genes)
                if gene_strand == strand
            ]
            
            if overlapping_same_strand:
                # denovo gene overlaps with known gene(s) on same strand
                # ONLY keep if ALL transcripts are putative_novel_isoform
                transcript_classes = set(
                    attrs.get('transcript_class', '')
                    for _, _, attrs in tx_list
                )
                
                # Keep ONLY if ALL transcripts are novel isoforms with the SAME source gene
                all_novel_isoforms = all(
                    attrs.get('transcript_class') == 'putative_novel_isoform'
                    for _, _, attrs in tx_list
                )
                
                # Check if source_gene matches one of the overlapping known genes
                matches_known_gene = any(
                    source_gene == source_gene_known
                    for source_gene_known, _ in overlapping_same_strand
                )
                
                should_keep = all_novel_isoforms and matches_known_gene
                
                if not should_keep:
                    # Not a valid novel isoform - remove all transcripts from this denovo gene
                    for _, tx_id, _ in tx_list:
                        genes_to_remove.add(tx_id)
                    conflicts_resolved += 1
                    overlapping_names = ', '.join(set(src for src, _ in overlapping_same_strand))
                    tc_str = ', '.join(sorted(transcript_classes - {''}))
                    reason = "not all novel isoforms" if not all_novel_isoforms else "source_gene mismatch"
                    logger.info(f"  Removed denovo gene {source_gene} [{tc_str}] ({reason}, overlaps: {overlapping_names})")
        
        # Also check denovo genes WITHOUT source_gene (from genes_without_source)
        if chrom in genes_without_source:
            for assigned_gene_id, tx_list in genes_without_source[chrom].items():
                # Check if this is denovo
                first_tx_id = tx_list[0][1]
                if not _is_denovo(first_tx_id):
                    continue  # Not denovo, keep it
                
                # Get the span of this denovo gene
                starts = [tx_obj.start for tx_obj, _, _ in tx_list]
                ends = [tx_obj.stop for tx_obj, _, _ in tx_list]
                gene_start = min(starts)
                gene_end = max(ends)
                strand = tx_list[0][0].strand
                
                # Check if it overlaps with any known gene on the same strand
                overlapping_known_genes = known_gene_tree.overlap(gene_start, gene_end)
                overlapping_same_strand = [
                    (source_gene, gene_strand) for (source_gene, gene_strand) 
                    in (item.data for item in overlapping_known_genes)
                    if gene_strand == strand
                ]
                
                if overlapping_same_strand:
                    # augPB gene overlaps with known gene(s) on same strand
                    # Remove it (no source_gene means it shouldn't be there)
                    for _, tx_id, _ in tx_list:
                        genes_to_remove.add(tx_id)
                    conflicts_resolved += 1
                    overlapping_names = ', '.join(set(src for src, _ in overlapping_same_strand))
                    logger.info(f"  Removed augPB gene {assigned_gene_id} [no source_gene] (overlaps: {overlapping_names})")
    
    # SECOND: Resolve conflicts between known genes (original logic)
    for chrom in genes_by_chrom:
        # Build interval tree for this chromosome
        gene_intervals = []
        for source_gene, tx_list in genes_by_chrom[chrom].items():
            # IMPORTANT: Check if this source_gene has transcripts at multiple disjoint loci
            # Cluster by OVERLAPPING LOCUS (regardless of mode) to avoid fragmentation
            # Transcripts from different modes at the same locus should be in one gene interval
            # Only create separate intervals for true multi-locus paralogs
            
            # Cluster by overlapping loci (all modes together)
            locus_clusters = []
            for tx_obj, tx_id, attrs in tx_list:
                # Find which cluster this transcript belongs to
                found_cluster = False
                for cluster in locus_clusters:
                    # Check if this transcript overlaps with any transcript in the cluster
                    cluster_start = min(t[0].start for t in cluster)
                    cluster_end = max(t[0].stop for t in cluster)
                    if not (tx_obj.stop <= cluster_start or tx_obj.start >= cluster_end):
                        # Overlaps - add to this cluster
                        cluster.append((tx_obj, tx_id, attrs))
                        found_cluster = True
                        break
                if not found_cluster:
                    # Start a new cluster
                    locus_clusters.append([(tx_obj, tx_id, attrs)])
            
            # Create a gene interval for each locus cluster
            for cluster_idx, cluster in enumerate(locus_clusters):
                # For overlap analysis, use only known-gene transcripts (exclude augPB) so that
                # spurious connections (e.g. CD38 linked to BST1) aren't created by augPB extending spans
                known_txs = [(tx_obj, tx_id, attrs) for tx_obj, tx_id, attrs in cluster
                            if attrs.get('alignment_mode') != 'augPB']
                if known_txs:
                    # Guard against outlier modes with greatly inflated spans for this gene.
                    # These can spuriously bridge nearby genes and trigger \"spurious overlap\" removals.
                    #
                    # We only disregard the outlier mode for span computation (overlap detection),
                    # not for transcript selection/scoring.
                    ratio = float(disregard_long_mode_ratio)
                    min_extra = int(disregard_long_mode_min_bp)

                    by_mode = collections.defaultdict(list)
                    for tx_obj, tx_id, attrs in known_txs:
                        by_mode[attrs.get("alignment_mode", "unknown")].append((tx_obj, tx_id, attrs))

                    # Compute span length per mode for this gene at this locus.
                    mode_spans = {}
                    for mode, items in by_mode.items():
                        s = min(t[0].start for t in items)
                        e = max(t[0].stop for t in items)
                        mode_spans[mode] = (s, e, e - s)

                    def _median(vals):
                        v = sorted(vals)
                        n = len(v)
                        if n == 0:
                            return 0
                        mid = n // 2
                        if n % 2 == 1:
                            return v[mid]
                        return int((v[mid - 1] + v[mid]) / 2)

                    excluded_modes = set()
                    if len(mode_spans) >= 2:
                        for mode, (_s, _e, span_len) in mode_spans.items():
                            # Ignore tiny fragments (common in txTM) when computing "typical" spans.
                            other_lens = [
                                l
                                for m, (_s2, _e2, l) in mode_spans.items()
                                if m != mode and l >= 1000
                            ]
                            med_other = _median(other_lens)
                            if med_other <= 0:
                                continue
                            if span_len >= ratio * med_other and (span_len - med_other) >= min_extra:
                                excluded_modes.add(mode)

                    filtered_cluster = list(cluster)
                    if excluded_modes:
                        span_txs = [
                            (tx_obj, tx_id, attrs)
                            for tx_obj, tx_id, attrs in known_txs
                            if attrs.get("alignment_mode", "unknown") not in excluded_modes
                        ]
                        # If we excluded everything (shouldn't happen), fall back to original.
                        if span_txs:
                            starts = [tx_obj.start for tx_obj, _, _ in span_txs]
                            ends = [tx_obj.stop for tx_obj, _, _ in span_txs]
                            gene_start = min(starts)
                            gene_end = max(ends)
                            # Outlier modes are excluded from span/overlap detection only.
                            # Do not drop their transcripts from the final consensus.
                            logger.info(
                                f"  Disregarded long-span mode(s) for {source_gene} cluster {cluster_idx+1}: "
                                f"{sorted(excluded_modes)} (span/overlap only; transcripts retained)"
                            )
                        else:
                            starts = [tx_obj.start for tx_obj, _, _ in known_txs]
                            ends = [tx_obj.stop for tx_obj, _, _ in known_txs]
                            gene_start = min(starts)
                            gene_end = max(ends)
                    else:
                        starts = [tx_obj.start for tx_obj, _, _ in known_txs]
                        ends = [tx_obj.stop for tx_obj, _, _ in known_txs]
                        gene_start = min(starts)
                        gene_end = max(ends)
                else:
                    # Fallback: cluster has only augPB transcripts, use full cluster
                    starts = [tx_obj.start for tx_obj, _, _ in cluster]
                    ends = [tx_obj.stop for tx_obj, _, _ in cluster]
                    gene_start = min(starts)
                    gene_end = max(ends)
                    filtered_cluster = list(cluster)
                gene_score = sum(attrs.get('score', 0) for _, _, attrs in filtered_cluster)
                
                # Log multi-locus genes
                if len(locus_clusters) > 1:
                    modes = set(attrs.get('alignment_mode', 'unknown') for _, _, attrs in filtered_cluster)
                    logger.info(f"  Multi-locus gene: {source_gene} ({modes}) cluster {cluster_idx+1}/{len(locus_clusters)} at {chrom}:{gene_start}-{gene_end} with {len(cluster)} transcripts")
                
                gene_intervals.append((gene_start, gene_end, source_gene, gene_score, filtered_cluster))
        
        # Sort by start position (only compare start and end, not the transcript objects)
        gene_intervals.sort(key=lambda x: (x[0], x[1]))
        
        # STEP 1: Identify multi-way overlap groups (3+ genes)
        # Build overlap graph to find connected components
        from collections import defaultdict
        overlap_graph = defaultdict(set)
        
        for i in range(len(gene_intervals)):
            start_i, end_i, source_i, score_i, tx_list_i = gene_intervals[i]
            strand_i = tx_list_i[0][0].strand
            
            for j in range(i + 1, len(gene_intervals)):
                start_j, end_j, source_j, score_j, tx_list_j = gene_intervals[j]
                
                # Stop if no overlap possible (intervals sorted by start)
                if start_j >= end_i:
                    break
                
                # Require actual overlap: [start_i,end_i] and [start_j,end_j] overlap iff
                # start_i < end_j AND start_j < end_i. We have start_j < end_i; also need start_i < end_j.
                if end_j <= start_i:
                    continue
                
                strand_j = tx_list_j[0][0].strand
                
                # Only consider same-strand overlaps
                if strand_i != strand_j:
                    continue
                
                # Add edge to overlap graph
                overlap_graph[i].add(j)
                overlap_graph[j].add(i)
        
        # Find connected components (overlap groups) using DFS
        visited = set()
        overlap_groups = []
        
        def dfs(node, group):
            if node in visited:
                return
            visited.add(node)
            group.add(node)
            for neighbor in overlap_graph[node]:
                dfs(neighbor, group)
        
        for i in range(len(gene_intervals)):
            if i not in visited:
                group = set()
                dfs(i, group)
                if len(group) >= 3:  # Only care about 3+ overlaps
                    overlap_groups.append(group)
        
        # STEP 2: Process multi-way overlaps (3+ genes)
        # For each group, remove genes that don't overlap with ANY other in the group in reference
        genes_removed_in_multiway = set()
        
        for group in overlap_groups:
            group_list = sorted(list(group))
            gene_ids = [gene_intervals[idx][2] for idx in group_list]
            
            logger.info(f"  Found {len(group_list)}-way overlap: {', '.join(gene_ids)}")

            # PASS 1 (no state mutation): collect every gene in the cluster that is
            # "spurious" — i.e., it does not overlap ANY other gene in this cluster
            # in the reference (same-biotype only). The original implementation
            # removed every spurious gene one-by-one with order-dependent state, which
            # wiped out legitimate gene-family clusters (e.g. VCX, MAGED4, ARHGEF35)
            # where N ≥ 2 candidates were simultaneously "spurious". We now keep the
            # best of those candidates and only remove the rest.
            spurious_candidates = []  # list of dicts: idx, source_gene, modes, tx_list, score
            for idx in group_list:
                start, end, source_gene, score, tx_list = gene_intervals[idx]
                normalized_source = re.sub(r'_\d+$', '', source_gene)

                overlaps_any_in_ref = False
                if ref_gene_coords and normalized_source in ref_gene_coords:
                    ref_chrom, ref_start, ref_end = ref_gene_coords[normalized_source]
                    for other_idx in group_list:
                        if other_idx == idx:
                            continue
                        other_source = gene_intervals[other_idx][2]
                        other_tx_list = gene_intervals[other_idx][4]
                        if _different_biotypes_allow_coexistence(tx_list, other_tx_list):
                            overlaps_any_in_ref = True
                            break
                        normalized_other = re.sub(r'_\d+$', '', other_source)
                        if normalized_other in ref_gene_coords:
                            other_ref_chrom, other_ref_start, other_ref_end = ref_gene_coords[normalized_other]
                            if (ref_chrom == other_ref_chrom and
                                ref_start < other_ref_end and other_ref_start < ref_end):
                                overlaps_any_in_ref = True
                                break

                if overlaps_any_in_ref:
                    continue

                modes = set(attrs.get('alignment_mode', 'unknown') for _, _, attrs in tx_list)
                spurious_candidates.append({
                    'idx': idx,
                    'source_gene': source_gene,
                    'modes': modes,
                    'tx_list': tx_list,
                    'score': score,
                })

            if not spurious_candidates:
                continue

            if skip_spurious_reference_overlap_removals:
                for cand in spurious_candidates:
                    logger.info(
                        f"  [skip spurious filter] Keeping {cand['source_gene']} ({cand['modes']}): would remove as spurious in "
                        f"{len(group_list)}-way overlap (doesn't overlap any other in reference)"
                    )
                continue

            # PASS 2: one protein-coding gene per assembly overlap cluster when genes do not
            # overlap each other in the reference (collapsed paralog / array projections).
            _MODE_PRIORITY = [
                'transMap', 'transMap_pairwise', 'txTM', 'augTM',
                'augTM_pairwise', 'augMP', 'augPB', 'strg',
            ]  # earlier = preferred; augMP last among projection modes
            def _mode_rank(modes_set):
                if not modes_set:
                    return len(_MODE_PRIORITY)
                return min(
                    _MODE_PRIORITY.index(m) if m in _MODE_PRIORITY else len(_MODE_PRIORITY)
                    for m in modes_set
                )

            def _candidate_sort_key(cand):
                return (
                    len(cand['modes']),
                    -_mode_rank(cand['modes']),
                    float(cand['score'] or 0),
                )

            spurious_candidates.sort(key=_candidate_sort_key, reverse=True)
            if len(spurious_candidates) >= 2:
                winner = spurious_candidates[0]
                logger.info(
                    f"  Kept {winner['source_gene']} ({winner['modes']}): single winner among "
                    f"{len(spurious_candidates)} genes in {len(group_list)}-way overlap "
                    f"(no reference overlap between candidates)"
                )
                for cand in spurious_candidates[1:]:
                    logger.info(
                        f"  Removed {cand['source_gene']} ({cand['modes']}): co-located gene "
                        f"in {len(group_list)}-way overlap (kept {winner['source_gene']})"
                    )
                    for _, tx_id, _ in cand['tx_list']:
                        genes_to_remove.add(tx_id)
                    genes_removed_in_multiway.add(cand['idx'])
                    conflicts_resolved += 1
            elif len(spurious_candidates) == 1:
                cand = spurious_candidates[0]
                logger.info(
                    f"  Kept {cand['source_gene']} ({cand['modes']}): sole spurious candidate "
                    f"in {len(group_list)}-way overlap"
                )
        
        # STEP 3: Find overlapping gene pairs (standard pairwise resolution)
        for i in range(len(gene_intervals)):
            # Skip if already removed in multi-way overlap
            if i in genes_removed_in_multiway:
                continue
            
            start_i, end_i, source_i, score_i, tx_list_i = gene_intervals[i]
            
            for j in range(i + 1, len(gene_intervals)):
                # Skip if already removed in multi-way overlap
                if j in genes_removed_in_multiway:
                    continue
                
                start_j, end_j, source_j, score_j, tx_list_j = gene_intervals[j]
                
                # If genes don't overlap, skip (intervals sorted by start)
                if start_j >= end_i:
                    break
                # Require actual overlap: need start_i < end_j as well
                if end_j <= start_i:
                    continue
                
                # Get strands for both genes
                strand_i = tx_list_i[0][0].strand  # Get strand from first transcript
                strand_j = tx_list_j[0][0].strand  # Get strand from first transcript
                
                # CRITICAL: Genes on opposite strands are NOT in conflict - skip them!
                if strand_i != strand_j:
                    continue
                
                # Get alignment modes for each gene (needed for logging)
                modes_i = set(attrs.get('alignment_mode', 'unknown') for _, _, attrs in tx_list_i)
                modes_j = set(attrs.get('alignment_mode', 'unknown') for _, _, attrs in tx_list_j)
                # Normalize gene IDs for reference lookups (strip txTM CNV suffixes like _1, _2)
                normalized_i = re.sub(r'_\d+$', '', source_i)
                normalized_j = re.sub(r'_\d+$', '', source_j)
                
                # Check if either gene is a readthrough gene
                # Readthrough genes overlap their parent genes by design - they should always coexist
                is_readthrough_i = False
                is_readthrough_j = False
                
                if readthrough_gene_set:
                    # First check if gene ID is in the readthrough set (from reference database)
                    is_readthrough_i = normalized_i in readthrough_gene_set
                    is_readthrough_j = normalized_j in readthrough_gene_set
                
                # Also check gene names - readthrough genes often have hyphenated names
                # But only if not already identified by database
                if not is_readthrough_i and not is_readthrough_j:
                    name_i = tx_list_i[0][2].get('source_gene_common_name', '')
                    name_j = tx_list_j[0][2].get('source_gene_common_name', '')
                    # Readthrough genes like GENE1-GENE2 have hyphens (but not pseudogenes like GENE1P-GENE2P)
                    is_readthrough_i = '-' in str(name_i) and 'P-' not in str(name_i) if name_i else False
                    is_readthrough_j = '-' in str(name_j) and 'P-' not in str(name_j) if name_j else False
                
                # If one is a readthrough gene, skip conflict resolution - they should coexist
                if is_readthrough_i or is_readthrough_j:
                    logger.info(f"  Kept readthrough gene: {source_i if is_readthrough_i else source_j} overlaps with {source_j if is_readthrough_i else source_i}")
                    continue
                
                # CRITICAL: Check for CNV genes at different loci
                # If two genes have the SAME normalized gene ID (e.g., both are ENSG00000215354.10)
                # but don't overlap in genomic coordinates, they are CNV copies at different loci
                # → Keep both, don't resolve conflict
                same_gene_id = (normalized_i == normalized_j)
                genes_overlap_in_assembly = not (end_j <= start_i or end_i <= start_j)
                
                if same_gene_id and not genes_overlap_in_assembly:
                    # Same gene, different loci (no overlap) → CNV copies, keep both
                    distance = min(abs(start_j - end_i), abs(start_i - end_j))
                    logger.info(f"  Kept CNV copies of {source_i}: non-overlapping loci separated by {distance:,}bp")
                    continue
                
                # Check if these genes should be allowed to overlap
                # ONLY allow if they actually overlap in the reference
                # This prevents spurious overlaps from mis-mapped genes
                should_keep_both = False
                reason = None
                
                # Check if the specific pair overlaps in reference
                # BUT: only apply this exemption for DIFFERENT genes (source_i != source_j)
                # Same gene at overlapping loci should be resolved, not exempted
                if ref_gene_coords and source_i != source_j and normalized_i in ref_gene_coords and normalized_j in ref_gene_coords:
                    ref_chrom_i, ref_start_i, ref_end_i = ref_gene_coords[normalized_i]
                    ref_chrom_j, ref_start_j, ref_end_j = ref_gene_coords[normalized_j]
                    
                    # Check if genes are on same chromosome and overlap in reference
                    if ref_chrom_i == ref_chrom_j:
                        # Check for overlap
                        if not (ref_end_i <= ref_start_j or ref_end_j <= ref_start_i):
                            should_keep_both = True
                            reason = "they also overlap in reference"
                
                # If genes should coexist, skip conflict resolution
                if should_keep_both:
                    logger.info(f"  Kept both genes: {source_i} ({modes_i}) and {source_j} ({modes_j}) ({reason})")
                    continue
                
                # transMap/txTM outrank augMP at co-located loci; augMP is a rescue mode.
                source_priority = OVERLAP_RESOLUTION_SOURCE_PRIORITY
                priority_i = max([source_priority.get(m, 0) for m in modes_i], default=0)
                priority_j = max([source_priority.get(m, 0) for m in modes_j], default=0)

                def _span_excluding_mode(tx_list, exclude_mode):
                    kept = [
                        (tx_obj, tx_id, attrs)
                        for tx_obj, tx_id, attrs in tx_list
                        if attrs.get("alignment_mode", "unknown") != exclude_mode
                    ]
                    if not kept:
                        return None
                    s = min(t[0].start for t in kept)
                    e = max(t[0].stop for t in kept)
                    return s, e

                def _overlaps(a, b):
                    return not (a[1] <= b[0] or b[1] <= a[0])

                # If this overlap is caused ONLY by one source mode expanding one of the genes,
                # ignore that mode for this overlap decision and keep both genes.
                if genes_overlap_in_assembly:
                    culprit_modes = []
                    # Try excluding each mode from gene i
                    for m in sorted(modes_i):
                        if len(modes_i) < 2:
                            break
                        span_i2 = _span_excluding_mode(tx_list_i, m)
                        if span_i2 is None:
                            continue
                        if not _overlaps(span_i2, (start_j, end_j)):
                            culprit_modes.append(("i", m))
                    # Try excluding each mode from gene j
                    for m in sorted(modes_j):
                        if len(modes_j) < 2:
                            break
                        span_j2 = _span_excluding_mode(tx_list_j, m)
                        if span_j2 is None:
                            continue
                        if not _overlaps((start_i, end_i), span_j2):
                            culprit_modes.append(("j", m))

                    # If exactly one mode exclusion resolves the overlap, treat this as boundary noise.
                    if len(culprit_modes) == 1:
                        side, m = culprit_modes[0]
                        logger.info(
                            f"  Kept both genes: {source_i} ({modes_i}) and {source_j} ({modes_j}) "
                            f"(overlap only due to {m} span in gene {side})"
                        )
                        continue
                
                # 2-WAY SPURIOUS OVERLAP FILTER
                # For pairs that don't overlap in reference, apply priority order
                # (Multi-way overlaps already handled above)
                # SKIP this check if either gene is a CNV copy (has _N suffix) - CNV copies are at different
                # loci than their reference position, so reference coordinate comparison is meaningless
                is_cnv_i = re.search(r'_\d+$', source_i) is not None
                is_cnv_j = re.search(r'_\d+$', source_j) is not None
                
                if (ref_gene_coords and source_i != source_j and
                    not is_cnv_i and not is_cnv_j and  # Skip if either is a CNV copy
                    normalized_i in ref_gene_coords and normalized_j in ref_gene_coords and
                    not should_keep_both):  # Already passed overlapping family check

                    ref_i = ref_gene_coords[normalized_i]
                    ref_j = ref_gene_coords[normalized_j]

                    # Check overlap in reference, including pairs on different chromosomes.
                    same_ref_chrom = ref_i[0] == ref_j[0]
                    ref_overlap = (
                        same_ref_chrom and (ref_i[1] < ref_j[2] and ref_j[1] < ref_i[2])
                    )

                    if not ref_overlap:
                        # Require meaningful assembly overlap before treating this as a collapsed-paralog
                        # / spurious overlap. Tiny overlaps are often just boundary noise from one mode.
                        overlap_start = max(start_i, start_j)
                        overlap_end = min(end_i, end_j)
                        overlap_len = max(0, overlap_end - overlap_start)
                        len_i = max(1, end_i - start_i)
                        len_j = max(1, end_j - start_j)
                        recip = min(overlap_len / len_i, overlap_len / len_j)
                        if (
                            overlap_len < int(spurious_overlap_min_assembly_overlap_bp)
                            and recip < float(spurious_overlap_min_reciprocal)
                        ):
                            # Keep both; skip the spurious removal logic for this pair.
                            continue

                        # Assembly overlap without reference overlap: keep one gene at this locus.
                        if same_ref_chrom:
                            ref_relation = f"{abs(ref_i[1] - ref_j[1])//1000}kb apart"
                        else:
                            ref_relation = f"different reference chromosomes ({ref_i[0]} vs {ref_j[0]})"

                        pc_i = ref_pc_ensg and normalized_i in ref_pc_ensg
                        pc_j = ref_pc_ensg and normalized_j in ref_pc_ensg
                        if pc_i and pc_j:
                            if priority_i > priority_j:
                                for _, tx_id, _ in tx_list_j:
                                    genes_to_remove.add(tx_id)
                                logger.info(
                                    f"  Removed {source_j} ({modes_j}): assembly overlap without ref overlap "
                                    f"with {source_i} ({ref_relation}, priority {priority_j} < {priority_i})"
                                )
                            elif priority_i < priority_j:
                                for _, tx_id, _ in tx_list_i:
                                    genes_to_remove.add(tx_id)
                                logger.info(
                                    f"  Removed {source_i} ({modes_i}): assembly overlap without ref overlap "
                                    f"with {source_j} ({ref_relation}, priority {priority_i} < {priority_j})"
                                )
                            elif score_i > score_j:
                                for _, tx_id, _ in tx_list_j:
                                    genes_to_remove.add(tx_id)
                                logger.info(
                                    f"  Removed {source_j} ({modes_j}): assembly overlap without ref overlap "
                                    f"with {source_i} ({ref_relation}, score {score_j:.1f} < {score_i:.1f})"
                                )
                            else:
                                for _, tx_id, _ in tx_list_i:
                                    genes_to_remove.add(tx_id)
                                logger.info(
                                    f"  Removed {source_i} ({modes_i}): assembly overlap without ref overlap "
                                    f"with {source_j} ({ref_relation}, score {score_i:.1f} < {score_j:.1f})"
                                )
                            conflicts_resolved += 1
                            continue
                        if pc_i and not pc_j:
                            for _, tx_id, _ in tx_list_j:
                                genes_to_remove.add(tx_id)
                            logger.info(
                                f"  Removed {source_j} ({modes_j}): spurious overlap with ref PC "
                                f"{source_i} ({ref_relation})"
                            )
                            conflicts_resolved += 1
                            continue
                        if pc_j and not pc_i:
                            for _, tx_id, _ in tx_list_i:
                                genes_to_remove.add(tx_id)
                            logger.info(
                                f"  Removed {source_i} ({modes_i}): spurious overlap with ref PC "
                                f"{source_j} ({ref_relation})"
                            )
                            conflicts_resolved += 1
                            continue

                        if skip_spurious_reference_overlap_removals:
                            logger.info(
                                f"  [skip spurious filter] Keeping both {source_i} and {source_j}: would remove one for "
                                f"spurious 2-way overlap (not in reference, {ref_relation})"
                            )
                        elif priority_i > priority_j:
                            for _, tx_id, _ in tx_list_j:
                                genes_to_remove.add(tx_id)
                            logger.info(f"  Removed {source_j} ({modes_j}): spurious 2-way overlap with {source_i} (not in reference, {ref_relation}, priority {priority_j} < {priority_i})")
                            conflicts_resolved += 1
                            continue
                        elif priority_i < priority_j:
                            for _, tx_id, _ in tx_list_i:
                                genes_to_remove.add(tx_id)
                            logger.info(f"  Removed {source_i} ({modes_i}): spurious 2-way overlap with {source_j} (not in reference, {ref_relation}, priority {priority_i} < {priority_j})")
                            conflicts_resolved += 1
                            continue
                        else:
                            # Same priority - use score to break tie
                            if score_i > score_j:
                                for _, tx_id, _ in tx_list_j:
                                    genes_to_remove.add(tx_id)
                                logger.info(f"  Removed {source_j} ({modes_j}): spurious 2-way overlap with {source_i} (not in reference, {ref_relation}, same priority, score {score_j:.1f} < {score_i:.1f})")
                                conflicts_resolved += 1
                                continue
                            else:
                                for _, tx_id, _ in tx_list_i:
                                    genes_to_remove.add(tx_id)
                                logger.info(f"  Removed {source_i} ({modes_i}): spurious 2-way overlap with {source_j} (not in reference, {ref_relation}, same priority, score {score_i:.1f} < {score_j:.1f})")
                                conflicts_resolved += 1
                                continue
                
                # Calculate reciprocal overlap
                overlap_start = max(start_i, start_j)
                overlap_end = min(end_i, end_j)
                overlap_len = overlap_end - overlap_start
                
                len_i = end_i - start_i
                len_j = end_j - start_j
                
                overlap_i = overlap_len / len_i if len_i > 0 else 0
                overlap_j = overlap_len / len_j if len_j > 0 else 0
                
                # Check if one gene is >90% contained within another (for protein_coding genes)
                one_gene_mostly_inside_other = (overlap_i > 0.9 or overlap_j > 0.9)
                
                both_protein_coding = (
                    _interval_is_protein_coding(tx_list_i)
                    and _interval_is_protein_coding(tx_list_j)
                )
                
                # If significant reciprocal overlap OR (both protein_coding AND one mostly inside other)
                if min(overlap_i, overlap_j) >= overlap_threshold or (both_protein_coding and one_gene_mostly_inside_other):
                    # Check if these are the same gene from different sources
                    # (modes_i and modes_j already defined earlier)
                    
                    # Only resolve conflict if genes are truly problematic:
                    # 1. Same source_gene but different loci (multi-locus from same method)
                    # 2. Very low-scoring gene overlapping much higher-scoring gene (likely error)
                    # 3. Different biotypes are ALWAYS okay to overlap (e.g., pseudogene + lncRNA)
                    
                    # Check if this is the same gene from the same method appearing twice
                    # Use normalized gene IDs (already computed above at line 3323)
                    # same_gene_id = (normalized_i == normalized_j)  # Already defined above
                    same_method = len(modes_i & modes_j) > 0  # overlapping methods
                    
                    different_biotypes = _different_biotypes_allow_coexistence(tx_list_i, tx_list_j)
                    
                    # Calculate score ratio to identify clearly inferior predictions
                    score_ratio = min(score_i, score_j) / max(score_i, score_j) if max(score_i, score_j) > 0 else 1.0
                    
                    # Check if one is denovo and the other is a known gene
                    is_augpb_i = any(attrs.get('alignment_mode') in ('augPB', 'strg') for _, _, attrs in tx_list_i)
                    is_augpb_j = any(attrs.get('alignment_mode') in ('augPB', 'strg') for _, _, attrs in tx_list_j)
                    
                    # Calculate overlap percentage (what % of the smaller gene is overlapped)
                    gene_i_len = end_i - start_i
                    gene_j_len = end_j - start_j
                    overlap_len = min(end_i, end_j) - max(start_i, start_j)
                    overlap_pct = overlap_len / min(gene_i_len, gene_j_len) if min(gene_i_len, gene_j_len) > 0 else 0
                    
                    source_priority = OVERLAP_RESOLUTION_SOURCE_PRIORITY
                    priority_i = max([source_priority.get(m, 0) for m in modes_i], default=0)
                    priority_j = max([source_priority.get(m, 0) for m in modes_j], default=0)
                    
                    # Check overlap patterns
                    nearly_identical = overlap_pct > 0.9  # >90% reciprocal overlap
                    one_gene_mostly_inside_other = (overlap_i > 0.9 or overlap_j > 0.9)  # one gene >90% contained in other
                    
                    # CRITICAL RULE: For protein_coding genes that are the SAME gene, 
                    # if they have >90% overlap (reciprocal OR one-sided), resolve by priority
                    # (only one copy of a protein_coding gene should exist at a locus)
                    high_overlap_protein_coding = both_protein_coding and same_gene_id and (nearly_identical or one_gene_mostly_inside_other)
                    
                    # Only remove if it's a clear case of redundancy or error:
                    # - Same gene from SAME method appearing twice WITH poor score ratio
                    # - OR SAME gene from DIFFERENT methods with poor score ratio (<30%)
                    # - OR SAME gene overlaps >90% reciprocally (nearly identical loci, duplicate)
                    # - OR SAME gene is protein_coding AND has >90% overlap (reciprocal OR one-sided)
                    # BUT never remove if genes have different biotypes
                    # ALSO never remove known gene in favor of augPB
                    # Different protein-coding genes at the same assembly locus: resolve unless they
                    # overlap in the reference (nested genes / gene families).
                    if (
                        both_protein_coding
                        and not same_gene_id
                        and not different_biotypes
                        and not should_keep_both
                    ):
                        if priority_i > priority_j or (priority_i == priority_j and score_i > score_j):
                            loser, winner = source_j, source_i
                            loser_tx, winner_tx = tx_list_j, tx_list_i
                            loser_modes, winner_modes = modes_j, modes_i
                        elif priority_j > priority_i or score_j > score_i:
                            loser, winner = source_i, source_j
                            loser_tx, winner_tx = tx_list_i, tx_list_j
                            loser_modes, winner_modes = modes_i, modes_j
                        else:
                            if len(tx_list_i) >= len(tx_list_j):
                                loser, winner = source_j, source_i
                                loser_tx, winner_tx = tx_list_j, tx_list_i
                                loser_modes, winner_modes = modes_j, modes_i
                            else:
                                loser, winner = source_i, source_j
                                loser_tx, winner_tx = tx_list_i, tx_list_j
                                loser_modes, winner_modes = modes_i, modes_j
                        for _, tx_id, _ in loser_tx:
                            genes_to_remove.add(tx_id)
                        conflicts_resolved += 1
                        logger.info(
                            f"  Resolved conflict: {loser} ({list(loser_modes)}) removed in favor of {winner} "
                            f"({list(winner_modes)}): co-located protein_coding without reference overlap"
                        )
                        continue

                    should_resolve = (((same_gene_id and same_method and score_ratio < 0.4) or 
                                      (same_gene_id and not same_method and score_ratio < 0.3) or  # Same gene, different methods
                                      (same_gene_id and nearly_identical) or  # CRITICAL: only if SAME gene
                                      high_overlap_protein_coding) 
                                     and not different_biotypes)
                    
                    # If one is augPB and other is known gene, prefer the known gene regardless of score
                    if (is_augpb_i or is_augpb_j) and is_augpb_i != is_augpb_j:
                        # Keep the known gene, remove augPB
                        if is_augpb_j:
                            # Gene j is augPB, remove it
                            for _, tx_id, _ in tx_list_j:
                                genes_to_remove.add(tx_id)
                            # Track competing source modes for the winner
                            for _, tx_id_i, _ in tx_list_i:
                                competing_sources[tx_id_i].update(modes_j)
                            conflicts_resolved += 1
                            logger.info(f"  Resolved conflict: {source_j} (augPB, score {score_j:.1f}) removed in favor of {source_i} (known gene, score {score_i:.1f})")
                        else:
                            # Gene i is augPB, remove it
                            for _, tx_id, _ in tx_list_i:
                                genes_to_remove.add(tx_id)
                            # Track competing source modes for the winner
                            for _, tx_id_j, _ in tx_list_j:
                                competing_sources[tx_id_j].update(modes_i)
                            conflicts_resolved += 1
                            logger.info(f"  Resolved conflict: {source_i} (augPB, score {score_i:.1f}) removed in favor of {source_j} (known gene, score {score_j:.1f})")
                    elif should_resolve:
                        # For protein_coding genes with high overlap OR nearly-identical overlaps, prefer by source priority
                        if (high_overlap_protein_coding and priority_i != priority_j) or (nearly_identical and priority_i != priority_j):
                            if priority_i > priority_j:
                                for _, tx_id, _ in tx_list_j:
                                    genes_to_remove.add(tx_id)
                                # Track competing source modes for the winner
                                for _, tx_id_i, _ in tx_list_i:
                                    competing_sources[tx_id_i].update(modes_j)
                                conflicts_resolved += 1
                                overlap_desc = f"{overlap_i*100:.0f}%/{overlap_j*100:.0f}%" if not nearly_identical else f"{overlap_pct*100:.0f}%"
                                logger.info(f"  Resolved conflict (protein_coding, {overlap_desc} overlap): {source_j} (priority {priority_j}, {list(modes_j)}) removed in favor of {source_i} (priority {priority_i}, {list(modes_i)})")
                            else:
                                for _, tx_id, _ in tx_list_i:
                                    genes_to_remove.add(tx_id)
                                # Track competing source modes for the winner
                                for _, tx_id_j, _ in tx_list_j:
                                    competing_sources[tx_id_j].update(modes_i)
                                conflicts_resolved += 1
                                overlap_desc = f"{overlap_i*100:.0f}%/{overlap_j*100:.0f}%" if not nearly_identical else f"{overlap_pct*100:.0f}%"
                                logger.info(f"  Resolved conflict (protein_coding, {overlap_desc} overlap): {source_i} (priority {priority_i}, {list(modes_i)}) removed in favor of {source_j} (priority {priority_j}, {list(modes_j)})")
                        # Otherwise, keep the gene with higher score
                        elif score_i > score_j:
                            for _, tx_id, _ in tx_list_j:
                                genes_to_remove.add(tx_id)
                            # Track competing source modes for the winner
                            for _, tx_id_i, _ in tx_list_i:
                                competing_sources[tx_id_i].update(modes_j)
                            conflicts_resolved += 1
                            logger.info(f"  Resolved conflict: {source_j} (score {score_j:.1f}) removed in favor of {source_i} (score {score_i:.1f}) [ratio: {score_ratio:.2f}]")
                        elif score_j > score_i:
                            for _, tx_id, _ in tx_list_i:
                                genes_to_remove.add(tx_id)
                            # Track competing source modes for the winner
                            for _, tx_id_j, _ in tx_list_j:
                                competing_sources[tx_id_j].update(modes_i)
                            conflicts_resolved += 1
                            logger.info(f"  Resolved conflict: {source_i} (score {score_i:.1f}) removed in favor of {source_j} (score {score_j:.1f}) [ratio: {score_ratio:.2f}]")
                        else:
                            # Equal scores - keep the one with more transcripts
                            # Sanity check: log intervals so we can verify they actually overlap
                            overlap_ok = (start_i < end_j and start_j < end_i)
                            if not overlap_ok:
                                logger.warning(f"  [BUG] Resolving tie but intervals do not overlap: {source_i} [{start_i}-{end_i}] vs {source_j} [{start_j}-{end_j}] - skipping removal")
                                continue
                            if len(tx_list_i) >= len(tx_list_j):
                                for _, tx_id, _ in tx_list_j:
                                    genes_to_remove.add(tx_id)
                                # Track competing source modes for the winner
                                for _, tx_id_i, _ in tx_list_i:
                                    competing_sources[tx_id_i].update(modes_j)
                                conflicts_resolved += 1
                                logger.info(f"  Resolved conflict (tie): {source_j} removed in favor of {source_i}")
                            else:
                                for _, tx_id, _ in tx_list_i:
                                    genes_to_remove.add(tx_id)
                                # Track competing source modes for the winner
                                for _, tx_id_j, _ in tx_list_j:
                                    competing_sources[tx_id_j].update(modes_i)
                                conflicts_resolved += 1
                                logger.info(f"  Resolved conflict (tie): {source_i} removed in favor of {source_j}")
                    else:
                        # Non-protein-coding or same-gene edge cases not handled above
                        reason = "different biotypes" if different_biotypes else f"similar scores (ratio: {score_ratio:.2f})"
                        logger.info(f"  Kept both overlapping genes: {source_i} and {source_j} [{reason}]")
    
    metrics['Discarded by source gene conflicts'] = len(genes_to_remove)
    
    # Add competing sources to transcript attributes before returning
    for tx_id in competing_sources:
        if tx_id in deduplicated_consensus and tx_id not in genes_to_remove:
            sources_list = sorted(competing_sources[tx_id])
            if sources_list:
                deduplicated_consensus[tx_id]['AdditionalSources'] = ','.join(sources_list)
    
    # Return filtered consensus
    return {tx_id: attrs for tx_id, attrs in deduplicated_consensus.items() if tx_id not in genes_to_remove}


def resolve_opposite_strand(deduplicated_consensus, tx_dict, metrics):
    """Resolve transcripts on opposite strands"""
    gene_dict = collections.defaultdict(list)
    deduplicated_strand_resolved_consensus = []
    
    for tx_id, attrs in deduplicated_consensus.items():
        tx_obj = tx_dict[tx_id]
        # For denovo putative_novel_isoform, use source_gene to group with the existing gene
        # For other transcripts, use original_gene_id to keep multi-locus copies separate
        is_augpb_isoform = (attrs.get('alignment_mode') in ('augPB', 'strg') and 
                           attrs.get('transcript_class') == 'putative_novel_isoform')
        if is_augpb_isoform and attrs.get('source_gene'):
            gene_key = attrs.get('source_gene')
        else:
            gene_key = attrs.get('original_gene_id', attrs.get('source_gene'))
        if gene_key and gene_key != 'N/A':
            gene_dict[gene_key].append([tx_obj, attrs])
        else:
            deduplicated_strand_resolved_consensus.append([tx_obj.name, attrs])
    
    # CRITICAL: For each gene_key, check if transcripts are at non-overlapping loci
    # If so, split them into separate genes (paralogs at different loci)
    final_gene_groups = []
    for gene_key, tx_list in gene_dict.items():
        # Sort by chromosome and start position
        tx_list_sorted = sorted(tx_list, key=lambda x: (x[0].chromosome, x[0].start))
        
        # Group transcripts by overlapping loci (paralogs at different loci should be separate)
        locus_groups = []
        current_group = []
        current_chrom = None
        current_end = 0
        
        for tx_obj, attrs in tx_list_sorted:
            if current_chrom is None or tx_obj.chromosome != current_chrom or tx_obj.start > current_end:
                # Start new locus group (different chromosome or no overlap)
                if current_group:
                    locus_groups.append(current_group)
                current_group = [[tx_obj, attrs]]
                current_chrom = tx_obj.chromosome
                current_end = tx_obj.stop
            else:
                # Extend current locus group (overlapping or adjacent)
                current_group.append([tx_obj, attrs])
                current_end = max(current_end, tx_obj.stop)
        
        if current_group:
            locus_groups.append(current_group)
        
        final_gene_groups.extend(locus_groups)
    
    # Now resolve strand conflicts within each locus group
    for gene_group in final_gene_groups:
        tx_objs, attrs_list = list(zip(*gene_group))

        if len(set(tx_obj.strand for tx_obj in tx_objs)) > 1:
            # Multiple strands - pick best
            strand_scores = collections.Counter()
            for tx_obj, attrs in gene_group:
                strand_scores[tx_obj.strand] += attrs.get('score', 0)
            best_strand = max(strand_scores.items(), key=lambda x: x[1])[0]

            for tx_obj, attrs in gene_group:
                if tx_obj.strand == best_strand:
                    deduplicated_strand_resolved_consensus.append([tx_obj.name, attrs])
                else:
                    metrics['Discarded by strand resolution'] += 1
        else:
            deduplicated_strand_resolved_consensus.extend([[tx_obj.name, attrs] for tx_obj, attrs in gene_group])
    
    return deduplicated_strand_resolved_consensus


def resolve_overlapping_different_genes(deduplicated_strand_resolved_consensus, tx_dict):
    """
    Resolve genes with different source_gene IDs that map to overlapping loci.
    When two different genes overlap, keep the one with better support.
    """
    import logging
    logger = logging.getLogger(__name__)
    
    # Group transcripts by source_gene and compute gene boundaries
    gene_groups = collections.defaultdict(list)
    for tx_id, attrs in deduplicated_strand_resolved_consensus:
        source_gene = attrs.get('source_gene')
        if source_gene and source_gene != 'N/A':
            gene_groups[source_gene].append((tx_id, attrs))
    
    # Compute gene boundaries (min start, max end) for each gene
    gene_boundaries = {}  # source_gene -> (chrom, strand, min_start, max_end, total_score, tx_count)
    for source_gene, tx_list in gene_groups.items():
        if not tx_list:
            continue
        chrom = tx_dict[tx_list[0][0]].chromosome
        strand = tx_dict[tx_list[0][0]].strand
        min_start = min(tx_dict[tx_id].start for tx_id, _ in tx_list)
        max_end = max(tx_dict[tx_id].stop for tx_id, _ in tx_list)
        total_score = sum(attrs.get('score', 0) for _, attrs in tx_list)
        tx_count = len(tx_list)
        gene_boundaries[source_gene] = (chrom, strand, min_start, max_end, total_score, tx_count)
    
    # Find overlapping gene pairs with different source_gene IDs on the same strand
    genes_to_remove = set()
    checked_pairs = set()
    
    for gene1, (chrom1, strand1, start1, end1, score1, count1) in gene_boundaries.items():
        if gene1 in genes_to_remove:
            continue
        for gene2, (chrom2, strand2, start2, end2, score2, count2) in gene_boundaries.items():
            if gene1 == gene2 or gene2 in genes_to_remove:
                continue
            if (gene1, gene2) in checked_pairs or (gene2, gene1) in checked_pairs:
                continue
            checked_pairs.add((gene1, gene2))
            
            # Check if genes overlap (same chromosome, same strand, overlapping positions)
            if chrom1 == chrom2 and strand1 == strand2:
                overlap = max(0, min(end1, end2) - max(start1, start2))
                if overlap > 0:
                    # Genes overlap - decide which one to keep
                    # Priority: 1) Total score, 2) Number of transcripts, 3) Gene span
                    avg_score1 = score1 / count1 if count1 > 0 else 0
                    avg_score2 = score2 / count2 if count2 > 0 else 0
                    
                    # Keep the gene with higher average score, or more transcripts if tied
                    if avg_score1 > avg_score2 or (avg_score1 == avg_score2 and count1 > count2):
                        genes_to_remove.add(gene2)
                        logger.info(f"  Removing overlapping gene '{gene2}' (score={avg_score1:.1f} vs {avg_score2:.1f}, "
                                   f"overlap={overlap}bp at {chrom1}:{start1}-{end1}) - keeping '{gene1}'")
                    else:
                        genes_to_remove.add(gene1)
                        logger.info(f"  Removing overlapping gene '{gene1}' (score={avg_score2:.1f} vs {avg_score1:.1f}, "
                                   f"overlap={overlap}bp at {chrom2}:{start2}-{end2}) - keeping '{gene2}'")
                        break  # gene1 is removed, no need to check further
    
    # Filter out transcripts from removed genes
    filtered_consensus = []
    removed_count = 0
    for tx_id, attrs in deduplicated_strand_resolved_consensus:
        source_gene = attrs.get('source_gene')
        if source_gene not in genes_to_remove:
            filtered_consensus.append([tx_id, attrs])
        else:
            removed_count += 1
    
    if removed_count > 0:
        logger.info(f"  Removed {removed_count} transcripts from {len(genes_to_remove)} overlapping genes")
    
    return filtered_consensus


def finalize_consensus_after_source_gene_resolution(
    consensus_seed,
    tx_dict,
    metrics,
    args,
    readthrough_gene_set,
    ref_gene_coords,
    genes_with_overlaps_in_ref,
    run_resolve_overlapping_different_genes,
    gene_biotype_map=None,
):
    """
    Shared tail of the consensus pipeline after per-chromosome merge: source-gene conflict resolution,
    invalid source_gene filter, strand resolution, optional overlapping-gene resolution, optional CDS overlap,
    then sort by position.
    """
    ref_pc_ensg = getattr(args, 'ref_pc_ensg', None)
    if not ref_pc_ensg and gene_biotype_map:
        ref_pc_ensg = {norm_ensg(g) for g, b in gene_biotype_map.items() if b == 'protein_coding'}

    logger.info("\nResolving conflicting source genes at overlapping loci...")
    logger.info("  Priority: Known genes (TxTM/transMap/augTM) > augPB (unless novel isoform or paralog)")
    deduplicated_consensus = resolve_conflicting_source_genes(
        consensus_seed,
        tx_dict,
        metrics,
        readthrough_gene_set=readthrough_gene_set,
        ref_gene_coords=ref_gene_coords,
        genes_with_overlaps_in_ref=genes_with_overlaps_in_ref,
        disregard_long_mode_ratio=float(getattr(args, "disregard_long_mode_ratio", 2.0)),
        disregard_long_mode_min_bp=int(getattr(args, "disregard_long_mode_min_bp", 50_000)),
        ref_pc_ensg=ref_pc_ensg,
    )
    logger.info(f"✓ After resolving source gene conflicts: {len(deduplicated_consensus)} transcripts")
    if metrics.get('Discarded by source gene conflicts', 0) > 0:
        logger.info(f"  Discarded {metrics['Discarded by source gene conflicts']} transcripts due to overlapping source genes")

    min_len_ratio = float(getattr(args, "min_pc_len_ratio_vs_reference", 0.0) or 0.0)
    if min_len_ratio > 0 and ref_gene_coords:
        logger.info("\nFiltering short transcripts vs reference gene span...")
        logger.info(
            f"  Removing individual transcripts with (tx_span/reference_gene_span) "
            f"ratio < {min_len_ratio:.3f} (genes kept if any isoform passes; "
            f"applies to all gene biotypes)"
        )
        deduplicated_consensus = filter_short_transcripts_vs_reference_gene_span(
            deduplicated_consensus,
            tx_dict,
            ref_gene_coords,
            metrics,
            min_length_ratio=min_len_ratio,
        )
        logger.info(f"✓ After short-transcript filter: {len(deduplicated_consensus)} transcripts")
        if metrics.get("Discarded short transcripts vs reference span", 0) > 0:
            logger.info(
                f"  Discarded {metrics['Discarded short transcripts vs reference span']} transcripts "
                f"short vs reference gene span"
            )

    logger.info("\nFiltering genes with invalid source_gene IDs...")
    logger.info("  Removing genes where source_gene is not a valid gene ID (augMP and augPB exempted)")
    deduplicated_consensus = filter_gene_fragments(deduplicated_consensus, tx_dict, ref_gene_coords, metrics, min_length_ratio=0.7)
    logger.info(f"✓ After filtering invalid genes: {len(deduplicated_consensus)} transcripts")
    if metrics.get('Discarded gene fragments', 0) > 0:
        logger.info(f"  Discarded {metrics['Discarded gene fragments']} transcripts with invalid source_gene IDs")

    logger.info("\nResolving opposite strand conflicts...")
    deduplicated_strand_resolved_consensus = resolve_opposite_strand(deduplicated_consensus, tx_dict, metrics)
    logger.info(f"✓ After strand resolution: {len(deduplicated_strand_resolved_consensus)} transcripts")
    if metrics.get('Discarded by strand resolution', 0) > 0:
        logger.info(f"  Discarded {metrics['Discarded by strand resolution']} transcripts due to strand conflicts")

    # Optional post-strand cleanup: remove overlapping protein-coding gene loci that do not overlap
    # in the reference. This is stricter than the default conflict resolver and is intended to
    # suppress spurious gene inflation from mis-mapped projections.
    filter_spurious_pc = bool(getattr(args, "filter_spurious_pc_overlaps_not_in_reference", False))
    if filter_spurious_pc and ref_gene_coords:
        logger.info("\nFiltering spurious protein-coding overlaps not supported by reference...")
        deduplicated_strand_resolved_consensus = filter_spurious_pc_overlaps_not_in_reference(
            deduplicated_strand_resolved_consensus,
            tx_dict,
            metrics,
            ref_gene_coords,
            spurious_overlap_min_assembly_overlap_bp=int(getattr(args, "spurious_overlap_min_assembly_overlap_bp", 2000)),
            spurious_overlap_min_reciprocal=float(getattr(args, "spurious_overlap_min_reciprocal", 0.02)),
        )
        logger.info(f"✓ After spurious overlap filter: {len(deduplicated_strand_resolved_consensus)} transcripts")

    if run_resolve_overlapping_different_genes:
        logger.info("\nResolving overlapping genes with different source_gene IDs...")
        deduplicated_strand_resolved_consensus = resolve_overlapping_different_genes(
            deduplicated_strand_resolved_consensus, tx_dict
        )
        logger.info(f"✓ After resolving overlapping different genes: {len(deduplicated_strand_resolved_consensus)} transcripts")

    if args.filter_overlapping_genes:
        logger.info("\nResolving overlapping CDS intervals...")
        gene_resolved_consensus = resolve_overlapping_cds_intervals(
            args.overlapping_ignore_bases,
            deduplicated_strand_resolved_consensus,
            tx_dict,
        )
        logger.info(f"✓ After overlap resolution: {len(gene_resolved_consensus)} transcripts")
    else:
        logger.info("\nSkipping overlap resolution (--filter-overlapping-genes not set)")
        gene_resolved_consensus = deduplicated_strand_resolved_consensus

    logger.info("\nSorting final consensus by genomic position...")
    final_consensus = sorted(
        gene_resolved_consensus,
        key=lambda tx_attrs: (tx_dict[tx_attrs[0]].chromosome, tx_dict[tx_attrs[0]].start),
    )
    logger.info(f"\n✓ Final consensus set: {len(final_consensus)} transcripts")
    return final_consensus


def filter_spurious_pc_overlaps_not_in_reference(
    deduplicated_strand_resolved_consensus,
    tx_dict,
    metrics,
    ref_gene_coords,
    spurious_overlap_min_assembly_overlap_bp=2000,
    spurious_overlap_min_reciprocal=0.02,
):
    """
    Remove overlapping protein_coding gene loci that do NOT overlap in the reference.

    Operates on the strand-resolved transcript list. Groups transcripts into per-locus gene intervals
    (by source_gene), then removes the lower-priority (or lower-score) gene interval when two
    different protein-coding genes overlap in the target but do not overlap in reference.

    This is a strict cleanup pass intended to reduce spurious inflation.
    """
    import re
    from collections import defaultdict

    def _parse_modes(tx_list):
        return set(attrs.get("alignment_mode", "unknown") for _, _, attrs in tx_list)

    def _priority(modes):
        return max((OVERLAP_RESOLUTION_SOURCE_PRIORITY.get(m, 0) for m in modes), default=0)

    def _overlap_len(a0, a1, b0, b1):
        s = max(a0, b0)
        e = min(a1, b1)
        return max(0, e - s)

    # Group transcripts by chromosome and source_gene.
    by_chrom = defaultdict(lambda: defaultdict(list))
    for tx_id, attrs in deduplicated_strand_resolved_consensus:
        tx_obj = tx_dict[tx_id]
        sg = attrs.get("source_gene")
        if not sg or sg == "N/A":
            continue
        by_chrom[tx_obj.chromosome][sg].append((tx_obj, tx_id, attrs))

    to_remove = set()
    removed_pairs = 0

    for chrom, gene_map in by_chrom.items():
        # Build per-locus intervals per source_gene (split multi-locus).
        intervals = []
        for sg, tx_list in gene_map.items():
            # Split into overlapping locus clusters.
            clusters = []
            for tx_obj, tx_id, attrs in tx_list:
                found = False
                for cl in clusters:
                    cs = min(t[0].start for t in cl)
                    ce = max(t[0].stop for t in cl)
                    if not (tx_obj.stop <= cs or tx_obj.start >= ce):
                        cl.append((tx_obj, tx_id, attrs))
                        found = True
                        break
                if not found:
                    clusters.append([(tx_obj, tx_id, attrs)])

            for cl in clusters:
                if not _interval_is_protein_coding(cl):
                    continue
                s = min(t[0].start for t in cl)
                e = max(t[0].stop for t in cl)
                score = sum(a.get("score", 0) for _, _, a in cl)
                modes = _parse_modes(cl)
                intervals.append((s, e, sg, score, modes, cl))

        intervals.sort(key=lambda x: (x[0], x[1]))

        for i in range(len(intervals)):
            s1, e1, g1, sc1, m1, cl1 = intervals[i]
            if any(tx_id in to_remove for _, tx_id, _ in cl1):
                continue
            for j in range(i + 1, len(intervals)):
                s2, e2, g2, sc2, m2, cl2 = intervals[j]
                if s2 >= e1:
                    break
                if g1 == g2:
                    continue
                if any(tx_id in to_remove for _, tx_id, _ in cl2):
                    continue

                ov = _overlap_len(s1, e1, s2, e2)
                if ov <= 0:
                    continue
                len1 = max(1, e1 - s1)
                len2 = max(1, e2 - s2)
                recip = min(ov / len1, ov / len2)
                if ov < int(spurious_overlap_min_assembly_overlap_bp) and recip < float(spurious_overlap_min_reciprocal):
                    continue

                n1 = re.sub(r"_\d+$", "", g1)
                n2 = re.sub(r"_\d+$", "", g2)
                if n1 not in ref_gene_coords or n2 not in ref_gene_coords:
                    continue
                ref1 = ref_gene_coords[n1]
                ref2 = ref_gene_coords[n2]
                same_ref_chrom = ref1[0] == ref2[0]
                ref_overlap = same_ref_chrom and (ref1[1] < ref2[2] and ref2[1] < ref1[2])
                if ref_overlap:
                    continue

                # Remove the lower-priority/score gene interval.
                p1 = _priority(m1)
                p2 = _priority(m2)
                remove_g = None
                if p1 > p2:
                    remove_g = (g2, cl2, p2, sc2, g1, p1, sc1)
                elif p2 > p1:
                    remove_g = (g1, cl1, p1, sc1, g2, p2, sc2)
                else:
                    # same priority: drop lower score
                    if sc1 >= sc2:
                        remove_g = (g2, cl2, p2, sc2, g1, p1, sc1)
                    else:
                        remove_g = (g1, cl1, p1, sc1, g2, p2, sc2)

                loser, loser_cl, loser_p, loser_sc, winner, winner_p, winner_sc = remove_g
                for _, tx_id, _ in loser_cl:
                    to_remove.add(tx_id)
                removed_pairs += 1
                ref_relation = (
                    f"{abs(ref1[1] - ref2[1])//1000}kb apart" if same_ref_chrom
                    else f"different reference chromosomes ({ref1[0]} vs {ref2[0]})"
                )
                logger.info(
                    f"  Removed {loser} ({sorted(list(_parse_modes(loser_cl)))}): spurious overlap with {winner} "
                    f"(not in reference, {ref_relation}, overlap={ov}bp, recip={recip:.3f}, "
                    f"priority {loser_p} vs {winner_p}, score {loser_sc:.1f} vs {winner_sc:.1f})"
                )

    if to_remove:
        metrics["Discarded spurious pc overlaps not in reference"] = len(to_remove)
        metrics["Discarded spurious pc overlaps not in reference (pairs)"] = removed_pairs
        return [[tx_id, attrs] for tx_id, attrs in deduplicated_strand_resolved_consensus if tx_id not in to_remove]
    return deduplicated_strand_resolved_consensus


def filter_short_transcripts_vs_reference_gene_span(
    deduplicated_consensus,
    tx_dict,
    ref_gene_coords,
    metrics,
    min_length_ratio=0.4,
):
    """
    Transcript-level length filter against the reference gene span for the same source_gene.

    Removes individual transcripts whose span is too short relative to the reference
    gene span. Genes are not removed as a unit; a gene remains represented if at
    least one of its isoforms passes (e.g. a full transMap model survives while
    fragmented txTM isoforms for the same ENSG are dropped).

    Normalization:
    - Strips trailing txTM/transMap CNV suffixes like _N for reference lookup.
    """
    import re

    if min_length_ratio <= 0:
        return deduplicated_consensus

    to_remove = set()

    for tx_id, attrs in deduplicated_consensus.items():
        source_gene = attrs.get("source_gene", "N/A")
        if not source_gene or source_gene == "N/A":
            continue
        normalized = re.sub(r"_\d+$", "", str(source_gene))
        if normalized not in ref_gene_coords:
            continue

        ref_chrom, ref_start, ref_end = ref_gene_coords[normalized]
        ref_len = max(1, int(ref_end) - int(ref_start))

        tx_obj = tx_dict[tx_id]
        tx_len = max(1, int(tx_obj.stop) - int(tx_obj.start))
        ratio = tx_len / ref_len

        if ratio < min_length_ratio:
            to_remove.add(tx_id)

    if to_remove:
        metrics["Discarded short transcripts vs reference span"] = len(to_remove)
        return {tx_id: attrs for tx_id, attrs in deduplicated_consensus.items() if tx_id not in to_remove}
    return deduplicated_consensus


def resolve_overlapping_cds_intervals(overlapping_ignore_bases, deduplicated_strand_resolved_consensus, tx_dict):
    """Resolve overlapping CDS from different genes"""
    attr_df = []
    with tools.fileOps.TemporaryFilePath() as tmp_gp, tools.fileOps.TemporaryFilePath() as tmp_clustered:
        with open(tmp_gp, 'w') as outf:
            for tx_id, attrs in deduplicated_strand_resolved_consensus:
                tx_obj = tx_dict[tx_id]
                tools.fileOps.print_row(outf, tx_obj.get_gene_pred())
                attr_df.append([tx_id, attrs.get('transcript_class', 'ortholog'), 
                               attrs.get('gene_biotype', 'unknown'),
                               attrs.get('source_gene', tx_obj.name2), 
                               attrs.get('score', 0)])
        
        cmd = ['clusterGenes', '-cds', f'-ignoreBases={overlapping_ignore_bases}',
               tmp_clustered, 'no', tmp_gp]
        tools.procOps.run_proc(cmd)
        cluster_df = pd.read_csv(tmp_clustered, sep='\t')
    
    attr_df = pd.DataFrame(attr_df, columns=['transcript_id', 'transcript_class', 'gene_biotype', 'gene_id', 'score'])
    m = attr_df.merge(cluster_df, left_on='transcript_id', right_on='gene')
    
    to_remove = set()
    for cluster_id, group in m.groupby('#cluster'):
        if len(set(group['gene_id'])) > 1:
            # Pick gene with highest average score
            avg_scores = group[['gene_id', 'score']].groupby('gene_id', as_index=False).mean()
            best_gene = avg_scores.sort_values('score', ascending=False).iloc[0]['gene_id']
            to_remove.update(set(group[group.gene_id != best_gene].transcript_id))
    
    return [[tx_id, attrs] for tx_id, attrs in deduplicated_strand_resolved_consensus if tx_id not in to_remove]


def calculate_completeness(final_consensus, metrics, gene_biotype_map, transcript_biotype_map):
    """Calculate final gene/transcript completeness and collect per-transcript metrics.
    
    This is the single point where Transcript Modes, Coverage, Identity, and support
    metrics are calculated — always from the FINAL consensus after all filtering.
    
    gene_biotype_map: {gene_id: biotype} from sqlInterface.get_gene_biotype_map
    transcript_biotype_map: {transcript_id: biotype} from sqlInterface.get_transcript_biotype_map
    """
    found_genes = collections.defaultdict(set)
    found_txs = collections.Counter()
    
    # Reset accumulation metrics so they reflect only the final consensus
    metrics['Transcript Modes'] = collections.Counter()
    metrics['Coverage'] = collections.defaultdict(list)
    metrics['Identity'] = collections.defaultdict(list)
    for key in ['Splice Support', 'Exon Support', 'Original Introns',
                'Splice Annotation Support', 'Exon Annotation Support']:
        metrics[key] = collections.defaultdict(list)
    
    for aln_id, c in final_consensus:
        source_gene = c.get('source_gene')
        if source_gene and source_gene != 'N/A':
            found_genes[c.get('gene_biotype', 'unknown')].add(source_gene)
        
        biotype = c.get('transcript_biotype', 'unknown')
        found_txs[biotype] += 1
        
        mode = c.get('alignment_mode', 'unknown')
        if biotype in ('protein_coding', 'unknown_likely_coding'):
            metrics['Transcript Modes'][mode] += 1
        
        # Collect per-transcript plotting metrics stored during selection
        tx_metrics = c.get('_metrics')
        if tx_metrics:
            metrics['Coverage'][biotype].append(tx_metrics.get('Coverage', 0))
            metrics['Identity'][biotype].append(tx_metrics.get('Identity', 0))
            for key in ['Splice Support', 'Exon Support', 'Original Introns',
                        'Splice Annotation Support', 'Exon Annotation Support']:
                metrics[key][biotype].append(tx_metrics.get(key, 0))
    
    found_genes_counts = {biotype: len(gene_set) for biotype, gene_set in found_genes.items()}
    metrics['Completeness'] = {'Gene': found_genes_counts, 'Transcript': dict(found_txs)}
    
    ref_genes_by_biotype = collections.Counter()
    for gene_id, biotype in gene_biotype_map.items():
        ref_genes_by_biotype[biotype] += 1
    
    ref_txs_by_biotype = collections.Counter()
    for tx_id, biotype in transcript_biotype_map.items():
        ref_txs_by_biotype[biotype] += 1
    
    for biotype, total in ref_genes_by_biotype.items():
        found = found_genes_counts.get(biotype, 0)
        metrics['Gene Missing'][biotype] = max(0, total - found)
    
    for biotype, total in ref_txs_by_biotype.items():
        found = found_txs.get(biotype, 0)
        metrics['Transcript Missing'][biotype] = max(0, total - found)


def write_consensus_gps(consensus_gp, consensus_gp_info, final_consensus, tx_dict, genome):
    """Write consensus genePred and info files"""
    import logging
    logger = logging.getLogger(__name__)
    
    genes_seen = collections.defaultdict(dict)
    gene_count = 0
    consensus_gene_dict = DefaultOrderedDict(lambda: DefaultOrderedDict(list))
    gp_infos = []
    
    # Track gene loci to detect non-overlapping paralogs
    # Format: gene_loci[chromosome][source_gene] = [(gene_key, min_start, max_end), ...]
    # min_start and max_end represent the span of all transcripts at this locus
    gene_loci = collections.defaultdict(lambda: collections.defaultdict(list))
    
    # Track boundaries per gene_key for isoform grouping
    # Format: gene_key_boundaries[chromosome][gene_key] = (min_start, max_end)
    gene_key_boundaries = collections.defaultdict(dict)
    
    with open(consensus_gp, 'w') as out_gp:
        for tx_count, (tx, attrs) in enumerate(final_consensus, 1):
            tx_obj = tx_dict[tx]
            name = ID_TEMPLATE.format(genome=genome, tag_type='T', unique_id=tx_count)
            score = int(round(attrs.get('score', 0)))
            
            # Gene grouping logic:
            # - For putative_novel_isoform: group by source_gene (the known gene they're an isoform of)
            # - For CNV copies: extract copy number from source_gene or alignment_id and group by (base_gene, copy_num)
            # - For others: use source_gene
            source_gene = attrs.get('source_gene')
            if source_gene is None or source_gene == 'None' or source_gene == 'N/A':
                source_gene = tx_obj.name2
            
            alignment_id = attrs.get('alignment_id', '')
            transcript_class = attrs.get('transcript_class', 'ortholog')
            
            import re
            
            # Check if this is a CNV copy
            # CNV copies can have _N suffix in either source_gene or alignment_id
            copy_num = None
            base_source_gene = source_gene
            
            # First check source_gene for _N suffix (transMap style: ENSG00000236424.7_1)
            # Only treat as CNV copy if the number is small (≤20) - larger numbers are likely
            # part of the gene ID (e.g., hg002_chrY_paternal_166), not a copy number
            if source_gene:
                source_match = re.search(r'_(\d+)$', source_gene)
                if source_match and int(source_match.group(1)) <= 20:
                    copy_num = source_match.group(1)
                    base_source_gene = source_gene[:source_match.start()]
            
            # If not found in source_gene, check alignment_id (txTM style: ENST00000428845.6_1)
            if copy_num is None and alignment_id:
                align_match = re.search(r'_(\d+)$', alignment_id)
                if align_match and int(align_match.group(1)) > 0:
                    copy_num = align_match.group(1)
            
            # Debug logging for TSPY10
            if 'TSPY10' in str(attrs.get('source_gene_common_name', '')):
                import logging
                logger = logging.getLogger(__name__)
                logger.info(f"  TSPY10 transcript: source_gene={source_gene}, copy_num={copy_num}, class={transcript_class}, at {tx_obj.start}-{tx_obj.stop}")
            
            if transcript_class == 'putative_novel_isoform' and attrs.get('source_gene'):
                # Novel isoforms: For multi-locus genes, assign to the nearest copy based on genomic distance
                # First, check if there are existing genes with this source_gene on this chromosome
                existing_genes = []
                for existing_key, existing_gene_id in genes_seen[tx_obj.chromosome].items():
                    # Check if this existing gene is for the same source gene
                    if isinstance(existing_key, tuple):
                        # CNV copy: (base_gene, copy_num)
                        if existing_key[0] == base_source_gene:
                            existing_genes.append((existing_key, existing_gene_id))
                    elif existing_key == attrs['source_gene'] or existing_key == base_source_gene:
                        # Single-copy gene with same source
                        existing_genes.append((existing_key, existing_gene_id))
                
                if existing_genes:
                    # Find the nearest existing gene by calculating distance to gene boundaries
                    min_distance = float('inf')
                    nearest_gene_key = None
                    
                    for gene_key_candidate, gene_id_candidate in existing_genes:
                        # Get all transcripts for this gene to find its boundaries
                        gene_name = ID_TEMPLATE.format(genome=genome, tag_type='G', unique_id=gene_id_candidate)
                        if gene_name in consensus_gene_dict[tx_obj.chromosome]:
                            gene_txs = consensus_gene_dict[tx_obj.chromosome][gene_name]
                            gene_starts = [t.start for t, _ in gene_txs]
                            gene_ends = [t.stop for t, _ in gene_txs]
                            gene_start = min(gene_starts)
                            gene_end = max(gene_ends)
                            
                            # Calculate distance: 0 if overlapping, otherwise minimum gap
                            if tx_obj.start <= gene_end and tx_obj.stop >= gene_start:
                                # Overlaps
                                distance = 0
                            else:
                                # Gap distance
                                distance = min(abs(tx_obj.start - gene_end), abs(gene_start - tx_obj.stop))
                            
                            if distance < min_distance:
                                min_distance = distance
                                nearest_gene_key = gene_key_candidate
                    
                    if nearest_gene_key:
                        gene_key = nearest_gene_key
                    else:
                        # Fallback to source_gene if we couldn't find distances
                        gene_key = attrs['source_gene']
                else:
                    # No existing genes yet, use source_gene
                    gene_key = attrs['source_gene']
            elif copy_num is not None:
                # CNV copy: use (base_gene, copy_number) so all transcripts from the same copy are grouped
                gene_key = (base_source_gene, copy_num)
            else:
                # Regular gene: use source_gene
                # BUT: if there's already a gene with this source_gene at a different location,
                # create a separate gene entry (paralogs/multi-locus genes)
                gene_key = source_gene
                
                # Check if a gene with this source_gene exists at a non-overlapping location (paralog)
                # Use gene_loci tracker to find all existing loci for this source_gene
                existing_loci = gene_loci[tx_obj.chromosome].get(source_gene, [])
                
                # Check if current transcript meaningfully overlaps an existing locus.
                # Use reciprocal overlap (not just any 1bp touch) to avoid chaining distant loci
                # into one giant gene span via long/bridge transcripts.
                #
                # IMPORTANT: Use the FIXED locus boundaries stored in gene_loci for overlap detection.
                # Using expanding boundaries here can accidentally "chain" disjoint loci if a transcript
                # was previously mis-assigned to the wrong locus.
                matching_locus = None
                for locus_gene_key, locus_min_start, locus_max_end in existing_loci:
                    overlap_start = max(tx_obj.start, locus_min_start)
                    overlap_end = min(tx_obj.stop, locus_max_end)
                    overlap_len = max(0, overlap_end - overlap_start)
                    if overlap_len <= 0:
                        continue

                    tx_len = max(1, tx_obj.stop - tx_obj.start)
                    locus_len = max(1, locus_max_end - locus_min_start)
                    reciprocal = min(overlap_len / tx_len, overlap_len / locus_len)

                    # Require meaningful overlap to join loci.
                    if overlap_len >= 2000 or reciprocal >= 0.02:
                        matching_locus = (locus_gene_key, locus_min_start, locus_max_end)
                        break
                
                if matching_locus:
                    # Use the existing gene_key for this locus
                    gene_key = matching_locus[0]
                else:
                    # Non-overlapping paralog at a different locus
                    # Create a unique gene_key using (source_gene, locus_index)
                    locus_index = len(existing_loci)
                    gene_key = (source_gene, locus_index) if locus_index > 0 else source_gene
                    
                    # Debug logging for TSPY10
                    if 'TSPY10' in str(attrs.get('source_gene_common_name', '')):
                        import logging
                        logger = logging.getLogger(__name__)
                        logger.info(f"  TSPY10 paralog detected: {source_gene} at {tx_obj.start}-{tx_obj.stop}, locus_index={locus_index}, gene_key={gene_key}")
            
            if gene_key not in genes_seen[tx_obj.chromosome]:
                gene_count += 1
                genes_seen[tx_obj.chromosome][gene_key] = gene_count
                
                # Track this gene locus for future paralog detection with FIXED boundaries
                # These boundaries are used ONLY for detecting non-overlapping paralogs
                gene_loci[tx_obj.chromosome][source_gene].append((gene_key, tx_obj.start, tx_obj.stop))
                
                # Also track boundaries per gene_key for isoform grouping
                gene_key_boundaries[tx_obj.chromosome][gene_key] = (tx_obj.start, tx_obj.stop)
            else:
                # Update the gene_key boundaries for isoform grouping
                # This allows isoforms at the same locus to be grouped correctly
                if gene_key in gene_key_boundaries[tx_obj.chromosome]:
                    curr_start, curr_stop = gene_key_boundaries[tx_obj.chromosome][gene_key]
                    gene_key_boundaries[tx_obj.chromosome][gene_key] = (
                        min(curr_start, tx_obj.start),
                        max(curr_stop, tx_obj.stop)
                    )
            
            gene_id = genes_seen[tx_obj.chromosome][gene_key]
            name2 = ID_TEMPLATE.format(genome=genome, tag_type='G', unique_id=gene_id)
            out_gp.write('\t'.join(tx_obj.get_gene_pred(name=name, name2=name2, score=score)) + '\n')
            
            attrs['transcript_id'] = name
            attrs['gene_id'] = name2
            gp_info_attrs = {k: v for k, v in attrs.items() if not k.startswith('_')}
            gp_infos.append(gp_info_attrs)
            consensus_gene_dict[tx_obj.chromosome][name2].append([tx_obj, attrs])
    
    gp_info_df = pd.DataFrame(gp_infos).set_index(['gene_id', 'transcript_id'])
    if 'alternative_source_transcripts' not in gp_info_df.columns:
        gp_info_df['alternative_source_transcripts'] = 'N/A'
    
    with open(consensus_gp_info, 'w') as outf:
        gp_info_df.to_csv(outf, sep='\t', na_rep='N/A')
    
    return consensus_gene_dict


def write_consensus_gff3(consensus_gene_dict, consensus_gff3):
    """Write consensus GFF3 file"""
    from cat.consensus import write_consensus_gff3 as original_write_gff3
    original_write_gff3(consensus_gene_dict, consensus_gff3)


def write_consensus_fastas(consensus_gene_dict, consensus_fasta, consensus_protein_fasta, fasta):
    """Write consensus FASTA files"""
    import logging
    logger = logging.getLogger(__name__)
    
    seq_dict = tools.bio.get_sequence_dict(fasta)
    skipped_count = 0
    
    with open(consensus_fasta, 'w') as cfa, open(consensus_protein_fasta, 'w') as cpfa:
        for chrom in sorted(consensus_gene_dict):
            for gene_id, tx_list in consensus_gene_dict[chrom].items():
                for tx_obj, attrs in tx_list:
                    try:
                        mrna_seq = tx_obj.get_mrna(seq_dict)
                        tools.bio.write_fasta(cfa, attrs['transcript_id'], mrna_seq)
                        if tx_obj.cds_size > 0:
                            protein_seq = tx_obj.get_protein_sequence(seq_dict)
                        elif str(attrs.get('transcript_class', '')).startswith('putative_novel'):
                            # For novel predictions without CDS annotations (e.g., StringTie),
                            # treat the full spliced transcript as coding sequence.
                            protein_seq = tools.bio.translate_sequence(mrna_seq)
                        else:
                            protein_seq = ''
                        if protein_seq:
                            tools.bio.write_fasta(cpfa, attrs['transcript_id'], protein_seq)
                    except (ValueError, KeyError) as e:
                        # Skip transcripts with coordinate issues (e.g., coordinates exceed chromosome length)
                        logger.warning(f"  Skipping transcript {attrs['transcript_id']}: {e}")
                        skipped_count += 1
                        continue
    
    if skipped_count > 0:
        logger.warning(f"⚠ Skipped {skipped_count} transcripts due to coordinate issues")


if __name__ == "__main__":
    main()

