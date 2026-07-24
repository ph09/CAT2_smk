#!/usr/bin/env python3
"""
CLI entry point for transcript-level minimap2 mapping (txTM).

Invoked from Snakemake (directly in local mode, via SLURM in cluster mode).
"""

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path

import yaml

from cat.transcript_map import run_transcript_map


def transcript_map_kwargs(cfg: dict) -> dict:
    """Build run_transcript_map keyword args from a config dict."""
    cfg_tm = cfg.get('transcript_map', {})
    keep_ge = cfg_tm.get('gene_locus_keep_if_match_frac_ge', 0.80)
    return dict(
        min_coverage=cfg_tm.get('min_coverage', 0.50),
        min_identity=cfg_tm.get('min_identity', 0.50),
        secondary_ratio=cfg_tm.get('secondary_ratio', 0.50),
        copy_min_identity=cfg_tm.get('copy_min_identity', cfg.get('txTM_sc', 0.80)),
        max_span_ratio=cfg_tm.get('max_span_ratio', 3.60),
        max_span_extra_bp=cfg_tm.get('max_span_extra_bp', 28000),
        copy_min_coverage=cfg_tm.get('copy_min_coverage', 0.78),
        max_blocks_factor=cfg_tm.get('max_blocks_factor', 0.0),
        max_blocks_extra=cfg_tm.get('max_blocks_extra', 28),
        max_blocks_abs=cfg_tm.get('max_blocks_abs', 220),
        primary_max_blocks_factor=cfg_tm.get('primary_max_blocks_factor', 0.0),
        primary_max_blocks_extra=cfg_tm.get('primary_max_blocks_extra', 40),
        primary_max_blocks_abs=cfg_tm.get('primary_max_blocks_abs', 350),
        merged_exon_max_factor=cfg_tm.get('merged_exon_max_factor', 3.5),
        merged_exon_max_extra=cfg_tm.get('merged_exon_max_extra', 15),
        merged_exon_max_abs=cfg_tm.get('merged_exon_max_abs', 200),
        copy_merged_exon_max_factor=cfg_tm.get('copy_merged_exon_max_factor', 2.5),
        copy_merged_exon_max_extra=cfg_tm.get('copy_merged_exon_max_extra', 8),
        copy_merged_exon_max_abs=cfg_tm.get('copy_merged_exon_max_abs', 80),
        min_target_span_ratio=cfg_tm.get('min_target_span_ratio'),
        min_ref_span_for_span_ratio_bp=cfg_tm.get(
            'min_ref_span_for_span_ratio_bp', 15000
        ),
        end_block_min_bp=cfg_tm.get('end_block_min_bp', 14),
        end_block_max_gap_bp=cfg_tm.get('end_block_max_gap_bp', 2600),
        filter_weak_gene_loci=cfg_tm.get('filter_weak_gene_loci', True),
        gene_locus_min_transcript_frac=cfg_tm.get('gene_locus_min_transcript_frac', 0.45),
        gene_locus_min_match_frac=cfg_tm.get('gene_locus_min_match_frac', 0.38),
        gene_locus_keep_if_match_frac_ge=(
            None if keep_ge is False else keep_ge
        ),
        gene_locus_escape_min_transcript_frac=cfg_tm.get(
            'gene_locus_escape_min_transcript_frac', 0.15
        ),
        chimeric_intron_ratio=cfg_tm.get('chimeric_intron_ratio', 8.0),
        chimeric_intron_pad_bp=cfg_tm.get('chimeric_intron_pad_bp', 8000),
        chimeric_intron_floor_bp=cfg_tm.get('chimeric_intron_floor_bp', 35000),
        gene_region_rescue=cfg_tm.get('gene_region_rescue', True),
        gene_region_flank_bp=cfg_tm.get('gene_region_flank_bp', 2000),
        gene_region_minimap2_preset=cfg_tm.get('gene_region_minimap2_preset', 'asm10'),
        gene_region_min_coverage=cfg_tm.get('gene_region_min_coverage', 0.50),
        gene_region_min_identity=cfg_tm.get('gene_region_min_identity', 0.50),
        gene_region_max_secondary=cfg_tm.get('gene_region_max_secondary', 10),
        gene_region_secondary_ratio=cfg_tm.get('gene_region_secondary_ratio', 0.5),
        gene_region_snap_max_bp=cfg_tm.get('gene_region_snap_max_bp', 50),
        gene_region_rescue_deep=cfg_tm.get('gene_region_rescue_deep', True),
        gene_region_deep_flank_bp=cfg_tm.get('gene_region_deep_flank_bp', 5000),
        gene_region_deep_minimap2_preset=cfg_tm.get(
            'gene_region_deep_minimap2_preset', 'asm20'
        ),
        gene_region_deep_min_coverage=cfg_tm.get('gene_region_deep_min_coverage', 0.50),
        gene_region_deep_min_identity=cfg_tm.get('gene_region_deep_min_identity', 0.50),
        gene_region_deep_max_secondary=cfg_tm.get('gene_region_deep_max_secondary', 50),
        gene_region_deep_secondary_ratio=cfg_tm.get(
            'gene_region_deep_secondary_ratio', 0.5
        ),
        gene_region_deep_extra_minimap2=cfg_tm.get(
            'gene_region_deep_extra_minimap2', '--end-bonus 5'
        ),
        gene_region_deep_snap_max_bp=cfg_tm.get('gene_region_deep_snap_max_bp', 50),
        extra_minimap2=cfg_tm.get('extra_minimap2', ''),
        two_pass=cfg_tm.get('two_pass', True),
    )


def run_for_genome(cfg: dict, genome: str, threads: int, log_path: str,
                   output_gp: str, output_psl: str, output_gtf: str,
                   output_attrs: str, output_dups: str):
    """Run transcript_map for one genome and finalize Snakemake outputs."""
    work_dir = cfg['work_dir']
    ref_genome = cfg['ref_genome']
    out_dir = f"{work_dir}/txTM"
    ref_db = f"{work_dir}/databases/{ref_genome}.db"

    # Prefer TMPDIR when writable; otherwise fall back to work_dir/tmp then
    # the system temp dir. Do not hard-require site-specific paths like /data/tmp.
    candidates = []
    if os.environ.get('TMPDIR'):
        candidates.append(os.environ['TMPDIR'])
    candidates.append(f"{work_dir}/tmp")
    candidates.append(tempfile.gettempdir())
    tmp_base = None
    for cand in candidates:
        try:
            Path(cand).mkdir(parents=True, exist_ok=True)
            probe = Path(cand) / f".cat2_write_probe_{os.getpid()}"
            probe.touch()
            probe.unlink()
            tmp_base = cand
            break
        except OSError:
            continue
    if tmp_base is None:
        raise PermissionError(
            "No writable temp directory found. Set TMPDIR to a writable path "
            f"(tried: {', '.join(candidates)})."
        )

    with tempfile.TemporaryDirectory(prefix='txmap_', dir=tmp_base) as tmp_dir:
        run_transcript_map(
            work_dir=work_dir,
            ref_genome=ref_genome,
            genome=genome,
            threads=threads,
            out_dir=out_dir,
            tmp_dir=tmp_dir,
            log_path=log_path,
            ref_db_path=ref_db,
            **transcript_map_kwargs(cfg),
        )

    raw_gp = Path(out_dir) / f"{genome}.gp"
    if raw_gp.exists() and str(raw_gp) != output_gp:
        shutil.move(str(raw_gp), output_gp)

    raw_psl = Path(out_dir) / f"{genome}_filtered.psl"
    if raw_psl.exists() and str(raw_psl) != output_psl:
        shutil.move(str(raw_psl), output_psl)

    # gp_attrs is written by run_transcript_map; attrs path must match output.
    written_attrs = Path(out_dir) / f"{genome}_txTM.gp_attrs"
    if written_attrs.exists() and str(written_attrs) != output_attrs:
        shutil.move(str(written_attrs), output_attrs)

    for path in (output_gtf, output_dups):
        Path(path).touch()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Run transcript-level minimap2 mapping (txTM) for one genome.'
    )
    parser.add_argument('--config', required=True,
                        help='YAML config snapshot (work_dir, ref_genome, transcript_map)')
    parser.add_argument('--genome', required=True)
    parser.add_argument('--threads', type=int, required=True)
    parser.add_argument('--log', required=True)
    parser.add_argument('--output-gp', required=True)
    parser.add_argument('--output-psl', required=True)
    parser.add_argument('--output-gtf', required=True)
    parser.add_argument('--output-attrs', required=True)
    parser.add_argument('--output-dups', required=True)
    args = parser.parse_args(argv)

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    for key in ('work_dir', 'ref_genome'):
        if key not in cfg:
            sys.exit(f"Config snapshot missing required key: {key}")

    run_for_genome(
        cfg=cfg,
        genome=args.genome,
        threads=args.threads,
        log_path=args.log,
        output_gp=args.output_gp,
        output_psl=args.output_psl,
        output_gtf=args.output_gtf,
        output_attrs=args.output_attrs,
        output_dups=args.output_dups,
    )


if __name__ == '__main__':
    main()
