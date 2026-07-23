#!/usr/bin/env python3
"""
CAT (Comparative Annotation Toolkit) Snakemake Pipeline
A comprehensive pipeline for comparative genomic annotation using multiple approaches.
"""

import os
import sys
from pathlib import Path

# Load configuration. No default configfile is hardcoded on purpose: pass your own
# with `--configfile <file>` (e.g. input.yaml or panprimates.yaml). Hardcoding a
# default caused CLI configs to be deep-MERGED on top of it, leaking stale keys
# (e.g. another run's stringtie_genomes / transcriptomic_data) into your run.
if not config:
    raise SystemExit(
        "No configuration loaded. Pass one with:  "
        "snakemake --configfile <your_config.yaml> ...\n"
        "(e.g. --configfile input.yaml or --configfile panprimates.yaml)"
    )

# Validate configuration parameters up-front and fail fast with a single,
# aggregated, human-readable error listing everything that is wrong, rather
# than crashing deep inside a rule hours into a run.
def validate_config():
    """Validate the loaded config; raise ValueError listing all problems."""
    errors = []
    warnings = []

    # 1. Required keys must be present before any other check can run.
    required_config = ["work_dir", "hal", "annotation", "ref_genome", "genomes"]
    missing = [p for p in required_config if p not in config]
    if missing:
        raise ValueError(
            "CONFIG ERROR: missing required key(s): " + ", ".join(missing)
        )

    genomes = config["genomes"]
    if not isinstance(genomes, (list, tuple)) or len(genomes) == 0:
        raise ValueError("CONFIG ERROR: 'genomes' must be a non-empty list")
    gset = set(genomes)
    if config["ref_genome"] not in gset:
        warnings.append(
            f"ref_genome '{config['ref_genome']}' is not listed in 'genomes'; "
            "it will be added automatically"
        )

    # 2. Input files that the pipeline actually consumes must exist.
    file_checks = [("hal", config.get("hal")), ("annotation", config.get("annotation"))]
    # Protein set for miniprot/augMP: either a static protein_fasta, or a
    # protein_db block that builds one from UniProt reference proteomes. augMP
    # (part of the augustus path) needs one of the two.
    pdb = config.get("protein_db")
    if pdb is not None:
        if not isinstance(pdb, dict):
            errors.append("'protein_db' must be a mapping (species / taxa / clades / base_fasta / out)")
        else:
            if not (pdb.get("species") or pdb.get("taxa") or pdb.get("clades") or pdb.get("base_fasta")):
                errors.append("'protein_db' must set at least one of species, taxa, clades, or base_fasta")
            if pdb.get("base_fasta"):
                file_checks.append(("protein_db.base_fasta", pdb["base_fasta"]))
    elif config.get("protein_fasta"):
        file_checks.append(("protein_fasta", config["protein_fasta"]))
    elif config.get("augustus", False):
        errors.append("augustus is enabled (augMP aligns proteins with miniprot) but "
                      "neither 'protein_fasta' nor 'protein_db' is set")
    if config.get("augustus", False):
        for key in ("tm_cfg_path", "tmr_cfg_path", "minisplice_model", "minisplice_calibration"):
            if config.get(key):
                file_checks.append((key, config[key]))
        if config.get("augustus_pb", False) and config.get("pb_cfg_path"):
            file_checks.append(("pb_cfg_path", config["pb_cfg_path"]))
    for name, path in file_checks:
        if path and not Path(path).exists():
            errors.append(f"file for '{name}' not found: {path}")

    # 3. Genome subsets referenced elsewhere must be declared in 'genomes'.
    for key in ("rnaseq_genomes", "isoseq_genomes", "stringtie_genomes"):
        for g in (config.get(key) or []):
            if g not in gset:
                errors.append(f"'{key}' lists '{g}', which is not in 'genomes'")

    # 4. transcriptomic_data: entries should map to known genomes and its
    #    referenced BAM files must exist (they are required inputs to rules).
    for g, streams in (config.get("transcriptomic_data") or {}).items():
        if g not in gset:
            warnings.append(f"transcriptomic_data has entry for unknown genome '{g}' (ignored)")
        for btype in ("bam", "intronbam", "isoseq_bam"):
            for bam in ((streams or {}).get(btype) or []):
                if not Path(bam).exists():
                    errors.append(f"transcriptomic_data[{g}].{btype} file not found: {bam}")

    # 5. execution_mode must be one we support ('auto' is resolved later).
    mode = config.get("execution_mode", "auto")
    if mode not in ("slurm", "sge", "local", "auto"):
        errors.append(f"execution_mode must be one of auto/slurm/sge/local, got '{mode}'")

    # 6. Ancestor annotation sanity.
    if config.get("annotate_ancestors", False):
        for m in (config.get("ancestor_modes") or []):
            if m not in ("transMap", "transMap_pairwise", "txTM"):
                errors.append(
                    f"ancestor_modes contains unsupported mode '{m}' "
                    "(allowed: transMap, transMap_pairwise, txTM)"
                )

    for w in warnings:
        print(f"CONFIG WARNING: {w}", file=sys.stderr)
    if errors:
        raise ValueError(
            "CONFIG ERROR: the configuration has "
            f"{len(errors)} problem(s):\n  - " + "\n  - ".join(errors)
        )

# Run validation
validate_config()

# Global variables - use consistent naming
WORK_DIR = Path(config["work_dir"])
TOIL_BASE_JOB_STORE_DIR = WORK_DIR / "toil_job_stores"

# Vendored UCSC/kent binaries (pslMapPostChain, transMapPslToGenePred, etc.) for cluster nodes
# where the submitter's conda PATH is minimal or incomplete.
try:
    CAT2_ROOT = Path(workflow.basedir)
except NameError:
    CAT2_ROOT = Path(config.get("cat2_root", ".")).resolve()
CAT2_STANDALONES = CAT2_ROOT / "standalones"

# Genome sets
ALL_GENOMES = config['genomes'].copy()
REF_GENOME = config['ref_genome']
if REF_GENOME not in ALL_GENOMES:
    ALL_GENOMES.append(REF_GENOME)
TARGET_GENOMES = [g for g in ALL_GENOMES if g != REF_GENOME]
RNASEQ_GENOMES = config.get("rnaseq_genomes", [])
ISOSEQ_GENOMES = config.get("isoseq_genomes", [])
NON_RNASEQ_GENOMES = [g for g in TARGET_GENOMES if g not in RNASEQ_GENOMES]

# Conditional genome sets based on enabled modes
AUGUSTUS_GENOMES = TARGET_GENOMES if config.get("augustus", False) else []
AUGUSTUS_PB_GENOMES = [g for g in ISOSEQ_GENOMES if config.get("augustus_pb", False)]
STRG_GENOMES = [g for g in ISOSEQ_GENOMES if config.get("stringtie", False) and g in config.get("stringtie_genomes", TARGET_GENOMES)]

# --- Protein reference set for miniprot / augMP --------------------------------
# augMP aligns a protein database to every target genome with miniprot. Using only
# reference (e.g. human) proteins means augMP can never find a gene that has no
# reference ortholog. A `protein_db:` block builds a broader, multi-species protein
# set from UniProt reference proteomes (scripts/build_protein_db.py); its output
# then becomes the miniprot input. Without that block the static `protein_fasta:`
# file is used directly. PROTEIN_FASTA is the single source of truth downstream.
_PROTEIN_DB_CFG = config.get("protein_db")
if _PROTEIN_DB_CFG:
    PROTEIN_FASTA = _PROTEIN_DB_CFG.get("out") or f"{config['work_dir']}/protein_db/protein_db.fa"
    BUILD_PROTEIN_DB = True
else:
    PROTEIN_FASTA = config.get("protein_fasta")
    BUILD_PROTEIN_DB = False
# Fall back to a placeholder path so rule definitions never receive a None input.
# This is only ever consumed by the augMP path, which validate_config guarantees
# has a real protein source; otherwise run_miniprot is not part of the DAG.
if not PROTEIN_FASTA:
    PROTEIN_FASTA = f"{config['work_dir']}/protein_db/protein_db.fa"

# Pipeline modes and constraints
VALID_MODES = ["transMap", "transMap_pairwise", "augTM", "augTMR", "augTM_pairwise", "augTMR_pairwise", "augMP", "txTM", "augPB", "strg"]
VALID_ALIGNMENT_MODES = ["transMap", "transMap_pairwise", "augTM", "augTMR", "augTM_pairwise", "augTMR_pairwise", "augMP", "txTM"]

# Active alignment modes (filtered based on config)
ACTIVE_ALIGNMENT_MODES = ["transMap", "transMap_pairwise"]  # transMap and transMap_pairwise are always active
if config.get("augustus", False):
    ACTIVE_ALIGNMENT_MODES.append("augTM")
    ACTIVE_ALIGNMENT_MODES.append("augTM_pairwise")
    if RNASEQ_GENOMES:
        ACTIVE_ALIGNMENT_MODES.append("augTMR")
        ACTIVE_ALIGNMENT_MODES.append("augTMR_pairwise")
    ACTIVE_ALIGNMENT_MODES.append("augMP")  # Augustus MP mode
if config.get("txTM", False):
    ACTIVE_ALIGNMENT_MODES.append("txTM")

# Ancestor genome annotation (internal Cactus HAL nodes)
ANCESTOR_MODES_DEFAULT = ["transMap", "transMap_pairwise", "txTM"]
VALID_ANCESTOR_MODES = ["transMap", "transMap_pairwise", "txTM"]

def _genome_wc(genomes):
    return "|".join(genomes) if genomes else "never_match"

def _resolve_ancestor_genomes():
    if not config.get("annotate_ancestors", False):
        return []
    explicit = config.get("ancestor_genomes")
    if explicit:
        return list(explicit)
    from tools.hal import extract_ancestor_genomes
    target_set = list(set(config["genomes"] + [REF_GENOME]))
    return list(extract_ancestor_genomes(config["hal"], target_set))

ANCESTOR_GENOMES = _resolve_ancestor_genomes()
ANCESTOR_MODES = config.get("ancestor_modes", ANCESTOR_MODES_DEFAULT)
if ANCESTOR_GENOMES:
    invalid = [m for m in ANCESTOR_MODES if m not in VALID_ANCESTOR_MODES]
    if invalid:
        raise ValueError(
            f"Invalid ancestor_modes {invalid}; allowed: {', '.join(VALID_ANCESTOR_MODES)}"
        )
    from tools.hal import list_hal_genomes
    hal_genomes = set(list_hal_genomes(config["hal"]))
    missing = [g for g in ANCESTOR_GENOMES if g not in hal_genomes]
    if missing:
        raise ValueError(f"ancestor_genomes not found in HAL: {', '.join(missing)}")
    for g in ANCESTOR_GENOMES:
        if g not in ALL_GENOMES:
            ALL_GENOMES.append(g)

ANNOTATION_GENOMES = TARGET_GENOMES + ANCESTOR_GENOMES
ANCESTOR_ALIGNMENT_MODES = [m for m in ANCESTOR_MODES if m in VALID_ANCESTOR_MODES]
ALL_ALIGNMENT_MODES = list(dict.fromkeys(ACTIVE_ALIGNMENT_MODES + ANCESTOR_ALIGNMENT_MODES))
TXTM_GENOMES = ([g for g in TARGET_GENOMES if config.get("txTM", False)] +
                [g for g in ANCESTOR_GENOMES if "txTM" in ANCESTOR_MODES])
ANNOTATION_GENOME_WC = _genome_wc(ANNOTATION_GENOMES)
TXTM_GENOME_WC = _genome_wc(TXTM_GENOMES)
# Per-mode genome wildcard patterns (reused across the Augustus rules).
AUGUSTUS_GENOME_WC = _genome_wc(AUGUSTUS_GENOMES)
RNASEQ_GENOME_WC = _genome_wc(RNASEQ_GENOMES)
AUGUSTUS_RNASEQ_GENOME_WC = _genome_wc([g for g in RNASEQ_GENOMES if g in AUGUSTUS_GENOMES])
AUGUSTUS_NON_RNASEQ_GENOME_WC = _genome_wc([g for g in NON_RNASEQ_GENOMES if g in AUGUSTUS_GENOMES])
AUGUSTUS_PB_GENOME_WC = _genome_wc(AUGUSTUS_PB_GENOMES)
STRG_GENOME_WC = _genome_wc(STRG_GENOMES)

# Wildcard constraints for validation
wildcard_constraints:
    genome = f"({'|'.join(TARGET_GENOMES)})",
    mode = f"({'|'.join(VALID_MODES)})",
    alignment_mode = f"({'|'.join([mode for mode in VALID_ALIGNMENT_MODES if mode != 'augTMR' or RNASEQ_GENOMES])})"

# Output file collections
consensus_files = []
for genome in ANNOTATION_GENOMES:
    consensus_files.extend([
        str(WORK_DIR / f"consensus_gene_set/{genome}_consensus.gp"),
        str(WORK_DIR / f"consensus_gene_set/{genome}_consensus.gp_info"), 
        str(WORK_DIR / f"consensus_gene_set/{genome}_consensus.gff3"),
        str(WORK_DIR / f"consensus_gene_set/{genome}_consensus.fasta"),
        str(WORK_DIR / f"consensus_gene_set/{genome}_consensus_protein.fasta")
    ])

# Helper functions for dynamic file generation
def mode_gp_paths(genome, work_dir=None):
    """Canonical genePred output path for every annotation mode.

    Single source of truth for the per-mode ``*.gp`` paths so the various
    consumers (alignment, consensus, single-mode lookup) stay in sync. Note
    that consensus deliberately overrides ``transMap_pairwise`` to the
    *unfiltered* GP; that divergence is applied explicitly at the call site.
    """
    wd = work_dir if work_dir is not None else config['work_dir']
    return {
        'transMap':           f"{wd}/transMap/{genome}_filtered.gp",
        'transMap_pairwise':  f"{wd}/transMap_pairwise/{genome}_filtered.gp",
        'augTM':              f"{wd}/augustus/{genome}_augTM.gp",
        'augTM_pairwise':     f"{wd}/augustus/{genome}_augTM_pairwise.gp",
        'augTMR':             f"{wd}/augustus/{genome}_augTMR.gp",
        'augTMR_pairwise':    f"{wd}/augustus/{genome}_augTMR_pairwise.gp",
        'augMP':              f"{wd}/augustus/{genome}_augMP.gp",
        'txTM':               f"{wd}/txTM/{genome}_txTM.gp",
        'augPB':              f"{wd}/augustus_pb/{genome}_augPB.gp",
        'strg':               f"{wd}/stringtie/{genome}_strg.gp",
    }

def get_stringtie_bams(wildcards):
    """Get BAM files for StringTie processing."""
    genome_config = config.get("transcriptomic_data", {}).get(wildcards.genome, {})
    sr_bams = genome_config.get("intronbam", []) or genome_config.get("bam", [])
    lr_bams = genome_config.get("isoseq_bam", [])
    return {"sr": sr_bams, "lr": lr_bams}

def get_gp_path_for_mode(mode, genome):
    """Get the genePred file path for a specific mode and genome."""
    path_map = {
        'transMap': WORK_DIR / f"transMap/{genome}_filtered.gp",
        'transMap_pairwise': WORK_DIR / f"transMap_pairwise/{genome}_filtered.gp",
        'augTM': WORK_DIR / f"augustus/{genome}_augTM.gp", 
        'augTMR': WORK_DIR / f"augustus/{genome}_augTMR.gp",
        'augTM_pairwise': WORK_DIR / f"augustus/{genome}_augTM_pairwise.gp",
        'augTMR_pairwise': WORK_DIR / f"augustus/{genome}_augTMR_pairwise.gp",
        'augMP': WORK_DIR / f"augustus/{genome}_augMP.gp",
        'txTM': WORK_DIR / f"txTM/{genome}_txTM.gp",
        'augPB': WORK_DIR / f"augustus_pb/{genome}_augPB.gp",
        'strg': WORK_DIR / f"stringtie/{genome}_strg.gp"
    }
    return str(path_map[mode])

# Rule execution order
ruleorder: prepare_reference_files > prepare_genome_files
ruleorder: augustus_run_tm_and_tmr > augustus_run_tm_only
ruleorder: augustus_run_tm_pairwise_and_tmr_pairwise > augustus_run_tm_pairwise_only
ruleorder: generate_hints > run_transcript_map

# ─── Execution mode ───────────────────────────────────────────────────────────
# `auto` (the default) detects the scheduler from the environment (sbatch ->
# slurm, qsub+$SGE_ROOT -> sge, otherwise local) so the pipeline runs on any
# system without editing config. Resolve it once here to a concrete backend so
# every rule and child script sees the same value.
from cat.scheduler import get_scheduler, resolve_execution_mode
EXECUTION_MODE = resolve_execution_mode(config.get("execution_mode", "auto"))
if EXECUTION_MODE not in ("slurm", "sge", "local"):
    raise ValueError(f"execution_mode must be 'slurm', 'sge', 'local', or 'auto', got '{EXECUTION_MODE}'")
IS_CLUSTER = (EXECUTION_MODE in ("slurm", "sge"))

SCHEDULER = get_scheduler(EXECUTION_MODE, config)

# ─── Resource resolution ──────────────────────────────────────────────────────
# Fallback defaults when a key is absent from config['slurm']['rules'][rule].
_RULE_DEFAULTS = {
    "prepare_genome_files":      {"mem": "128G", "cpus": 64,  "time": "01:00:00"},
    "prepare_reference_files":   {"mem": "64G",  "cpus": 2,   "time": "01:00:00"},
    "transmap_map_psl":          {"mem": "16G",  "cpus": 4,   "time": "01:00:00"},
    "transmap_unfiltered_gtf":   {"mem": "16G",  "cpus": 1,   "time": "01:00:00"},
    "minimap2_bam":              {"mem": "64G",  "cpus": 32,  "time": "12:00:00"},
    "bam_to_chain":              {"mem": "16G",  "cpus": 4,   "time": "01:00:00"},
    "transmap_pairwise_map_psl": {"mem": "16G",  "cpus": 4,   "time": "01:00:00"},
    "run_chaining_per_genome":   {"mem": "128G", "cpus": 64,  "time": "12:00:00",  "timeout_hours": 24},
    "run_miniprot":              {"mem": "128G", "cpus": 64,  "time": "04:00:00",  "timeout_hours": 12,
                                 "minisplice_parallel": True, "minisplice_cpus": 4, "minisplice_mem": "8G",
                                 "minisplice_time": "01:00:00", "minisplice_max_concurrent": 25},
    "run_txTM":               {"mem": "128G", "cpus": 64,  "time": "04:00:00", "timeout_hours": 12},
    "stringtie_merge":           {"mem": "16G",  "cpus": 8,   "time": "04:00:00"},
    "stringtie_sort":            {"mem": "128G", "cpus": 64,  "time": "06:00:00"},
    "stringtie_run":             {"mem": "128G", "cpus": 64,  "time": "12:00:00"},
    "stringtie_convert":         {"mem": "8G",   "cpus": 2,   "time": "02:00:00"},
    "stringtie_gp":              {"mem": "8G",   "cpus": 2,   "time": "02:00:00"},
    "generate_consensus":        {"mem": "256G", "cpus": 32,  "time": "12:00:00"},
    "annotate_novel_genes":      {"mem": "32G",  "cpus": 16,  "time": "02:00:00"},
    "generate_hints":            {"mem": "256G", "cpus": 128, "time": "12:00:00",  "max_concurrent_jobs": 50},
    "align_transcripts":         {"mem": "64G",  "cpus": 64,  "time": "02:00:00",  "max_concurrent_jobs": 200, "timeout_hours": 12, "chunk_size": 500},
    "evaluate_transcripts":      {"mem": "16G",  "cpus": 1,   "time": "01:00:00",  "max_concurrent_jobs": 20,  "chunk_size": 500},
    "find_denovo_parents":       {"mem": "32G",  "cpus": 1,   "time": "01:00:00",  "max_concurrent_jobs": 20},
}

_LOCAL_DEFAULTS = {
    "prepare_genome_files":    {"threads": 64,  "mem_gb": 128, "time_h": 4},
    "prepare_reference_files": {"threads": 2,   "mem_gb": 64,  "time_h": 2},
    "transmap_map_psl":        {"threads": 4,   "mem_gb": 16,  "time_h": 2},
    "transmap_unfiltered_gtf": {"threads": 1,   "mem_gb": 16,  "time_h": 1},
    "minimap2_bam":            {"threads": 32,  "mem_gb": 64,  "time_h": 12},
    "bam_to_chain":            {"threads": 4,   "mem_gb": 16,  "time_h": 2},
    "transmap_pairwise_map_psl": {"threads": 4, "mem_gb": 16,  "time_h": 2},
    "run_chaining_per_genome": {"threads": 64,  "mem_gb": 128, "time_h": 24},
    "run_miniprot":            {"threads": 64,  "mem_gb": 128, "time_h": 4,
                                 "minisplice_parallel": True, "minisplice_cpus": 4, "minisplice_max_jobs": 8},
    "run_txTM":             {"threads": 64,  "mem_gb": 128, "time_h": 8},
    "stringtie_merge":         {"threads": 8,   "mem_gb": 16,  "time_h": 4},
    "stringtie_sort":          {"threads": 64,  "mem_gb": 128, "time_h": 6},
    "stringtie_run":           {"threads": 64,  "mem_gb": 128, "time_h": 12},
    "stringtie_convert":       {"threads": 2,   "mem_gb": 8,   "time_h": 2},
    "stringtie_gp":            {"threads": 2,   "mem_gb": 8,   "time_h": 2},
    "generate_consensus":      {"threads": 32,  "mem_gb": 128, "time_h": 12},
    "annotate_novel_genes":    {"threads": 16,  "mem_gb": 32,  "time_h": 2},
    "generate_hints":          {"threads": 128, "mem_gb": 256, "time_h": 12},
    "align_transcripts":       {"threads": 64,  "mem_gb": 128, "time_h": 24},
    "evaluate_transcripts":    {"threads": 4,   "mem_gb": 16,  "time_h": 24},
    "find_denovo_parents":               {"threads": 4,  "mem_gb": 32,  "time_h": 24},
    "augustus_run_tm_and_tmr":           {"threads": 64, "mem_gb": 128, "time_h": 48},
    "augustus_run_tm_only":              {"threads": 64, "mem_gb": 128, "time_h": 24},
    "augustus_run_tm_pairwise_and_tmr_pairwise": {"threads": 64, "mem_gb": 128, "time_h": 48},
    "augustus_run_tm_pairwise_only":     {"threads": 64, "mem_gb": 128, "time_h": 24},
    "augustus_run_mp":                   {"threads": 64, "mem_gb": 128, "time_h": 24},
    "run_augustus_pb":                   {"threads": 64, "mem_gb": 128, "time_h": 24},
}

def get_res(rule_name, key):
    """Return a resource value for a rule. Config overrides built-in defaults."""
    cfg = config.get("slurm", {}).get("rules", {}).get(rule_name, {})
    if key in cfg:
        return cfg[key]
    return _RULE_DEFAULTS[rule_name][key]

def get_local_res(rule_name, key):
    """Return local-mode resource for Snakemake scheduling (threads/mem_gb/time_h)."""
    cfg = config.get("local", {}).get("rules", {}).get(rule_name, {})
    if key in cfg:
        return cfg[key]
    return _LOCAL_DEFAULTS.get(rule_name, {}).get(key, 1)

# ─── miniprot mapping sensitivity ─────────────────────────────────────────────
# Defaults are deliberately more permissive than miniprot's own defaults so the
# augMP path recovers more paralogous / divergent gene copies. Every knob can be
# overridden under config['miniprot']; the converter (convert_miniprot_to_genepred.py)
# keeps each PAF row with its real metrics so downstream consensus filtering stays
# honest. miniprot defaults shown in brackets for reference.
_MINIPROT_MAP_DEFAULTS = {
    "splice_model":        2,        # -j  vertebrate/insect splice model            [1]
    "max_intron":          "auto",   # "auto" => -I (3.6*sqrt(refLen)); or e.g. "500k" => -G500k
    "min_secondary_ratio": 0.2,      # -p  min secondary-to-primary score ratio      [0.7]
    "max_secondary":       100,      # -N  consider at most N secondary alignments   [30]
    "out_n":               100,      # --outn  max alignments emitted per query       [1000]
    "out_s":               0.5,      # --outs  emit if score >= this * bestScore      [0.99]
    "out_c":               0.1,      # --outc  emit if this fraction of query aligned [0.1]
    "extra_flags":         "",       # escape hatch for any additional miniprot flags
}

# ─── high-recall preset ───────────────────────────────────────────────────────
# Opt-in (config['high_recall']: true). When enabled, the recall-limiting gates
# across every mode are loosened in one place so fewer genes are dropped, at the
# expense of precision / more false positives. When disabled (default) behaviour
# is unchanged. These overrides take precedence over per-key config values, since
# the whole point of the preset is to relax them; leave high_recall off and set
# individual keys if you want fine-grained control instead.
HIGH_RECALL = bool(config.get("high_recall", False))
_HIGH_RECALL_OVERRIDES = {
    # transMap filtering (filter_transmap.py)
    "tm_global_near_best":              0.30,   # keep alignments within 30% of best (vs 0.1)
    "tm_filter_overlapping":            False,  # do not collapse overlapping paralogous genes
    "tm_min_cover":                     0.0,    # pslCDnaFilter -minCover
    "tm_min_span":                      0.0,    # pslCDnaFilter -minSpan
    "tm_max_ref_span":                  20,     # allow larger target spans vs reference
    "tm_paralog_rescue_min_coverage":  0.20,    # rescue paralog-collapse victims down to 20% cov
    # txTM consensus gating
    "txTM_min_coverage":                0.0,
    "txTM_strict_metrics":              False,
    "txTM_min_coverage_no_transmap":    0,
    "txTM_strict_metrics_no_transmap":  False,
    "txTM_transmap_anchor_overlap":     0.0,
    "rescue_min_txTM_coverage":         0,
    "rescue_min_txTM_coverage_noncoding": 0,
    "min_nc_len_ratio_txTM_only_rescue": 0.0,
    # augMP consensus gating
    "augMP_min_coverage_no_anchor":     0,
    "augMP_strict_metrics_no_anchor":   False,
    # augMP pre-consensus redundancy filter (filter_augMP.py): loosen the
    # per-locus backstop so recall-mode keeps more paralogs in dense arrays.
    # (Structure/locus de-dup itself never drops a distinct copy.)
    "augMP_filter_max_models_per_locus": 60,
    # generic consensus gates
    "min_pc_len_ratio_vs_reference":    0.0,
    "min_pc_len_ratio_txTM_only_rescue": 0.0,
    "cnv_score_similarity":             0.0,
    "consensus_fragment_max_coverage":  0.0,    # 0 disables fragment reclassification
    "consensus_fragment_max_identity":  0.0,
    # protein-only novel genes (augMP models with no reference transcript)
    "keep_protein_only_novel":          True,   # retain + cleanly label lineage-specific protein-only genes
    "protein_novel_min_exons":          1,      # keep single-exon novel too (recall over precision)
    "protein_novel_min_cds_aa":         0,      # no minimum ORF length
    "protein_novel_min_coverage":       0.0,    # no extra coverage gate
    "protein_novel_min_identity":       0.0,    # no extra identity gate
    "protein_novel_keep_overlapping":   True,   # keep novel overlapping projected features too
    # expressed non-coding -> protein_coding rescue (recall). On by default; in
    # high_recall also allow single-exon and shorter ORFs to maximise recovery.
    "rescue_expressed_noncoding_to_pc": True,
    "rescue_expressed_allow_single_exon": True,
    "rescue_expressed_min_cds_aa":      50,
    # consensus postprocess (kept on, but its low-support drop rule is disabled)
    "postprocess_low_support_fraction": 0.0,
    # denovo (augPB / strg)
    "denovo_allow_unsupported":         True,
}

def rcfg(key, default):
    """Recall-aware ``config.get``: honours the high_recall preset overrides."""
    if HIGH_RECALL and key in _HIGH_RECALL_OVERRIDES:
        return _HIGH_RECALL_OVERRIDES[key]
    return config.get(key, default)

def get_miniprot_opt(key):
    """Return a miniprot mapping option, config override first then default."""
    return config.get("miniprot", {}).get(key, _MINIPROT_MAP_DEFAULTS[key])

def build_miniprot_map_flags(cpus):
    """Assemble the shared miniprot mapping flags (used for both PAF and GTF passes).

    Excludes ``--spsc``/``--gtf``/index/query, which the caller appends.
    """
    max_intron = str(get_miniprot_opt("max_intron")).lower()
    intron_flag = "-I" if max_intron in ("auto", "i", "") else f"-G{max_intron}"
    flags = [
        intron_flag,
        "-u",
        f"-t{cpus}",
        f"-j{get_miniprot_opt('splice_model')}",
        f"-p{get_miniprot_opt('min_secondary_ratio')}",
        f"-N{get_miniprot_opt('max_secondary')}",
        f"--outn={get_miniprot_opt('out_n')}",
        f"--outs={get_miniprot_opt('out_s')}",
        f"--outc={get_miniprot_opt('out_c')}",
    ]
    extra = str(get_miniprot_opt("extra_flags")).strip()
    if extra:
        flags.append(extra)
    return " ".join(flags)

def _slurm_partition(rule_name=None):
    if rule_name:
        r = config.get("slurm", {}).get("rules", {}).get(rule_name, {})
        if "partition" in r:
            return r["partition"]
    return config.get("slurm", {}).get("partition", "high_priority")

def _slurm_exclude():
    return config.get("slurm", {}).get("exclude_nodes", "")

def _aug_slurm_args(rule_key):
    """Build cluster CLI args for augustus_parallel.py / augustus_pb_parallel.py.
    """
    cfg = config.get("slurm", {}).get("rules", {}).get(rule_key, {})
    exclude = _slurm_exclude()
    arg_map = {
        "preprocessing_partition": "--slurm_partition",
        "jobs_partition":          "--slurm_jobs_partition",
        "transcripts_partition":   "--slurm_transcripts_partition",
        "hints_mem":               "--slurm_hints_mem",
        "jobs_mem":                "--slurm_jobs_mem",
        "transcripts_mem":         "--slurm_transcripts_mem",
        "setup_time":              "--slurm_setup_time",
        "hints_time":              "--slurm_hints_time",
        "jobs_time":               "--slurm_jobs_time",
        "transcripts_time":        "--slurm_transcripts_time",
        "hints_concurrency":       "--slurm_hints_concurrency",
        "jobs_concurrency":        "--slurm_jobs_concurrency",
        "transcripts_concurrency": "--slurm_transcripts_concurrency",
    }
    # augustus_pb_parallel.py doesn't have a transcripts stage; strip those keys.
    if rule_key == "augustus_pb":
        arg_map = {k: v for k, v in arg_map.items() if not k.startswith("transcripts_")}
    args = []
    for cfg_key, arg_name in arg_map.items():
        if cfg_key in cfg:
            args += [arg_name, str(cfg[cfg_key])]
    if exclude:
        args += ["--slurm_exclude_nodes", exclude]
    args += ["--execution_mode", EXECUTION_MODE]
    module_load = config.get("slurm", {}).get("module_load", "") or ""
    if module_load:
        args += ["--module_load", module_load]
    sge_cfg = (config.get("cluster", {}) or {}).get("sge", {}) or {}
    if sge_cfg.get("parallel_env"):
        args += ["--sge_parallel_env", str(sge_cfg["parallel_env"])]
    if sge_cfg.get("memory_flag"):
        args += ["--sge_memory_flag", str(sge_cfg["memory_flag"])]
    return args

def _slurm_module():
    return config.get("slurm", {}).get("module_load", "")

def run_augustus_parallel(*, input, output, params, wildcards, log_path,
                          genome_work_dir, with_tmr):
    """Shared driver for the augTM / augTMR (and their pairwise) ``run:`` blocks.

    Builds and executes the ``augustus_parallel.py`` command, cleans up the
    per-genome temp dir, and touches the done markers. ``with_tmr`` toggles the
    RNA-seq (TMR) arguments and the extra done file.
    """
    os.makedirs(genome_work_dir, exist_ok=True)

    cmd = [
        "python", input.script,
        "--genome_fasta", input.fasta,
        "--coding_gp", input.coding_gp,
        "--filtered_tm_psl", input.filtered_tm_psl,
        "--ref_psl", input.ref_psl,
        "--annotation_gp", input.annotation_gp,
        "--tm_cfg", params.tm_cfg,
        "--genome", wildcards.genome,
        "--augustus_species", params.species,
        "--utr", str(params.utr),
        "--augustus_tm_gtf", output.tm_gtf,
    ]
    if with_tmr:
        cmd += [
            "--augustus_tmr_gtf", output.tmr_gtf,
            "--augustus_hints_db", input.hints_db,
            "--tmr_cfg", params.tmr_cfg,
        ]
    cmd += ["--work_dir", genome_work_dir]
    cmd += (["--no_slurm_preprocessing", "--no_slurm_transcripts"]
            if not IS_CLUSTER else _aug_slurm_args("augustus_tm"))

    # try/finally so the (potentially large) per-genome temp dir is always removed,
    # even when Augustus fails — otherwise failed/retried jobs leak temp dirs.
    try:
        shell(" ".join(cmd) + f" > {log_path} 2>&1")
    finally:
        if os.path.exists(genome_work_dir):
            shell(f"rm -rf {genome_work_dir}")

    done_files = [output.tm_gtf_done] + ([output.tmr_gtf_done] if with_tmr else [])
    shell("touch " + " ".join(str(d) for d in done_files))

# Shared body for the Augustus GTF->GenePred conversion rules. An empty GTF
# (a genome with no Augustus predictions for that mode) yields an empty GenePred.
AUGUSTUS_GTF_TO_GP_SHELL = r"""
if [ -s {input.gtf} ]; then
    gtfToGenePred -genePredExt {input.gtf} {output.gp} &> {log}
else
    echo "Empty Augustus GTF - creating empty GenePred" > {log}
    touch {output.gp}
fi
"""

def build_sbatch_header(rule_name, job_name, log_out, log_err):
    """Backend-agnostic scheduler header for a job script.
    """
    # Append standalones after the submitter PATH (same precedence as cat/transcript_map.py).
    path_line = f'export PATH="$PATH:{CAT2_STANDALONES}"\n'
    return SCHEDULER.header(
        job_name=job_name,
        cpus=get_res(rule_name, 'cpus'),
        mem=get_res(rule_name, 'mem'),
        walltime=get_res(rule_name, 'time'),
        log_out=log_out,
        log_err=log_err,
        partition=_slurm_partition(rule_name),
        queue=_slurm_partition(rule_name),
    ) + path_line + f"cd {CAT2_ROOT}\n" + "set -euo pipefail\n"

def build_minisplice_step(
    *,
    work_dir,
    genome,
    genome_fasta,
    chrom_sizes,
    minisplice_model,
    minisplice_calibration,
    splice_scores_out,
    log_path,
    use_slurm_array,
):
    """Return bash for minisplice predict, optionally one SLURM task per chromosome."""
    parallel = get_res("run_miniprot", "minisplice_parallel")
    chrom_list = f"{work_dir}/miniprot/{genome}_minisplice_chroms.txt"
    chrom_fa_dir = f"{work_dir}/miniprot/{genome}_minisplice_chrom_fasta"
    chrom_scores_dir = f"{work_dir}/miniprot/{genome}_minisplice_by_chrom"
    array_script = f"{work_dir}/miniprot/{genome}_minisplice_array.sh"

    with open(chrom_sizes) as chrom_sizes_f:
        chroms = [line.split()[0] for line in chrom_sizes_f if line.strip()]
    if not chroms:
        raise ValueError(f"No chromosomes found in {chrom_sizes}")

    with open(chrom_list, "w") as chrom_list_f:
        chrom_list_f.write("\n".join(chroms) + "\n")

    num_chroms = len(chroms)
    merge_scores = f"""
# Merge per-chromosome splice scores (order matches chrom.sizes)
: > {splice_scores_out}
while read -r chrom _rest; do
    [[ -z "$chrom" ]] && continue
    chrom_tsv="{chrom_scores_dir}/${{chrom}}.tsv"
    if [[ ! -s "$chrom_tsv" ]]; then
        echo "ERROR: missing minisplice scores for $chrom" >> {log_path}
        exit 1
    fi
    cat "$chrom_tsv" >> {splice_scores_out}
done < {chrom_sizes}
echo "minisplice predict completed ({num_chroms} chromosomes merged)" >> {log_path}
"""

    if not parallel:
        return f"""
# Step 0: Run minisplice on the full genome
echo "Running minisplice predict..." >> {log_path}
minisplice predict -t {get_res('run_miniprot', 'cpus')} -c {minisplice_calibration} {minisplice_model} {genome_fasta} > {splice_scores_out} 2>> {log_path}
echo "minisplice predict completed" >> {log_path}
"""

    minisplice_cpus = get_res("run_miniprot", "minisplice_cpus")
    per_chrom_body = f"""
set -euo pipefail
CHROM=$(sed -n "${{SLURM_ARRAY_TASK_ID:-$1}}p" {chrom_list})
if [[ -z "$CHROM" ]]; then
    echo "ERROR: no chromosome for task id ${{SLURM_ARRAY_TASK_ID:-$1}}" >&2
    exit 1
fi
mkdir -p {chrom_fa_dir} {chrom_scores_dir}
CHROM_FA="{chrom_fa_dir}/${{CHROM}}.fa"
CHROM_TSV="{chrom_scores_dir}/${{CHROM}}.tsv"
samtools faidx {genome_fasta} "$CHROM" > "$CHROM_FA"
minisplice predict -t {minisplice_cpus} -c {minisplice_calibration} {minisplice_model} "$CHROM_FA" > "$CHROM_TSV"
"""

    if use_slurm_array:
        array_header = SCHEDULER.header(
            job_name=f"minisplice-{genome}",
            cpus=minisplice_cpus,
            mem=get_res("run_miniprot", "minisplice_mem"),
            walltime=get_res("run_miniprot", "minisplice_time"),
            log_out=f"{work_dir}/logs/miniprot/{genome}_minisplice_%a_slurm.out",
            log_err=f"{work_dir}/logs/miniprot/{genome}_minisplice_%a_slurm.err",
            partition=_slurm_partition("run_miniprot"),
            array=(1, num_chroms),
            max_concurrent=get_res("run_miniprot", "minisplice_max_concurrent"),
        )
        with open(array_script, "w") as array_script_f:
            array_script_f.write(array_header + per_chrom_body)

        return f"""
# Step 0: Run minisplice predict per chromosome (SLURM array)
echo "Running minisplice predict on {num_chroms} chromosomes..." >> {log_path}
mkdir -p {chrom_fa_dir} {chrom_scores_dir}
ARRAY_JOB=$(sbatch --parsable {array_script})
echo "Submitted minisplice array job $ARRAY_JOB ({num_chroms} tasks)" >> {log_path}
while squeue -h -j "$ARRAY_JOB" 2>/dev/null | grep -q .; do
    sleep 30
done
{merge_scores}
"""

    local_jobs = get_local_res("run_miniprot", "minisplice_max_jobs")
    return f"""
# Step 0: Run minisplice predict per chromosome (local parallel)
echo "Running minisplice predict on {num_chroms} chromosomes..." >> {log_path}
mkdir -p {chrom_fa_dir} {chrom_scores_dir}
run_minisplice_chrom() {{
    local chrom="$1"
    local chrom_fa="{chrom_fa_dir}/${{chrom}}.fa"
    local chrom_tsv="{chrom_scores_dir}/${{chrom}}.tsv"
    samtools faidx {genome_fasta} "$chrom" > "$chrom_fa"
    minisplice predict -t {minisplice_cpus} -c {minisplice_calibration} {minisplice_model} "$chrom_fa" > "$chrom_tsv"
}}
export -f run_minisplice_chrom
if command -v parallel >/dev/null 2>&1; then
    parallel -j {local_jobs} run_minisplice_chrom :::: {chrom_list} >> {log_path} 2>&1
else
    while read -r chrom; do
        [[ -z "$chrom" ]] && continue
        run_minisplice_chrom "$chrom"
    done < {chrom_list}
fi
{merge_scores}
"""

def run_or_submit(script_body, job_script_path, outputs_to_check,
                  log_handle, rule_name, max_wait_s=14400, check_interval_s=30):
    """Write a job script and run it via the active backend.
    """
    import subprocess as _sp
    import time as _t
    import os as _os

    if IS_CLUSTER:
        SCHEDULER.write_script(script_body, job_script_path)
        job_id = SCHEDULER.submit(job_script_path)
        log_handle.write(f"Submitted: {SCHEDULER.name} job {job_id}\n")
        elapsed = 0
        while elapsed < max_wait_s:
            if not SCHEDULER.job_present(job_id):
                result = SCHEDULER.verify_completed(job_id)
                if not result.ok:
                    raise RuntimeError(
                        f"[{rule_name}] cluster job {job_id} failed: {result.detail}"
                    )
                missing = [p for p in outputs_to_check if not _os.path.exists(str(p))]
                if missing:
                    raise RuntimeError(
                        f"[{rule_name}] cluster job {job_id} finished but outputs "
                        f"missing: {', '.join(str(p) for p in missing)}"
                    )
                log_handle.write(f"[{rule_name}] completed after {elapsed}s\n")
                return
            _t.sleep(check_interval_s)
            elapsed += check_interval_s
        raise RuntimeError(f"[{rule_name}] timed out after {max_wait_s}s")
    else:
        clean = "\n".join(
            line for line in script_body.splitlines()
            if not line.startswith("#SBATCH") and not line.startswith("#$") and line != "#!/bin/bash"
        )
        log_handle.write(f"[{rule_name}] running locally\n")
        result = _sp.run(['bash', '-euo', 'pipefail', '-c', clean])
        if result.returncode != 0:
            raise RuntimeError(f"[{rule_name}] local run failed (rc={result.returncode})")

# ─── localrules ───────────────────────────────────────────────────────────────
# In SLURM mode, the sbatch-submitting rules are lightweight wrappers (use
# minimal controller resources). In local mode they do real work but are still
# declared local so Snakemake doesn't double-count their threads — actual
# concurrency is gated by threads: / mem_gb: in each rule's resources block.
localrules: setup_pipeline_directories, build_db, prepare_genome_files, \
    prepare_reference_files, transmap_map_psl, run_miniprot, run_transcript_map, \
    run_chaining_per_genome, transmap_unfiltered_gtf, minimap2_bam, \
    bam_to_chain, transmap_pairwise_map_psl, augustus_run_tm_and_tmr, \
    augustus_run_tm_only, augustus_run_tm_pairwise_and_tmr_pairwise, \
    augustus_run_tm_pairwise_only, run_augustus_pb, \
    stringtie_merge_bams, stringtie_sort_bams, stringtie_run, \
    convert_stringtie_to_strg, convert_strg_gtf_to_gp

rule setup_pipeline_directories:
    """Create all necessary directories for the pipeline."""
    output:
        setup_done = str(WORK_DIR / ".setup_done")
    params:
        work_dir = WORK_DIR,
        toil_dir = TOIL_BASE_JOB_STORE_DIR
    resources:
        mem_gb=2,
        time_h=1,
        job_id=lambda wildcards, attempt: f"setup-{attempt}"
    run:
        # Define directory structure
        directories = [
            "genome_files", "reference", "hints_database", "chaining",
            "miniprot", "stringtie", "txTM", "transMap", "augustus",
            "augustus_pb", "databases", "transcript_alignment", 
            "consensus_gene_set", "chaining_bam", "transMap_pairwise",
            "plots"
        ]
        
        log_directories = [
            "prepare_genome_files", "prepare_reference_files", "generate_hints",
            "build_db", "chaining", "stringtie_merge", "stringtie_sort", 
            "stringtie_run",
            "miniprot", "txTM", "transcript_map", "transmap_map",
            "transmap_filter", "transmap_evaluate", "transmap_gtf",
            "augustus_extract_coding_gp", "augustus_run", "augustus_convert",
            "augustus_pb", "find_denovo_parents", "fix_augmp_gene_names",
            "isoseq_structures", "align_transcripts", "evaluate_transcripts",
            "consensus", "minimap2_bam", "bam_to_chain", "transmap_pairwise_map",
            "transmap_pairwise_filter"
        ]
        
        # Create main directories
        for directory in directories:
            Path(params.work_dir, directory).mkdir(parents=True, exist_ok=True)
            
        # Create log directories  
        for log_dir in log_directories:
            Path(params.work_dir, "logs", log_dir).mkdir(parents=True, exist_ok=True)
            
        # Create toil job store directory
        Path(params.toil_dir).mkdir(parents=True, exist_ok=True)
        
        # Touch completion marker
        Path(output.setup_done).touch()

rule all:
    """Main target rule that triggers the complete pipeline.

    Marked default_target so a bare `snakemake` invocation runs the whole pipeline
    even though `setup_pipeline_directories` is defined earlier in the file (in
    Snakemake the first-defined rule is otherwise the default target).
    """
    default_target: True
    input:
        pipeline_complete = str(WORK_DIR / "pipeline.complete.done"),
        setup_done = rules.setup_pipeline_directories.output.setup_done,
        consensus_files = consensus_files

rule all_with_cleanup:
    """Complete pipeline with cleanup of log files and done markers."""
    input:
        cleanup_done = str(WORK_DIR / "cleanup.complete")

rule prepare_genome_files:
    """
    Submits genome file preparation as a SLURM job to extract genome sequences 
    from HAL and create required indices.
    """
    input:
        hal = config['hal']
    output:
        fasta = WORK_DIR / "genome_files/{genome}.fa",
        two_bit = WORK_DIR / "genome_files/{genome}.2bit", 
        sizes = WORK_DIR / "genome_files/{genome}.chrom.sizes",
        fasta_index = WORK_DIR / "genome_files/{genome}.fa.fai",
        # The miniprot index is only needed by the augMP path; skip it otherwise.
        **({"protein_index": WORK_DIR / "genome_files/{genome}.mpi"}
           if config.get("augustus", False) else {}),
    wildcard_constraints:
        genome = f"({'|'.join(ALL_GENOMES)})"
    threads: 1 if IS_CLUSTER else get_local_res("prepare_genome_files", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("prepare_genome_files", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("prepare_genome_files", "time_h"),
        job_id=lambda wildcards, attempt: f"prep-genome-submit-{wildcards.genome}-{attempt}"
    log:
        WORK_DIR / "logs/prepare_genome_files/{genome}.log"
    run:
        import os
        import subprocess
        import time

        # Get paths and variables
        work_dir = config['work_dir']
        genome = wildcards.genome
        job_script = f"{work_dir}/genome_files/{genome}_prepare_job.sh"

        # The miniprot index is only needed by the augMP path; skip it otherwise.
        if config.get("augustus", False):
            miniprot_index_cmd = (
                "\n# Create miniprot index (augMP path only)\n"
                f"miniprot -t{get_res('prepare_genome_files', 'cpus')} "
                f"-d {output.protein_index} {output.fasta} 2>> {log[0]}\n"
            )
        else:
            miniprot_index_cmd = ""

        script_content = build_sbatch_header(
            "prepare_genome_files",
            f"prep-genome-{genome}",
            f"{work_dir}/logs/prepare_genome_files/{genome}_slurm.out",
            f"{work_dir}/logs/prepare_genome_files/{genome}_slurm.err"
        ) + f"""
# Log start time
echo "Starting genome preparation for: {genome}"
echo "Job ID: ${{SLURM_JOB_ID}}"
echo "Node: $(hostname)"
echo "Start time: $(date)"

# Extract FASTA from HAL
hal2fasta {input.hal} {genome} > {output.fasta} 2>> {log[0]}

# Create 2bit file
faToTwoBit {output.fasta} {output.two_bit} 2>> {log[0]}

# Generate chromosome sizes
faSize -detailed {output.fasta} > {output.sizes} 2>> {log[0]}

# Index FASTA for pysam
samtools faidx {output.fasta} 2>> {log[0]}
{miniprot_index_cmd}
# Log completion
echo "Completed genome preparation for: {genome}"
echo "End time: $(date)"
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting genome preparation job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.fasta, output.two_bit, output.sizes, output.fasta_index]
                + ([output.protein_index] if config.get("augustus", False) else []),
                log_file,
                "prepare_genome_files",
                max_wait_s=14400
            )

rule init_target_genome_database:
    """Create empty SQLite placeholder for target genomes (ref DB comes from prepare_reference_files)."""
    input:
        fasta=WORK_DIR / "genome_files/{genome}.fa",
    output:
        genome_db=WORK_DIR / "databases/{genome}.db",
    wildcard_constraints:
        genome=f"({ANNOTATION_GENOME_WC})",
    resources:
        mem_gb=1,
        time_h=1,
        job_id=lambda wildcards, attempt: f"init-db-{wildcards.genome}-{attempt}",
    shell:
        "touch {output.genome_db}"

rule prepare_reference_files:
    """
    Submits reference file preparation as a SLURM job to process reference 
    annotation into various formats and create databases.
    """
    input:
        gff3 = config["annotation"],
        ref_fasta = WORK_DIR / f"genome_files/{REF_GENOME}.fa",
        ref_sizes = WORK_DIR / f"genome_files/{REF_GENOME}.chrom.sizes"
    output:
        gp = WORK_DIR / f"reference/{REF_GENOME}.gp",
        attrs = WORK_DIR / f"reference/{REF_GENOME}.gp_attrs", 
        gtf = WORK_DIR / f"reference/{REF_GENOME}.gtf",
        bed = WORK_DIR / f"reference/{REF_GENOME}.bed",
        transcript_fasta = WORK_DIR / f"reference/{REF_GENOME}.fa",
        transcript_fasta_index = WORK_DIR / f"reference/{REF_GENOME}.fa.fai",
        psl = WORK_DIR / f"reference/{REF_GENOME}.psl",
        duplicates = WORK_DIR / f"reference/{REF_GENOME}.duplicates.txt",
        db = WORK_DIR / f"databases/{REF_GENOME}.db",
        gff3_db = WORK_DIR / f"reference/{REF_GENOME}.gff3_db",
        db_ready = WORK_DIR / f"databases/{REF_GENOME}.db.ready"
    threads: 1 if IS_CLUSTER else get_local_res("prepare_reference_files", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("prepare_reference_files", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("prepare_reference_files", "time_h"),
        job_id=lambda wildcards, attempt: f"prep-ref-submit-{attempt}"
    log:
        WORK_DIR / f"logs/prepare_reference_files/{REF_GENOME}.log"
    run:
        import os
        import subprocess
        import time

        # Get paths and variables
        work_dir = config['work_dir']
        job_script = f"{work_dir}/reference/{REF_GENOME}_prepare_job.sh"

        script_content = build_sbatch_header(
            "prepare_reference_files",
            f"prep-ref-{REF_GENOME}",
            f"{work_dir}/logs/prepare_reference_files/{REF_GENOME}_slurm.out",
            f"{work_dir}/logs/prepare_reference_files/{REF_GENOME}_slurm.err"
        ) + f"""
# Log start time
echo "Starting reference preparation for: {REF_GENOME}"
echo "Job ID: ${{SLURM_JOB_ID}}"
echo "Node: $(hostname)"
echo "Start time: $(date)"

# Convert GFF3 to GenePred with attributes
gff3ToGenePred -rnaNameAttr=transcript_id -geneNameAttr=gene_id \\
    -attrsOut={output.attrs} {input.gff3} {output.gp} 2>> {log[0]}

# Generate GTF from GenePred
genePredToGtf file {output.gp} -utr -honorCdsStat -source=CAT {output.gtf} 2>> {log[0]}

# Generate BED from GenePred
genePredToBed {output.gp} {output.bed} 2>> {log[0]}

# Extract transcript sequences
gffread {output.gtf} -g {input.ref_fasta} -w {output.transcript_fasta} 2>> {log[0]}

# Index transcript FASTA
samtools faidx {output.transcript_fasta} 2>> {log[0]}

# Generate PSL alignment
genePredToFakePsl -chromSize={input.ref_sizes} noDB {output.gp} {output.psl} /dev/null 2>> {log[0]}

# Find duplicate transcript IDs
( gff3ToGenePred -rnaNameAttr=transcript_id -geneNameAttr=gene_id \\
    -honorStartStopCodons -refseqHacks \\
    -attrsOut=/dev/null {input.gff3} /dev/stdout 2>/dev/null || true ) \\
| awk '{{print $1}}' | sort | uniq -d > {output.duplicates} 2>> {log[0]}

# Build the GFF3 sqlite db FIRST so the augmentation step can read it
python3 - <<EOF 2>> {log[0]}
import gffutils
gffutils.create_db("{input.gff3}", "{output.gff3_db}",
                   merge_strategy="create_unique", force=True)
EOF

python3 cat/augment_reference_for_lifting.py \\
    --gff3-db   {output.gff3_db} \\
    --ref-gp    {output.gp} \\
    --ref-fa    {output.transcript_fasta} \\
    --genome-fa {input.ref_fasta} \\
    --gp-attrs  {output.attrs} 2>> {log[0]}

# Regenerate GTF / BED / PSL so they include the synthetic entries
genePredToGtf file {output.gp} -utr -honorCdsStat -source=CAT {output.gtf} 2>> {log[0]}
genePredToBed {output.gp} {output.bed} 2>> {log[0]}
genePredToFakePsl -chromSize={input.ref_sizes} noDB {output.gp} {output.psl} /dev/null 2>> {log[0]}

# Create annotation database from the augmented GP + gp_attrs
python3 - <<EOF 2>> {log[0]}
from tools.sqlInterface import Annotation
from tools.sqlite import ExclusiveSqlConnection
from tools.gff3 import parse_gff3

# Parse and load annotation data
df = parse_gff3("{output.attrs}", "{output.gp}")
table = Annotation.__tablename__
with ExclusiveSqlConnection("{output.db}") as engine:
    df.to_sql(table, engine, if_exists="replace", index=True)
EOF

# Create completion marker to signal database is ready
touch {output.db_ready}

# Log completion
echo "Completed reference preparation for: {REF_GENOME}"
echo "End time: $(date)"
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting reference preparation job for {REF_GENOME}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.gp, output.attrs, output.gtf, output.bed, output.transcript_fasta,
                 output.transcript_fasta_index, output.psl, output.duplicates, output.db,
                 output.gff3_db, output.db_ready],
                log_file,
                "prepare_reference_files",
                max_wait_s=10800
            )

rule transmap_map_psl:
    """
    Submits transMap PSL mapping as a SLURM job to map reference PSL alignments 
    to target genome using chain files.
    """
    input:
        target_2bit=f"{config['work_dir']}/genome_files/{{genome}}.2bit",
        chain_file=f"{config['work_dir']}/chaining/{{genome}}/{config['ref_genome']}-{{genome}}.chain",
        ref_psl=f"{config['work_dir']}/reference/{config['ref_genome']}.psl",
        ref_fa=f"{config['work_dir']}/reference/{config['ref_genome']}.fa",
        ref_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp"
    wildcard_constraints:
        genome = ANNOTATION_GENOME_WC
    output:
        tm_psl=temp(f"{config['work_dir']}/transMap/{{genome}}.psl"),
        tm_gp=f"{config['work_dir']}/transMap/{{genome}}.gp"
    threads: 1 if IS_CLUSTER else get_local_res("transmap_map_psl", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("transmap_map_psl", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("transmap_map_psl", "time_h"),
        job_id=lambda wildcards, attempt: f"transmap-map-submit-{wildcards.genome}-{attempt}"
    log:
        f"{config['work_dir']}/logs/transmap_map/{{genome}}.log"
    run:
        import os
        import subprocess
        import time

        # Get paths and variables
        work_dir = config['work_dir']
        genome = wildcards.genome
        job_script = f"{work_dir}/transMap/{genome}_map_job.sh"

        script_content = build_sbatch_header(
            "transmap_map_psl",
            f"transmap-{genome}",
            f"{work_dir}/logs/transmap_map/{genome}_slurm.out",
            f"{work_dir}/logs/transmap_map/{genome}_slurm.err"
        ) + f"""
# Log start time
echo "Starting transMap PSL mapping for: {genome}"
echo "Job ID: ${{SLURM_JOB_ID}}"
echo "Node: $(hostname)"
echo "Start time: $(date)"

# Run pslMap, trying both normal and -swapMap in case chain is reversed
(pslMap -chainMapFile {input.ref_psl} {input.chain_file} stdout || \\
 pslMap -swapMap -chainMapFile {input.ref_psl} {input.chain_file} stdout) | \\
pslMapPostChain stdin stdout | \\
sort --parallel={get_res('transmap_map_psl', 'cpus')} -k14,14 -k16,16n | \\
pslRecalcMatch stdin {input.target_2bit} {input.ref_fa} stdout | \\
sort --parallel={get_res('transmap_map_psl', 'cpus')} -k10,10 > {output.tm_psl} 2> {log[0]}

# Convert the raw PSL to GenePred
transMapPslToGenePred -nonCodingGapFillMax=80 -codingGapFillMax=50 \\
    {input.ref_gp} {output.tm_psl} {output.tm_gp} 2>> {log[0]}

# Log completion
echo "Completed transMap PSL mapping for: {genome}"
echo "End time: $(date)"
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting transMap PSL mapping job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.tm_psl, output.tm_gp],
                log_file,
                "transmap_map_psl",
                max_wait_s=7200
            )

rule transmap_filter:
    input:
        tm_psl=f"{config['work_dir']}/transMap/{{genome}}.psl",
        tm_gp=f"{config['work_dir']}/transMap/{{genome}}.gp",
        ref_psl=f"{config['work_dir']}/reference/{config['ref_genome']}.psl",
        ref_db=f"{config['work_dir']}/databases/{config['ref_genome']}.db",
        annotation_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp"
    output:
        filtered_psl=f"{config['work_dir']}/transMap/{{genome}}_filtered.psl",
        filtered_gp=f"{config['work_dir']}/transMap/{{genome}}_filtered.gp",
        metrics_json=f"{config['work_dir']}/transMap/{{genome}}_filter_tm_metrics.json",
        resolved_df=temp(f"{config['work_dir']}/transMap/{{genome}}_resolved_df.pkl")
    wildcard_constraints:
        genome = ANNOTATION_GENOME_WC
    params:
        db_path=f"{config['work_dir']}/databases/{{genome}}.db",
        # Honour input.yaml's `global_near_best` as a fallback (historical key name).
        global_near_best=rcfg("tm_global_near_best", config.get("global_near_best", 0.1)),
        filter_overlapping_genes="--filter-overlapping-genes" if rcfg("tm_filter_overlapping", True) else "",
        overlapping_ignore_bases=config.get("tm_overlapping_ignore_bases", 0),
        min_cover=rcfg("tm_min_cover", 0.1),
        min_span=rcfg("tm_min_span", 0.2),
        max_ref_span=rcfg("tm_max_ref_span", 5),
        paralog_rescue_min_coverage=rcfg("tm_paralog_rescue_min_coverage", 0.5),
        script="cat/filter_transmap.py"
    log:
        f"{config['work_dir']}/logs/transmap_filter/{{genome}}.log"
    resources:
        mem_gb=8,
        time_h=4,
        job_id=lambda wildcards, attempt: f"tm-filter-{wildcards.genome}-{attempt}"
    shell:
        """
        set -euo pipefail
        python {params.script} \
            --tm-psl {input.tm_psl} \
            --ref-psl {input.ref_psl} \
            --tm-gp {input.tm_gp} \
            --db-path {input.ref_db} \
            --filtered-psl {output.filtered_psl} \
            --filtered-gp {output.filtered_gp} \
            --metrics-json {output.metrics_json} \
            --resolved-df {output.resolved_df} \
            --global-near-best {params.global_near_best} \
            --overlapping-ignore-bases {params.overlapping_ignore_bases} \
            --min-cover {params.min_cover} \
            --min-span {params.min_span} \
            --max-ref-span {params.max_ref_span} \
            --paralog-rescue-min-coverage {params.paralog_rescue_min_coverage} \
            {params.filter_overlapping_genes}

        python3 - <<PY
import pandas as pd
from tools.sqlInterface import TmFilterEval
from tools.sqlite import ExclusiveSqlConnection
table = TmFilterEval.__tablename__
genome_db="{params.db_path}"
res_df = pd.read_pickle("{output.resolved_df}")
with ExclusiveSqlConnection(genome_db) as engine:
    res_df.to_sql(table, engine, if_exists="replace")
PY

        transMapPslToGenePred -nonCodingGapFillMax=80 -codingGapFillMax=50 {input.annotation_gp} {output.filtered_psl} {output.filtered_gp}
        """

rule transmap_evaluate:
    """
    Evaluates filtered transMap alignments on various quality metrics
    and loads the results into the genome database.
    """
    input:
        filtered_psl=f"{config['work_dir']}/transMap/{{genome}}_filtered.psl",
        filtered_gp=f"{config['work_dir']}/transMap/{{genome}}_filtered.gp",
        ref_psl=f"{config['work_dir']}/reference/{config['ref_genome']}.psl",
        ref_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp",
        fasta=f"{config['work_dir']}/genome_files/{{genome}}.fa"
    output:
        db_done=temp(f"{config['work_dir']}/databases/{{genome}}.tm_eval.done"),
        resolved_df=temp(f"{config['work_dir']}/transMap/{{genome}}_resolved_df_2.pkl")
    wildcard_constraints:
        genome = ANNOTATION_GENOME_WC
    params:
        db_path=f"{config['work_dir']}/databases/{{genome}}.db",
        script="cat/transmap_classify.py"
    log:
        f"{config['work_dir']}/logs/transmap_evaluate/{{genome}}.log"
    resources:
        mem_gb=8,
        time_h=4,
        job_id=lambda wildcards, attempt: f"tm-eval-{wildcards.genome}-{attempt}"
    shell:
        """
        python {params.script} \
            --filtered-tm-psl {input.filtered_psl} \
            --filtered-tm-gp {input.filtered_gp} \
            --ref-psl {input.ref_psl} \
            --annotation-gp {input.ref_gp} \
            --fasta {input.fasta} \
            --db-path {params.db_path} \
            --resolved-df {output.resolved_df}
        touch {output.db_done}
        python3 - <<PY
import pandas as pd
from tools.sqlInterface import TmEval
from tools.sqlite import ExclusiveSqlConnection
table = TmEval.__tablename__
genome_db="{params.db_path}"
res_df = pd.read_pickle("{output.resolved_df}")
with ExclusiveSqlConnection(genome_db) as engine:
    res_df.to_sql(table, engine, if_exists="replace", index=True)
PY
        """

rule transmap_unfiltered_gtf:
    """
    Submits transMap unfiltered GTF creation as a SLURM job to create a GTF 
    from the original, unfiltered transmap results for visualization or analysis.
    """
    input:
        tm_psl=f"{config['work_dir']}/transMap/{{genome}}.psl",
        ref_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp"
    output:
        gtf=f"{config['work_dir']}/transMap/{{genome}}.gtf"
    wildcard_constraints:
        genome = ANNOTATION_GENOME_WC
    threads: 1 if IS_CLUSTER else get_local_res("transmap_unfiltered_gtf", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("transmap_unfiltered_gtf", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("transmap_unfiltered_gtf", "time_h"),
        job_id=lambda wildcards, attempt: f"tm-gtf-submit-{wildcards.genome}-{attempt}"
    log:
        f"{config['work_dir']}/logs/transmap_gtf/{{genome}}.log"
    run:
        import os
        import subprocess
        import time

        # Get paths and variables
        work_dir = config['work_dir']
        genome = wildcards.genome
        job_script = f"{work_dir}/transMap/{genome}_unfiltered_gtf_job.sh"

        script_content = build_sbatch_header(
            "transmap_unfiltered_gtf",
            f"transmap-gtf-{genome}",
            f"{work_dir}/logs/transmap_gtf/{genome}_slurm.out",
            f"{work_dir}/logs/transmap_gtf/{genome}_slurm.err"
        ) + f"""
# Log start time
echo "Starting transMap unfiltered GTF creation for: {genome}"
echo "Job ID: ${{SLURM_JOB_ID}}"
echo "Node: $(hostname)"
echo "Start time: $(date)"

# Convert PSL to GenePred and then to GTF
transMapPslToGenePred -nonCodingGapFillMax=80 -codingGapFillMax=50 \\
    {input.ref_gp} {input.tm_psl} stdout | \\
genePredToGtf file stdin {output.gtf} >> {log[0]}

# Log completion
echo "Completed transMap unfiltered GTF creation for: {genome}"
echo "End time: $(date)"
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting transMap unfiltered GTF job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.gtf],
                log_file,
                "transmap_unfiltered_gtf",
                max_wait_s=7200
            )

rule minimap2_bam:
    """
    Generates BAM files using minimap2 alignment between target and query genomes.
    This creates alignments that will be used to generate chains.
    """
    input:
        target_fa=f"{config['work_dir']}/genome_files/{{genome}}.fa",
        query_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp",
        query_fa=f"{config['work_dir']}/genome_files/{config['ref_genome']}.fa"
    output:
        bam=f"{config['work_dir']}/chaining_bam/{{genome}}/{config['ref_genome']}-{{genome}}.bam",
        bam_bai=f"{config['work_dir']}/chaining_bam/{{genome}}/{config['ref_genome']}-{{genome}}.bam.bai",
    wildcard_constraints:
        genome = ANNOTATION_GENOME_WC
    params:
        work_dir=config['work_dir'],
        ref_genome=config['ref_genome']
    threads: 1 if IS_CLUSTER else get_local_res("minimap2_bam", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("minimap2_bam", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("minimap2_bam", "time_h"),
        job_id=lambda wildcards, attempt: f"minimap2-bam-{wildcards.genome}-{attempt}"
    log:
        f"{config['work_dir']}/logs/minimap2_bam/{{genome}}.log"
    run:
        import os
        import subprocess
        import time

        work_dir = config['work_dir']
        genome = wildcards.genome
        job_script = f"{work_dir}/chaining_bam/{genome}_minimap2_job.sh"

        # Create output directory
        os.makedirs(f"{work_dir}/chaining_bam/{genome}", exist_ok=True)

        script_content = build_sbatch_header(
            "minimap2_bam",
            f"minimap2-{genome}",
            f"{work_dir}/logs/minimap2_bam/{genome}_slurm.out",
            f"{work_dir}/logs/minimap2_bam/{genome}_slurm.err"
        ) + f"""
echo "Starting minimap2 alignment for: {genome}"
echo "Job ID: ${{SLURM_JOB_ID}}"
echo "Start time: $(date)"

# Run minimap2 to align query (reference) to target genome
minimap2 -ax asm5 -t {get_res('minimap2_bam', 'cpus')} \\
    {input.target_fa} {input.query_fa} | \\
samtools view -bS - | \\
samtools sort -o {output.bam} -

# Index the BAM file
samtools index {output.bam}

echo "Completed minimap2 alignment for: {genome}"
echo "End time: $(date)"
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting minimap2 BAM generation job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.bam, output.bam_bai],
                log_file,
                "minimap2_bam",
                max_wait_s=43200
            )

rule bam_to_chain:
    """
    Converts BAM files to chain files using bamToPsl and axtChain.
    Also generates net files using chainSort and chainNet.
    """
    input:
        bam=f"{config['work_dir']}/chaining_bam/{{genome}}/{config['ref_genome']}-{{genome}}.bam",
        bam_bai=f"{config['work_dir']}/chaining_bam/{{genome}}/{config['ref_genome']}-{{genome}}.bam.bai",
        target_2bit=f"{config['work_dir']}/genome_files/{{genome}}.2bit",
        query_2bit=f"{config['work_dir']}/genome_files/{config['ref_genome']}.2bit",
        target_sizes=f"{config['work_dir']}/genome_files/{{genome}}.chrom.sizes",
        query_sizes=f"{config['work_dir']}/genome_files/{config['ref_genome']}.chrom.sizes"
    output:
        chain=f"{config['work_dir']}/chaining_bam/{{genome}}/{config['ref_genome']}-{{genome}}.chain",
        net=f"{config['work_dir']}/chaining_bam/{{genome}}/{config['ref_genome']}-{{genome}}.net"
    wildcard_constraints:
        genome = ANNOTATION_GENOME_WC
    params:
        work_dir=config['work_dir'],
        ref_genome=config['ref_genome']
    threads: 1 if IS_CLUSTER else get_local_res("bam_to_chain", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("bam_to_chain", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("bam_to_chain", "time_h"),
        job_id=lambda wildcards, attempt: f"bam-to-chain-{wildcards.genome}-{attempt}"
    log:
        f"{config['work_dir']}/logs/bam_to_chain/{{genome}}.log"
    run:
        import os
        import subprocess
        import time

        work_dir = config['work_dir']
        genome = wildcards.genome
        job_script = f"{work_dir}/chaining_bam/{genome}_bam_to_chain_job.sh"

        script_content = build_sbatch_header(
            "bam_to_chain",
            f"bam-to-chain-{genome}",
            f"{work_dir}/logs/bam_to_chain/{genome}_slurm.out",
            f"{work_dir}/logs/bam_to_chain/{genome}_slurm.err"
        ) + f"""
echo "Starting BAM to chain conversion for: {genome}"
echo "Job ID: ${{SLURM_JOB_ID}}"
echo "Start time: $(date)"

# Check if BAM file exists and has alignments
if [ ! -f {input.bam} ]; then
    echo "ERROR: BAM file does not exist: {input.bam}"
    exit 1
fi

BAM_COUNT=$(samtools view -c {input.bam} 2>/dev/null || echo "0")
echo "BAM file has $BAM_COUNT alignments"

if [ "$BAM_COUNT" -eq "0" ]; then
    echo "ERROR: BAM file is empty or has no alignments!"
    exit 1
fi

# Convert BAM to PSL and then to chain
# Note: target_asm = {genome} (target), query_asm = {config['ref_genome']} (reference/query)
# Pass BAM file directly to bamToPsl (matching original bam-to-bigchain script)
bamToPsl -nohead {input.bam} stdout | \\
axtChain -linearGap=medium -psl stdin {input.target_2bit} {input.query_2bit} {output.chain}

# Sort chains and create net (note: -wholeChains removed as per user request)
chainSort {output.chain} stdout | \\
chainNet -inclHap stdin \\
    {input.target_sizes} {input.query_sizes} \\
    {output.net} /dev/null

echo "Completed BAM to chain conversion for: {genome}"
echo "End time: $(date)"
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting BAM to chain conversion job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.chain, output.net],
                log_file,
                "bam_to_chain",
                max_wait_s=21600
            )

rule transmap_pairwise_map_psl:
    """
    Runs transMap on BAM-based chains to map reference PSL alignments to target genome.
    Similar to transmap_map_psl but uses chains generated from BAM files.
    """
    input:
        target_2bit=f"{config['work_dir']}/genome_files/{{genome}}.2bit",
        chain_file=f"{config['work_dir']}/chaining_bam/{{genome}}/{config['ref_genome']}-{{genome}}.chain",
        ref_psl=f"{config['work_dir']}/reference/{config['ref_genome']}.psl",
        ref_fa=f"{config['work_dir']}/reference/{config['ref_genome']}.fa",
        ref_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp"
    wildcard_constraints:
        genome = ANNOTATION_GENOME_WC
    output:
        tm_psl=temp(f"{config['work_dir']}/transMap_pairwise/{{genome}}.psl"),
        tm_gp=f"{config['work_dir']}/transMap_pairwise/{{genome}}.gp"
    threads: 1 if IS_CLUSTER else get_local_res("transmap_pairwise_map_psl", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("transmap_pairwise_map_psl", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("transmap_pairwise_map_psl", "time_h"),
        job_id=lambda wildcards, attempt: f"transmap-bam-map-{wildcards.genome}-{attempt}"
    log:
        f"{config['work_dir']}/logs/transmap_pairwise_map/{{genome}}.log"
    run:
        import os
        import subprocess
        import time

        work_dir = config['work_dir']
        genome = wildcards.genome
        job_script = f"{work_dir}/transMap_pairwise/{genome}_map_job.sh"

        # Create output directory
        os.makedirs(f"{work_dir}/transMap_pairwise", exist_ok=True)

        script_content = build_sbatch_header(
            "transmap_pairwise_map_psl",
            f"transmap-{wildcards.genome}",
            f"{work_dir}/logs/transmap_pairwise_map/{wildcards.genome}_slurm.out",
            f"{work_dir}/logs/transmap_pairwise_map/{wildcards.genome}_slurm.err"
        ) + f"""
echo "Starting transMap BAM-based PSL mapping for: {genome}"
echo "Job ID: ${{SLURM_JOB_ID}}"
echo "Start time: $(date)"

# Run pslMap, trying both normal and -swapMap in case chain is reversed
(pslMap -chainMapFile {input.ref_psl} {input.chain_file} stdout || \\
 pslMap -swapMap -chainMapFile {input.ref_psl} {input.chain_file} stdout) | \\
pslMapPostChain stdin stdout | \\
sort --parallel={get_res('transmap_pairwise_map_psl', 'cpus')} -k14,14 -k16,16n | \\
pslRecalcMatch stdin {input.target_2bit} {input.ref_fa} stdout | \\
sort --parallel={get_res('transmap_pairwise_map_psl', 'cpus')} -k10,10 > {output.tm_psl} 2> {log[0]}

# Convert the raw PSL to GenePred
transMapPslToGenePred -nonCodingGapFillMax=80 -codingGapFillMax=50 \\
    {input.ref_gp} {output.tm_psl} {output.tm_gp} 2>> {log[0]}

echo "Completed transMap BAM-based PSL mapping for: {genome}"
echo "End time: $(date)"
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting transMap BAM-based PSL mapping job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.tm_psl, output.tm_gp],
                log_file,
                "transmap_pairwise_map_psl",
                max_wait_s=7200
            )

rule transmap_pairwise_filter:
    """
    Filters transMap results from BAM-based chains using the same filtering logic
    as the standard transMap filtering.
    """
    input:
        tm_psl=f"{config['work_dir']}/transMap_pairwise/{{genome}}.psl",
        tm_gp=f"{config['work_dir']}/transMap_pairwise/{{genome}}.gp",
        ref_psl=f"{config['work_dir']}/reference/{config['ref_genome']}.psl",
        ref_db=f"{config['work_dir']}/databases/{config['ref_genome']}.db",
        annotation_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp"
    output:
        filtered_psl=f"{config['work_dir']}/transMap_pairwise/{{genome}}_filtered.psl",
        filtered_gp=f"{config['work_dir']}/transMap_pairwise/{{genome}}_filtered.gp",
        metrics_json=f"{config['work_dir']}/transMap_pairwise/{{genome}}_filter_tm_metrics.json",
        resolved_df=temp(f"{config['work_dir']}/transMap_pairwise/{{genome}}_resolved_df.pkl")
    wildcard_constraints:
        genome = ANNOTATION_GENOME_WC
    params:
        db_path=f"{config['work_dir']}/databases/{{genome}}.db",
        global_near_best=config.get("tm_global_near_best", 0.1),
        filter_overlapping_genes="--filter-overlapping-genes" if config.get("tm_filter_overlapping", True) else "",
        overlapping_ignore_bases=config.get("tm_overlapping_ignore_bases", 0),
        script="cat/filter_transmap.py"
    log:
        f"{config['work_dir']}/logs/transmap_pairwise_filter/{{genome}}.log"
    resources:
        mem_gb=8,
        time_h=4,
        job_id=lambda wildcards, attempt: f"tm-bam-filter-{wildcards.genome}-{attempt}"
    shell:
        """
        python {params.script} \
            --tm-psl {input.tm_psl} \
            --ref-psl {input.ref_psl} \
            --tm-gp {input.tm_gp} \
            --db-path {input.ref_db} \
            --filtered-psl {output.filtered_psl} \
            --filtered-gp {output.filtered_gp} \
            --metrics-json {output.metrics_json} \
            --resolved-df {output.resolved_df} \
            --global-near-best {params.global_near_best} \
            {params.filter_overlapping_genes} \
            --overlapping-ignore-bases {params.overlapping_ignore_bases} \
            > {log} 2>&1
        
        # Recalculate GenePred after filtering
        transMapPslToGenePred -nonCodingGapFillMax=80 -codingGapFillMax=50 {input.annotation_gp} {output.filtered_psl} {output.filtered_gp}
        
        # Store resolved dataframe in database
        python3 - <<PY
import pandas as pd
from tools.sqlInterface import TmPwFilterEval
from tools.sqlite import ExclusiveSqlConnection
table = TmPwFilterEval.__tablename__
genome_db="{params.db_path}"
res_df = pd.read_pickle("{output.resolved_df}")
with ExclusiveSqlConnection(genome_db) as engine:
    res_df.to_sql(table, engine, if_exists="replace", index=True)
PY
        """

rule generate_hints:
    input:
        fasta = f"{WORK_DIR}/genome_files/{{genome}}.fa",
        bams = lambda wc: config.get("transcriptomic_data", {}).get(wc.genome, {}).get("bam", []),
        intron_bams = lambda wc: config.get("transcriptomic_data", {}).get(wc.genome, {}).get("intronbam", []),
        iso_bams = lambda wc: config.get("transcriptomic_data", {}).get(wc.genome, {}).get("isoseq_bam", []),
        annotation_gp = lambda wc: [p] if (p := config.get("transcriptomic_data", {}).get(wc.genome, {}).get("annotation")) else [],
        protein_fasta = lambda wc: [p] if (p := config.get("transcriptomic_data", {}).get(wc.genome, {}).get("protein_fasta")) else []
        #setup_done = rules.setup_pipeline_directories.output.setup_done
    output:
        hints = f"{WORK_DIR}/hints_database/{{genome}}_extrinsic_hints.gff"
    wildcard_constraints:
        # Only run if there are rnaseq_genomes or isoseq_genomes provided in the input YAML
        genome = f"({'|'.join(RNASEQ_GENOMES + ISOSEQ_GENOMES)})" if (RNASEQ_GENOMES + ISOSEQ_GENOMES) else "never_match"
    params:
        bams_arg = lambda wildcards, input: f"--bams {' '.join(input.bams)}" if input.bams else "",
        intron_bams_arg = lambda wildcards, input: f"--intron_bams {' '.join(input.intron_bams)}" if input.intron_bams else "",
        iso_bams_arg = lambda wildcards, input: f"--iso_bams {' '.join(input.iso_bams)}" if input.iso_bams else "",
        annotation_arg = lambda wildcards, input: f"--annotation_gp {input.annotation_gp[0]}" if input.annotation_gp else "",
        protein_arg = lambda wildcards, input: f"--protein_fasta {input.protein_fasta[0]}" if input.protein_fasta else "",
        slurm_memory = lambda wc: int(get_res("generate_hints", "mem").rstrip("G")),
        slurm_cpus = lambda wc: get_res("generate_hints", "cpus"),
        slurm_time = lambda wc: get_res("generate_hints", "time"),
        slurm_partition = lambda wc: _slurm_partition("generate_hints"),
        slurm_max_jobs = lambda wc: get_res("generate_hints", "max_concurrent_jobs"),
        jobStore = f"{WORK_DIR}/toil_job_stores/{{genome}}_hints",
        workDir = WORK_DIR,
        use_cluster = IS_CLUSTER,
        execution_mode = EXECUTION_MODE,
        exclude_nodes = _slurm_exclude(),
        module_load = _slurm_module(),
        sge_parallel_env = config.get('cluster', {}).get('sge', {}).get('parallel_env', 'smp'),
        sge_memory_flag = config.get('cluster', {}).get('sge', {}).get('memory_flag', 'h_vmem')
    threads: 1 if IS_CLUSTER else get_local_res("generate_hints", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("generate_hints", "mem_gb"),
        time_h=4 if IS_CLUSTER else get_local_res("generate_hints", "time_h"),
        job_id=lambda wildcards, attempt: f"hints-{wildcards.genome}-{attempt}"
    log:
        f"{WORK_DIR}/logs/generate_hints/{{genome}}.log"
    shell:
        r"""
        if [ "{params.use_cluster}" = "True" ]; then
            echo "Using {params.execution_mode}-based parallel hints generation for {wildcards.genome}" > {log}
            python3 cat/hints_db.py \
              --mode cluster \
              --execution_mode {params.execution_mode} \
              --exclude_nodes "{params.exclude_nodes}" \
              --module_load "{params.module_load}" \
              --sge_parallel_env {params.sge_parallel_env} \
              --sge_memory_flag {params.sge_memory_flag} \
              --genome {wildcards.genome} \
              --fasta {input.fasta} \
              {params.bams_arg} \
              {params.intron_bams_arg} \
              {params.iso_bams_arg} \
              {params.annotation_arg} \
              {params.protein_arg} \
              --hints_out {output.hints} \
              --slurm_memory {params.slurm_memory} \
              --slurm_cpus {params.slurm_cpus} \
              --slurm_time {params.slurm_time} \
              --slurm_partition {params.slurm_partition} \
              --slurm_max_jobs {params.slurm_max_jobs} >> {log} 2>&1
        else
            echo "Using local hints generation for {wildcards.genome}" > {log}
            rm -rf {params.jobStore}
            python3 cat/hints_db.py \
              --mode toil \
              --genome {wildcards.genome} \
              --fasta {input.fasta} \
              {params.bams_arg} \
              {params.intron_bams_arg} \
              {params.iso_bams_arg} \
              {params.annotation_arg} \
              {params.protein_arg} \
              --hints_out {output.hints} \
              --batchSystem single_machine \
              --maxCores {threads} \
              --maxMemory {resources.mem_gb}G \
              --workDir {params.workDir} \
              {params.jobStore} >> {log} 2>&1
        fi
        """

rule build_db:
    input:
        fastas=expand(f"{config['work_dir']}/genome_files/{{genome}}.fa", genome=TARGET_GENOMES),
        # Only include hints if there are rnaseq_genomes or isoseq_genomes provided in the input YAML
        hints=expand(f"{config['work_dir']}/hints_database/{{genome}}_extrinsic_hints.gff", genome=RNASEQ_GENOMES + ISOSEQ_GENOMES) if (RNASEQ_GENOMES + ISOSEQ_GENOMES) else []
    output:
        done=touch(f"{config['work_dir']}/hints_database/hints.db.done"),
        hints_db=f"{config['work_dir']}/hints_database/hints.db"
    priority: 100  # High priority to unblock augustus
    resources:
        mem_gb=16,
        time_h=2,
        job_id="build-hints-db"
    log:
        f"{config['work_dir']}/logs/build_db/hints.log"
    run:
        base_cmd = f"load2sqlitedb --noIdx --clean --dbaccess={output.hints_db}"

        # Load each genome's data
        for i, genome in enumerate(TARGET_GENOMES):
            fasta = input.fastas[i]
            
            shell(f"echo 'Loading sequence for {genome}' >> {log}; "
                  f"{base_cmd} --species={genome} {fasta} 2>> {log}")

        # Only load hints if there are rnaseq or isoseq genomes
        if RNASEQ_GENOMES + ISOSEQ_GENOMES:
            for i, genome in enumerate(RNASEQ_GENOMES + ISOSEQ_GENOMES):
                hints_file = input.hints[i]
                
                if os.path.getsize(hints_file) > 0:
                    shell(f"echo 'Loading hints for {genome}' >> {log}; "
                          f"{base_cmd} --species={genome} {hints_file} 2>> {log}")
                else:
                    shell(f"echo 'No hints to load for {genome}' >> {log}")
        else:
            shell("echo 'No rnaseq or isoseq genomes - skipping hints loading' >> {log}")

        shell("echo 'Indexing database...' >> {log}; "
              "load2sqlitedb --makeIdx --clean --dbaccess={output.hints_db} 2>> {log}")

rule run_chaining_per_genome:
    """
    Submits separate SLURM jobs for each target genome for parallel chaining.
    Each genome gets its own Slurm job, enabling true parallel execution.
    """
    input:
        hal=config['hal'],
        ref_two_bit=f"{config['work_dir']}/genome_files/{config['ref_genome']}.2bit",
        ref_sizes=f"{config['work_dir']}/genome_files/{config['ref_genome']}.chrom.sizes",
        target_two_bits=expand(f"{config['work_dir']}/genome_files/{{genome}}.2bit", genome=ANNOTATION_GENOMES),
        target_sizes=expand(f"{config['work_dir']}/genome_files/{{genome}}.chrom.sizes", genome=ANNOTATION_GENOMES),
        setup_done=rules.setup_pipeline_directories.output.setup_done
    output:
        chains=expand(f"{config['work_dir']}/chaining/{{genome}}/{config['ref_genome']}-{{genome}}.chain", genome=ANNOTATION_GENOMES)
    params:
        ref_genome=config['ref_genome'],
        work_dir=config['work_dir'],
        target_genomes=ANNOTATION_GENOMES,
        target_two_bit_args=lambda wildcards, input: " ".join([f"--target_genome {genome} --target_two_bit {input.target_two_bits[i]} --chain_file {config['work_dir']}/chaining/{genome}/{config['ref_genome']}-{genome}.chain" for i, genome in enumerate(ANNOTATION_GENOMES)])
    log:
        f"{config['work_dir']}/logs/chaining_parallel.log"
    threads: 1 if IS_CLUSTER else get_local_res("run_chaining_per_genome", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("run_chaining_per_genome", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("run_chaining_per_genome", "time_h"),
        job_id=lambda wildcards, attempt: f"chaining-parallel-submit-{attempt}"
    run:
        import os
        import subprocess
        import time
        import re

        # Get paths and variables
        work_dir = config['work_dir']
        timeout_hours = get_res('run_chaining_per_genome', 'timeout_hours')

        # Create logs directory for per-genome jobs
        os.makedirs(f"{work_dir}/logs/chaining", exist_ok=True)

        def _make_chaining_body(genome, target_two_bit, chain_file):
            """Return the bash body (without SLURM header) for a single genome."""
            return f"""
# Log start time
echo "Starting chaining for genome: {genome}"
echo "Job ID: ${{SLURM_JOB_ID}}"
echo "Node: $(hostname)"
echo "Start time: $(date)"

# Create output directory and temporary working directory
mkdir -p "{os.path.abspath(os.path.dirname(chain_file))}"
TEMP_DIR="{os.path.abspath(os.path.dirname(chain_file))}/temp_${{SLURM_JOB_ID}}"
mkdir -p "$TEMP_DIR"
cd "$TEMP_DIR"

# Link/copy input files to temp directory (using absolute paths)
# Use symbolic links for large files to avoid slow copies
echo "Setting up input files..."
ln -s "{os.path.abspath(input.hal)}" hal_file.hal
cp "{os.path.abspath(input.ref_sizes)}" query_sizes.txt
ln -s "{os.path.abspath(input.ref_two_bit)}" query.2bit
ln -s "{os.path.abspath(target_two_bit)}" target.2bit

# Create error log file
touch chaining_errors.log

echo "Processing genome: {genome}"

# Process each chromosome from the reference genome
while read -r chrom size; do
    if [[ -z "$chrom" || -z "$size" ]]; then
        continue
    fi

    echo "Processing chromosome: $chrom ($size bp)"

    # Create BED file for this chromosome
    echo -e "$chrom\\t0\\t$size" > "${{chrom}}.bed"

    # Run halLiftover -> pslPosTarget -> axtChain pipeline
    if halLiftover --outPSL hal_file.hal "{config['ref_genome']}" "${{chrom}}.bed" "{genome}" /dev/stdout 2>> chaining_errors.log | \\
       pslPosTarget /dev/stdin /dev/stdout 2>> chaining_errors.log | \\
       axtChain -psl -verbose=0 -linearGap=medium /dev/stdin target.2bit query.2bit "${{chrom}}.chain" 2>> chaining_errors.log; then
        echo "Successfully processed chromosome $chrom"
    else
        echo "Warning: Failed to process chromosome $chrom, creating empty chain"
        echo "  Check chaining_errors.log for details"
        touch "${{chrom}}.chain"
    fi

    # Clean up BED file
    rm -f "${{chrom}}.bed"

done < query_sizes.txt

# Merge all chromosome chains
echo "Merging chain files for {genome}..."

# Count total and non-empty chain files
total_files=$(find . -maxdepth 1 -name "*.chain" | wc -l)
non_empty_files=$(find . -maxdepth 1 -name "*.chain" -size +0 | wc -l)
empty_files=$((total_files - non_empty_files))

echo "Found $total_files chain files ($non_empty_files non-empty, $empty_files empty)"

# Create file list for chainMergeSort (only non-empty files)
find . -maxdepth 1 -name "*.chain" -size +0 | sort > chain_files.lst

if [[ -s chain_files.lst ]]; then
    echo "Merging $non_empty_files non-empty chain files..."

    # Create temporary directory for chainMergeSort
    mkdir -p temp_merge

    if chainMergeSort -inputList=chain_files.lst -tempDir=temp_merge/ > "{os.path.abspath(chain_file)}"; then
        if [[ -f "{os.path.abspath(chain_file)}" && -s "{os.path.abspath(chain_file)}" ]]; then
            chain_lines=$(wc -l < "{os.path.abspath(chain_file)}")
            echo "Successfully created chain file: {os.path.abspath(chain_file)} ($chain_lines lines)"
        else
            echo "Warning: chainMergeSort produced empty output, creating empty file"
            touch "{os.path.abspath(chain_file)}"
        fi
    else
        echo "Error: chainMergeSort failed, creating empty file"
        touch "{os.path.abspath(chain_file)}"
    fi

    # Clean up temporary directory
    rm -rf temp_merge
else
    echo "Warning: No valid chain files found, creating empty chain file"
    touch "{os.path.abspath(chain_file)}"
fi

# Copy error log to output directory if it has content
if [[ -s chaining_errors.log ]]; then
    cp chaining_errors.log "{os.path.abspath(os.path.dirname(chain_file))}/chaining_errors_{genome}.log"
    echo "Errors logged to chaining_errors_{genome}.log"
fi

# Clean up temporary directory
cd /
rm -rf "$TEMP_DIR"

# Log completion
echo "Completed chaining for genome: {genome}"
echo "End time: $(date)"
"""

        with open(log[0], 'w') as log_file:
            if IS_CLUSTER:
                submitted_jobs = {}
                log_file.write(f"Submitting {len(ANNOTATION_GENOMES)} separate {SCHEDULER.name} jobs for parallel chaining...\n\n")

                for i, genome in enumerate(ANNOTATION_GENOMES):
                    target_two_bit = input.target_two_bits[i]
                    chain_file = f"{work_dir}/chaining/{genome}/{config['ref_genome']}-{genome}.chain"
                    job_script = f"{work_dir}/chaining/chaining_{genome}_job.sh"

                    script_content = build_sbatch_header(
                        "run_chaining_per_genome",
                        f"chain-{genome}",
                        f"{work_dir}/logs/chaining/{genome}_cluster.out",
                        f"{work_dir}/logs/chaining/{genome}_cluster.err"
                    ) + _make_chaining_body(genome, target_two_bit, chain_file)

                    SCHEDULER.write_script(script_content, job_script)
                    try:
                        job_id = SCHEDULER.submit(job_script)
                    except Exception as e:
                        log_file.write(f"Error submitting job for {genome}: {e}\n")
                        raise

                    submitted_jobs[genome] = job_id
                    log_file.write(f"[{i+1}/{len(ANNOTATION_GENOMES)}] Submitted job {job_id} for genome: {genome}\n")

                log_file.write(f"\nAll {len(submitted_jobs)} jobs submitted successfully!\n")
                log_file.write(f"Job IDs: {', '.join(submitted_jobs.values())}\n\n")
                log_file.write(f"Waiting for all jobs to complete (timeout: {timeout_hours} hours)...\n")
                log_file.flush()

                max_wait = timeout_hours * 3600
                check_interval = 30
                elapsed = 0

                while elapsed < max_wait:
                    all_files_exist = all(os.path.exists(cf) for cf in output.chains)

                    if all_files_exist:
                        log_file.write(f"All chaining jobs completed successfully after {elapsed//60} minutes!\n")
                        log_file.write("\nChain file summary:\n")
                        for genome, cf in zip(ANNOTATION_GENOMES, output.chains):
                            if os.path.exists(cf):
                                size = os.path.getsize(cf)
                                log_file.write(f"  {genome}: {size:,} bytes\n")
                        log_file.flush()
                        break

                    completed_files = sum(1 for cf in output.chains if os.path.exists(cf))
                    if elapsed % 300 == 0 and elapsed > 0:
                        log_file.write(f"  Progress: {completed_files}/{len(output.chains)} genomes completed ({elapsed//60} min elapsed)\n")
                        log_file.flush()

                    time.sleep(check_interval)
                    elapsed += check_interval
                else:
                    failed = [g for g, cf in zip(ANNOTATION_GENOMES, output.chains) if not os.path.exists(cf)]
                    log_file.write(f"\nTimeout reached after {timeout_hours} hours!\n")
                    if failed:
                        log_file.write(f"Failed/Incomplete: {', '.join(failed)}\n")
                    log_file.flush()
                    raise Exception(f"Chaining jobs timed out. Failed genomes: {', '.join(failed)}")

            else:
                # Local mode: run each genome sequentially
                log_file.write(f"Running chaining locally for {len(ANNOTATION_GENOMES)} genomes...\n\n")

                for i, genome in enumerate(ANNOTATION_GENOMES):
                    target_two_bit = input.target_two_bits[i]
                    chain_file = f"{work_dir}/chaining/{genome}/{config['ref_genome']}-{genome}.chain"

                    body = _make_chaining_body(genome, target_two_bit, chain_file)
                    # Strip scheduler directive lines (SLURM #SBATCH and SGE #$)
                    # so the body is safe to run via plain bash locally.
                    clean = "\n".join(
                        line for line in body.splitlines()
                        if not line.startswith("#SBATCH") and not line.startswith("#$") and line != "#!/bin/bash"
                    )
                    log_file.write(f"[{i+1}/{len(ANNOTATION_GENOMES)}] Running chaining locally for {genome}...\n")
                    log_file.flush()
                    result = subprocess.run(['bash', '-euo', 'pipefail', '-c', clean])
                    if result.returncode != 0:
                        raise Exception(f"Local chaining failed for {genome} (rc={result.returncode})")
                    log_file.write(f"  Completed {genome}\n")

                log_file.write(f"\nAll {len(ANNOTATION_GENOMES)} genomes chained successfully.\n")

def get_stringtie_bams(wildcards):
    genome_config = config.get("transcriptomic_data", {}).get(wildcards.genome, {})
    sr_bams = genome_config.get("intronbam", []) or genome_config.get("bam", [])
    lr_bams = genome_config.get("isoseq_bam", [])
    return {"sr": sr_bams, "lr": lr_bams}

rule stringtie_merge_bams:
    """
    Merges multiple short-read and long-read BAMs into one file each.
    If only one BAM exists, it is copied. If none exist, an empty file is touched.
    """
    input:
        # Use the helper function to dynamically get lists of BAMs
        unpack(get_stringtie_bams)
    output:
        sr_merged=f"{config['work_dir']}/stringtie/{{genome}}_short_read_merged.bam",
        lr_merged=f"{config['work_dir']}/stringtie/{{genome}}_long_read_merged.bam"
    log:
        f"{config['work_dir']}/logs/stringtie_merge/{{genome}}.log"
    threads: 1 if IS_CLUSTER else get_local_res("stringtie_merge", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("stringtie_merge", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("stringtie_merge", "time_h"),
        job_id=lambda wildcards, attempt: f"stringtie-merge-{wildcards.genome}-{attempt}"
    run:
        work_dir = config['work_dir']
        genome = wildcards.genome
        job_script = f"{work_dir}/stringtie/{genome}_merge_job.sh"

        def _merge_cmd(bams, out):
            bams = list(bams)
            if len(bams) > 1:
                return f"samtools merge -f -@ {get_res('stringtie_merge', 'cpus')} {out} {' '.join(bams)} >> {log[0]} 2>&1"
            elif len(bams) == 1:
                return f"cp {bams[0]} {out} >> {log[0]} 2>&1"
            return f"touch {out}"

        script_content = build_sbatch_header(
            "stringtie_merge",
            f"strg-merge-{genome}",
            f"{work_dir}/logs/stringtie_merge/{genome}_slurm.out",
            f"{work_dir}/logs/stringtie_merge/{genome}_slurm.err",
        ) + f"""
echo "Merging StringTie BAMs for {genome}" >> {log[0]}
{_merge_cmd(input.sr, output.sr_merged)}
{_merge_cmd(input.lr, output.lr_merged)}
echo "Done merging BAMs for {genome}" >> {log[0]}
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting StringTie BAM merge job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.sr_merged, output.lr_merged],
                log_file,
                "stringtie_merge",
                max_wait_s=6 * 3600,
            )

rule stringtie_sort_bams:
    """Sorts the merged BAMs in preparation for StringTie."""
    input:
        # Use a wildcard to handle both short and long read types
        bam=f"{config['work_dir']}/stringtie/{{genome}}_{{type}}_read_merged.bam"
    output:
        bam=f"{config['work_dir']}/stringtie/{{genome}}_{{type}}_read_merged_sorted.bam"
    log:
        f"{config['work_dir']}/logs/stringtie_sort/{{genome}}_{{type}}.log"
    threads: 1 if IS_CLUSTER else get_local_res("stringtie_sort", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("stringtie_sort", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("stringtie_sort", "time_h"),
        job_id=lambda wildcards, attempt: f"stringtie-sort-{wildcards.genome}-{wildcards.type}-{attempt}"
    run:
        work_dir = config['work_dir']
        genome = wildcards.genome
        rtype = wildcards.type
        job_script = f"{work_dir}/stringtie/{genome}_{rtype}_sort_job.sh"
        cpus = get_res("stringtie_sort", "cpus")

        script_content = build_sbatch_header(
            "stringtie_sort",
            f"strg-sort-{genome}-{rtype}",
            f"{work_dir}/logs/stringtie_sort/{genome}_{rtype}_slurm.out",
            f"{work_dir}/logs/stringtie_sort/{genome}_{rtype}_slurm.err",
        ) + f"""
# Only sort if the input file is not empty
if [ -s {input.bam} ]; then
    samtools sort -@ {cpus} {input.bam} -o {output.bam} >> {log[0]} 2>&1
else
    touch {output.bam}
fi
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting StringTie sort job for {genome} ({rtype})...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.bam],
                log_file,
                "stringtie_sort",
                max_wait_s=8 * 3600,
            )

rule stringtie_run:
    """
    Runs StringTie itself. The command is chosen based on whether
    short-read, long-read, or both types of evidence are available.
    Uses --mix when both data types are present.
    """
    input:
        sr_sorted=f"{config['work_dir']}/stringtie/{{genome}}_short_read_merged_sorted.bam",
        lr_sorted=f"{config['work_dir']}/stringtie/{{genome}}_long_read_merged_sorted.bam"
    output:
        temp_gtf=f"{config['work_dir']}/stringtie/{{genome}}_temp_stringtie.gtf"
    log:
        f"{config['work_dir']}/logs/stringtie_run/{{genome}}.log"
    threads: 1 if IS_CLUSTER else get_local_res("stringtie_run", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("stringtie_run", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("stringtie_run", "time_h"),
        job_id=lambda wildcards, attempt: f"stringtie-run-{wildcards.genome}-{attempt}"
    run:
        work_dir = config['work_dir']
        genome = wildcards.genome
        job_script = f"{work_dir}/stringtie/{genome}_run_job.sh"
        cpus = get_res("stringtie_run", "cpus")

        sr_exists = os.path.getsize(input.sr_sorted) > 0
        lr_exists = os.path.getsize(input.lr_sorted) > 0

        if sr_exists and lr_exists:
            strg_cmd = f"stringtie --mix -o {output.temp_gtf} -p {cpus} {input.sr_sorted} {input.lr_sorted}"
        elif lr_exists:
            strg_cmd = f"stringtie -L -o {output.temp_gtf} -p {cpus} {input.lr_sorted}"
        elif sr_exists:
            strg_cmd = f"stringtie -o {output.temp_gtf} -p {cpus} {input.sr_sorted}"
        else:
            strg_cmd = f"touch {output.temp_gtf}"

        script_content = build_sbatch_header(
            "stringtie_run",
            f"strg-run-{genome}",
            f"{work_dir}/logs/stringtie_run/{genome}_slurm.out",
            f"{work_dir}/logs/stringtie_run/{genome}_slurm.err",
        ) + f"""
echo "Running StringTie for {genome}" >> {log[0]}
{strg_cmd} >> {log[0]} 2>&1
echo "StringTie completed for {genome}" >> {log[0]}
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting StringTie run job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.temp_gtf],
                log_file,
                "stringtie_run",
                max_wait_s=14 * 3600,
            )

if BUILD_PROTEIN_DB:
    rule build_protein_db:
        """
        Assemble the multi-species protein database used by miniprot/augMP from
        UniProt reference proteomes named in config['protein_db'] -- by species,
        by taxon id, and/or by clade (genus/family; one proteome per member
        species, best BUSCO first) -- plus an optional local base_fasta. Runs once
        on the controller node (it needs outbound internet); per-proteome downloads
        are cached next to the output so reruns are fast. Its output is the miniprot
        protein input for every genome.
        """
        output:
            fasta=PROTEIN_FASTA
        params:
            species=_PROTEIN_DB_CFG.get("species", []),
            taxa=_PROTEIN_DB_CFG.get("taxa", []),
            clades=_PROTEIN_DB_CFG.get("clades", []),
            max_per_clade=_PROTEIN_DB_CFG.get("max_per_clade", 25),
            clade_include_other=_PROTEIN_DB_CFG.get("clade_include_other", False),
            base_fasta=_PROTEIN_DB_CFG.get("base_fasta", ""),
            min_len=_PROTEIN_DB_CFG.get("min_len", 20),
            dedup=_PROTEIN_DB_CFG.get("dedup", True),
            script=str(CAT2_ROOT / "scripts" / "build_protein_db.py"),
        log:
            f"{config['work_dir']}/logs/build_protein_db.log"
        resources:
            mem_gb=4,
            time_h=6,
            job_id=lambda wildcards, attempt: f"build-protein-db-{attempt}"
        run:
            import shlex
            args = ["python3", params.script, "--out", output.fasta]
            for s in (params.species or []):
                args += ["--species", str(s)]
            for t in (params.taxa or []):
                args += ["--taxa", str(t)]
            for c in (params.clades or []):
                args += ["--clades", str(c)]
            args += ["--max-per-clade", str(params.max_per_clade)]
            if params.clade_include_other:
                args += ["--clade-include-other"]
            if params.base_fasta:
                args += ["--base-fasta", str(params.base_fasta)]
            args += ["--min-len", str(params.min_len)]
            if not params.dedup:
                args += ["--no-dedup"]
            cmd = " ".join(shlex.quote(a) for a in args)
            shell("mkdir -p $(dirname {log[0]})")
            shell(cmd + " > {log[0]} 2>&1")

rule run_miniprot:
    """
    Submits miniprot as a SLURM job to the short queue to align a protein reference 
    set against a target genome, then converts the alignments to a GFF hints file.
    Runs minisplice predict first to generate splice-site scores for miniprot.
    When minisplice_parallel is enabled (default), minisplice runs one task per
    chromosome and merges the TSV scores before miniprot.
    """
    input:
        genome_fasta=f"{config['work_dir']}/genome_files/{{genome}}.fa",
        genome_index=f"{config['work_dir']}/genome_files/{{genome}}.mpi",
        chrom_sizes=f"{config['work_dir']}/genome_files/{{genome}}.chrom.sizes",
        protein_ref=PROTEIN_FASTA,
        minisplice_model=config["minisplice_model"],
        minisplice_calibration=config["minisplice_calibration"]
    output:
        hints=f"{config['work_dir']}/miniprot/{{genome}}_miniprot_hints.gff",
        paf=f"{config['work_dir']}/miniprot/{{genome}}_miniprot.paf",
        splice_scores=f"{config['work_dir']}/miniprot/{{genome}}_minisplice_scores.tsv"
    wildcard_constraints:
        genome = "|".join(TARGET_GENOMES)
    priority: 90  # High priority to unblock augustus
    log:
        f"{config['work_dir']}/logs/miniprot/{{genome}}.log"
    threads: 1 if IS_CLUSTER else get_local_res("run_miniprot", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("run_miniprot", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("run_miniprot", "time_h"),
        job_id=lambda wildcards, attempt: f"miniprot-submit-{wildcards.genome}-{attempt}"
    run:
        import os
        import subprocess
        import time

        # Get paths and variables
        work_dir = config['work_dir']
        genome = wildcards.genome
        job_script = f"{work_dir}/miniprot/{genome}_miniprot_job.sh"
        minisplice_step = build_minisplice_step(
            work_dir=work_dir,
            genome=genome,
            genome_fasta=input.genome_fasta,
            chrom_sizes=input.chrom_sizes,
            minisplice_model=input.minisplice_model,
            minisplice_calibration=input.minisplice_calibration,
            splice_scores_out=output.splice_scores,
            log_path=log[0],
            use_slurm_array=IS_CLUSTER,
        )
        timeout_s = int(get_res("run_miniprot", "timeout_hours") * 3600)
        map_flags = build_miniprot_map_flags(get_res('run_miniprot', 'cpus'))

        script_content = build_sbatch_header(
            "run_miniprot",
            f"miniprot-{genome}",
            f"{work_dir}/logs/miniprot/{genome}_slurm.out",
            f"{work_dir}/logs/miniprot/{genome}_slurm.err"
        ) + minisplice_step + f"""

# Step 1a: Run miniprot in PAF mode with splice scores (templates for augMP)
echo "Running miniprot for PAF (templates)..." >> {log[0]}
echo "miniprot map flags: {map_flags}" >> {log[0]}
miniprot {map_flags} --spsc={output.splice_scores} {input.genome_index} {input.protein_ref} > {output.paf} 2>> {log[0]}

# Step 1b: Run miniprot with --gtf and splice scores for hints generation
echo "Running miniprot for GTF (hints)..." >> {log[0]}
TEMP_GTF={work_dir}/miniprot/{genome}_miniprot_temp.gtf
miniprot {map_flags} --gtf --spsc={output.splice_scores} {input.genome_index} {input.protein_ref} > $TEMP_GTF 2>> {log[0]}

# Step 2: Convert the GTF output to an Augustus hints file
echo "Converting GTF to hints..." >> {log[0]}
aln2hints.pl --genome_file={input.genome_fasta} \\
             --in=$TEMP_GTF \\
             --out={output.hints} \\
             --prg=miniprot >> {log[0]} 2>&1

# Clean up temporary GTF
rm -f $TEMP_GTF
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting miniprot job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.hints, output.paf, output.splice_scores],
                log_file,
                "run_miniprot",
                max_wait_s=timeout_s
            )

rule run_transcript_map:
    """
    Transcript-level minimap2 mapping (txTM).

    Submits a SLURM job (cluster mode) or runs locally with full resources.
    Config keys under ``config['transcript_map']`` apply **only** to this rule.
    (minimap2 → filtered PSL → genePred here). They do not change transMap,
    miniprot, or other annotation paths.

    Maps reference transcripts directly to each target genome using minimap2
    (splice-aware), converts the alignments to genePredExt using reference-guided
    exon merging, and transfers CDS annotations from the reference GenePred.
    Supports CNV gene annotation via secondary alignments.

    Outputs go to the ``txTM/`` work-dir subdirectory; downstream consensus rules
    pick the same paths up by mode name.
    """
    input:
        ref_tx_fa=f"{config['work_dir']}/reference/{config['ref_genome']}.fa",
        ref_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp",
        ref_db=f"{config['work_dir']}/databases/{config['ref_genome']}.db",
        target_fa=f"{config['work_dir']}/genome_files/{{genome}}.fa"
    output:
        gp=f"{config['work_dir']}/txTM/{{genome}}_txTM.gp",
        psl=f"{config['work_dir']}/txTM/{{genome}}_filtered.psl",
        gtf=f"{config['work_dir']}/txTM/{{genome}}_txTM.gtf",
        attrs=f"{config['work_dir']}/txTM/{{genome}}_txTM.gp_attrs",
        dups=f"{config['work_dir']}/txTM/{{genome}}_txTM.duplicates.txt"
    wildcard_constraints:
        genome = TXTM_GENOME_WC
    log:
        f"{config['work_dir']}/logs/transcript_map/{{genome}}.log"
    threads: 1 if IS_CLUSTER else get_local_res("run_txTM", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("run_txTM", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("run_txTM", "time_h"),
        job_id=lambda wildcards, attempt: f"transcript-map-submit-{wildcards.genome}-{attempt}"
    run:
        import yaml
        from pathlib import Path

        work_dir = config['work_dir']
        genome   = wildcards.genome
        threads  = get_res('run_txTM', 'cpus') if IS_CLUSTER else snakemake.threads
        job_script = f"{work_dir}/txTM/{genome}_txTM_job.sh"

        Path(f"{work_dir}/txTM").mkdir(parents=True, exist_ok=True)
        cfg_snapshot = f"{work_dir}/txTM/{genome}_txmap_config.yaml"
        with open(cfg_snapshot, 'w') as fh:
            yaml.dump({
                'work_dir': work_dir,
                'ref_genome': config['ref_genome'],
                'txTM_sc': config.get('txTM_sc', 0.80),
                'transcript_map': config.get('transcript_map', {}),
            }, fh)

        script_content = build_sbatch_header(
            "run_txTM",
            f"txTM-{genome}",
            f"{work_dir}/logs/transcript_map/{genome}_slurm.out",
            f"{work_dir}/logs/transcript_map/{genome}_slurm.err"
        ) + f"""
python3 cat/transcript_map_runner.py \\
    --config {cfg_snapshot} \\
    --genome {genome} \\
    --threads {threads} \\
    --log {log[0]} \\
    --output-gp {output.gp} \\
    --output-psl {output.psl} \\
    --output-gtf {output.gtf} \\
    --output-attrs {output.attrs} \\
    --output-dups {output.dups}
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting transcript_map job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.gp, output.psl, output.attrs],
                log_file,
                "run_txTM",
                max_wait_s=int(get_res("run_txTM", "timeout_hours") * 3600),
            )


rule augustus_extract_coding_gp:
    """
    Extracts only coding transcripts from the filtered transMap genePred.
    This is a pre-processing step for the main Augustus pipeline.
    """
    input:
        filtered_gp=f"{config['work_dir']}/transMap/{{genome}}_filtered.gp"
    output:
        coding_gp=f"{config['work_dir']}/transMap/{{genome}}_coding.gp"
    wildcard_constraints:
        genome = AUGUSTUS_GENOME_WC
    resources:
        mem_gb=2,
        time_h=1,
        job_id=lambda wildcards, attempt: f"extract-coding-{wildcards.genome}-{attempt}"
    shell:
        # A transcript is coding if the CDS start (col 6) != CDS end (col 7)
        "awk '$6 != $7' {input} > {output}"

rule augustus_run_tm_and_tmr:
    """
    Launches the parallel augustus_parallel.py pipeline for genomes with RNA-seq,
    generating both augTM and augTMR outputs using chromosome-based parallelization.
    """
    input:
        script="cat/augustus_parallel.py",
        fasta=f"{config['work_dir']}/genome_files/{{genome}}.fa",
        coding_gp=f"{config['work_dir']}/transMap/{{genome}}_coding.gp",
        filtered_tm_psl=f"{config['work_dir']}/transMap/{{genome}}_filtered.psl",
        ref_psl=f"{config['work_dir']}/reference/{config['ref_genome']}.psl",
        annotation_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp",
        hints_db=f"{config['work_dir']}/hints_database/hints.db"
    output:
        tm_gtf=f"{config['work_dir']}/augustus/{{genome}}_augTM.gtf",
        tmr_gtf=f"{config['work_dir']}/augustus/{{genome}}_augTMR.gtf",
        tmr_gtf_done=f"{config['work_dir']}/augustus/{{genome}}_augTMR.gtf.done",
        tm_gtf_done=f"{config['work_dir']}/augustus/{{genome}}_augTM.gtf.done"
    wildcard_constraints:
        # Only run for genomes that have RNA-seq data AND when augustus is enabled
        genome = AUGUSTUS_RNASEQ_GENOME_WC
    params:
        work_dir=config['work_dir'],
        species=config['augustus_species'],
        tm_cfg=config['tm_cfg_path'],
        tmr_cfg=config['tmr_cfg_path'],
        utr=1 if config.get('predict_utr') else 0
    threads: 1 if IS_CLUSTER else get_local_res("augustus_run_tm_and_tmr", "threads")
    resources:
        mem_gb=4 if IS_CLUSTER else get_local_res("augustus_run_tm_and_tmr", "mem_gb"),
        time_h=24 if IS_CLUSTER else get_local_res("augustus_run_tm_and_tmr", "time_h"),
        job_id=lambda wildcards, attempt: f"augustus-parallel-tmr-{wildcards.genome}-{attempt}"
    log:
        f"{config['work_dir']}/logs/augustus_run/{{genome}}.log"
    run:
        run_augustus_parallel(
            input=input, output=output, params=params, wildcards=wildcards,
            log_path=log[0],
            genome_work_dir=f"{params.work_dir}/augustus_parallel_temp_{wildcards.genome}",
            with_tmr=True,
        )

rule augustus_run_tm_only:
    """
    Launches the parallel augustus_parallel.py pipeline for genomes without RNA-seq,
    generating only the augTM output using chromosome-based parallelization.
    """
    input:
        script="cat/augustus_parallel.py",
        fasta=f"{config['work_dir']}/genome_files/{{genome}}.fa",
        coding_gp=f"{config['work_dir']}/transMap/{{genome}}_coding.gp",
        filtered_tm_psl=f"{config['work_dir']}/transMap/{{genome}}_filtered.psl",
        ref_psl=f"{config['work_dir']}/reference/{config['ref_genome']}.psl",
        annotation_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp"
    output:
        tm_gtf=f"{config['work_dir']}/augustus/{{genome}}_augTM.gtf",
        tm_gtf_done=f"{config['work_dir']}/augustus/{{genome}}_augTM.gtf.done"
    wildcard_constraints:
        # Only run for genomes that do NOT have RNA-seq data AND when augustus is enabled
        genome = AUGUSTUS_NON_RNASEQ_GENOME_WC
    params:
        work_dir=config['work_dir'],
        species=config['augustus_species'],
        tm_cfg=config['tm_cfg_path'],
        utr=1 if config.get('predict_utr') else 0
    log:
        f"{config['work_dir']}/logs/augustus_run/{{genome}}.log"
    threads: 1 if IS_CLUSTER else get_local_res("augustus_run_tm_only", "threads")
    resources:
        mem_gb=4 if IS_CLUSTER else get_local_res("augustus_run_tm_only", "mem_gb"),
        time_h=12 if IS_CLUSTER else get_local_res("augustus_run_tm_only", "time_h"),
        job_id=lambda wildcards, attempt: f"augustus-parallel-tm-{wildcards.genome}-{attempt}"
    run:
        run_augustus_parallel(
            input=input, output=output, params=params, wildcards=wildcards,
            log_path=log[0],
            genome_work_dir=f"{params.work_dir}/augustus_parallel_temp_{wildcards.genome}",
            with_tmr=False,
        )

rule augustus_convert_tm_gtf_to_gp:
    """
    Post-processing step to convert Augustus TM GTF files to GenePred.
    Runs for each genome with Augustus TM output.
    """
    input:
        gtf=f"{config['work_dir']}/augustus/{{genome}}_augTM.gtf",
        gtf_done=f"{config['work_dir']}/augustus/{{genome}}_augTM.gtf.done"
    output:
        gp=f"{config['work_dir']}/augustus/{{genome}}_augTM.gp"
    wildcard_constraints:
        # Only run for genomes that have Augustus enabled
        genome = AUGUSTUS_GENOME_WC
    resources:
        mem_gb=16,
        time_h=2,
        job_id=lambda wildcards, attempt: f"aug-convert-tm-{wildcards.genome}-{attempt}"
    log:
        f"{config['work_dir']}/logs/augustus_convert/{{genome}}_augTM.log"
    shell:
        """
        if [ -s {input.gtf} ]; then
            gtfToGenePred -genePredExt {input.gtf} {output.gp} &> {log}
        else
            echo "Empty Augustus TM GTF - creating empty GenePred" > {log}
            touch {output.gp}
        fi
        """

rule augustus_convert_tmr_gtf_to_gp:
    """
    Post-processing step to convert Augustus TMR GTF files to GenePred.
    Only runs for genomes with RNA-seq data (rnaseq_genomes).
    """
    input:
        gtf=f"{config['work_dir']}/augustus/{{genome}}_augTMR.gtf",
        gtf_done=f"{config['work_dir']}/augustus/{{genome}}_augTMR.gtf.done"
    output:
        gp=f"{config['work_dir']}/augustus/{{genome}}_augTMR.gp"
    wildcard_constraints:
        # Only run for genomes that have RNA-seq data
        genome = RNASEQ_GENOME_WC
    resources:
        mem_gb=16,
        time_h=2,
        job_id=lambda wildcards, attempt: f"aug-convert-tmr-{wildcards.genome}-{attempt}"
    log:
        f"{config['work_dir']}/logs/augustus_convert/{{genome}}_augTMR.log"
    shell:
        """
        if [ -s {input.gtf} ]; then
            gtfToGenePred -genePredExt {input.gtf} {output.gp} &> {log}
        else
            echo "Empty Augustus TMR GTF - creating empty GenePred" > {log}
            touch {output.gp}
        fi
        """

# ═══════════════════════════════════════════════════════════════════════════════
# AUGUSTUS PAIRWISE MODE: Uses transMap_pairwise instead of transMap
# ═══════════════════════════════════════════════════════════════════════════════

rule augustus_extract_coding_gp_pairwise:
    """
    Extracts only coding transcripts from the filtered transMap_pairwise genePred.
    This is a pre-processing step for the main Augustus pipeline using pairwise chains.
    """
    input:
        filtered_gp=f"{config['work_dir']}/transMap_pairwise/{{genome}}_filtered.gp"
    output:
        coding_gp=f"{config['work_dir']}/transMap_pairwise/{{genome}}_coding.gp"
    wildcard_constraints:
        genome = AUGUSTUS_GENOME_WC
    resources:
        mem_gb=2,
        time_h=1,
        job_id=lambda wildcards, attempt: f"extract-coding-pairwise-{wildcards.genome}-{attempt}"
    shell:
        # A transcript is coding if the CDS start (col 6) != CDS end (col 7)
        "awk '$6 != $7' {input} > {output}"

rule augustus_run_tm_pairwise_and_tmr_pairwise:
    """
    Launches the parallel augustus_parallel.py pipeline for genomes with RNA-seq,
    generating both augTM_pairwise and augTMR_pairwise outputs using chromosome-based parallelization.
    Uses transMap_pairwise instead of transMap.
    """
    input:
        script="cat/augustus_parallel.py",
        fasta=f"{config['work_dir']}/genome_files/{{genome}}.fa",
        coding_gp=f"{config['work_dir']}/transMap_pairwise/{{genome}}_coding.gp",
        filtered_tm_psl=f"{config['work_dir']}/transMap_pairwise/{{genome}}_filtered.psl",
        ref_psl=f"{config['work_dir']}/reference/{config['ref_genome']}.psl",
        annotation_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp",
        hints_db=f"{config['work_dir']}/hints_database/hints.db"
    output:
        tm_gtf=f"{config['work_dir']}/augustus/{{genome}}_augTM_pairwise.gtf",
        tmr_gtf=f"{config['work_dir']}/augustus/{{genome}}_augTMR_pairwise.gtf",
        tmr_gtf_done=f"{config['work_dir']}/augustus/{{genome}}_augTMR_pairwise.gtf.done",
        tm_gtf_done=f"{config['work_dir']}/augustus/{{genome}}_augTM_pairwise.gtf.done"
    wildcard_constraints:
        # Only run for genomes that have RNA-seq data AND when augustus is enabled
        genome = AUGUSTUS_RNASEQ_GENOME_WC
    params:
        work_dir=config['work_dir'],
        species=config['augustus_species'],
        tm_cfg=config['tm_cfg_path'],
        tmr_cfg=config['tmr_cfg_path'],
        utr=1 if config.get('predict_utr') else 0
    threads: 1 if IS_CLUSTER else get_local_res("augustus_run_tm_pairwise_and_tmr_pairwise", "threads")
    resources:
        mem_gb=4 if IS_CLUSTER else get_local_res("augustus_run_tm_pairwise_and_tmr_pairwise", "mem_gb"),
        time_h=24 if IS_CLUSTER else get_local_res("augustus_run_tm_pairwise_and_tmr_pairwise", "time_h"),
        job_id=lambda wildcards, attempt: f"augustus-parallel-tmr-pairwise-{wildcards.genome}-{attempt}"
    log:
        f"{config['work_dir']}/logs/augustus_run/{{genome}}_pairwise.log"
    run:
        run_augustus_parallel(
            input=input, output=output, params=params, wildcards=wildcards,
            log_path=log[0],
            genome_work_dir=f"{params.work_dir}/augustus_parallel_temp_{wildcards.genome}_pairwise",
            with_tmr=True,
        )

rule augustus_run_tm_pairwise_only:
    """
    Launches the parallel augustus_parallel.py pipeline for genomes without RNA-seq,
    generating only the augTM_pairwise output using chromosome-based parallelization.
    Uses transMap_pairwise instead of transMap.
    """
    input:
        script="cat/augustus_parallel.py",
        fasta=f"{config['work_dir']}/genome_files/{{genome}}.fa",
        coding_gp=f"{config['work_dir']}/transMap_pairwise/{{genome}}_coding.gp",
        filtered_tm_psl=f"{config['work_dir']}/transMap_pairwise/{{genome}}_filtered.psl",
        ref_psl=f"{config['work_dir']}/reference/{config['ref_genome']}.psl",
        annotation_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp"
    output:
        tm_gtf=f"{config['work_dir']}/augustus/{{genome}}_augTM_pairwise.gtf",
        tm_gtf_done=f"{config['work_dir']}/augustus/{{genome}}_augTM_pairwise.gtf.done"
    wildcard_constraints:
        # Only run for genomes that do NOT have RNA-seq data AND when augustus is enabled
        genome = AUGUSTUS_NON_RNASEQ_GENOME_WC
    params:
        work_dir=config['work_dir'],
        species=config['augustus_species'],
        tm_cfg=config['tm_cfg_path'],
        utr=1 if config.get('predict_utr') else 0
    log:
        f"{config['work_dir']}/logs/augustus_run/{{genome}}_pairwise.log"
    threads: 1 if IS_CLUSTER else get_local_res("augustus_run_tm_pairwise_only", "threads")
    resources:
        mem_gb=4 if IS_CLUSTER else get_local_res("augustus_run_tm_pairwise_only", "mem_gb"),
        time_h=12 if IS_CLUSTER else get_local_res("augustus_run_tm_pairwise_only", "time_h"),
        job_id=lambda wildcards, attempt: f"augustus-parallel-tm-pairwise-{wildcards.genome}-{attempt}"
    run:
        run_augustus_parallel(
            input=input, output=output, params=params, wildcards=wildcards,
            log_path=log[0],
            genome_work_dir=f"{params.work_dir}/augustus_parallel_temp_{wildcards.genome}_pairwise",
            with_tmr=False,
        )

rule augustus_convert_tm_pairwise_gtf_to_gp:
    """
    Post-processing step to convert Augustus TM_pairwise GTF files to GenePred.
    Runs for each genome with Augustus TM_pairwise output.
    """
    input:
        gtf=f"{config['work_dir']}/augustus/{{genome}}_augTM_pairwise.gtf",
        gtf_done=f"{config['work_dir']}/augustus/{{genome}}_augTM_pairwise.gtf.done"
    output:
        gp=f"{config['work_dir']}/augustus/{{genome}}_augTM_pairwise.gp"
    wildcard_constraints:
        # Only run for genomes that have Augustus enabled
        genome = AUGUSTUS_GENOME_WC
    resources:
        mem_gb=16,
        time_h=2,
        job_id=lambda wildcards, attempt: f"aug-convert-tm-pairwise-{wildcards.genome}-{attempt}"
    log:
        f"{config['work_dir']}/logs/augustus_convert/{{genome}}_augTM_pairwise.log"
    shell:
        """
        if [ -s {input.gtf} ]; then
            gtfToGenePred -genePredExt {input.gtf} {output.gp} &> {log}
        else
            echo "Empty Augustus TM_pairwise GTF - creating empty GenePred" > {log}
            touch {output.gp}
        fi
        """

rule augustus_convert_tmr_pairwise_gtf_to_gp:
    """
    Post-processing step to convert Augustus TMR_pairwise GTF files to GenePred.
    Only runs for genomes with RNA-seq data (rnaseq_genomes).
    """
    input:
        gtf=f"{config['work_dir']}/augustus/{{genome}}_augTMR_pairwise.gtf",
        gtf_done=f"{config['work_dir']}/augustus/{{genome}}_augTMR_pairwise.gtf.done"
    output:
        gp=f"{config['work_dir']}/augustus/{{genome}}_augTMR_pairwise.gp"
    wildcard_constraints:
        # Only run for genomes that have RNA-seq data
        genome = RNASEQ_GENOME_WC
    resources:
        mem_gb=16,
        time_h=2,
        job_id=lambda wildcards, attempt: f"aug-convert-tmr-pairwise-{wildcards.genome}-{attempt}"
    log:
        f"{config['work_dir']}/logs/augustus_convert/{{genome}}_augTMR_pairwise.log"
    shell:
        """
        if [ -s {input.gtf} ]; then
            gtfToGenePred -genePredExt {input.gtf} {output.gp} &> {log}
        else
            echo "Empty Augustus TMR_pairwise GTF - creating empty GenePred" > {log}
            touch {output.gp}
        fi
        """

# ═══════════════════════════════════════════════════════════════════════════════
# AUGUSTUS MP MODE: Uses miniprot predictions as templates to recover missing genes
# ═══════════════════════════════════════════════════════════════════════════════

rule miniprot_paf_to_genepred:
    """
    Convert miniprot PAF output to GenePred + real protein-space PSL.

    Each PAF row becomes one record. Multi-exon CIGAR (cg:Z:) is parsed
    into proper exon boundaries; paralogs keep distinct ``_2``, ``_3``, ...
    copy suffixes. The PSL holds real miniprot-derived alignment statistics
    that downstream rules turn into truthful AlnCoverage / AlnIdentity for
    augMP transcripts (no more fake 100/100 metrics).

    Filtering knobs are exposed under ``config['miniprot']`` but default to
    "off" so every alignment carries its real metrics into consensus, where
    they can be filtered honestly alongside transMap / txTM hits:
        min_coverage  default 0.0   (no filter)
        min_identity  default 0.0   (no filter)
        min_mapq      default 0     (no filter)
        min_score     default 0     (no filter)
    """
    input:
        paf=f"{config['work_dir']}/miniprot/{{genome}}_miniprot.paf"
    output:
        gp=f"{config['work_dir']}/miniprot/{{genome}}_miniprot.gp",
        psl=f"{config['work_dir']}/miniprot/{{genome}}_miniprot.psl"
    wildcard_constraints:
        genome = AUGUSTUS_GENOME_WC
    log:
        f"{config['work_dir']}/logs/miniprot_to_gp/{{genome}}.log"
    resources:
        mem_gb=8,
        time_h=2,
        job_id=lambda wildcards, attempt: f"mp-to-gp-{wildcards.genome}-{attempt}"
    params:
        min_coverage=config.get('miniprot', {}).get('min_coverage', 0.0),
        min_identity=config.get('miniprot', {}).get('min_identity', 0.0),
        min_mapq=config.get('miniprot', {}).get('min_mapq', 0),
        min_score=config.get('miniprot', {}).get('min_score', 0),
    shell:
        """
        python3 cat/convert_miniprot_to_genepred.py {input.paf} {output.gp} \
            --psl {output.psl} \
            --min-coverage {params.min_coverage} \
            --min-identity {params.min_identity} \
            --min-mapq {params.min_mapq} \
            --min-score {params.min_score} \
            > {log} 2>&1
        """

rule augustus_run_mp:
    """
    Augustus MP mode: Uses all miniprot transcripts as templates (like augTM uses TransMap).
    Runs Augustus on every miniprot prediction to refine gene structures.
    """
    input:
        script="cat/augustus_parallel.py",
        fasta=f"{config['work_dir']}/genome_files/{{genome}}.fa",
        coding_gp=f"{config['work_dir']}/miniprot/{{genome}}_miniprot.gp",
        filtered_tm_psl=f"{config['work_dir']}/transMap/{{genome}}_filtered.psl",
        ref_psl=f"{config['work_dir']}/reference/{config['ref_genome']}.psl",
        annotation_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp",
        miniprot_hints=f"{config['work_dir']}/miniprot/{{genome}}_miniprot_hints.gff"
    output:
        mp_gtf=f"{config['work_dir']}/augustus/{{genome}}_augMP.gtf",
        mp_gtf_done=f"{config['work_dir']}/augustus/{{genome}}_augMP.gtf.done"
    wildcard_constraints:
        genome = AUGUSTUS_GENOME_WC
    params:
        work_dir=config['work_dir'],
        species=config['augustus_species'],
        tm_cfg=config['tm_cfg_path'],
        utr=1 if config.get('predict_utr') else 0
    threads: 1 if IS_CLUSTER else get_local_res("augustus_run_mp", "threads")
    resources:
        mem_gb=4 if IS_CLUSTER else get_local_res("augustus_run_mp", "mem_gb"),
        time_h=12 if IS_CLUSTER else get_local_res("augustus_run_mp", "time_h"),
        job_id=lambda wildcards, attempt: f"augustus-mp-{wildcards.genome}-{attempt}"
    log:
        f"{config['work_dir']}/logs/augustus_mp/{{genome}}.log"
    run:
        import os
        import subprocess

        import shutil

        genome_work_dir = f"{params.work_dir}/augustus_mp_temp_{wildcards.genome}"
        os.makedirs(genome_work_dir, exist_ok=True)

        # On SUCCESS the (large) per-genome temp dir is removed to avoid leaking disk.
        # On FAILURE it is intentionally kept: the per-task SLURM .err/.out logs live
        # inside it and are the only way to diagnose why augMP produced nothing.
        mp_succeeded = False
        try:
            # Check if there are any miniprot-only transcripts
            with open(input.coding_gp) as f:
                num_transcripts = sum(1 for line in f if line.strip())

            shell(f"echo 'Genome: {wildcards.genome}' > {log[0]}")
            shell(f"echo 'Miniprot-only transcripts: {num_transcripts}' >> {log[0]}")

            if num_transcripts == 0:
                shell(f"touch {output.mp_gtf}")
                shell(f"echo 'No miniprot-only transcripts for {wildcards.genome}' > {output.mp_gtf_done}")
                shell(f"echo 'No miniprot-only transcripts - created empty output' >> {log[0]}")
                mp_succeeded = True
            else:
                shell(f"echo 'Running Augustus MP for {num_transcripts} miniprot-only transcripts' >> {log[0]}")

                cmd = [
                    "python3", input.script,
                    "--genome_fasta", input.fasta,
                    "--coding_gp", input.coding_gp,
                    "--filtered_tm_psl", input.filtered_tm_psl,
                    "--ref_psl", input.ref_psl,
                    "--annotation_gp", input.annotation_gp,
                    "--tm_cfg", params.tm_cfg,
                    "--genome", wildcards.genome,
                    "--augustus_species", params.species,
                    "--utr", str(params.utr),
                    "--augustus_tm_gtf", output.mp_gtf,
                    "--miniprot_hints_gff", input.miniprot_hints,
                    "--work_dir", genome_work_dir
                ] + (["--no_slurm_preprocessing", "--no_slurm_transcripts"] if not IS_CLUSTER else _aug_slurm_args("augustus_tm"))

                result = subprocess.run(cmd, capture_output=True, text=True)

                with open(log[0], 'a') as log_file:
                    log_file.write("\n=== COMMAND ===\n")
                    log_file.write(" ".join(cmd) + "\n\n")
                    log_file.write("=== STDOUT ===\n")
                    log_file.write(result.stdout + "\n\n")
                    log_file.write("=== STDERR ===\n")
                    log_file.write(result.stderr + "\n")

                if result.returncode != 0:
                    raise Exception(f"Augustus MP failed for {wildcards.genome}")

                with open(output.mp_gtf_done, 'w') as f:
                    f.write(f"Completed Augustus MP for {wildcards.genome}\n")
                    f.write(f"Processed {num_transcripts} miniprot-only transcripts\n")
            mp_succeeded = True
        finally:
            if mp_succeeded and os.path.isdir(genome_work_dir):
                # Best-effort cleanup: the outputs are already written, so a flaky
                # rmtree (common on NFS: lingering .nfs* handles or lag raising
                # "Directory not empty") must NEVER fail the rule and discard the
                # good augMP results. Retry a few times, then give up quietly.
                import time as _time
                cleaned = False
                last_err = None
                for _attempt in range(5):
                    try:
                        shutil.rmtree(genome_work_dir)
                        cleaned = True
                        break
                    except OSError as _e:
                        last_err = _e
                        _time.sleep(3)
                if not cleaned:
                    shutil.rmtree(genome_work_dir, ignore_errors=True)
                with open(log[0], 'a') as log_file:
                    if cleaned:
                        log_file.write(f"\nCleaned up temp directory: {genome_work_dir}\n")
                    else:
                        log_file.write(
                            f"\nWARNING: could not fully remove temp directory "
                            f"{genome_work_dir} ({last_err}); left best-effort cleanup. "
                            f"augMP outputs are valid and complete.\n"
                        )
            elif os.path.isdir(genome_work_dir):
                with open(log[0], 'a') as log_file:
                    log_file.write(
                        f"\nAugustus MP FAILED — keeping temp directory for diagnosis: "
                        f"{genome_work_dir}\n"
                        f"Inspect per-task logs under "
                        f"{genome_work_dir}/slurm_transcript_temp_{wildcards.genome}/ "
                        f"(and augustus_parallel_temp/) for the batch failure reason.\n"
                    )

rule augustus_convert_mp_gtf_to_gp:
    """
    Convert Augustus MP GTF to GenePred for consensus building.
    """
    input:
        gtf=f"{config['work_dir']}/augustus/{{genome}}_augMP.gtf",
        gtf_done=f"{config['work_dir']}/augustus/{{genome}}_augMP.gtf.done"
    output:
        gp=f"{config['work_dir']}/augustus/{{genome}}_augMP.raw.gp"
    wildcard_constraints:
        genome = AUGUSTUS_GENOME_WC
    log:
        f"{config['work_dir']}/logs/augustus_convert/{{genome}}_augMP.log"
    resources:
        mem_gb=16,
        time_h=2,
        job_id=lambda wildcards, attempt: f"aug-convert-mp-{wildcards.genome}-{attempt}"
    shell:
        """
        if [ -s {input.gtf} ]; then
            gtfToGenePred -genePredExt {input.gtf} {output.gp} &> {log}
        else
            echo "Empty Augustus MP GTF - creating empty GenePred" > {log}
            touch {output.gp}
        fi
        """

rule fix_augmp_gene_names:
    """
    Fix augMP genePred files by replacing transcript IDs in name2 field with gene names.
    Uses reference database to map transcript IDs to gene names.
    """
    input:
        gp=f"{config['work_dir']}/augustus/{{genome}}_augMP.raw.gp",
        ref_db=f"{config['work_dir']}/databases/{config['ref_genome']}.db"
    output:
        done=touch(f"{config['work_dir']}/databases/{{genome}}_augMP_gene_names_fixed.done")
    wildcard_constraints:
        genome = AUGUSTUS_GENOME_WC
    log:
        f"{config['work_dir']}/logs/fix_augmp_gene_names/{{genome}}.log"
    resources:
        mem_gb=4,
        time_h=1,
        job_id=lambda wildcards, attempt: f"fix-augmp-{wildcards.genome}-{attempt}"
    shell:
        """
        python cat/fix_augmp_gene_names.py \
            --ref-db {input.ref_db} \
            --augmp-files {input.gp} \
            >> {log} 2>&1 && touch {output.done}
        """

# ═══════════════════════════════════════════════════════════════════════════════

rule run_augustus_pb:
    input:
        script="cat/augustus_pb_parallel.py",
        fasta=f"{config['work_dir']}/genome_files/{{genome}}.fa",
        sizes=f"{config['work_dir']}/genome_files/{{genome}}.chrom.sizes",
        hints_gff=f"{config['work_dir']}/hints_database/{{genome}}_extrinsic_hints.gff"
    output:
        raw_gtf=f"{config['work_dir']}/augustus_pb/{{genome}}_raw_augPB.gtf",
        gtf=f"{config['work_dir']}/augustus_pb/{{genome}}_augPB.gtf",
        gp=f"{config['work_dir']}/augustus_pb/{{genome}}_augPB.gp"
    wildcard_constraints:
        genome = AUGUSTUS_PB_GENOME_WC
    params:
        work_dir=config['work_dir'],
        species=config['augustus_species'],
        pb_cfg=config.get('pb_cfg_path', ''),
        utr=1 if config.get('predict_utr') else 0,
        chunksize=config.get('pb_genome_chunksize', 11000000),
        overlap=config.get('pb_genome_overlap', 1000000)
    log:
        f"{config['work_dir']}/logs/augustus_pb/{{genome}}.log"
    threads: 1 if IS_CLUSTER else get_local_res("run_augustus_pb", "threads")
    resources:
        mem_gb=4 if IS_CLUSTER else get_local_res("run_augustus_pb", "mem_gb"),
        time_h=12 if IS_CLUSTER else get_local_res("run_augustus_pb", "time_h"),
        job_id=lambda wildcards, attempt: f"augustus-pb-parallel-{wildcards.genome}-{attempt}"
    run:
        # Create temporary work directory for this genome
        genome_work_dir = f"{params.work_dir}/augustus_pb_parallel_temp_{wildcards.genome}"
        os.makedirs(genome_work_dir, exist_ok=True)

        # Run the parallel Augustus PB pipeline
        cmd = [
            "python", input.script,
            "--genome_fasta", input.fasta,
            "--chrom_sizes", input.sizes,
            "--hints_gff", input.hints_gff,
            "--pb_cfg", params.pb_cfg,
            "--raw_gtf", output.raw_gtf,
            "--gtf", output.gtf,
            "--gp", output.gp,
            "--species", params.species,
            "--utr", str(params.utr),
            "--chunksize", str(params.chunksize),
            "--overlap", str(params.overlap),
            "--work_dir", genome_work_dir
        ] + (["--no_slurm_preprocessing", "--no_slurm_jobs"] if not IS_CLUSTER else _aug_slurm_args("augustus_pb"))

        # try/finally so the per-genome temp dir is always cleaned up, even on
        # failure — otherwise failed/retried jobs leak (large) temp dirs.
        try:
            shell(" ".join(cmd) + f" > {log} 2>&1")
        finally:
            if os.path.exists(genome_work_dir):
                shell(f"rm -rf {genome_work_dir}")

rule convert_stringtie_to_strg:
    """
    Convert StringTie GTF to a clean strg-prefixed GTF. This produces
    a standalone strg mode that feeds into the consensus pipeline.
    Only runs for genomes with isoseq data and stringtie enabled.
    """
    input:
        script="cat/convert_stringtie_to_augpb_format.py",
        temp_gtf=f"{config['work_dir']}/stringtie/{{genome}}_temp_stringtie.gtf"
    output:
        gtf=f"{config['work_dir']}/stringtie/{{genome}}_strg.gtf"
    wildcard_constraints:
        genome = STRG_GENOME_WC
    log:
        f"{config['work_dir']}/logs/stringtie_run/{{genome}}_strg_convert.log"
    threads: 1 if IS_CLUSTER else get_local_res("stringtie_convert", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("stringtie_convert", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("stringtie_convert", "time_h"),
        job_id=lambda wildcards, attempt: f"strg-convert-{wildcards.genome}-{attempt}"
    run:
        work_dir = config['work_dir']
        genome = wildcards.genome
        job_script = f"{work_dir}/stringtie/{genome}_strg_convert_job.sh"

        script_content = build_sbatch_header(
            "stringtie_convert",
            f"strg-convert-{genome}",
            f"{work_dir}/logs/stringtie_run/{genome}_strg_convert_slurm.out",
            f"{work_dir}/logs/stringtie_run/{genome}_strg_convert_slurm.err",
        ) + f"""
python3 {input.script} {input.temp_gtf} {output.gtf} >> {log[0]} 2>&1
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting StringTie->strg conversion job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.gtf],
                log_file,
                "stringtie_convert",
                max_wait_s=3 * 3600,
            )

rule convert_strg_gtf_to_gp:
    """Convert the strg GTF to GenePred for consensus building."""
    input:
        gtf=f"{config['work_dir']}/stringtie/{{genome}}_strg.gtf"
    output:
        gp=f"{config['work_dir']}/stringtie/{{genome}}_strg.gp"
    wildcard_constraints:
        genome = STRG_GENOME_WC
    log:
        f"{config['work_dir']}/logs/stringtie_run/{{genome}}_strg_gp.log"
    threads: 1 if IS_CLUSTER else get_local_res("stringtie_gp", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("stringtie_gp", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("stringtie_gp", "time_h"),
        job_id=lambda wildcards, attempt: f"strg-gp-{wildcards.genome}-{attempt}"
    run:
        work_dir = config['work_dir']
        genome = wildcards.genome
        job_script = f"{work_dir}/stringtie/{genome}_strg_gp_job.sh"

        script_content = build_sbatch_header(
            "stringtie_gp",
            f"strg-gp-{genome}",
            f"{work_dir}/logs/stringtie_run/{genome}_strg_gp_slurm.out",
            f"{work_dir}/logs/stringtie_run/{genome}_strg_gp_slurm.err",
        ) + f"""
if [ -s {input.gtf} ]; then
    gtfToGenePred -genePredExt {input.gtf} {output.gp} >> {log[0]} 2>&1
else
    touch {output.gp}
    echo "Empty input GTF, created empty GP" >> {log[0]}
fi
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting strg GTF->GP conversion job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.gp],
                log_file,
                "stringtie_gp",
                max_wait_s=3 * 3600,
            )

def get_denovo_gp(wildcards):
    """Helper function to find the correct de novo genePred based on the mode."""
    gp_map = {
        "augPB": f"{config['work_dir']}/augustus_pb/{wildcards.genome}_augPB.gp",
        "strg": f"{config['work_dir']}/stringtie/{wildcards.genome}_strg.gp",
    }
    if wildcards.mode not in gp_map:
        raise ValueError(f"Invalid mode for parent gene assignment: {wildcards.mode}")
    return gp_map[wildcards.mode]

def get_denovo_parent_tablename(wildcards):
    """Returns the correct SQL table name based on the mode."""
    table_map = {
        "augPB": "AugPbAlternativeGenes",
        "strg": "StrgAlternativeGenes",
    }
    if wildcards.mode not in table_map:
        raise ValueError(f"Invalid mode for parent gene assignment: {wildcards.mode}")
    return table_map[wildcards.mode]

rule find_denovo_parents:
    """
    Assigns parental genes to de novo transcripts from AugustusPB or StringTie
    by comparing them to TransMap projections.
    """
    input:
        denovo_gp=get_denovo_gp,
        filtered_tm_gp=f"{config['work_dir']}/transMap/{{genome}}_filtered.gp",
        unfiltered_tm_gp=f"{config['work_dir']}/transMap/{{genome}}.gp",
        sizes=f"{config['work_dir']}/genome_files/{{genome}}.chrom.sizes",
    output:
        done_file=temp(f"{config['work_dir']}/databases/{{genome}}_{{mode}}_parents.done")
    params:
        script="cat/parent_gene_assignment_cluster.py" if IS_CLUSTER else "cat/parent_gene_assignment.py",
        db_path=f"{config['work_dir']}/databases/{{genome}}.db",
        table_name=get_denovo_parent_tablename,
        cluster_args=(
            f"--execution-mode {EXECUTION_MODE} "
            f"--partition {_slurm_partition('find_denovo_parents')} "
            f"--exclude-nodes \"{_slurm_exclude()}\" "
            f"--module-load \"{_slurm_module()}\" "
            f"--sge-parallel-env {config.get('cluster', {}).get('sge', {}).get('parallel_env', 'smp')} "
            f"--sge-memory-flag {config.get('cluster', {}).get('sge', {}).get('memory_flag', 'h_vmem')} "
            f"--memory {int(get_res('find_denovo_parents', 'mem').rstrip('G'))} "
            f"--cpus {get_res('find_denovo_parents', 'cpus')} "
            f"--time {get_res('find_denovo_parents', 'time')} "
            f"--max-jobs {get_res('find_denovo_parents', 'max_concurrent_jobs')} "
            f"--cleanup"
        ) if IS_CLUSTER else ""
    log:
        f"{config['work_dir']}/logs/find_denovo_parents/{{genome}}_{{mode}}.log"
    threads: 1 if IS_CLUSTER else get_local_res("find_denovo_parents", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("find_denovo_parents", "mem_gb"),
        time_h=24,
        job_id=lambda wildcards, attempt: f"parents-{wildcards.genome}-{wildcards.mode}-{attempt}"
    shell:
        """
        python {params.script} \
            --filtered-tm-gp {input.filtered_tm_gp} \
            --unfiltered-tm-gp {input.unfiltered_tm_gp} \
            --chrom-sizes {input.sizes} \
            --denovo-gp {input.denovo_gp} \
            --db-path {params.db_path} \
            --table-name {params.table_name} \
            {params.cluster_args} > {log} 2>&1 && \
        touch {output.done_file}
        """

rule isoseq_construct_structures:
    input:
        # Depends on the per-genome hints file from the generate_extrinsic_hints rule
        hints_gff=f"{config['work_dir']}/hints_database/{{genome}}_extrinsic_hints.gff"
    output:
        done_file=touch(f"{config['work_dir']}/databases/{{genome}}_isoseq_structures.done")
    params:
        db_path=f"{config['work_dir']}/databases/{{genome}}.db"
    log:
        f"{config['work_dir']}/logs/isoseq_structures/{{genome}}.log"
    script:
        "cat/isoseq_transcripts_wrapper.py"

def get_active_modes_for_wildcards(wildcards):
    genome = wildcards.genome
    if genome in ANCESTOR_GENOMES:
        return [m for m in ANCESTOR_MODES if m in VALID_ANCESTOR_MODES]

    modes = ['transMap', 'transMap_pairwise'] # transMap and transMap_pairwise are always active for target genomes

    # Add modes based on the main config flags
    if config.get("txTM", False):
        modes.append('txTM')
    if config.get("augustus", False):
        modes.append('augTM')
        modes.append('augTM_pairwise')
        if genome in config.get("rnaseq_genomes", []):
            modes.append('augTMR')
            modes.append('augTMR_pairwise')
        modes.append('augMP')
    if config.get("augustus_pb", False) and genome in config.get("isoseq_genomes", []):
        modes.append('augPB')
    if config.get("stringtie", False) and genome in config.get("isoseq_genomes", []) and genome in config.get("stringtie_genomes", TARGET_GENOMES):
        modes.append('strg')
        
    return modes


def get_alignment_modes_for_genome(wildcards):
    """Per-genome alignment modes that require align_transcripts/evaluate_transcripts."""
    return [m for m in get_active_modes_for_wildcards(wildcards) if m in ALL_ALIGNMENT_MODES]


def get_align_transcript_outputs(wildcards):
    """Determines the expected output PSL files for the align_transcripts rule."""
    work_dir = config["work_dir"]
    output_dir=f"{work_dir}/transcript_alignment"
    active_modes = get_active_modes_for_wildcards(wildcards)
    
    psl_files = []
    for mode in active_modes:
        psl_files.append(f"{output_dir}/{wildcards.genome}_{mode}_mRNA.psl")
        psl_files.append(f"{output_dir}/{wildcards.genome}_{mode}_CDS.psl")
        
    return psl_files

def get_gps_for_alignment(wildcards):
    """Helper to get all GP paths for a genome's active modes."""
    # This dictionary provides a clear mapping from mode to file path
    gp_path_map = mode_gp_paths(wildcards.genome)
    return [gp_path_map[mode] for mode in get_active_modes_for_wildcards(wildcards)]

def get_gp_path_for_mode(wildcards):
    """Helper to get the genePred path for a single, specific transcript mode."""
    return mode_gp_paths(wildcards.genome)[wildcards.alignment_mode]


rule align_transcripts:
    input:
        script="cat/align_transcripts.py",
        ref_db=f"{config['work_dir']}/databases/{config['ref_genome']}.db",
        ref_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp",
        ref_fasta=f"{config['work_dir']}/genome_files/{config['ref_genome']}.fa",
        genome_fasta=f"{config['work_dir']}/genome_files/{{genome}}.fa",
        # Get the specific genePred for this single mode
        gp=get_gp_path_for_mode,
        setup_done=rules.setup_pipeline_directories.output.setup_done,
        # Ensure augMP gene names are fixed before alignment when using augMP mode
        augmp_fix=lambda w: f"{config['work_dir']}/databases/{w.genome}_augMP_gene_names_fixed.done" if w.alignment_mode == 'augMP' else config['work_dir'] + "/.setup_done"
    output:
        mrna_psl=f"{config['work_dir']}/transcript_alignment/{{genome}}_{{alignment_mode}}_mRNA.psl",
        cds_psl=f"{config['work_dir']}/transcript_alignment/{{genome}}_{{alignment_mode}}_CDS.psl"
    wildcard_constraints:
        # Only include modes that are actually enabled in the config
        genome=ANNOTATION_GENOME_WC,
        alignment_mode=f"({'|'.join(ALL_ALIGNMENT_MODES)})"
    params:
        mode_file_args=lambda w, input, output: (
            f"--mode-files {w.alignment_mode} {input.gp} {output.mrna_psl} {output.cds_psl}"
        ),
        job_store=f"file:{config['work_dir']}/toil_job_stores/align_transcripts_{{genome}}_{{alignment_mode}}",
        job_store_dir=f"{config['work_dir']}/toil_job_stores/align_transcripts_{{genome}}_{{alignment_mode}}",
        work_dir=config['work_dir'],
        genome_name=f"{{genome}}",
        use_cluster=IS_CLUSTER,
        execution_mode=EXECUTION_MODE,
        exclude_nodes=_slurm_exclude(),
        module_load=_slurm_module(),
        sge_parallel_env=config.get('cluster', {}).get('sge', {}).get('parallel_env', 'smp'),
        sge_memory_flag=config.get('cluster', {}).get('sge', {}).get('memory_flag', 'h_vmem'),
        memory=int(get_res("align_transcripts", "mem").rstrip("G")),
        cpus=get_res("align_transcripts", "cpus"),
        walltime=get_res("align_transcripts", "time"),
        partition=_slurm_partition("align_transcripts"),
        max_jobs=get_res("align_transcripts", "max_concurrent_jobs"),
        timeout_hours=get_res("align_transcripts", "timeout_hours"),
        chunk_size=get_res("align_transcripts", "chunk_size")
    log:
        f"{config['work_dir']}/logs/align_transcripts/{{genome}}_{{alignment_mode}}.log"
    threads: 1 if IS_CLUSTER else get_local_res("align_transcripts", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("align_transcripts", "mem_gb"),
        time_h=24,
        job_id=lambda wildcards, attempt: f"align-{wildcards.genome}-{wildcards.alignment_mode}-{attempt}"
    shell:
        r"""
        if [ "{params.use_cluster}" = "True" ]; then
            echo "Using {params.execution_mode}-based transcript alignment for {wildcards.genome} {wildcards.alignment_mode}" > {log}
            python {input.script} \
              --mode cluster \
              --execution-mode {params.execution_mode} \
              --genome {params.genome_name} \
              --ref-genome-fasta {input.ref_fasta} \
              --genome-fasta {input.genome_fasta} \
              --annotation-gp {input.ref_gp} \
              --ref-db-path {input.ref_db} \
              {params.mode_file_args} \
              --partition {params.partition} \
              --exclude-nodes "{params.exclude_nodes}" \
              --module-load "{params.module_load}" \
              --sge-parallel-env {params.sge_parallel_env} \
              --sge-memory-flag {params.sge_memory_flag} \
              --memory {params.memory} \
              --cpus {params.cpus} \
              --time {params.walltime} \
              --max-jobs {params.max_jobs} \
              --timeout-hours {params.timeout_hours} \
              --chunk-size {params.chunk_size} \
              --cleanup >> {log} 2>&1
        else
            echo "Using local transcript alignment for {wildcards.genome} {wildcards.alignment_mode}" > {log}
            rm -rf {params.job_store_dir}
            python {input.script} \
              --mode toil \
              {params.job_store} \
              --workDir {params.work_dir} \
              --batchSystem single_machine \
              --maxCores {threads} \
              --logFile {log} \
              --ref-genome-fasta {input.ref_fasta} \
              --genome-fasta {input.genome_fasta} \
              --annotation-gp {input.ref_gp} \
              --ref-db-path {input.ref_db} \
              --genome {params.genome_name} \
              {params.mode_file_args} >> {log} 2>&1
        fi
        """

rule evaluate_transcripts:
    """
    Evaluates transcript predictions for each mode (transMap, augTM, etc.)
    by comparing alignments to the reference annotation. This rule runs
    once per-genome, per-mode.
    """
    input:
        ref_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp",
        target_fasta=f"{config['work_dir']}/genome_files/{{genome}}.fa",
        gp=get_gp_path_for_mode,
        mrna_psl=f"{config['work_dir']}/transcript_alignment/{{genome}}_{{alignment_mode}}_mRNA.psl",
        cds_psl=f"{config['work_dir']}/transcript_alignment/{{genome}}_{{alignment_mode}}_CDS.psl"
    output:
        # Create a done file for each mode to signal completion
        done_file=temp(f"{config['work_dir']}/databases/{{genome}}_{{alignment_mode}}_evaluation.done"),
        resolved_df=temp(f"{config['work_dir']}/databases/{{genome}}_{{alignment_mode}}_resolved_df.pkl")
    wildcard_constraints:
        genome=ANNOTATION_GENOME_WC,
        alignment_mode=f"({'|'.join(ALL_ALIGNMENT_MODES)})"
    params:
        script="cat/classify_cluster.py" if IS_CLUSTER else "cat/classify.py",
        db_path=f"{config['work_dir']}/databases/{{genome}}.db",
        ref_db_path=f"{config['work_dir']}/databases/{config['ref_genome']}.db",
        mode_file_args=lambda w, input: f"{w.alignment_mode} {input.gp} {input.mrna_psl} {input.cds_psl}",
        cluster_args=(
            f"--execution-mode {EXECUTION_MODE} "
            f"--partition {_slurm_partition('evaluate_transcripts')} "
            f"--exclude-nodes \"{_slurm_exclude()}\" "
            f"--module-load \"{_slurm_module()}\" "
            f"--sge-parallel-env {config.get('cluster', {}).get('sge', {}).get('parallel_env', 'smp')} "
            f"--sge-memory-flag {config.get('cluster', {}).get('sge', {}).get('memory_flag', 'h_vmem')} "
            f"--memory {int(get_res('evaluate_transcripts', 'mem').rstrip('G'))} "
            f"--cpus {get_res('evaluate_transcripts', 'cpus')} "
            f"--time {get_res('evaluate_transcripts', 'time')} "
            f"--max-jobs {get_res('evaluate_transcripts', 'max_concurrent_jobs')} "
            f"--chunk-size {get_res('evaluate_transcripts', 'chunk_size')} "
            f"--cleanup"
        ) if IS_CLUSTER else ""
    log:
        f"{config['work_dir']}/logs/evaluate_transcripts/{{genome}}_{{alignment_mode}}.log"
    threads: 1 if IS_CLUSTER else get_local_res("evaluate_transcripts", "threads")
    resources:
        mem_gb=16 if IS_CLUSTER else get_local_res("evaluate_transcripts", "mem_gb"),
        time_h=24,
        job_id=lambda wildcards, attempt: f"eval-{wildcards.genome}-{wildcards.alignment_mode}-{attempt}"
    shell:
        """
        python {params.script} \\
            --annotation-gp {input.ref_gp} \\
            --ref-db-path {params.ref_db_path} \\
            --fasta {input.target_fasta} \\
            --db-path {params.db_path} \\
            --resolved-df {output.resolved_df} \\
            --mode-files {params.mode_file_args} \\
            {params.cluster_args}

        python3 - <<PY
import pandas as pd
import os
from tools.sqlite import ExclusiveSqlConnection
genome_db = "{params.db_path}"
pickle_file = "{output.resolved_df}"

# Check if pickle file exists before trying to read it
if not os.path.exists(pickle_file):
    print(f"Error: Pickle file {{pickle_file}} does not exist")
    exit(1)

try:
    results = pd.read_pickle(pickle_file)
    for table_name, df in results:
        # Always write the table so it exists even when empty (consensus expects all evaluation tables).
        with ExclusiveSqlConnection(genome_db) as engine:
            df.to_sql(table_name, engine, if_exists="replace", index=True)
except Exception as e:
    print(f"Error processing pickle file {{pickle_file}}: {{e}}")
    exit(1)
PY
        touch {output.done_file}
        """

def get_evaluation_done_inputs(wildcards):
    """
    Collect only the evaluation.done files for the supported alignment modes.
    Uses per-genome active modes so e.g. augTMR is not required for genomes
    without RNA-seq (bosTau8), while ancestors only get ancestor_modes.
    """
    alignment_modes = get_alignment_modes_for_genome(wildcards)
    return expand(
        f"{config['work_dir']}/databases/{{genome}}_{{alignment_mode}}_evaluation.done",
        work_dir=config["work_dir"],
        genome=wildcards.genome,
        alignment_mode=alignment_modes
    )

rule aggregate_evaluations:
    input:
        get_evaluation_done_inputs
    output:
        # This is the exact file that the generate_consensus rule is missing
        touch(f"{config['work_dir']}/logs/{{genome}}_evaluation.done")
    wildcard_constraints:
        genome=ANNOTATION_GENOME_WC,
    resources:
        mem_gb=128,
        time_h=1,
        job_id=lambda wildcards, attempt: f"agg-eval-{wildcards.genome}-{attempt}"
    shell:
        # The presence of the input files is enough, we just need to create the
        # output file to resolve the dependency for the next step.
        "echo 'All evaluation modes for {wildcards.genome} are complete.'"


def get_consensus_inputs(wildcards):
    """Dynamically determines all input files required for consensus generation."""
    work_dir = config["work_dir"]
    genome = wildcards.genome
    active_modes = get_active_modes_for_wildcards(wildcards) # Use the single source of truth

    # Base inputs that are always required
    inputs = {
        # ancient(): the per-genome DB is a shared mutable container; consensus
        # freshness is signalled by the explicit .done markers below (tm_eval_done,
        # parent_done, psl_metrics_done), so the DB's bumped mtime must not force a
        # spurious consensus rerun.
        'db_path': ancient(f"{work_dir}/databases/{genome}.db"),
        'ref_db_path': f"{work_dir}/databases/{config['ref_genome']}.db",
        'fasta': f"{work_dir}/genome_files/{genome}.fa",
        'ref_gp': f"{work_dir}/reference/{config['ref_genome']}.gp",
        'tm_eval_done': f"{work_dir}/databases/{genome}.tm_eval.done",
        'eval_done': f"{work_dir}/logs/{genome}_evaluation.done",
        # Pairwise map + filter must finish before consensus; consensus uses unfiltered .gp
        # (filter --filter-overlapping-genes drops many ref PC genes e.g. at collapsed loci).
        'transmap_pairwise_gp': f"{work_dir}/transMap_pairwise/{genome}.gp",
        'transmap_pairwise_filtered': f"{work_dir}/transMap_pairwise/{genome}_filtered.gp",
    }

    # Map each mode to its specific genePred file path
    gp_path_map = mode_gp_paths(genome, work_dir)
    # Consensus intentionally consumes the *unfiltered* pairwise GP (not _filtered.gp).
    gp_path_map['transMap_pairwise'] = f"{work_dir}/transMap_pairwise/{genome}.gp"
    gp_list = [gp_path_map[mode] for mode in active_modes]
    transmap_pairwise_gp = gp_path_map['transMap_pairwise']
    if transmap_pairwise_gp not in gp_list:
        gp_list.append(transmap_pairwise_gp)
    inputs['gp_list'] = gp_list

    parent_modes = {'augPB', 'strg'}
    inputs['parent_done'] = [f"{work_dir}/databases/{genome}_{mode}_parents.done"
                             for mode in active_modes if mode in parent_modes]
    
    # augMP requires fix_augmp_gene_names to run before consensus uses the GP file
    inputs['augmp_fix_done'] = [f"{work_dir}/databases/{genome}_augMP_gene_names_fixed.done"] if 'augMP' in active_modes else [f"{work_dir}/.setup_done"]
    
    # BAM files for computing real RNA-seq support
    genome_data = config.get("transcriptomic_data", {}).get(genome, {})
    bam_files = genome_data.get("bam", []) + genome_data.get("intronbam", [])
    isoseq_bam_files = genome_data.get("isoseq_bam", [])
    if bam_files:
        inputs['rnaseq_bams'] = bam_files
    if isoseq_bam_files:
        inputs['isoseq_bams'] = isoseq_bam_files

    # Ensure PSL-derived metrics are present for consensus (improves fragment filtering
    # and makes behavior closer to (or better than) txTM across modes).
    work_dir = config["work_dir"]
    genome = wildcards.genome
    done_files = []
    if genome in ANCESTOR_GENOMES:
        if "txTM" in ANCESTOR_MODES:
            done_files.append(f"{work_dir}/databases/{genome}_txTM_psl_metrics.done")
        if "transMap" in ANCESTOR_MODES:
            done_files.append(f"{work_dir}/databases/{genome}_transMap_psl_metrics.done")
        if "transMap_pairwise" in ANCESTOR_MODES:
            done_files.append(f"{work_dir}/databases/{genome}_transMap_pairwise_psl_metrics.done")
    else:
        if config.get("txTM", False):
            done_files.append(f"{work_dir}/databases/{genome}_txTM_psl_metrics.done")
        done_files.append(f"{work_dir}/databases/{genome}_transMap_psl_metrics.done")
        done_files.append(f"{work_dir}/databases/{genome}_transMap_pairwise_psl_metrics.done")
        # augMP has no native DB metrics tables; add the real miniprot-PSL-derived
        # coverage/identity (see generate_augMP_psl) so consensus filtering treats
        # augMP like other modes.
        if "augMP" in active_modes:
            done_files.append(f"{work_dir}/databases/{genome}_augMP_psl_metrics.done")
    inputs["psl_metrics_done"] = done_files
    
    return inputs


def psl_metrics_prior_after_txTM(wildcards):
    """Order PSL metric DB writes: transMap runs after txTM (or setup) to avoid sqlite lock errors."""
    if wildcards.genome in ANCESTOR_GENOMES:
        if "txTM" in ANCESTOR_MODES:
            return f"{config['work_dir']}/databases/{wildcards.genome}_txTM_psl_metrics.done"
        return f"{config['work_dir']}/.setup_done"
    if config.get("txTM", False):
        return f"{config['work_dir']}/databases/{wildcards.genome}_txTM_psl_metrics.done"
    return f"{config['work_dir']}/.setup_done"


rule store_psl_metrics_txTM:
    """
    Store PSL-derived AlnCoverage/AlnIdentity for txTM-mode transcripts.
    This avoids sparse DB metrics when sequence-pair evaluation is skipped.
    """
    input:
        psl=f"{config['work_dir']}/txTM/{{genome}}_filtered.psl",
        # The per-genome DB is a shared, mutable container written by many later
        # stages; freshness is tracked via the per-stage *_psl_metrics.done markers,
        # not the DB mtime. ancient() stops its bumped mtime from falsely marking
        # this rule stale and cascading a consensus rerun.
        db_path=ancient(f"{config['work_dir']}/databases/{{genome}}.db"),
        ref_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp"
    output:
        done=touch(f"{config['work_dir']}/databases/{{genome}}_txTM_psl_metrics.done")
    wildcard_constraints:
        genome = TXTM_GENOME_WC
    run:
        shell(
            "python3 cat/store_psl_metrics.py "
            f"--db-path {input.db_path} "
            f"--psl {input.psl} "
            "--mode txTM "
            f"--ref-gp {input.ref_gp} "
            f"&& touch {output.done}"
        )


rule store_psl_metrics_transMap:
    """Store PSL-derived AlnCoverage/AlnIdentity for transMap-mode transcripts."""
    input:
        psl=f"{config['work_dir']}/transMap/{{genome}}_filtered.psl",
        # ancient(): DB mtime is not a valid rerun trigger (shared mutable container);
        # freshness is tracked by the *_psl_metrics.done markers.
        db_path=ancient(f"{config['work_dir']}/databases/{{genome}}.db"),
        ref_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp",
        prior_txTM=psl_metrics_prior_after_txTM,
    output:
        done=touch(f"{config['work_dir']}/databases/{{genome}}_transMap_psl_metrics.done")
    wildcard_constraints:
        genome = ANNOTATION_GENOME_WC
    run:
        shell(
            "python3 cat/store_psl_metrics.py "
            f"--db-path {input.db_path} "
            f"--psl {input.psl} "
            "--mode transMap "
            f"--ref-gp {input.ref_gp} "
            f"&& touch {output.done}"
        )


rule store_psl_metrics_transMap_pairwise:
    """Store PSL-derived AlnCoverage/AlnIdentity for transMap_pairwise transcripts."""
    input:
        psl=f"{config['work_dir']}/transMap_pairwise/{{genome}}_filtered.psl",
        # ancient(): DB mtime is not a valid rerun trigger (shared mutable container);
        # freshness is tracked by the *_psl_metrics.done markers.
        db_path=ancient(f"{config['work_dir']}/databases/{{genome}}.db"),
        ref_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp",
        prior_transMap=f"{config['work_dir']}/databases/{{genome}}_transMap_psl_metrics.done",
    output:
        done=touch(f"{config['work_dir']}/databases/{{genome}}_transMap_pairwise_psl_metrics.done")
    wildcard_constraints:
        genome = ANNOTATION_GENOME_WC
    run:
        shell(
            "python3 cat/store_psl_metrics.py "
            f"--db-path {input.db_path} "
            f"--psl {input.psl} "
            "--mode transMap_pairwise "
            f"--ref-gp {input.ref_gp} "
            f"&& touch {output.done}"
        )

rule generate_augMP_psl:
    """
    Generate a REAL PSL for augMP transcripts by mapping each
    ``augMP-<source>[_<copy>]`` AlignmentId back to its miniprot source PSL row
    (real alignment statistics from the cg:Z: CIGAR + AS/np/ms tags).

    Replaces the previous ``generate_augMP_fake_psl`` rule, which used
    ``genePredToFakePsl`` to assign every augMP record 100 % cov / 100 % id.
    augMP records that have no matching miniprot PSL row are simply omitted
    from the output — they will stay NaN in the DB and be filtered honestly by
    the consensus stage (no synthetic perfect metrics anywhere).

    Must run after fix_augmp_gene_names so the augMP genePred is populated
    (an early run against an empty GP produced a 0-byte PSL and stale metrics).
    """
    input:
        gp=f"{config['work_dir']}/augustus/{{genome}}_augMP.raw.gp",
        miniprot_psl=f"{config['work_dir']}/miniprot/{{genome}}_miniprot.psl",
        augmp_fix=f"{config['work_dir']}/databases/{{genome}}_augMP_gene_names_fixed.done",
    output:
        psl=f"{config['work_dir']}/augustus/{{genome}}_augMP.raw.psl",
        # Sentinel so store_psl_metrics_augMP cannot run on a stale 0-byte PSL file.
        generated=f"{config['work_dir']}/databases/{{genome}}_augMP_psl_generated.done",
    wildcard_constraints:
        genome = "|".join(TARGET_GENOMES)
    log:
        f"{config['work_dir']}/logs/generate_augMP_psl/{{genome}}.log"
    shell:
        """
        rm -f {output.psl}
        python3 cat/generate_augMP_psl.py \
            --augmp-gp {input.gp} \
            --miniprot-psl {input.miniprot_psl} \
            --out-psl {output.psl} > {log} 2>&1
        n_gp=$(grep -c '^augMP-' {input.gp} || true)
        n_psl=$(wc -l < {output.psl} | tr -d ' ')
        if [ "$n_gp" -gt 0 ] && [ "$n_psl" -eq 0 ]; then
            echo "ERROR: {output.psl} is empty but {input.gp} has $n_gp augMP records" >> {log}
            exit 1
        fi
        echo "gp_augMP=$n_gp psl_lines=$n_psl" > {output.generated}
        """


rule filter_augMP:
    """Collapse within-locus augMP redundancy BEFORE consensus (filter_augMP.py).

    miniprot is run permissively so every paralog / lineage-specific copy is
    found; that means one genomic copy is hit by many homologous DB proteins, each
    producing its own near-identical augMP model. Those pileups (up to ~500 models
    on a single locus, 826k genome-wide) make consensus intractable. This rule
    keeps the best model per distinct gene structure per locus, so distinct copies
    (distinct coordinates/structures) are all preserved while pure redundancy is
    removed. Produces the canonical ``_augMP.gp``/``_augMP.psl`` consumed by
    align_transcripts, evaluate_transcripts, store_psl_metrics_augMP and consensus.
    """
    input:
        raw_gp=f"{config['work_dir']}/augustus/{{genome}}_augMP.raw.gp",
        raw_psl=f"{config['work_dir']}/augustus/{{genome}}_augMP.raw.psl",
        generated=f"{config['work_dir']}/databases/{{genome}}_augMP_psl_generated.done",
        augmp_fix=f"{config['work_dir']}/databases/{{genome}}_augMP_gene_names_fixed.done",
    output:
        gp=f"{config['work_dir']}/augustus/{{genome}}_augMP.gp",
        psl=f"{config['work_dir']}/augustus/{{genome}}_augMP.psl",
    wildcard_constraints:
        genome = "|".join(TARGET_GENOMES)
    params:
        disabled_flag="" if rcfg("augMP_filter_enabled", True) else "--disabled",
        max_per_locus=rcfg("augMP_filter_max_models_per_locus", 25),
        structure_round=config.get("augMP_filter_structure_round_bp", 10),
        single_round=config.get("augMP_filter_single_exon_round_bp", 50),
        single_min_cov=config.get("augMP_filter_single_exon_min_coverage", 0.0),
    log:
        f"{config['work_dir']}/logs/filter_augMP/{{genome}}.log"
    shell:
        """
        python3 cat/filter_augMP.py \
            --in-gp {input.raw_gp} \
            --in-psl {input.raw_psl} \
            --out-gp {output.gp} \
            --out-psl {output.psl} \
            --max-models-per-locus {params.max_per_locus} \
            --structure-round-bp {params.structure_round} \
            --single-exon-round-bp {params.single_round} \
            --single-exon-min-coverage {params.single_min_cov} \
            {params.disabled_flag} > {log} 2>&1
        """


rule store_psl_metrics_augMP:
    """Store REAL PSL-derived AlnCoverage/AlnIdentity for augMP-mode transcripts."""
    input:
        psl=f"{config['work_dir']}/augustus/{{genome}}_augMP.psl",
        generated=f"{config['work_dir']}/databases/{{genome}}_augMP_psl_generated.done",
        augmp_gp=f"{config['work_dir']}/augustus/{{genome}}_augMP.gp",
        # ancient(): DB mtime is not a valid rerun trigger (shared mutable container);
        # freshness is tracked by the *_psl_metrics.done markers.
        db_path=ancient(f"{config['work_dir']}/databases/{{genome}}.db"),
        ref_gp=f"{config['work_dir']}/reference/{config['ref_genome']}.gp",
        prior_pairwise=f"{config['work_dir']}/databases/{{genome}}_transMap_pairwise_psl_metrics.done",
    output:
        done=f"{config['work_dir']}/databases/{{genome}}_augMP_psl_metrics.done"
    wildcard_constraints:
        genome = "|".join(TARGET_GENOMES)
    shell:
        """
        set -euo pipefail
        n_gp=$(grep -c '^augMP-' {input.augmp_gp} || true)
        n_psl=$(wc -l < {input.psl} | tr -d ' ')
        if [ "$n_gp" -gt 0 ] && [ "$n_psl" -eq 0 ]; then
            echo "ERROR: augMP PSL empty — rerun generate_augMP_psl (missing {input.generated})" >&2
            exit 1
        fi
        python3 cat/store_psl_metrics.py \
            --db-path {input.db_path} \
            --psl {input.psl} \
            --mode augMP \
            --ref-gp {input.ref_gp}
        echo "stored_augMP_psl_metrics n_psl=$n_psl" > {output.done}
        """

rule generate_consensus:
    """
    Submits consensus generation as a SLURM job to process each genome.
    Similar to how the run_transcript_map (txTM) rule is run, this creates and submits a SLURM job script.
    """
    input:
        unpack(get_consensus_inputs)
    output:
        gp=f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus.gp",
        gp_info=f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus.gp_info",
        gff3=f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus.gff3",
        fasta=f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus.fasta",
        protein_fasta=f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus_protein.fasta",
        metrics_json=f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus.json",
        done=touch(f"{config['work_dir']}/{{genome}}_consensus.done")
    wildcard_constraints:
        genome = ANNOTATION_GENOME_WC
    log:
        f"{config['work_dir']}/logs/consensus/{{genome}}.log"
    threads: 1 if IS_CLUSTER else get_local_res("generate_consensus", "threads")
    resources:
        mem_gb=1 if IS_CLUSTER else get_local_res("generate_consensus", "mem_gb"),
        time_h=1 if IS_CLUSTER else get_local_res("generate_consensus", "time_h"),
        job_id=lambda wildcards, attempt: f"consensus-submit-{wildcards.genome}-{attempt}"
    run:
        import os
        import subprocess
        import time

        # Get paths and variables
        work_dir = config['work_dir']
        genome = wildcards.genome
        job_script = f"{work_dir}/consensus_gene_set/{genome}_consensus_job.sh"

        cpus = get_res('generate_consensus', 'cpus')

        # Get active modes and build denovo_tx_modes argument
        active_modes = get_active_modes_for_wildcards(wildcards)
        denovo_modes = [m for m in active_modes if m in ['augPB', 'strg']]
        denovo_tx_modes_arg = f"--denovo-tx-modes {' '.join(denovo_modes)}" if denovo_modes else ""

        # Build BAM file arguments for real RNA-seq support computation
        genome_data = config.get("transcriptomic_data", {}).get(genome, {})
        rnaseq_bams = genome_data.get("bam", []) + genome_data.get("intronbam", [])
        isoseq_bams = genome_data.get("isoseq_bam", [])
        bam_args = f"--bam-files {' '.join(rnaseq_bams)}" if rnaseq_bams else ""
        isoseq_bam_args = f"--isoseq-bam-files {' '.join(isoseq_bams)}" if isoseq_bams else ""
        ref_gp_arg = f"--ref-gp {input.ref_gp}"

        # Build gp-list argument
        gp_list_str = ' '.join(input.gp_list)

        # Build optional flag arguments
        optional_flags = []
        if config.get("denovo_ignore_novel_genes", False):
            optional_flags.append("--denovo-ignore-novel-genes")
        if config.get("denovo_only_novel_genes", False):
            optional_flags.append("--denovo-only-novel-genes")
        if rcfg("denovo_allow_unsupported", False):
            optional_flags.append("--denovo-allow-unsupported")
        if config.get("denovo_allow_bad_annot_or_tm", False):
            optional_flags.append("--denovo-allow-bad-annot-or-tm")
        if config.get("denovo_allow_novel_ends", False):
            optional_flags.append("--denovo-allow-novel-ends")
        if config.get("require_pacbio_support", False):
            optional_flags.append("--require-pacbio-support")
        if len(RNASEQ_GENOMES) > 0:
            optional_flags.append("--hints-db-has-rnaseq")
        # rebuild_consensus in config is documentation only; regen via rm outputs + --forcerun
        if config.get("filter_overlapping_genes", False):
            optional_flags.append("--filter-overlapping-genes")
        if config.get("filter_spurious_pc_overlaps_not_in_reference", False):
            optional_flags.append("--filter-spurious-pc-overlaps-not-in-reference")
        if config.get("in_species_rna_support_only", False):
            optional_flags.append("--in-species-rna-support-only")
        if not config.get("consensus_postprocess", True):
            optional_flags.append("--no-consensus-postprocess")

        optional_flags_str = ' '.join(optional_flags) if optional_flags else ""

        script_content = build_sbatch_header(
            "generate_consensus",
            f"consensus-{genome}",
            f"{work_dir}/logs/consensus/{genome}_slurm.out",
            f"{work_dir}/logs/consensus/{genome}_slurm.err"
        ) + f"""
# Run consensus_runner.py with the specified parameters
python cat/consensus_runner.py \\
    --gp-list {gp_list_str} \\
    --db-path {input.db_path} \\
    --ref-db-path {input.ref_db_path} \\
    --fasta {input.fasta} \\
    --genome {genome} \\
    --consensus-gp {output.gp} \\
    --consensus-gp-info {output.gp_info} \\
    --consensus-gff3 {output.gff3} \\
    --consensus-fasta {output.fasta} \\
    --protein-fasta {output.protein_fasta} \\
    --metrics-json {output.metrics_json} \\
    --intron-rnaseq-support {config.get("intron_rnaseq_support", 0)} \\
    --exon-rnaseq-support {config.get("exon_rnaseq_support", 0)} \\
    --intron-annot-support {config.get("intron_annot_support", 0)} \\
    --exon-annot-support {config.get("exon_annot_support", 0)} \\
    --original-intron-support {config.get("original_intron_support", 0)} \\
    --cnv-score-similarity {rcfg("cnv_score_similarity", 0.80)} \\
    --fragment-max-coverage {rcfg("consensus_fragment_max_coverage", 30.0)} \\
    --fragment-max-identity {rcfg("consensus_fragment_max_identity", 30.0)} \\
    {"--keep-protein-only-novel " if rcfg("keep_protein_only_novel", False) else ""}\\
    {f"--protein-novel-min-coverage {rcfg('protein_novel_min_coverage', 0.0)} --protein-novel-min-identity {rcfg('protein_novel_min_identity', 0.0)} --protein-novel-min-exons {rcfg('protein_novel_min_exons', 2)} --protein-novel-min-cds-aa {rcfg('protein_novel_min_cds_aa', 100)} " if rcfg("keep_protein_only_novel", False) else ""}\\
    {"--protein-novel-keep-overlapping " if rcfg("protein_novel_keep_overlapping", False) else ""}\\
    {f"--rescue-expressed-noncoding-to-pc --rescue-expressed-min-cds-aa {rcfg('rescue_expressed_min_cds_aa', 100)} " if rcfg("rescue_expressed_noncoding_to_pc", True) else ""}\\
    {"--rescue-expressed-allow-single-exon " if rcfg("rescue_expressed_allow_single_exon", False) else ""}\\
    {"" if rcfg("rescue_noncoding_require_protein_evidence", True) else "--no-rescue-noncoding-require-protein-evidence "}\\
    {f"--rescue-dropped-augMP --rescue-augMP-min-exons {rcfg('rescue_augMP_min_exons', 2)} --rescue-augMP-min-cds-aa {rcfg('rescue_augMP_min_cds_aa', 100)} --rescue-augMP-single-exon-min-cds-aa {rcfg('rescue_augMP_single_exon_min_cds_aa', 300)} --rescue-augMP-min-coverage {rcfg('rescue_augMP_min_coverage', 0.0)} --rescue-augMP-min-identity {rcfg('rescue_augMP_min_identity', 0.0)} " if rcfg("rescue_dropped_augMP", True) else ""}\\
    --txTM-min-coverage {rcfg("txTM_min_coverage", 0.0)} \\
    {"--txTM-strict-metrics " if rcfg("txTM_strict_metrics", False) else ""}\\
    --txTM-min-coverage-no-transmap {rcfg("txTM_min_coverage_no_transmap", 80)} \\
    {"--txTM-strict-metrics-no-transmap " if rcfg("txTM_strict_metrics_no_transmap", True) else "--no-txTM-strict-metrics-no-transmap "}\\
    {f"--txTM-min-coverage-noncoding {config['txTM_min_coverage_noncoding']} " if config.get("txTM_min_coverage_noncoding") is not None else ""}\\
    {f"--txTM-min-coverage-no-transmap-noncoding {config['txTM_min_coverage_no_transmap_noncoding']} " if config.get("txTM_min_coverage_no_transmap_noncoding") is not None else ""}\\
    {"--txTM-strict-metrics-no-transmap-noncoding " if config.get("txTM_strict_metrics_no_transmap_noncoding") else ""}\\
    {"--no-txTM-strict-metrics-no-transmap-noncoding " if config.get("txTM_strict_metrics_no_transmap_noncoding") is False else ""}\\
    --augMP-min-coverage-no-anchor {rcfg("augMP_min_coverage_no_anchor", 80)} \\
    {"--augMP-strict-metrics-no-anchor " if rcfg("augMP_strict_metrics_no_anchor", True) else "--no-augMP-strict-metrics-no-anchor "}\\
    --txTM-transmap-anchor-overlap {rcfg("txTM_transmap_anchor_overlap", 0.25)} \\
    --min-pc-len-ratio-txTM-only-rescue {rcfg("min_pc_len_ratio_txTM_only_rescue", rcfg("min_pc_len_ratio_vs_reference", 0.4))} \\
    {"--no-rescue-reference-isoforms " if not config.get("rescue_reference_isoforms", True) else ""}\\
    {"--no-rescue-reference-noncoding-genes " if not config.get("rescue_reference_noncoding_genes", True) else ""}\\
    {"--no-rescue-alternative-isoforms " if not config.get("rescue_alternative_isoforms", True) else ""}\\
    --rescue-min-txTM-coverage {rcfg("rescue_min_txTM_coverage", 80)} \\
    --rescue-min-txTM-coverage-noncoding {rcfg("rescue_min_txTM_coverage_noncoding", 50)} \\
    --min-nc-len-ratio-txTM-only-rescue {rcfg("min_nc_len_ratio_txTM_only_rescue", 0.25)} \\
    --denovo-num-introns {config.get("denovo_num_introns", 1)} \\
    --denovo-splice-support {config.get("denovo_splice_support", 1)} \\
    --denovo-exon-support {config.get("denovo_exon_support", 1)} \\
    --denovo-novel-end-distance {config.get("denovo_novel_end_distance", 0)} \\
    --strg-min-single-exon-len {config.get("strg_min_single_exon_len", 500)} \\
    --disregard-long-mode-ratio {config.get("disregard_long_mode_ratio", 2.0)} \\
    --disregard-long-mode-min-bp {config.get("disregard_long_mode_min_bp", 50000)} \\
    --min-pc-len-ratio-vs-reference {rcfg("min_pc_len_ratio_vs_reference", 0.4)} \\
    --postprocess-min-introns-low-support {config.get("postprocess_min_introns_low_support", 3)} \\
    --postprocess-low-support-fraction {rcfg("postprocess_low_support_fraction", 0.3)} \\
    --postprocess-augpb-chimera-exon-ratio {config.get("postprocess_augpb_chimera_exon_ratio", 1.5)} \\
    --num-workers {cpus} \\
    {denovo_tx_modes_arg} \\
    {bam_args} \\
    {isoseq_bam_args} \\
    {ref_gp_arg} \\
    {optional_flags_str}
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting consensus job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.gp, output.metrics_json],
                log_file,
                "generate_consensus",
                max_wait_s=43200
            )

            # Post-submission stability check (cluster mode only): wait for output
            # files to stabilise in size before Snakemake considers the rule done.
            if IS_CLUSTER:
                check_interval = 30
                last_size = -1
                stable_count = 0
                required_stable_checks = 2

                for _ in range(20):  # up to ~10 minutes of extra checks
                    if os.path.exists(output.gp) and os.path.exists(output.metrics_json):
                        gp_size = os.path.getsize(output.gp)
                        json_size = os.path.getsize(output.metrics_json)
                        if gp_size > 0 and json_size > 0:
                            current_size = gp_size + json_size
                            if current_size == last_size:
                                stable_count += 1
                                if stable_count >= required_stable_checks:
                                    log_file.write(f"Consensus output for {genome} validated (stable sizes)\n")
                                    break
                            else:
                                stable_count = 0
                                last_size = current_size
                    time.sleep(check_interval)


rule annotate_novel_genes:
    """
    For each novel gene (transcript_class == 'putative_novel') in the consensus
    gene set, runs DIAMOND blastp against reference proteins and assigns a
    'paralog of GENE_NAME' description. Produces annotated GFF3 and gp_info
    files alongside the standard consensus outputs.
    """
    input:
        protein_fasta   = f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus_protein.fasta",
        gp_info         = f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus.gp_info",
        gff3            = f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus.gff3",
        gp              = f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus.gp",
        metrics_json    = f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus.json",
        ref_gtf         = WORK_DIR / f"reference/{REF_GENOME}.gtf",
        ref_fasta       = WORK_DIR / f"genome_files/{REF_GENOME}.fa",
        ref_gp_attrs    = WORK_DIR / f"reference/{REF_GENOME}.gp_attrs",
    output:
        gff3    = f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus_novel_annotated.gff3",
        gp_info = f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus_novel_annotated.gp_info",
        gp      = f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus_novel_annotated.gp",
        metrics_json = f"{config['work_dir']}/consensus_gene_set/{{genome}}_consensus_novel_annotated.json",
        done    = touch(f"{config['work_dir']}/{{genome}}_novel_annotation.done"),
    wildcard_constraints:
        # Must cover ANNOTATION_GENOMES (targets + ancestors): gene_family_report
        # consumes *_consensus_novel_annotated.* for every annotated genome.
        genome = ANNOTATION_GENOME_WC
    log:
        f"{config['work_dir']}/logs/consensus/{{genome}}_novel_annotation.log"
    threads: 1 if IS_CLUSTER else get_local_res("annotate_novel_genes", "threads")
    resources:
        mem_gb  = 1 if IS_CLUSTER else get_local_res("annotate_novel_genes", "mem_gb"),
        time_h  = 1 if IS_CLUSTER else get_local_res("annotate_novel_genes", "time_h"),
        job_id  = lambda wildcards, attempt: f"novel-annot-submit-{wildcards.genome}-{attempt}"
    run:
        import os
        import subprocess
        import time

        work_dir = config['work_dir']
        genome   = wildcards.genome
        job_script = f"{work_dir}/consensus_gene_set/{genome}_novel_annotation_job.sh"
        cpus = get_res('annotate_novel_genes', 'cpus')

        novel_evalue     = config.get('novel_annotation_evalue', 1e-3)
        novel_min_id     = config.get('novel_annotation_min_identity', 20.0)
        novel_min_qcov   = config.get('novel_annotation_min_query_cover', 20.0)
        novel_paralog_cap = config.get('novel_paralog_max_copies', 20)

        script_content = build_sbatch_header(
            "annotate_novel_genes",
            f"novel-annot-{genome}",
            f"{work_dir}/logs/consensus/{genome}_novel_annotation_slurm.out",
            f"{work_dir}/logs/consensus/{genome}_novel_annotation_slurm.err"
        ) + f"""
echo "Annotating novel genes for: {genome}"
echo "Job ID: ${{SLURM_JOB_ID}}"
echo "Node: $(hostname)"
echo "Start time: $(date)"

python cat/annotate_novel_genes.py \\
    --consensus-protein-fasta {input.protein_fasta} \\
    --consensus-gp-info       {input.gp_info} \\
    --consensus-gff3          {input.gff3} \\
    --consensus-gp            {input.gp} \\
    --ref-gtf                 {input.ref_gtf} \\
    --ref-fasta               {input.ref_fasta} \\
    --ref-gp-attrs            {input.ref_gp_attrs} \\
    --output-gff3             {output.gff3} \\
    --output-gp-info          {output.gp_info} \\
    --output-gp               {output.gp} \\
    --consensus-metrics-json  {input.metrics_json} \\
    --output-metrics-json     {output.metrics_json} \\
    --novel-paralog-max-copies {novel_paralog_cap} \\
    --threads                 {cpus} \\
    --evalue                  {novel_evalue} \\
    --min-identity            {novel_min_id} \\
    --min-query-cover         {novel_min_qcov} \\
    2>&1 | tee -a {log[0]}

echo "Novel gene annotation complete for {genome}"
echo "End time: $(date)"
"""

        with open(log[0], 'a') as log_file:
            log_file.write(f"Submitting novel gene annotation job for {genome}...\n")
            run_or_submit(
                script_content,
                job_script,
                [output.gff3, output.gp_info, output.gp, output.metrics_json],
                log_file,
                "annotate_novel_genes",
                max_wait_s=7200
            )


rule generate_plots:
    """Generate summary PDF plots from pipeline metrics."""
    input:
        annotation_done = expand(str(WORK_DIR / "{genome}_novel_annotation.done"), genome=TARGET_GENOMES),
        consensus_jsons = expand(str(WORK_DIR / "consensus_gene_set/{genome}_consensus_novel_annotated.json"), genome=TARGET_GENOMES),
        tm_jsons = expand(str(WORK_DIR / "transMap/{genome}_filter_tm_metrics.json"), genome=TARGET_GENOMES),
        dbs = expand(str(WORK_DIR / "databases/{genome}.db"), genome=TARGET_GENOMES),
        annotation_db = str(WORK_DIR / f"databases/{config['ref_genome']}.db")
    output:
        done = touch(str(WORK_DIR / "plots.done"))
    log:
        str(WORK_DIR / "logs/consensus/generate_plots.log")
    resources:
        mem_gb=16,
        time_h=1,
        job_id="generate-plots"
    run:
        import sys
        sys.argv = ['plots']

        tm_json_args = []
        for genome in TARGET_GENOMES:
            tm_json_args.extend(['--tm-jsons', genome,
                                 str(WORK_DIR / f"transMap/{genome}_filter_tm_metrics.json")])
        metrics_json_args = []
        for genome in TARGET_GENOMES:
            metrics_json_args.extend(['--metrics-jsons', genome,
                                      str(WORK_DIR / f"consensus_gene_set/{genome}_consensus_novel_annotated.json")])
        db_args = []
        for genome in TARGET_GENOMES:
            db_args.extend(['--dbs', genome,
                            str(WORK_DIR / f"databases/{genome}.db")])

        pb_genome_args = []
        if ISOSEQ_GENOMES:
            pb_genome_args = ['--pb-genomes'] + list(ISOSEQ_GENOMES)

        sys.argv = ['plots'] + tm_json_args + metrics_json_args + db_args + [
            '--annotation-db', str(WORK_DIR / f"databases/{config['ref_genome']}.db"),
            '--ordered-genomes'] + list(TARGET_GENOMES) + pb_genome_args + [
            '--out-dir', str(WORK_DIR / "plots")]

        import cat.plots
        cat.plots.main()


rule gene_family_report:
    """
    Protein-coding ortholog / paralog / gene-family expansion-contraction report.

    Reuses CAT2's projection-based orthology (every consensus transcript carries its
    reference `source_gene`) to build copy-number matrices and expansion/contraction
    calls straight from the finished consensus gene sets. Cheap, read-only, opt-out
    via `gene_family_report: false` in the config.
    """
    input:
        annotation_done = expand(str(WORK_DIR / "{genome}_novel_annotation.done"), genome=ANNOTATION_GENOMES),
        gp_infos = expand(str(WORK_DIR / "consensus_gene_set/{genome}_consensus_novel_annotated.gp_info"), genome=ANNOTATION_GENOMES),
        gps = expand(str(WORK_DIR / "consensus_gene_set/{genome}_consensus_novel_annotated.gp"), genome=ANNOTATION_GENOMES),
        annotation_db = str(WORK_DIR / f"databases/{config['ref_genome']}.db")
    output:
        done = touch(str(WORK_DIR / "gene_family_report.done"))
    log:
        str(WORK_DIR / "logs/consensus/gene_family_report.log")
    resources:
        mem_gb=16,
        time_h=1,
        job_id="gene-family-report"
    params:
        script = str(CAT2_ROOT / "scripts/gene_family_analysis.py"),
        genomes = " ".join(ANNOTATION_GENOMES),
        hal = config["hal"],
        ref_genome = config["ref_genome"],
    shell:
        r"""
        MPLCONFIGDIR={WORK_DIR}/.mplcache \
        python {params.script} \
            --work-dir {WORK_DIR} \
            --ref-genome {params.ref_genome} \
            --genomes {params.genomes} \
            --hal {params.hal} > {log} 2>&1
        """


rule finish_pipeline:
    input:
        expand(str(WORK_DIR / "{genome}_consensus.done"), genome=ANNOTATION_GENOMES),
        expand(str(WORK_DIR / "{genome}_novel_annotation.done"), genome=ANNOTATION_GENOMES),
        str(WORK_DIR / "plots.done"),
        *( [str(WORK_DIR / "gene_family_report.done")]
           if config.get("gene_family_report", True) else [] )
    output:
        str(WORK_DIR / "pipeline.complete.done")
    resources:
        mem_gb=1,
        time_h=1,
        job_id="finish-pipeline"
    run:
        from snakemake.logging import logger
        logger.info("CAT comparative annotation pipeline has completed successfully!")
        logger.info("To view a detailed summary, run: snakemake --report report.html")
        
        with open(str(output[0]), 'w') as f:
            f.write("CAT pipeline completed successfully.\n")


rule cleanup_logs_and_done_files:
    input:
        str(WORK_DIR / "pipeline.complete.done")
    output:
        str(WORK_DIR / "cleanup.complete")
    resources:
        mem_gb=1,
        time_h=1,
        job_id="cleanup-logs-done"
    shell:
        """
        # Remove all .done files (including consensus.done files)
        find {WORK_DIR} -name "*.done" -type f -delete
        
        # Remove toil log files
        find {WORK_DIR} -name "toil_*.err.log" -type f -delete
        find {WORK_DIR} -name "toil_*.out.log" -type f -delete
        
        echo "Cleanup completed: All .done files and toil log files have been removed" > {output}
        """
