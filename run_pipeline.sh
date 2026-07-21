#!/usr/bin/env bash
#
# CAT2 production launcher.
#
# Runs the whole pipeline (or resumes it) on a SLURM/SGE submit node under the
# conda `cat` environment. The `cat` env is the single, fully-provisioned
# interpreter: it owns Python (>=3.11) and snakemake (>=9) plus every python
# dependency, while `setup_env.sh` also puts the cactus HAL tools and the repo
# `standalones/` on PATH. Nothing else needs to be activated.
#
# IMPORTANT: snakemake must be >=9. snakemake 8.x corrupts rule `run:`-block
# f-strings under Python 3.12; `run_pipeline.sh` refuses to start otherwise.
#
# Usage:
#   ./run_pipeline.sh [--work-dir DIR] [--cores N] [--configfile FILE]
#                     [--dry-run] [-- <extra snakemake args>]
#
# Long runs should be launched under tmux/nohup so they survive disconnects:
#   tmux new -s cat2 './run_pipeline.sh --work-dir my_run 2>&1 | tee my_run/run.log'
#   # or:
#   nohup ./run_pipeline.sh --work-dir my_run > my_run/run.log 2>&1 &
#
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

CONFIGFILE="input.yaml"
WORK_DIR=""
CORES=32
DRY_RUN=0
RESTART_TIMES="${RESTART_TIMES:-2}"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --work-dir)   WORK_DIR="$2"; shift 2 ;;
        --cores)      CORES="$2"; shift 2 ;;
        --configfile) CONFIGFILE="$2"; shift 2 ;;
        --dry-run|-n) DRY_RUN=1; shift ;;
        --)           shift; EXTRA_ARGS+=("$@"); break ;;
        *)            EXTRA_ARGS+=("$1"); shift ;;
    esac
done

# ── Environment ────────────────────────────────────────────────────────────
# Clear any stray PYTHONPATH from the caller's shell so it cannot shadow the
# cat env, then let setup_env.sh set the intended (cactus lib) PYTHONPATH.
unset PYTHONPATH
# shellcheck disable=SC1091
source "$REPO_DIR/setup_env.sh"

# ── Preflight checks ─────────────────────────────────────────────────────────
command -v snakemake >/dev/null 2>&1 || { echo "ERROR: snakemake not found on PATH (is the 'cat' env set up?)" >&2; exit 1; }
SMK_VER="$(snakemake --version)"
SMK_MAJOR="${SMK_VER%%.*}"
if [[ "$SMK_MAJOR" -lt 9 ]]; then
    echo "ERROR: snakemake ${SMK_VER} detected; this pipeline requires snakemake >=9 on Python 3.12" >&2
    echo "       (snakemake 8.x corrupts run:-block f-strings). Install with:" >&2
    echo "         conda env update -n cat -f environment.yaml   # or: pip install 'snakemake>=9'" >&2
    exit 1
fi
if [[ -n "$WORK_DIR" ]]; then
    mkdir -p "$WORK_DIR"
fi

CONFIG_ARGS=(--configfile "$CONFIGFILE")
if [[ -n "$WORK_DIR" ]]; then
    CONFIG_ARGS+=(--config "work_dir=$WORK_DIR")
fi

echo "CAT2 launcher: snakemake ${SMK_VER}, python $(python --version 2>&1), host $(hostname)"
echo "  configfile=${CONFIGFILE} work_dir=${WORK_DIR:-<from config>} cores=${CORES} dry_run=${DRY_RUN}"

# Clear a stale lock left by a previously killed controller (safe: we hold the
# only run here). Ignored if nothing is locked.
snakemake "${CONFIG_ARGS[@]}" --unlock >/dev/null 2>&1 || true

if [[ "$DRY_RUN" -eq 1 ]]; then
    exec snakemake "${CONFIG_ARGS[@]}" -n --rerun-triggers mtime "${EXTRA_ARGS[@]}"
fi

exec snakemake \
    "${CONFIG_ARGS[@]}" \
    --cores "$CORES" \
    --rerun-triggers mtime \
    --restart-times "$RESTART_TIMES" \
    --keep-going \
    --printshellcmds \
    "${EXTRA_ARGS[@]}"
