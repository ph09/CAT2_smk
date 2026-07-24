#!/bin/bash
# CAT2 runtime environment helper.
#
# Prerequisites:
#   1. A conda env with CAT2 dependencies (default name: cat2), OR an already
#      active conda env that has them.
#   2. CACTUS_BIN pointing at a Cactus binary install whose bin/ contains
#      HAL tools (halStats, hal2fasta, ...). Example:
#        export CACTUS_BIN=/path/to/cactus-bin-v3.2.1
#      On sites that provide Cactus via environment modules, load the module
#      first (or set CACTUS_BIN to that module's prefix), then source this file
#      so conda Python stays ahead of any Cactus-bundled Python.
#
# Usage:
#   export CACTUS_BIN=/path/to/cactus-bin-v3.2.1
#   source setup_env.sh

# Activate conda only if nothing is active yet. Prefer CAT2_CONDA_ENV (default cat2).
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
        # shellcheck source=/dev/null
        source "${HOME}/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "${HOME}/mambaforge/etc/profile.d/conda.sh" ]]; then
        # shellcheck source=/dev/null
        source "${HOME}/mambaforge/etc/profile.d/conda.sh"
    elif ! command -v conda >/dev/null 2>&1; then
        echo "ERROR: conda not found; activate your CAT2 env before sourcing setup_env.sh" >&2
        return 1 2>/dev/null || exit 1
    fi
    conda activate "${CAT2_CONDA_ENV:-cat2}"
fi

export PYTHONNOUSERSITE=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "${CACTUS_BIN:-}" ]]; then
    echo "ERROR: CACTUS_BIN is not set." >&2
    echo "  export CACTUS_BIN=/path/to/cactus-bin-vX.Y.Z   # directory with bin/halStats" >&2
    return 1 2>/dev/null || exit 1
fi
if [[ ! -d "${CACTUS_BIN}/bin" ]]; then
    echo "ERROR: CACTUS_BIN=${CACTUS_BIN} has no bin/ directory" >&2
    return 1 2>/dev/null || exit 1
fi

# Optional: some cactus binary releases ship a sibling venv with extra CLIs.
CACTUS_VENV="${CACTUS_VENV:-}"
if [[ -z "${CACTUS_VENV}" ]]; then
    for candidate in "${CACTUS_BIN}/venv-cactus-v3.2.1" "${CACTUS_BIN}/venv"; do
        if [[ -d "${candidate}/bin" ]]; then
            CACTUS_VENV="${candidate}"
            break
        fi
    done
fi

if [[ -d "${CACTUS_BIN}/lib" ]]; then
    export LD_LIBRARY_PATH="${CACTUS_BIN}/lib:${LD_LIBRARY_PATH:-}"
fi

# Order matters: conda → repo standalones → cactus HAL tools → optional cactus venv CLIs.
# Conda must stay first so Cactus-bundled Python never shadows the CAT2 env.
PATH_ADDITIONS="${CONDA_PREFIX}/bin:${SCRIPT_DIR}/standalones:${CACTUS_BIN}/bin"
if [[ -n "${CACTUS_VENV}" && -d "${CACTUS_VENV}/bin" ]]; then
    PATH_ADDITIONS="${PATH_ADDITIONS}:${CACTUS_VENV}/bin"
fi
export PATH="${PATH_ADDITIONS}:${PATH}"
