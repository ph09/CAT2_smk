#!/bin/bash
# CAT2 runtime environment.
# - conda "cat" owns Python and Snakemake (must stay first on PATH)
# - cactus HAL binaries/libs are added without sourcing the cactus venv activate
#   script, which would prepend venv Python and break Snakemake

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cat

export PYTHONNOUSERSITE=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Override if cactus is installed elsewhere: export CACTUS_BIN=/path/to/cactus-bin-v3.2.1
CACTUS_BIN="${CACTUS_BIN:-/private/groups/cgl/cactus/cactus-bin-v3.2.1}"
CACTUS_VENV="${CACTUS_BIN}/venv-cactus-v3.2.1"

export LD_LIBRARY_PATH="${CACTUS_BIN}/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${CACTUS_BIN}/lib:${PYTHONPATH:-}"

# Order matters: conda → repo standalones → cactus HAL tools → cactus venv CLIs
export PATH="${CONDA_PREFIX}/bin:${SCRIPT_DIR}/standalones:${CACTUS_BIN}/bin:${CACTUS_VENV}/bin:${PATH}"
