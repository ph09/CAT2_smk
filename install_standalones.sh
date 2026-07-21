#!/usr/bin/env bash
#
# Fetch the third-party command-line binaries that CAT2 expects on PATH.
#
# These tools are NOT committed to the repository (they are large,
# platform-specific binaries). This script downloads the UCSC "Kent" utilities
# for linux.x86_64 into ./standalones/ and marks them executable. `setup_env.sh`
# then puts ./standalones/ on PATH.
#
# Usage:
#   chmod +x install_standalones.sh
#   ./install_standalones.sh            # linux.x86_64 (default)
#   UCSC_ARCH=macOSX.x86_64 ./install_standalones.sh   # override arch
#
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="$REPO_DIR/standalones"
ARCH="${UCSC_ARCH:-linux.x86_64}"
BASE="http://hgdownload.soe.ucsc.edu/admin/exe/${ARCH}"

# UCSC Kent command-line utilities required by the pipeline.
UCSC_TOOLS=(
  axtChain
  bamToPsl
  bedSort
  bedToBigBed
  chainMergeSort
  chainNet
  chainSort
  chainToBigChain
  clusterGenes
  faSize
  faToTwoBit
  genePredToBed
  genePredToBigGenePred
  genePredToFakePsl
  genePredToGtf
  gff3ToGenePred
  gtfToGenePred
  hgLoadChain
  liftOver
  netChainSubset
  netFilter
  netSyntenic
  pslCDnaFilter
  pslCheck
  pslMap
  pslMapPostChain
  pslPosTarget
  pslRecalcMatch
  pslToBigPsl
  transMapPslToGenePred
  wigToBigWig
)

DL="curl -fSL"
if ! command -v curl >/dev/null 2>&1; then
  if command -v wget >/dev/null 2>&1; then
    DL="wget -O -"
  else
    echo "ERROR: need curl or wget on PATH." >&2
    exit 1
  fi
fi

echo "Downloading ${#UCSC_TOOLS[@]} UCSC utilities (${ARCH}) into $DEST ..."
fail=0
for tool in "${UCSC_TOOLS[@]}"; do
  out="$DEST/$tool"
  if [[ -x "$out" ]]; then
    echo "  [skip] $tool (already present)"
    continue
  fi
  printf '  [get ] %s ... ' "$tool"
  if $DL "$BASE/$tool" > "$out" 2>/dev/null && [[ -s "$out" ]]; then
    chmod +x "$out"
    echo "ok"
  else
    rm -f "$out"
    echo "FAILED"
    fail=1
  fi
done

echo
echo "Done. UCSC tools are in $DEST"
if [[ "$fail" -ne 0 ]]; then
  echo "WARNING: one or more UCSC downloads failed; re-run or fetch manually from $BASE" >&2
fi

cat <<'EOF'
EOF
