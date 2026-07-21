#!/bin/bash
set -euo pipefail

# Install the Python package
"${PYTHON}" -m pip install . --no-deps --no-build-isolation -v

# Install data files to $PREFIX/share/cat/
SHARE_DIR="${PREFIX}/share/cat"
mkdir -p "${SHARE_DIR}"

# Minisplice model files (used by miniprot/augustus_pb rules)
for f in standalones/vi2-7k.kan standalones/vi2-7k.kan.cali; do
    [ -f "$f" ] && install -m 644 "$f" "${SHARE_DIR}/$(basename $f)"
done

# Augustus extrinsic config files
for cfg in augustus_cfgs/*.cfg; do
    [ -f "$cfg" ] && install -m 644 "$cfg" "${SHARE_DIR}/$(basename $cfg)"
done

# Install custom standalone tools that do not have individual bioconda packages.
# These are pre-compiled Linux x86_64 binaries from the UCSC Kent source tree.
# Tools that DO have bioconda packages (ucsc-axtchain, ucsc-bamtopsl, etc.)
# are declared as conda dependencies in meta.yaml and are NOT installed here.
CUSTOM_TOOLS=(
    bam-to-bigchain
    chainToBigChain
    clusterGenes
    genePredToBigGenePred
    genePredToFakePsl
    hgLoadChain
    pslMapPostChain
    pslPosTarget
    pslToBigPsl
    transMapPslToGenePred
)

for tool in "${CUSTOM_TOOLS[@]}"; do
    src="standalones/${tool}"
    if [ -f "${src}" ]; then
        install -m 755 "${src}" "${PREFIX}/bin/${tool}"
    fi
done
