# CAT2

This is the current version of CAT2 and is under active development.

## Installation

```bash
cd /path/to/CAT2_smk
conda env create -f environment.yaml
conda activate cat2
pip install --no-deps -e .
chmod +x install_standalones.sh
./install_standalones.sh
```

## Setup

```bash
export CACTUS_BIN=/your/path/to/cactus-bin-v3.2.1
source setup_env.sh
```

## Usage

Run the pipeline with:

```bash
snakemake --configfile input.yaml -n all --cores 4
```
The hal associated with the test_data is here: https://public.gi.ucsc.edu/~pnhebbar/share/testData/vertebrates.hal
The bams associated with the test_data are here: https://public.gi.ucsc.edu/~pnhebbar/share/testData/bams/
The input.yaml is the config file that you will need to modify.  bam means short read bams, isoseq_bam means long read bams, and intron_bam means noisy short read bams.
