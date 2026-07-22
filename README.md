# CAT2

CAT2 (Comparative Annotation Toolkit, v2) is an improved version of CAT, a pipeline that annotates
one or more target genomes from a [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus)
HAL alignment plus a reference annotation. It combines several complementary
gene-finding strategies and reconciles them into a single consensus gene set per
genome:

- **transMap** — projects the reference annotation across the HAL alignment
  (recovers orthologs of reference genes).
- **AugustusTM/TMR** — refines transMap projections with Augustus using RNA-seq
  and annotation hints.
- **AugustusPB** — Augustus with PacBio/IsoSeq hints for isoform discovery.
- **augMP (miniprot)** — aligns a **protein database** to each genome; See
  [Protein set](#protein-set-for-augmp--miniprot).
- **StringTie** — reference-guided transcript assembly from RNA-seq or Iso-seq.

> This is the current version of CAT2 and is under active development.

---

## Requirements

- Linux x86-64, `conda`/`mamba` (miniforge recommended)
- A [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) install (for
  the HAL tools: `halStats`, `hal2fasta`, `halLiftover`, …)
- A cluster (SLURM or SGE) is recommended for real genomes; small inputs run
  locally.

## Installation

```bash
git clone https://github.com/ph09/CAT2_smk.git
cd CAT2_smk
conda env create -f environment.yaml
conda activate cat2
pip install --no-deps -e .        # installs the cat2 package itself
chmod +x install_standalones.sh
./install_standalones.sh          # fetches UCSC Kent binaries into ./standalones/
```

`install_standalones.sh` downloads the UCSC utilities the pipeline calls. Two
non-UCSC helpers (`aln2hints.pl` from AUGUSTUS, `pal2nal.pl`) must be placed in
`./standalones/` by hand — the script prints a reminder.

## Setup

`setup_env.sh` activates the conda env and puts Cactus's HAL tools and the repo
`standalones/` on `PATH`:

```bash
export CACTUS_BIN=/your/path/to/cactus-bin-v3.2.1
source setup_env.sh
```

## Quick start (bundled test data)

The repo ships small parser/test fixtures, but the test HAL and BAMs are hosted
separately. Download them into `test_data/` first:

- HAL: <https://public.gi.ucsc.edu/~pnhebbar/share/testData/vertebrates.hal>
- BAMs: <https://public.gi.ucsc.edu/~pnhebbar/share/testData/bams/>

Then dry-run, then run:

```bash
# dry-run: parse Snakefile, validate config, build the DAG
snakemake --configfile input.yaml -n all --cores 4

# run locally
snakemake --configfile input.yaml all --cores 4
```

`input.yaml` is the template config you will copy and modify for your own data.

---

## Configuring `input.yaml`

`input.yaml` is heavily commented; this section highlights the parts you will
most often change. Every key not shown here has a sensible built-in default.

### Core inputs (required)

| Key | Meaning |
| --- | --- |
| `work_dir` | Output directory for this run (created if missing). |
| `hal` | Path to the Cactus HAL alignment. |
| `annotation` | Reference annotation (GFF3) for `ref_genome`. |
| `ref_genome` | Name of the reference genome **as it appears in the HAL**. |
| `genomes` | List of target genome names (HAL names) to annotate. |

```yaml
work_dir: "my_run"
hal: "/data/primates.hal"
annotation: "/data/GRCh38.gff3"
ref_genome: "hg38"
genomes: ["hg38", "rheMac10", "calJac4"]
```

### Protein set (for augMP / miniprot)

augMP aligns a protein database to every target genome. With a **single-species**
(e.g. human-only) protein set it can only re-find genes that already exist in the
reference. To discover lineage-specific / non-reference genes you must give it
proteins from those lineages. There are two ways to provide proteins.

**(A) A protein FASTA you already have:**

```yaml
protein_fasta: "/data/my_proteins.fa"
```

**(B) Auto-built multi-species DB** — list species, NCBI taxon IDs, and/or clades,
and CAT2 downloads UniProt reference proteomes, merges and de-duplicates them, and
uses the result as the miniprot input. This runs **once on the controller machine**
(needs outbound internet) and is cached. If `protein_db` is set it takes
precedence over `protein_fasta`.

```yaml
protein_db:
  species:                     # names OR HAL genome names (auto-normalised)
    - "Homo sapiens"
    - "Macaca mulatta"
    - "Callithrix jacchus"
    - "Mus musculus"           # an outgroup broadens coverage
  # taxa: [9606, 9544]         # alternative/supplement: NCBI taxon IDs

  # clades: pull proteomes for members of a genus/family (or any higher rank).
  # Use this to cover poorly annotated species via their better-annotated
  # relatives. CAT2 finds clade members that have a proteome and keeps the best
  # (highest quality first, up to max_per_clade).
  clades:
    - "Cercopithecidae"        # Old World monkeys (family)
    - "Callitrichidae"         # marmosets/tamarins (family)
  max_per_clade: 25            # cap proteomes per clade (best quality first)
  clade_include_other: false   # also accept non-reference proteomes when expanding clades

  base_fasta: ""               # optionally fold in your own proteins too
  dedup: true                  # collapse identical sequences (keeps miniprot fast)
  min_len: 20                  # drop proteins shorter than this
  out: "my_run/protein_db/protein_db.fa"   # built DB path (this is the default location)
```

Tip: pick a handful of well-spread representatives rather than every species —
most genes are shared, so dumping hundreds of proteomes is slower with little
recall gain. Genus-/family-level `clades` are the best way to get coverage for
species that lack their own good proteome.

**Building the DB by hand** (same machinery, standalone) and pointing
`protein_fasta` at the result:

```bash
python scripts/build_protein_db.py \
    --out my_run/protein_db.fa \
    --clades "Cercopithecidae,Callitrichidae" \
    --species "Homo sapiens,Mus musculus" \
    --base-fasta /data/GRCh38.prot.fa
```

Useful `build_protein_db.py` flags: `--species`, `--taxa`, `--clades`,
`--max-per-clade` (default 25), `--clade-include-other`, `--base-fasta`,
`--min-len` (default 20), `--no-dedup`, `--cache-dir`, `--strict` (fail on any
species that will not resolve), `--summary` (per-species TSV). Species names are
normalised leniently, so HAL genome names paste in directly
(`PR00246~Eulemur_fulvus.pri` → `Eulemur fulvus`).

**miniprot sensitivity** is tunable under a `miniprot:` block (all keys optional;
defaults are already tuned more permissively than miniprot's own to recover more
divergent/paralogous copies). Raising sensitivity finds more candidates at the
cost of runtime and false positives (which consensus filtering then prunes). See
the `miniprot:` block in `input.yaml` for every knob (`splice_model`,
`max_intron`, `min_secondary_ratio`, `max_secondary`, `out_n/out_s/out_c`, …).

### Transcriptomic data (optional but recommended)

Per-genome RNA-seq / long-read BAMs, plus which genomes use them:

```yaml
transcriptomic_data:
  hg38:
    bam:                        # short-read RNA-seq
      - "/data/hg38/rnaseq1.bam"
    isoseq_bam:                 # long reads (PacBio/IsoSeq/ONT)
      - "/data/hg38/isoseq.bam"
    intronbam:                  # noisy short reads, used for intron hints only
      - "/data/hg38/noisy.bam"

rnaseq_genomes: ["hg38"]        # genomes with usable short-read RNA-seq
isoseq_genomes: ["hg38"]        # genomes with usable long reads
```

- `bam` — short-read RNA-seq (splice + coverage hints)
- `isoseq_bam` — long reads (isoform evidence, drives AugustusPB / StringTie)
- `intronbam` — noisy short reads used only for intron hints

### Mode toggles

`liftoff` is **off by default**. Turn it on only when annotating another genome
or haplotype of the **same species** (near-identical assemblies). For
cross-species projection, leave it off and use transMap / txTM / augMP.

```yaml
augustus: true            # AugustusTM/TMR refinement
augustus_pb: true         # AugustusPB (needs isoseq)
stringtie: true           # StringTie assembly
stringtie_genomes: ["hg38"]
txTM: true                # transcript-level minimap2 map (best on close relatives)
transmap_pairwise: true   # BAM/minimap2 pairwise transMap (+ augTM/TMR_pairwise).
                          # Set false for highly diverged genomes where pairwise
                          # chaining fails — skips minimap2_bam → bam_to_chain →
                          # the whole pairwise path.
liftoff: false            # External Liftoff; same-species haplotypes only (see above)
# liftoff_sc: 0.85        # optional Liftoff -sc sequence-identity cutoff
augustus_species: "human" # Augustus species parameter set
```

### Recall vs precision

`high_recall` is an opt-in master switch. When `true`, the recall-limiting gates
across **every** mode (transMap paralog/overlap filtering, coverage floors,
consensus length/CNV/fragment cutoffs, denovo support, postprocess low-support
drop) are loosened together so fewer genes are missed, at the cost of more false
positives. It overrides the individual knobs below.

```yaml
high_recall: false        # leave false and tune individual keys for fine control
```

Fine-grained knobs (all documented inline in `input.yaml`) include transMap
filtering (`global_near_best`, `tm_filter_overlapping`, `tm_min_cover`, …),
consensus fragment reclassification (`consensus_fragment_max_coverage/identity`),
and postprocess drop tuning (`postprocess_*`).

### Ancestral genomes

Annotate reconstructed ancestral (internal HAL tree) genomes using
alignment-only modes:

```yaml
annotate_ancestors: true
ancestor_genomes: ["Anc0", "Anc1"]         # optional; auto-detected if omitted
ancestor_modes: ["transMap", "transMap_pairwise", "txTM"]
```

### Execution mode and cluster resources

```yaml
execution_mode: "auto"    # auto | slurm | sge | local
```

- `auto` — detect the scheduler (`sbatch` → slurm; `qsub` + `$SGE_ROOT` → sge;
  otherwise local).
- `slurm` / `sge` — submit jobs to the cluster; per-rule `mem`/`cpus`/`time` come
  from the `slurm:` → `rules:` block (shared by both schedulers). SGE-specific
  knobs (queue, parallel environment, memory flag) live under `cluster.sge`.
- `local` — run everything on the current machine; per-rule hints come from the
  `local:` block.

The `slurm:` block also has site-specific fields to review before running:
`partition` and `exclude_nodes` (the shipped value is UCSC-specific — change or
clear it).

---

## Running the pipeline

```bash
source setup_env.sh                                    # every new shell

# recommended: use the launcher (enforces snakemake >= 9, resumable)
./run_pipeline.sh --work-dir my_run --configfile input.yaml --cores 32

# long runs: keep alive under tmux/nohup
tmux new -s cat2 './run_pipeline.sh --work-dir my_run 2>&1 | tee my_run/run.log'

# or drive snakemake directly
snakemake --configfile input.yaml all --cores 32 --keep-going --rerun-incomplete
```

> snakemake **must be ≥ 9**. Snakemake 8.x corrupts multi-line f-strings inside
> rule `run:` blocks under Python 3.12; the launcher refuses to start otherwise.

Final per-genome consensus annotations are written under `work_dir` (GFF3 / GenePred).

## Development / CI

A fast, cluster-free smoke test (used by CI) checks the snakemake version, a
`run:`-block f-string regression guard, and that the `cat2` package imports:

```bash
./scripts/smoke_test.sh
```

It skips the full DAG dry-run automatically when `halStats` / the test HAL are
not present (e.g. in CI).
