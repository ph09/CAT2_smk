#!/usr/bin/env bash
#
# CAT2 run status reporter.
#
# Summarises the progress of a pipeline run by inspecting the sentinel/output
# files in a work directory. Safe to run at any time (read-only), including
# while the pipeline is still running.
#
# Usage:
#   scripts/pipeline_status.sh [WORK_DIR] [CONFIGFILE]
#
# If WORK_DIR is omitted it is read from CONFIGFILE (default input.yaml).
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

CONFIGFILE="${2:-input.yaml}"
WORK_DIR="${1:-}"
if [[ -z "$WORK_DIR" ]]; then
    WORK_DIR="$(awk '/^work_dir:/{gsub(/["'\'' ]/,"",$2); print $2; exit}' "$CONFIGFILE" 2>/dev/null || true)"
fi
if [[ -z "$WORK_DIR" || ! -d "$WORK_DIR" ]]; then
    echo "ERROR: work_dir '${WORK_DIR:-<unset>}' not found. Pass it as the first argument." >&2
    exit 1
fi

# Genomes from config (best-effort: parse the genomes: [...] line).
mapfile -t GENOMES < <(awk -F'[][]' '/^genomes:/{print $2}' "$CONFIGFILE" 2>/dev/null \
    | tr ',' '\n' | sed 's/["'\'' ]//g' | grep -v '^$' || true)
REF_GENOME="$(awk '/^ref_genome:/{gsub(/["'\'' ]/,"",$2); print $2; exit}' "$CONFIGFILE" 2>/dev/null || true)"

mark() { [[ -e "$1" ]] && echo "  [x] $2" || echo "  [ ] $2"; }

echo "==================================================================="
echo "CAT2 status: work_dir=$WORK_DIR"
echo "==================================================================="

echo "Overall:"
mark "$WORK_DIR/pipeline.complete.done" "pipeline.complete.done (whole pipeline finished)"
mark "$WORK_DIR/plots.done"             "plots.done (QC plots)"
mark "$WORK_DIR/gene_family_report.done" "gene_family_report.done"

echo
echo "Per-genome consensus / novel annotation:"
printf "  %-14s %-10s %-14s %-8s\n" "GENOME" "CONSENSUS" "NOVEL_ANNOT" "GENES"
for g in "${GENOMES[@]}"; do
    if [[ "$g" == "$REF_GENOME" ]]; then
        printf "  %-14s %-10s %-14s %-8s\n" "$g" "(reference)" "-" "-"
        continue
    fi
    cons="pending"; novel="-"; genes="-"
    [[ -e "$WORK_DIR/${g}_consensus.done" ]] && cons="done"
    [[ -e "$WORK_DIR/${g}_novel_annotation.done" ]] && novel="done"
    gp="$WORK_DIR/consensus_gene_set/${g}_consensus.gp"
    [[ -s "$gp" ]] && genes="$(wc -l < "$gp" | tr -d ' ')"
    printf "  %-14s %-10s %-14s %-8s\n" "$g" "$cons" "$novel" "$genes"
done

# Ancestors (if annotated).
ANC_DONE=("$WORK_DIR"/Anc*_consensus.done)
if [[ -e "${ANC_DONE[0]}" ]]; then
    echo
    echo "Ancestor genomes annotated: $(ls "$WORK_DIR"/Anc*_consensus.done 2>/dev/null | wc -l | tr -d ' ')"
fi

echo
echo "Outputs:"
echo "  consensus gene sets : $(ls "$WORK_DIR"/consensus_gene_set/*_consensus.gff3 2>/dev/null | wc -l | tr -d ' ')"
echo "  QC plots            : $(ls "$WORK_DIR"/plots/*.pdf 2>/dev/null | wc -l | tr -d ' ')"
echo "  gene family report  : $WORK_DIR/gene_family_analysis/summary.md"

# Surface recent errors from the snakemake log if present.
if [[ -f "$WORK_DIR/run.log" ]]; then
    echo
    ERRC="$(grep -icE 'error|exception|failed' "$WORK_DIR/run.log" 2>/dev/null || echo 0)"
    echo "run.log: $ERRC line(s) mentioning error/exception/failed (tail -f $WORK_DIR/run.log to follow)"
fi
echo "==================================================================="
