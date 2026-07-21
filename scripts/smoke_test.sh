#!/usr/bin/env bash
#
# CAT2 smoke test — fast checks that require no cluster and no cactus install.
# Runs locally and in CI. Exits non-zero on the first failed check.
#
#   1. snakemake is >= 9 (older versions corrupt run:-block f-strings).
#   2. Regression guard: multi-line / backslash-continued f-strings inside a
#      run: block are emitted intact by the installed snakemake.
#   3. The real Snakefile parses, the config validates, and the full DAG builds
#      (dry-run) on the bundled test_data.
#
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

pass() { echo "PASS: $1"; }
fail() { echo "FAIL: $1" >&2; exit 1; }

command -v snakemake >/dev/null 2>&1 || fail "snakemake not on PATH"

# ── 1. version ───────────────────────────────────────────────────────────────
SMK_VER="$(snakemake --version)"
[[ "${SMK_VER%%.*}" -ge 9 ]] || fail "snakemake ${SMK_VER} < 9 (run:-block f-strings unsafe)"
pass "snakemake ${SMK_VER} (>= 9)"

# ── 2. run:-block f-string regression guard ──────────────────────────────────
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cat > "$TMP/Snakefile" <<'SMK'
rule all:
    input: "out.txt"
rule a:
    output: "out.txt"
    run:
        x = "TOKEN"
        open(output[0], "w").write(f"""#!/bin/bash
echo start
python cat/consensus_runner.py \\
    --gp-list {x} \\
    --db-path DB
echo done
""")
SMK
( cd "$TMP" && snakemake --cores 1 -q >/dev/null 2>&1 ) || fail "self-test Snakefile did not run"
for needle in "echo start" "consensus_runner.py" "--gp-list TOKEN" "--db-path DB" "echo done"; do
    grep -qF -- "$needle" "$TMP/out.txt" || fail "run:-block f-string corrupted (missing: '$needle') — snakemake regression"
done
pass "run:-block f-strings emitted intact"

# ── 3. real Snakefile: parse + config validation + DAG build ─────────────────
snakemake -n --configfile input.yaml --config work_dir=.smoke_work >/dev/null \
    || fail "dry-run failed (parse / config validation / DAG build)"
pass "Snakefile parses, config validates, full DAG builds (dry-run)"

echo "SMOKE TEST OK"
