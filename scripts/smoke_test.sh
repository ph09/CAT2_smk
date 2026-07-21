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

# ── 3. cat2 package imports (validates the editable install) ─────────────────
python -c "import cat, tools; import cat.scheduler; import tools.hal; import tools.sqlInterface" \
    || fail "cat2 package import failed (did 'pip install --no-deps -e .' run?)"
pass "cat2 package imports (cat, tools, cat.scheduler, tools.hal)"

# ── 4. real Snakefile: parse + config validation + DAG build ─────────────────
# This parses the whole Snakefile, which calls `halStats` on the input HAL at
# import time to enumerate genomes. It therefore needs the HAL tools on PATH
# (cactus install) and the bundled test HAL, neither of which exists in a bare
# CI runner. Run it only when both are available; otherwise skip (not fail).
if command -v halStats >/dev/null 2>&1 && [[ -f test_data/vertebrates.hal ]]; then
    snakemake -n --configfile input.yaml --config work_dir=.smoke_work >/dev/null \
        || fail "dry-run failed (parse / config validation / DAG build)"
    pass "Snakefile parses, config validates, full DAG builds (dry-run)"
else
    echo "SKIP: full DAG dry-run (needs halStats on PATH + test_data/vertebrates.hal)"
fi

echo "SMOKE TEST OK"
