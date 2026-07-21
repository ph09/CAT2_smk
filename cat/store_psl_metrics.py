#!/usr/bin/env python3
"""
store_psl_metrics.py
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
import time
from pathlib import Path

import tools.nameConversions


def _base_id(aln_id: str) -> str:
    """Normalize AlignmentId to a base transcript ID for ref mapping."""
    return tools.nameConversions.alignment_id_to_ref_transcript_id(aln_id)


def _load_ref_name2(ref_gp_path: Path) -> dict[str, str]:
    """Reference GenePredExt: name (tx) in col0, name2 (gene) in col11."""
    name2: dict[str, str] = {}
    with ref_gp_path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 12:
                name2[f[0]] = f[11]
    return name2


def _iter_psl_metrics(psl_path: Path):
    """Yield (aln_id, cov_pct, ident_pct)."""
    with psl_path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 21:
                continue
            try:
                matches = int(f[0])
                mismatches = int(f[1])
                rep_matches = int(f[2])
                q_size = int(f[10])
                block_sizes = [int(x) for x in f[18].rstrip(",").split(",") if x]
            except ValueError:
                continue

            if q_size <= 0:
                continue

            aligned = matches + mismatches + rep_matches
            cov = (sum(block_sizes) / q_size) * 100.0
            ident = (matches / aligned) * 100.0 if aligned > 0 else 0.0
            aln_id = f[9]
            yield aln_id, cov, ident


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db-path", required=True, help="Genome SQLite DB path")
    ap.add_argument("--psl", required=True, help="Filtered PSL file path")
    ap.add_argument(
        "--mode",
        required=True,
        choices=["txTM", "transMap", "transMap_pairwise", "augMP"],
        help="Alignment mode whose mRNA metrics table should be updated",
    )
    ap.add_argument("--ref-gp", required=True, help="Reference GenePredExt for tx→gene mapping")
    args = ap.parse_args()

    db_path = Path(args.db_path)
    psl_path = Path(args.psl)
    ref_gp_path = Path(args.ref_gp)

    table = f"mRNA_{args.mode}_Metrics"

    ref_name2 = _load_ref_name2(ref_gp_path)

    # Collect metrics rows and ids for cleanup.
    ids: list[str] = []
    rows: list[tuple[str, str, str, str, str]] = []

    for aln_id, cov, ident in _iter_psl_metrics(psl_path):
        base = _base_id(aln_id)
        gene_id = ref_name2.get(base, "")
        tx_id = base
        ids.append(aln_id)
        rows.append((aln_id, gene_id, tx_id, "AlnCoverage", f"{cov:.3f}"))
        rows.append((aln_id, gene_id, tx_id, "AlnIdentity", f"{ident:.3f}"))

    # An empty filtered PSL is legitimate: a small/divergent genome can yield zero
    # alignments for a mode (e.g. transMap_pairwise). This is not an error — fall
    # through so the (empty) metrics table is still created and the rule succeeds;
    # downstream consensus already handles modes with no records.
    if not rows:
        print(f"WARNING: no PSL records for mode {args.mode}: {psl_path} "
              f"(creating empty {table} table)", file=sys.stderr)

    # De-duplicate ids to keep temp table small.
    ids = sorted(set(ids))

    # Multiple Snakemake jobs used to hit the same genome DB in parallel → "database is locked".
    # Use a long busy timeout, WAL, and retries; Snakefile also serializes these rules per genome.
    max_attempts = 12
    for attempt in range(max_attempts):
        con = None
        try:
            con = sqlite3.connect(str(db_path), timeout=120.0)
            con.execute("PRAGMA journal_mode=WAL")
            cur = con.cursor()

            # Ensure table exists with expected columns (matches current DB schema).
            cur.execute(
                f"""
                CREATE TABLE IF NOT EXISTS {table} (
                    AlignmentId  TEXT,
                    GeneId       TEXT,
                    TranscriptId TEXT,
                    classifier   TEXT,
                    value        TEXT
                )
                """
            )
            cur.execute(f"CREATE INDEX IF NOT EXISTS ix_{table}_AlignmentId ON {table}(AlignmentId)")

            cur.execute("DROP TABLE IF EXISTS _tmp_aln_ids")
            cur.execute("CREATE TEMP TABLE _tmp_aln_ids(id TEXT PRIMARY KEY)")
            cur.executemany("INSERT INTO _tmp_aln_ids(id) VALUES (?)", [(i,) for i in ids])

            # Remove any existing PSL-derived coverage/identity rows for these AlignmentIds.
            cur.execute(
                f"""
                DELETE FROM {table}
                WHERE classifier IN ('AlnCoverage', 'AlnIdentity')
                  AND AlignmentId IN (SELECT id FROM _tmp_aln_ids)
                """
            )

            cur.executemany(
                f"INSERT INTO {table} (AlignmentId, GeneId, TranscriptId, classifier, value) VALUES (?, ?, ?, ?, ?)",
                rows,
            )
            con.commit()
            con.close()
            break
        except sqlite3.OperationalError as e:
            if con is not None:
                try:
                    con.close()
                except Exception:
                    pass
            if "locked" in str(e).lower() and attempt + 1 < max_attempts:
                time.sleep(2.0 * (attempt + 1))
                continue
            raise


if __name__ == "__main__":
    main()

