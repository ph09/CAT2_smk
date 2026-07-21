#!/usr/bin/env python3
"""
Generate a REAL augMP PSL by mapping the miniprot PSL onto augMP transcript IDs.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path


AUG_PREFIX = 'augMP-'
_COPY_PAT = re.compile(r'_(\d+)$')


def _strip_copy(name: str) -> tuple[str, str]:
    """Split ``XXX_2`` → (``XXX``, ``_2``). For ``XXX`` (no suffix) → (``XXX``, '')."""
    m = _COPY_PAT.search(name)
    if not m:
        return name, ''
    return name[:m.start()], name[m.start():]


def _read_augmp_ids(augmp_gp: Path) -> list[str]:
    ids: list[str] = []
    with augmp_gp.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            name = line.split('\t', 1)[0]
            if not name.startswith(AUG_PREFIX):
                continue
            ids.append(name)
    return ids


def _psl_lookup_keys(source_id: str) -> list[str]:
    """Keys to match augMP-<source> against miniprot PSL qName (column 10)."""
    keys: list[str] = []

    def add(s: str) -> None:
        if s and s not in keys:
            keys.append(s)

    add(source_id)
    if source_id.startswith('rna-'):
        add(source_id[4:])
    else:
        add(f'rna-{source_id}')
    base, copy_suffix = _strip_copy(source_id)
    if copy_suffix:
        add(base)
        if base.startswith('rna-'):
            add(base[4:])
        else:
            add(f'rna-{base}')
    return keys


def _index_miniprot_psl(miniprot_psl: Path) -> dict[str, list[str]]:
    """Map miniprot query name → list of PSL rows (raw, with original qName)."""
    idx: dict[str, list[str]] = defaultdict(list)
    with miniprot_psl.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 21:
                continue
            try:
                int(f[0])
            except ValueError:
                continue
            idx[f[9]].append(line.rstrip('\n'))
    return idx


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--augmp-gp', required=True, type=Path,
                    help='Augustus augMP GenePred (column 0 = augMP-<src>[_<copy>])')
    ap.add_argument('--miniprot-psl', required=True, type=Path,
                    help='Real miniprot PSL from convert_miniprot_to_genepred.py')
    ap.add_argument('--out-psl', required=True, type=Path,
                    help='Output augMP PSL (query names rewritten to augMP-*)')
    args = ap.parse_args()

    augmp_ids = _read_augmp_ids(args.augmp_gp)
    print(f"  read {len(augmp_ids):,} augMP records from {args.augmp_gp}",
          file=sys.stderr)

    mp_index = _index_miniprot_psl(args.miniprot_psl)
    print(f"  indexed {sum(len(v) for v in mp_index.values()):,} PSL rows "
          f"({len(mp_index):,} unique queries) from {args.miniprot_psl}",
          file=sys.stderr)

    n_matched = 0
    n_unmatched = 0
    with args.out_psl.open('w') as fout:
        for aug_id in augmp_ids:
            source_id = aug_id[len(AUG_PREFIX):]
            psl_rows = None
            for key in _psl_lookup_keys(source_id):
                psl_rows = mp_index.get(key)
                if psl_rows:
                    break
            if not psl_rows:
                n_unmatched += 1
                continue
            # If multiple copies (paralogs), there will be one PSL row per
            # source_id variant (_2, _3, ...). The augMP GP only tracks one
            # template-derived record per source_id, so take the first PSL
            # row (best alignment) for this source.
            row = psl_rows[0]
            f = row.split('\t')
            f[9] = aug_id
            fout.write('\t'.join(f) + '\n')
            n_matched += 1

    print(f"  wrote {n_matched:,} augMP PSL rows "
          f"({n_unmatched:,} augMP records had no miniprot PSL match → stay NaN)",
          file=sys.stderr)

    if len(augmp_ids) > 0 and n_matched == 0 and sum(len(v) for v in mp_index.values()) > 0:
        print(
            "  ERROR: augMP genePred has records but no PSL rows were written "
            "(check miniprot PSL query names vs augMP IDs)",
            file=sys.stderr,
        )
        sys.exit(1)
    if len(augmp_ids) > 0 and n_matched < len(augmp_ids) * 0.5:
        print(
            f"  WARNING: only {n_matched}/{len(augmp_ids)} augMP records matched miniprot PSL",
            file=sys.stderr,
        )


if __name__ == '__main__':
    main()
