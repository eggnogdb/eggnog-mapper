#!/usr/bin/env python3
"""Smoke-check a downloaded eggnog-mapper data dir.

Verifies that:
  - eggnog.db exists, is a v7 schema, has the expected core tables
  - prots and protein_names have plausible row counts
  - taxa.db is reachable and has a species table
  - eggnog_proteins.dmnd is present (basic byte-level sanity)

Exits non-zero if any check fails. Suitable for CI / post-download
smoke tests in release-db.sh's verification step.
"""
from __future__ import annotations
import argparse
import os
import sqlite3
import sys
from pathlib import Path


def fail(msg):
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def ok(msg):
    print(f"OK: {msg}")


def check_file(path, min_size=0):
    p = Path(path)
    if not p.is_file():
        fail(f"missing file: {p}")
    if p.stat().st_size < min_size:
        fail(f"{p}: size {p.stat().st_size} < {min_size}")
    ok(f"{p.name} present ({p.stat().st_size:,} bytes)")


def check_eggnog_db(path):
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    cur = con.cursor()

    # Tables we expect in v7+ schema
    for tbl in ("prots", "sp_events", "event_index", "protein_names", "ogs"):
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (tbl,))
        if cur.fetchone() is None:
            fail(f"eggnog.db missing table: {tbl}")
    ok("eggnog.db has all expected tables")

    # Plausible row counts
    n_prots = cur.execute("SELECT count(*) FROM prots").fetchone()[0]
    n_names = cur.execute("SELECT count(*) FROM protein_names").fetchone()[0]
    if n_prots < 1000:
        fail(f"prots has only {n_prots} rows; expected >= 1000")
    if n_names < n_prots:
        fail(f"protein_names ({n_names}) < prots ({n_prots})")
    ok(f"prots={n_prots:,}  protein_names={n_names:,}")

    # v7 schema marker (integer protein IDs)
    cur.execute("PRAGMA table_info(prots)")
    cols = [r[1] for r in cur.fetchall()]
    if "id" not in cols:
        fail("prots schema lacks `id` column (not v7?)")
    ok("prots.id present (v7 schema)")
    con.close()


def check_taxa_db(path):
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    cur = con.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='species'")
    if cur.fetchone() is None:
        fail("taxa.db missing `species` table")
    n = cur.execute("SELECT count(*) FROM species").fetchone()[0]
    if n < 1000:
        fail(f"species has only {n} rows")
    ok(f"taxa.db species={n:,}")
    con.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True,
                    help="Directory containing eggnog.db, eggnog.taxa.db, "
                         "eggnog_proteins.dmnd")
    args = ap.parse_args()

    d = Path(args.data_dir)
    if not d.is_dir():
        fail(f"--data-dir is not a directory: {d}")

    check_file(d / "eggnog.db", min_size=10_000_000)
    check_file(d / "eggnog.taxa.db", min_size=1_000_000)
    check_file(d / "eggnog_proteins.dmnd", min_size=10_000_000)

    check_eggnog_db(d / "eggnog.db")
    check_taxa_db(d / "eggnog.taxa.db")

    print()
    print("All checks passed.")


if __name__ == "__main__":
    main()
