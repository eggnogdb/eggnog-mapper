#!/usr/bin/env python3
"""Generate golden annotation outputs for the self-test dataset.

Runs emapper once per ``--itype`` against ``test_datasets/data`` and stores the
resulting ``.annotations`` as ``test_datasets/golden/<itype>.annotations``.
Re-run this only when the expected output legitimately changes; the pytest
self-test (``test_selftest.py``) diffs live runs against these goldens.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from os.path import join as pjoin

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DATA = pjoin(HERE, "data")
FIX = pjoin(HERE, "fixtures")
GOLDEN = pjoin(HERE, "golden")
EMAPPER = pjoin(REPO, "emapper.py")
DIAMOND_BIN_DIR = "/eggnog-eco/bin"

# itype -> fixture filename
CASES = {
    "proteins": "proteins.faa",
    "CDS": "cds.fna",
    "genome": "genome.fna",
    "metagenome": "metagenome.fna",
}


def run_case(itype: str, fixture: str, out_dir: str, name: str) -> str:
    """Run emapper for one itype; return path to the .annotations file."""
    env = dict(os.environ)
    env["PATH"] = DIAMOND_BIN_DIR + os.pathsep + env.get("PATH", "")
    cmd = [
        sys.executable, EMAPPER,
        "-i", pjoin(FIX, fixture),
        "--itype", itype,
        "-m", "diamond",
        "--data_dir", DATA,
        "--output_dir", out_dir,
        "-o", name,
        "--cpu", "1",
        "--override",
    ]
    subprocess.run(cmd, check=True, env=env, cwd=out_dir)
    ann = pjoin(out_dir, f"{name}.emapper.annotations")
    if not os.path.exists(ann):
        raise SystemExit(f"No annotations produced for itype={itype}: {ann}")
    return ann


def main():
    os.makedirs(GOLDEN, exist_ok=True)
    for itype, fixture in CASES.items():
        print(f"=== generating golden for --itype {itype} ===")
        with tempfile.TemporaryDirectory() as td:
            ann = run_case(itype, fixture, td, itype)
            dst = pjoin(GOLDEN, f"{itype}.annotations")
            shutil.copy2(ann, dst)
            print(f"    wrote {dst}")
    print("\nAll goldens generated.")


if __name__ == "__main__":
    main()
