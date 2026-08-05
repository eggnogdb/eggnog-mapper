#!/usr/bin/env python3
"""Offline self-test: emapper over every --itype against the downloadable
mini dataset, diffed against committed golden annotations.

Requires the ``test_datasets/data`` + ``test_datasets/fixtures`` +
``test_datasets/golden`` directories (produced by ``make_test_dataset.py`` and
``gen_golden.py``; the data itself is downloaded, not stored in the repo). If
those are absent, all tests skip. A DIAMOND binary must be importable via PATH
(the repo ships one under ``/eggnog-eco/bin`` for local runs).

Run:  pytest test_datasets/test_selftest.py -v
"""
from __future__ import annotations

import os
import subprocess
import sys
from os.path import join as pjoin

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DATA = pjoin(HERE, "data")
FIX = pjoin(HERE, "fixtures")
GOLDEN = pjoin(HERE, "golden")
EMAPPER = pjoin(REPO, "emapper.py")
DIAMOND_BIN_DIR = "/eggnog-eco/bin"

CASES = {
    "proteins": "proteins.faa",
    "CDS": "cds.fna",
    "genome": "genome.fna",
    "metagenome": "metagenome.fna",
}

_have_data = os.path.isdir(DATA) and os.path.isdir(FIX) and os.path.isdir(GOLDEN)
pytestmark = pytest.mark.skipif(
    not _have_data,
    reason="self-test dataset not present (download test_datasets/data etc.)",
)


def _data_rows(path: str) -> set[str]:
    """Return the set of annotation data rows, ignoring ``##`` metadata lines
    (command line, timestamps, versions, absolute paths)."""
    rows = set()
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            rows.add(line)
    return rows


def _run(itype: str, fixture: str, out_dir: str) -> str:
    env = dict(os.environ)
    env["PATH"] = DIAMOND_BIN_DIR + os.pathsep + env.get("PATH", "")
    cmd = [
        sys.executable, EMAPPER,
        "-i", pjoin(FIX, fixture),
        "--itype", itype,
        "-m", "diamond",
        "--data_dir", DATA,
        "--output_dir", out_dir,
        "-o", itype,
        "--cpu", "1",
        "--override",
    ]
    subprocess.run(cmd, check=True, env=env, cwd=out_dir)
    return pjoin(out_dir, f"{itype}.emapper.annotations")


@pytest.mark.parametrize("itype,fixture", list(CASES.items()))
def test_itype_annotations_match_golden(itype, fixture, tmp_path):
    golden = pjoin(GOLDEN, f"{itype}.annotations")
    if not os.path.exists(golden):
        pytest.skip(f"golden missing for {itype}")
    ann = _run(itype, fixture, str(tmp_path))
    assert os.path.exists(ann), f"no annotations produced for {itype}"
    got, want = _data_rows(ann), _data_rows(golden)
    missing = want - got
    extra = got - want
    assert not missing and not extra, (
        f"[{itype}] annotation mismatch vs golden\n"
        f"  missing ({len(missing)}): {sorted(missing)[:3]}\n"
        f"  extra   ({len(extra)}): {sorted(extra)[:3]}"
    )


@pytest.mark.parametrize("itype", ["proteins", "CDS"])
def test_nifh_seed_annotated_correctly(itype, tmp_path):
    """Biological sanity: the nifH query must recover the Fer4_NifH OG."""
    ann = _run(itype, CASES[itype], str(tmp_path))
    hits = [r for r in _data_rows(ann) if r.split("\t", 1)[0].endswith("_nifH")]
    assert hits, f"no nifH query row in {itype} output"
    assert any("Fer4_NifH" in r for r in hits), (
        f"nifH query did not recover Fer4_NifH OG in {itype}: {hits}"
    )
