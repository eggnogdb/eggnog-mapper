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


def _run(itype: str, fixture: str, out_dir: str, extra=None) -> str:
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
    if extra:
        cmd += list(extra)
    subprocess.run(cmd, check=True, env=env, cwd=out_dir)
    return pjoin(out_dir, f"{itype}.emapper.annotations")


def _read_fasta(path: str) -> dict:
    """Return {query_name: sequence} using the header's first whitespace token."""
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = ""
            elif name is not None:
                seqs[name] += line.strip()
    return seqs


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


def test_md5_column(tmp_path):
    """`--md5` appends a final `md5` column equal to md5(query sequence)."""
    import hashlib

    ann = _run("proteins", CASES["proteins"], str(tmp_path), extra=["--md5"])
    seqs = _read_fasta(pjoin(FIX, CASES["proteins"]))

    header, rows = None, []
    with open(ann) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#query"):
                header = line.split("\t")
                continue
            rows.append(line.split("\t"))

    assert header is not None and header[-1] == "md5", (
        f"--md5 did not add a trailing 'md5' header column: {header}"
    )
    assert rows, "no annotation rows produced with --md5"
    for cols in rows:
        query, got = cols[0], cols[-1]
        assert query in seqs, f"query {query!r} not found in the input fixture"
        want = hashlib.md5(seqs[query].encode()).hexdigest()
        assert got == want, (
            f"md5 mismatch for {query}: got {got}, expected md5(seq)={want}"
        )
