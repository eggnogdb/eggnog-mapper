##
## Self-test: verify an eggNOG-mapper installation reproduces known annotations.
##
"""Download the small reference self-test dataset and check that this build
reproduces the expected annotations across every input type.

The dataset is a self-consistent subsample of the reference eggNOG v7 mapper DB
(nif-operon query seeds plus their annotation-transfer closure) together with
per-``--itype`` query fixtures (proteins / CDS / genome / metagenome) and golden
outputs. Running emapper against it exercises the full pipeline — DIAMOND search,
auto-translation of CDS, prodigal gene prediction, and the annotation cascade —
entirely offline once downloaded.

Invoked via ``emapper.py --selftest``. Because it runs inside the same process
environment as a normal run (e.g. inside the Apptainer image), the bundled
DIAMOND / prodigal binaries are already on ``PATH`` — no external tools needed.
"""

import os
import sys
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from .version import __VERSION__
from .common import get_data_release
from .utils import colorify
from .emapperException import EmapperException

# itype -> (fixture filename, golden filename). The golden basename matches the
# ``-o`` project name used below, i.e. ``<itype>.emapper.annotations``.
SELFTEST_CASES = {
    "proteins":   "proteins.faa",
    "CDS":        "cds.fna",
    "genome":     "genome.fna",
    "metagenome": "metagenome.fna",
}

CHUNK = 1024 * 1024  # 1 MB stream chunks


class SelfTestError(EmapperException):
    """Raised when the self-test cannot be set up (download/extract failure)."""


def _selftest_url(release):
    """Return the bundle URL for a MAJOR.MINOR data release
    (``$EGGNOG_SELFTEST_URL`` overrides). The self-test data lives next to the
    DBs it was subsampled from, so it is pinned by the same ``emapper-<MAJOR.MINOR>/``
    folder and shared across a patch series."""
    override = os.environ.get("EGGNOG_SELFTEST_URL")
    if override:
        return override
    return (f"https://data.cgmlab.org/eggnog-mapper/emapper-{release}"
            f"/selftest/eggnog-mapper-selftest.tar.gz")


def _download(url, dest_path):
    """Stream ``url`` to ``dest_path``."""
    req = Request(url, headers={"User-Agent": f"eggnog-mapper/{__VERSION__} (selftest)"})
    try:
        with urlopen(req) as r, open(dest_path, "wb") as f:
            while True:
                buf = r.read(CHUNK)
                if not buf:
                    break
                f.write(buf)
    except HTTPError as e:
        raise SelfTestError(f"{url} -> HTTP {e.code} {e.reason}") from e
    except URLError as e:
        raise SelfTestError(f"{url} -> {e.reason}") from e


def _obtain_bundle(bundle_dir, release, workdir):
    """Resolve the dataset root holding ``data/``, ``fixtures/`` and ``golden/``.

    If ``bundle_dir`` already contains those, it is used as-is (no download).
    Otherwise the release bundle is downloaded into ``workdir`` and extracted.

    Args:
        bundle_dir: Optional pre-downloaded dataset dir (e.g. ``--data_dir``).
        release: Release string pinning the server subdirectory.
        workdir: Scratch directory for the download + extraction.

    Returns:
        Path to the dataset root.
    """
    def _has_all(root):
        return all((Path(root) / d).is_dir() for d in ("data", "fixtures", "golden"))

    if bundle_dir and _has_all(bundle_dir):
        print(colorify(f"Using local self-test dataset at {bundle_dir}", "green"))
        return Path(bundle_dir)

    url = _selftest_url(release)
    print(colorify(f"Downloading self-test dataset:\n  {url}", "cyan"))
    archive = Path(workdir) / "selftest.tar.gz"
    _download(url, archive)
    with tarfile.open(archive, "r:gz") as tf:
        tf.extractall(workdir)
    archive.unlink(missing_ok=True)

    root = Path(workdir)
    if not _has_all(root):
        # tolerate a single top-level wrapper directory inside the tarball
        subdirs = [p for p in root.iterdir() if p.is_dir()]
        for cand in subdirs:
            if _has_all(cand):
                root = cand
                break
    if not _has_all(root):
        raise SelfTestError(
            f"downloaded bundle is missing data/, fixtures/ or golden/ (looked in {root})")
    print(colorify(f"Dataset ready at {root}", "green"))
    return root


def _emapper_cmd():
    """Return ``[python, emapper.py]`` for re-invoking emapper as a subprocess.

    Works both from a source checkout (sibling ``emapper.py``) and from an
    installed image (``emapper.py`` on ``PATH``).
    """
    script = Path(__file__).resolve().parents[1] / "emapper.py"
    if not script.exists():
        found = shutil.which("emapper.py")
        if not found:
            raise SelfTestError("could not locate the emapper.py entrypoint")
        script = Path(found)
    return [sys.executable, str(script)]


def _data_rows(path):
    """Return the set of annotation data rows, ignoring ``##`` metadata lines
    (command line, timestamps, versions, absolute paths — all run-specific)."""
    rows = set()
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line and not line.startswith("##"):
                rows.add(line)
    return rows


def run_selftest(bundle_dir=None, release=None, cpu=1):
    """Run the full self-test and return a process exit code (0 = all passed).

    Args:
        bundle_dir: Optional dir already holding the dataset (skips download).
        release: Release to download; defaults to this emapper's version.
        cpu: CPUs to pass through to each emapper sub-run.

    Returns:
        ``0`` if every itype reproduces its golden annotations, else ``1``.
    """
    release = release or get_data_release()
    base_cmd = _emapper_cmd()

    workdir = tempfile.mkdtemp(prefix="emapper_selftest_")
    try:
        root = _obtain_bundle(bundle_dir, release, workdir)
        data = root / "data"
        fixtures = root / "fixtures"
        golden = root / "golden"
        out_dir = Path(workdir) / "out"
        out_dir.mkdir(exist_ok=True)

        print(colorify(f"\neggNOG-mapper self-test  (emapper {__VERSION__})", "lblue"))
        print(f"  dataset : {root}")
        print(f"  itypes  : {', '.join(SELFTEST_CASES)}\n")

        failures = []
        for itype, fixture in SELFTEST_CASES.items():
            print(colorify(f"[{itype}]", "yellow"), f"running {fixture} ...")
            cmd = base_cmd + [
                "-i", str(fixtures / fixture),
                "--itype", itype,
                "-m", "diamond",
                "--data_dir", str(data),
                "--output_dir", str(out_dir),
                "-o", itype,
                "--cpu", str(cpu),
                "--override",
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            got = out_dir / f"{itype}.emapper.annotations"
            want = golden / f"{itype}.annotations"

            if proc.returncode != 0 or not got.exists():
                failures.append(itype)
                print(colorify(f"  FAIL: emapper run errored (rc={proc.returncode})", "red"))
                tail = "\n".join(proc.stderr.strip().splitlines()[-8:])
                if tail:
                    print("\n".join("    " + l for l in tail.splitlines()))
                continue
            if not want.exists():
                failures.append(itype)
                print(colorify(f"  FAIL: golden missing ({want})", "red"))
                continue

            got_rows, want_rows = _data_rows(got), _data_rows(want)
            missing, extra = want_rows - got_rows, got_rows - want_rows
            if missing or extra:
                failures.append(itype)
                print(colorify("  FAIL: annotations differ from golden", "red"))
                for r in sorted(missing)[:2]:
                    print(colorify(f"    - expected: {r[:120]}", "red"))
                for r in sorted(extra)[:2]:
                    print(colorify(f"    + got:      {r[:120]}", "red"))
            else:
                n = sum(1 for r in got_rows if not r.startswith("#"))
                print(colorify(f"  PASS ({n} queries annotated)", "green"))

        print()
        if failures:
            print(colorify(
                f"SELF-TEST FAILED for: {', '.join(failures)}", "red"))
            return 1
        print(colorify(
            "SELF-TEST PASSED — all input types reproduce the reference "
            "annotations.", "green"))
        return 0
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
