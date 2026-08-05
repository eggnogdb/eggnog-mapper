#!/usr/bin/env python3
"""Download eggnog-mapper databases.

v3.0: refactored to fetch a `manifest.json` from the server first and
verify each artifact's sha256 against it before unpacking. Use
`--no-verify` to fall back to legacy unverified behaviour against
mirrors that don't carry a manifest (warns loudly).

The recursive HMMER-mode download path (-H / --hmmer) still uses
wget for its directory recursion; the per-clade hmm bundles are not
in the manifest scheme.

novel_fams (-F) was removed in v3.0 (the search mode is gone).
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from eggnogmapper.common import (
    DATA_PATH,
    existing_dir,
    get_data_path,
    get_eggnog_mmseqs_dbpath,
    get_eggnogdb_file,
    get_ncbitaxadb_file,
    get_pfam_dbpath,
    pexists,
    set_data_path,
)
from eggnogmapper.backends import DEFAULT_BACKEND, resolve_backend
from eggnogmapper.search.search_modes import (
    SEARCH_MODE_DIAMOND,
    get_eggnog_dmnd_db,
)
from eggnogmapper.utils import ask, colorify
from eggnogmapper.version import __VERSION__

if sys.version_info < (3, 9):
    sys.exit("Python < 3.9 is not supported")


# ---------------------------------------------------------------------------
# Manifest + checksum download primitives
# ---------------------------------------------------------------------------

USER_AGENT = f"eggnog-mapper/{__VERSION__} (download_eggnog_data.py)"
CHUNK = 1024 * 1024  # 1 MB stream chunks


class DownloadError(Exception):
    pass


def _http_get(url, dest_path=None):
    """GET a URL. Returns bytes if dest_path is None, else streams to disk
    and returns the sha256 hex digest of the downloaded bytes."""
    req = Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(req) as r:
            if dest_path is None:
                return r.read()
            h = hashlib.sha256()
            with open(dest_path, "wb") as f:
                while True:
                    buf = r.read(CHUNK)
                    if not buf:
                        break
                    h.update(buf)
                    f.write(buf)
            return h.hexdigest()
    except HTTPError as e:
        raise DownloadError(f"{url} -> HTTP {e.code} {e.reason}") from e
    except URLError as e:
        raise DownloadError(f"{url} -> {e.reason}") from e


def fetch_manifest(base_url):
    """GET base_url/manifest.json. Returns parsed dict or None on 404."""
    try:
        raw = _http_get(f"{base_url.rstrip('/')}/manifest.json")
        return json.loads(raw)
    except DownloadError as e:
        if "HTTP 404" in str(e):
            return None
        raise


def check_compat(manifest):
    """Refuse to install if mapper is older than manifest's minimum."""
    mn = manifest.get("min_mapper_version")
    if not mn:
        return
    cur = tuple(int(x) for x in __VERSION__.split("."))
    need = tuple(int(x) for x in mn.split("."))
    if cur < need:
        sys.exit(
            f"This DB (manifest min_mapper_version={mn}) is too new for "
            f"eggnog-mapper {__VERSION__}. Upgrade emapper first."
        )


def download_artifact(base_url, name, dest_dir, manifest, verify=True):
    """Stream-download one artifact, verifying sha256 against manifest.

    Returns the absolute path of the downloaded file.
    """
    url = f"{base_url.rstrip('/')}/{name}"
    dest = Path(dest_dir) / name
    print(colorify(f"  downloading {url}", "cyan"))
    got_sha = _http_get(url, dest_path=dest)

    if verify and manifest is not None:
        meta = (manifest.get("artifacts") or {}).get(name)
        if meta is None:
            raise DownloadError(f"manifest has no entry for {name!r}")
        want_sha = meta["sha256"]
        if got_sha != want_sha:
            dest.unlink(missing_ok=True)
            raise DownloadError(
                f"{name}: sha256 mismatch — expected {want_sha}, got {got_sha}; "
                "refusing to install. Re-run to retry."
            )
        size = meta.get("size")
        if size and dest.stat().st_size != size:
            raise DownloadError(
                f"{name}: size mismatch — expected {size}, got {dest.stat().st_size}"
            )
        print(colorify(f"  verified sha256 {got_sha[:12]}…", "green"))
    elif not verify:
        print(colorify(f"  WARNING: --no-verify set, sha256 not checked", "yellow"))

    return dest


def untar_gz(path, dest_dir):
    """Extract a .tar.gz into dest_dir; remove the archive."""
    with tarfile.open(path, "r:gz") as tf:
        tf.extractall(dest_dir)
    path.unlink()
    print(colorify(f"  extracted contents of {path.name} to {dest_dir}", "green"))


# ---------------------------------------------------------------------------
# Manifest-aware downloads (small, individually-named artifacts)
# ---------------------------------------------------------------------------

# The data server (data.cgmlab.org/eggnog-mapper/emapper-<release>/data/) serves
# every artifact UNCOMPRESSED, so downloads stream straight into the data dir —
# no gunzip/untar step. Multi-file databases (mmseqs, pfam) remain tarballs.

def download_annotations(base_url, data_path, manifest, verify, force):
    download_artifact(base_url, "eggnog.db", data_path, manifest, verify)
    # The prebuilt taxid cache is shipped next to the DB; grab it so first-run
    # doesn't have to rebuild it (~1 min scan). Optional: skip if absent.
    try:
        download_artifact(base_url, "eggnog.db.taxids.bin", data_path, manifest, verify)
    except DownloadError as exc:
        print(colorify(f"  (optional) eggnog.db.taxids.bin not available ({exc}); "
                       "it will be rebuilt on first run.", "yellow"))


# Fallback source when the data server does not ship the OBO. GO namespaces
# (MF/BP/CC) are stable across releases, so the current OBO Foundry release is a
# safe namespace map.
GO_OBO_FALLBACK_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"


def download_go_obo(base_url, data_path, manifest, verify, force):
    """Fetch the GO namespace OBO (uncompressed) used for the MF/BP/CC GO
    cascade. Without it, annotation silently falls back to legacy combined-GO.
    Falls back to the current OBO Foundry release if the server lacks it.
    Written to ``<data_path>/go-basic.obo`` where the annotator looks first.
    """
    dest = os.path.join(data_path, "go-basic.obo")
    try:
        download_artifact(base_url, "go-basic.obo", data_path, manifest, verify)
        return
    except DownloadError as exc:
        print(colorify(
            f"  go-basic.obo not on the data server ({exc}); "
            f"falling back to {GO_OBO_FALLBACK_URL}", "yellow"))
    _http_get(GO_OBO_FALLBACK_URL, dest_path=dest)
    print(colorify(f"  downloaded {dest}", "green"))


def download_taxa(base_url, data_path, manifest, verify, force):
    # Served as two plain files (no tarball).
    download_artifact(base_url, "eggnog.taxa.db", data_path, manifest, verify)
    download_artifact(base_url, "eggnog.taxa.db.traverse.pkl", data_path, manifest, verify)


def download_diamond_db(base_url, data_path, manifest, verify, force):
    download_artifact(base_url, "eggnog_proteins.dmnd", data_path, manifest, verify)


def download_mmseqs_db(base_url, data_path, manifest, verify, force):
    p = download_artifact(base_url, "mmseqs.tar.gz", data_path, manifest, verify)
    untar_gz(p, data_path)


def download_pfam_db(base_url, data_path, manifest, verify, force):
    p = download_artifact(base_url, "pfam.tar.gz", data_path, manifest, verify)
    untar_gz(p, data_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass


def main():
    parser = argparse.ArgumentParser(formatter_class=CustomFormatter,
                                     description=__doc__)

    parser.add_argument("-D", action="store_true", dest="skip_diamond",
                        help="Do not install the diamond database")
    parser.add_argument("-P", action="store_true", dest="pfam",
                        help="Install the Pfam database (de novo / realign)")
    parser.add_argument("-M", action="store_true", dest="mmseqs",
                        help="Install the MMseqs2 database (-m mmseqs)")

    parser.add_argument("-y", action="store_true", dest="allyes",
                        help='Assume "yes" to all questions')
    parser.add_argument("-f", action="store_true", dest="force",
                        help="Force download even if files exist")
    parser.add_argument("-s", action="store_true", dest="simulate",
                        help="Simulate; nothing downloaded")
    parser.add_argument("-q", action="store_true", dest="quiet",
                        help="Quiet mode")

    parser.add_argument("--data_dir", metavar="", type=existing_dir,
                        help="Directory to use for DATA_PATH.")

    # Release version pins the server subdirectory: emapper-<release>/data/.
    parser.add_argument("--release", type=str, default=__VERSION__.split("-")[0],
                        help=("eggNOG-mapper data release to download. Pins the server "
                              "subdirectory data.cgmlab.org/eggnog-mapper/emapper-<release>/data/. "
                              "Defaults to this emapper's release version."))
    parser.add_argument("--no-verify", action="store_true", dest="no_verify",
                        help=("Skip checksum verification when the server ships a "
                              "manifest.json. Currently a no-op when no manifest is present."))

    args = parser.parse_args()

    set_data_path(resolve_backend(DEFAULT_BACKEND))
    if "EGGNOG_DATA_DIR" in os.environ:
        set_data_path(os.environ["EGGNOG_DATA_DIR"])
    if args.data_dir:
        set_data_path(args.data_dir)
    data_path = get_data_path()

    base_url = f"https://data.cgmlab.org/eggnog-mapper/emapper-{args.release}/data"
    print(colorify(f"Release: {args.release}  base URL: {base_url}", "lblue"))

    # Optional integrity check: verify checksums only if the server ships a
    # manifest.json. The data server currently serves plain files without one,
    # so downloads proceed unverified (with a warning) rather than aborting.
    manifest = None
    verify = not args.no_verify
    if verify:
        try:
            manifest = fetch_manifest(base_url)
        except DownloadError as e:
            print(colorify(f"WARNING: could not fetch manifest.json ({e}); "
                           "downloads will NOT be checksum-verified.", "yellow"))
        if manifest is not None:
            check_compat(manifest)
            print(colorify(
                f"manifest OK: schema={manifest.get('schema_version')} "
                f"build_date={manifest.get('build_date')} "
                f"min_mapper={manifest.get('min_mapper_version')}",
                "green"))
        else:
            print(colorify("No manifest.json on the server; downloads will NOT be "
                           "checksum-verified.", "yellow"))
    else:
        print(colorify("WARNING: --no-verify; downloads will NOT be checksum-verified",
                       "yellow"))

    if args.simulate:
        print(colorify("simulate mode: would download artifacts under", "cyan"),
              data_path)
        return

    # ---- Annotation DB ----
    if args.force or not pexists(get_eggnogdb_file()):
        if args.allyes or ask("Download main annotation database?") == "y":
            print(colorify(f"Downloading eggnog.db at {data_path}…", "green"))
            download_annotations(base_url, data_path, manifest, verify, args.force)
        else:
            print("Skipping")
    elif not args.quiet:
        print(colorify("Skipping eggnog.db (already present). Use -f to force.", "lblue"))

    # ---- GO namespace OBO (needed for the MF/BP/CC GO cascade) ----
    obo_dest = os.path.join(data_path, "go-basic.obo")
    if args.force or not pexists(obo_dest):
        if args.allyes or ask("Download GO namespace OBO (go-basic.obo, ~35 MB)?") == "y":
            print(colorify(f"Downloading go-basic.obo at {data_path}…", "green"))
            download_go_obo(base_url, data_path, manifest, verify, args.force)
        else:
            print(colorify("Skipping go-basic.obo — GO output will use the legacy "
                           "combined cascade (no MF/BP/CC split).", "yellow"))
    elif not args.quiet:
        print(colorify("Skipping go-basic.obo (already present). Use -f to force.", "lblue"))

    # ---- NCBI taxa ----
    if args.force or not pexists(get_ncbitaxadb_file()):
        if args.allyes or ask("Download taxa database?") == "y":
            print(colorify(f"Downloading eggnog.taxa.db at {data_path}…", "green"))
            download_taxa(base_url, data_path, manifest, verify, args.force)
        else:
            print("Skipping")
    elif not args.quiet:
        print(colorify("Skipping eggnog.taxa.db (already present). Use -f to force.",
                       "lblue"))

    # ---- Diamond DB ----
    if not args.skip_diamond and (args.force or not pexists(
            get_eggnog_dmnd_db(None, SEARCH_MODE_DIAMOND, data_path))):
        if args.allyes or ask("Download diamond database (~22 GB)?") == "y":
            print(colorify(f"Downloading eggnog_proteins.dmnd at {data_path}…", "green"))
            download_diamond_db(base_url, data_path, manifest, verify, args.force)
        else:
            print("Skipping")
    elif not args.quiet:
        print(colorify("Skipping diamond DB (already present or -D). Use -f to force.",
                       "lblue"))

    # ---- PFAM (optional; may not be published for a given release) ----
    if args.pfam and (args.force or not pexists(get_pfam_dbpath())):
        if args.allyes or ask("Download Pfam database?") == "y":
            print(colorify(f"Downloading Pfam at {data_path}…", "green"))
            try:
                download_pfam_db(base_url, data_path, manifest, verify, args.force)
            except DownloadError as exc:
                print(colorify(f"  Pfam not available for release {args.release} "
                               f"({exc}); skipping.", "yellow"))
        else:
            print("Skipping")
    elif not args.quiet:
        print(colorify("Skipping Pfam DB (already present). Use -P -f to force.", "lblue"))

    # ---- MMseqs (optional; may not be published for a given release) ----
    if args.mmseqs and (args.force or not pexists(get_eggnog_mmseqs_dbpath())):
        if args.allyes or ask("Download MMseqs2 database?") == "y":
            print(colorify(f"Downloading MMseqs2 at {data_path}…", "green"))
            try:
                download_mmseqs_db(base_url, data_path, manifest, verify, args.force)
            except DownloadError as exc:
                print(colorify(f"  MMseqs2 not available for release {args.release} "
                               f"({exc}); skipping.", "yellow"))
        else:
            print("Skipping")
    elif not args.quiet:
        print(colorify("Skipping MMseqs2 DB (already present). Use -M -f to force.", "lblue"))

    print(colorify("Finished.", "green"))


if __name__ == "__main__":
    main()
