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
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from eggnogmapper.common import (
    DATA_PATH,
    HMMPRESS,
    existing_dir,
    get_data_path,
    get_eggnog_mmseqs_dbpath,
    get_eggnogdb_file,
    get_hmmer_base_dbpath,
    get_ncbitaxadb_file,
    get_pfam_dbpath,
    pexists,
    set_data_path,
)
from eggnogmapper.backends import DEFAULT_BACKEND, get_backend_names, resolve_backend
from eggnogmapper.search.search_modes import (
    SEARCH_MODE_DIAMOND,
    get_eggnog_dmnd_db,
)
from eggnogmapper.utils import ask, ask_name, colorify
from eggnogmapper.version import __DB_VERSION__, __VERSION__

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


def gunzip_in_place(path, force=False):
    """gunzip path → strip .gz suffix; remove the .gz."""
    out = path.with_suffix("")
    if out.exists() and not force:
        raise DownloadError(f"{out} already exists. Use -f to overwrite.")
    with gzip.open(path, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    path.unlink()
    print(colorify(f"  extracted {out}", "green"))


def untar_gz(path, dest_dir):
    """Extract a .tar.gz into dest_dir; remove the archive."""
    with tarfile.open(path, "r:gz") as tf:
        tf.extractall(dest_dir)
    path.unlink()
    print(colorify(f"  extracted contents of {path.name} to {dest_dir}", "green"))


# ---------------------------------------------------------------------------
# Manifest-aware downloads (small, individually-named artifacts)
# ---------------------------------------------------------------------------

def download_annotations(base_url, data_path, manifest, verify, force):
    p = download_artifact(base_url, "eggnog.db.gz", data_path, manifest, verify)
    gunzip_in_place(p, force=force)


# Fallback source when the DB server does not ship a build-time OBO. GO
# namespaces (MF/BP/CC) are stable across releases, so the current OBO Foundry
# release is a safe namespace map.
GO_OBO_FALLBACK_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"


def download_go_obo(base_url, data_path, manifest, verify, force):
    """Fetch the GO namespace OBO used for the MF/BP/CC GO cascade.

    Without it, annotation silently falls back to the legacy combined-GO
    output. Prefers the checksum-verified build-time artifact from the eggNOG
    download server (``go-basic.obo.gz`` in the manifest); falls back to the
    current OBO Foundry release when the server does not ship it. Written to
    ``<data_path>/go-basic.obo`` where the annotator looks for it first.
    """
    dest = os.path.join(data_path, "go-basic.obo")
    try:
        p = download_artifact(base_url, "go-basic.obo.gz", data_path, manifest, verify)
        gunzip_in_place(p, force=force)
        return
    except Exception as exc:
        print(colorify(
            f"  go-basic.obo.gz not available from the DB server ({exc}); "
            f"falling back to {GO_OBO_FALLBACK_URL}", "yellow"))
    _http_get(GO_OBO_FALLBACK_URL, dest_path=dest)
    print(colorify(f"  downloaded {dest}", "green"))


def download_taxa(base_url, data_path, manifest, verify, force):
    p = download_artifact(base_url, "eggnog.taxa.tar.gz", data_path, manifest, verify)
    untar_gz(p, data_path)


def download_diamond_db(base_url, data_path, manifest, verify, force):
    p = download_artifact(base_url, "eggnog_proteins.dmnd.gz", data_path, manifest, verify)
    gunzip_in_place(p, force=force)


def download_mmseqs_db(base_url, data_path, manifest, verify, force):
    p = download_artifact(base_url, "mmseqs.tar.gz", data_path, manifest, verify)
    untar_gz(p, data_path)


def download_pfam_db(base_url, data_path, manifest, verify, force):
    p = download_artifact(base_url, "pfam.tar.gz", data_path, manifest, verify)
    untar_gz(p, data_path)


# ---------------------------------------------------------------------------
# HMMER recursive download (kept on wget; not in manifest scheme)
# ---------------------------------------------------------------------------

def download_hmm_database(level, dbname, dbpath, force, simulate):
    """Per-clade HMMER DB build (recursive download + post-process).

    Unchanged from v2.x — the per-OG hmm bundles are too many to list
    individually in manifest.json. Verification happens via diamond
    search backend's manifest path; HMMER mode users opt in here.
    """
    if not os.path.exists(dbpath):
        os.makedirs(dbpath)

    EGGNOG_URL = "https://downloads.eggnogdb.org/emapper/eggnog_5.0/per_tax_level"
    baseurl = f"{EGGNOG_URL}/{level}/"
    hmmsurl = f"{baseurl}/{level}_hmms.tar.gz"
    seqsurl = f"{baseurl}/{level}_raw_algs.tar"
    flag = "" if force else "-N"

    cmd = (
        f"cd {dbpath}; "
        f"echo Downloading HMMs... && "
        f"wget {flag} -nH --user-agent=Mozilla/5.0 --relative -r --no-parent "
        f'--reject "index.html*" --cut-dirs=4 -e robots=off {hmmsurl} && '
        f"echo Decompressing HMMs... && "
        f"tar zxf {level}_hmms.tar.gz && "
        f"echo {level}/* | xargs mv -t ./ && rm -r {level} && "
        f"rm {level}_hmms.tar.gz; "
        "numf=$(find ./ | grep -c \".hmm$\"); "
        "curr=0; "
        f"cat /dev/null > {dbname}.hmm_tmp; "
        "for file in $(find ./ | grep \".hmm$\"); do "
        "curr=$((curr+1)); "
        'echo "merging HMMs... ${file} (${curr}/${numf})"; '
        f'cat "${{file}}" | sed -e "s/.faa.final_tree.fa//" '
        f'-e "s/.faa.final_tree//" >> {dbname}.hmm_tmp; '
        'rm "${file}"; done; '
        f"mv {dbname}.hmm_tmp {dbname}.hmm; "
        f"(if [ -f {dbname}.hmm.h3i ]; then rm {dbname}.hmm.h3*; fi) && "
        'echo "hmmpress-ing HMMs... " && '
        f"{HMMPRESS} {dbname}.hmm && "
        'echo "generating idmap file... " && '
        f"cat {dbname}.hmm | grep \"^NAME\" | sed -e \"s/^NAME *//\" "
        f'| awk \'{{print NR" "$0}}\' > {dbname}.hmm.idmap && '
        'echo "removing single OG hmm files... " && '
        "echo ./*hmm | xargs rm; "
    )
    if simulate:
        print(colorify(cmd, "cyan"))
    else:
        subprocess.run(cmd, shell=True, check=True)

    cmd = (
        f"cd {dbpath}; "
        f"echo Downloading FASTAs... && "
        f"wget {flag} -nH --user-agent=Mozilla/5.0 --relative -r --no-parent "
        f'--reject "index.html*" --cut-dirs=4 -e robots=off {seqsurl} && '
        f"echo Decompressing FASTAs... && "
        f"tar xf {level}_raw_algs.tar && "
        f"echo {level}/* | xargs mv -t ./ && rm -r {level} && "
        f"rm {level}_raw_algs.tar; "
        "for file in $(find ./ | grep \".faa.gz$\"); do "
        'outf=$(echo "$file" | sed "s/\\.raw_alg\\.faa\\.gz/\\.fa/"); '
        'zcat "$file" | awk \'/^>/{print; next}{gsub("-", ""); print}\' > "$outf" && '
        'rm "$file"; done'
    )
    if simulate:
        print(colorify(cmd, "cyan"))
    else:
        subprocess.run(cmd, shell=True, check=True)


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
    parser.add_argument("-H", action="store_true", dest="hmmer",
                        help='Install HMMER database for "-d TAXID"')
    parser.add_argument("-d", type=str, dest="hmmer_dbs",
                        help="HMMER tax ID (e.g. -H -d 2 for Bacteria)")
    parser.add_argument("--dbname", type=str, dest="dbname",
                        help="Name for the downloaded HMMER database dir")

    parser.add_argument("-y", action="store_true", dest="allyes",
                        help='Assume "yes" to all questions')
    parser.add_argument("-f", action="store_true", dest="force",
                        help="Force download even if files exist")
    parser.add_argument("-s", action="store_true", dest="simulate",
                        help="Simulate; nothing downloaded")
    parser.add_argument("-q", action="store_true", dest="quiet",
                        help="Quiet mode")

    parser.add_argument("--db", dest="db_backend", metavar="BACKEND", type=str,
                        default=DEFAULT_BACKEND, choices=get_backend_names(),
                        help=("Database backend. Overridden by --data_dir "
                              "and EGGNOG_DATA_DIR."))
    parser.add_argument("--data_dir", metavar="", type=existing_dir,
                        help="Directory to use for DATA_PATH.")

    # New in v3.0
    parser.add_argument("--db_version", type=str, default=__DB_VERSION__,
                        help=("DB version to download. Pins a specific server "
                              "subdirectory (downloads.eggnogdb.org/.../emapperdb-X.Y.Z/)."))
    parser.add_argument("--no-verify", action="store_true", dest="no_verify",
                        help=("Skip manifest+checksum verification. Use only "
                              "with mirrors that lack manifest.json. Will warn."))

    args = parser.parse_args()

    set_data_path(resolve_backend(args.db_backend))
    if "EGGNOG_DATA_DIR" in os.environ:
        set_data_path(os.environ["EGGNOG_DATA_DIR"])
    if args.data_dir:
        set_data_path(args.data_dir)
    data_path = get_data_path()

    base_url = f"https://downloads.eggnogdb.org/emapper/emapperdb-{args.db_version}"
    print(colorify(f"DB version: {args.db_version}  base URL: {base_url}", "lblue"))

    # Fetch manifest first (unless --no-verify)
    manifest = None
    verify = not args.no_verify
    if verify:
        try:
            manifest = fetch_manifest(base_url)
        except DownloadError as e:
            sys.exit(f"failed to fetch manifest.json: {e}")
        if manifest is None:
            sys.exit(
                f"server has no manifest.json at {base_url}. "
                "Pass --no-verify to bypass (you will not be able to "
                "detect corrupt downloads)."
            )
        check_compat(manifest)
        print(colorify(
            f"manifest OK: schema={manifest.get('schema_version')} "
            f"build_date={manifest.get('build_date')} "
            f"min_mapper={manifest.get('min_mapper_version')}",
            "green"))
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
        if args.allyes or ask("Download diamond database (~4 GB after decompression)?") == "y":
            print(colorify(f"Downloading eggnog_proteins.dmnd at {data_path}…", "green"))
            download_diamond_db(base_url, data_path, manifest, verify, args.force)
        else:
            print("Skipping")
    elif not args.quiet:
        print(colorify("Skipping diamond DB (already present or -D). Use -f to force.",
                       "lblue"))

    # ---- PFAM ----
    if args.pfam and (args.force or not pexists(get_pfam_dbpath())):
        if args.allyes or ask("Download Pfam database (~3 GB after decompression)?") == "y":
            print(colorify(f"Downloading Pfam at {data_path}…", "green"))
            download_pfam_db(base_url, data_path, manifest, verify, args.force)
        else:
            print("Skipping")
    elif not args.quiet:
        print(colorify("Skipping Pfam DB (already present). Use -P -f to force.", "lblue"))

    # ---- MMseqs ----
    if args.mmseqs and (args.force or not pexists(get_eggnog_mmseqs_dbpath())):
        if args.allyes or ask("Download MMseqs2 database (~10 GB after decompression)?") == "y":
            print(colorify(f"Downloading MMseqs2 at {data_path}…", "green"))
            download_mmseqs_db(base_url, data_path, manifest, verify, args.force)
        else:
            print("Skipping")
    elif not args.quiet:
        print(colorify("Skipping MMseqs2 DB (already present). Use -M -f to force.", "lblue"))

    # ---- HMMER (recursive wget; no manifest entry) ----
    if args.hmmer:
        if args.allyes or ask(f"Download HMMER database of tax ID {args.hmmer_dbs}?") == "y":
            if args.dbname:
                dbname = args.dbname
            elif args.allyes:
                dbname = args.hmmer_dbs
            else:
                dbname = ask_name(
                    "Specify a non-empty name for the database (e.g. Bacteria)",
                    args.hmmer_dbs,
                )
            dbspath = get_hmmer_base_dbpath(dbname)
            if args.force or not pexists(dbspath):
                print(colorify(
                    f"Downloading HMMER tax ID {args.hmmer_dbs} as {dbname!r} to {dbspath}",
                    "green"))
                print(colorify("This can take a long time for large clades", "red"))
                download_hmm_database(args.hmmer_dbs, dbname, dbspath,
                                      args.force, args.simulate)
            elif not args.quiet:
                print(colorify(
                    f"HMMER database {dbname} already present at {dbspath}. "
                    f"Use -f to force.", "lblue"))
        else:
            print(colorify("Skipping HMMER database", "lblue"))
    else:
        print(colorify(
            "No HMMER database requested. Use '-H -d taxid' to download one.",
            "lblue"))

    print(colorify("Finished.", "green"))


if __name__ == "__main__":
    main()
