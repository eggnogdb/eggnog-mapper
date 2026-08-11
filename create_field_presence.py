#!/usr/bin/env python3
"""Build the shipped field-presence mask (``eggnog.db.fieldpresence.bin``).

Maintainer / release tool. Produces the fixed-name lazy-cascade field-presence
mask that ships in the data bundle alongside ``eggnog.db``, so end users do not
pay the one-time build (a full scan of ``prots``) on their first annotation run.

The mask's ``gos_mf/bp/cc`` bits depend on the GO namespace map, so it is built
with the same ``go-basic.obo`` that resolves for the data dir (see
``resolve_go_obo_path``) — ship that OBO with the DB.

Run once per released DB (and per data dir the self-test ships), e.g.:

    python create_field_presence.py --data_dir data/
"""
from __future__ import annotations

import argparse
import os
import sys

from eggnogmapper.common import set_data_path, get_data_path, get_eggnogdb_file
from eggnogmapper.annotation.batch_annotate import resolve_go_obo_path
from eggnogmapper.annotator.e7.db import EggnogDB
from eggnogmapper.annotator.e7.annotate import AnnotationEngine


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data_dir", required=True,
                    help="Data dir containing eggnog.db (and go-basic.obo).")
    args = ap.parse_args()

    set_data_path(args.data_dir)
    db_path = get_eggnogdb_file()
    if not os.path.exists(db_path):
        sys.exit(f"eggnog.db not found at {db_path}")
    obo = resolve_go_obo_path()
    out = f"{db_path}.fieldpresence.bin"
    print(f"DB : {db_path}")
    print(f"OBO: {obo}  ({'found' if obo and os.path.exists(obo) else 'MISSING — GO bits would use legacy fallback'})")
    print(f"->  {out}")

    # taxid array not needed to build the mask (only a scan of `prots`).
    db = EggnogDB(db_path, load_taxids=False)
    engine = AnnotationEngine(db, go_obo_path=obo, lazy_cascade=True)
    ok = engine.load_field_presence(cache_path=out)
    if not ok or not os.path.exists(out):
        sys.exit("failed to build the field-presence mask")
    print(f"wrote {out} ({os.path.getsize(out):,} bytes)")


if __name__ == "__main__":
    main()
