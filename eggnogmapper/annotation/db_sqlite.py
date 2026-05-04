"""Thin sqlite3 connection holder for eggnog-mapper output and version display.

In v3 (Phase 2 commit 5) this module shrank to only what the live
v7-batch path needs:

- The shared singleton `AnnotDB` (open once, kept on a module global so
  output.py and the annotator share one connection).
- A 226 MB `_taxids` array preloaded once for O(1) species lookup; the
  v7 batch shim hands this array straight to `eggnogmapper.annotator.e7.EggnogDB`
  via `from_connection(...)` so memory is not duplicated.
- `decode_protein_ids` / `get_protein_name`: integer-id → display-name
  helpers used by `output.py` to write the .annotations TSV.
- `get_db_version`: read for the `--version` banner.

All v5-style query helpers (`get_member_ogs`, `get_ogs_description`,
`get_annotations`, `get_pfam_annotations`, `get_taxid`, `get_protein_id`,
`bulk_get_protein_ids`, `get_member_events`, the experimental
`bulk_get_*` family) and the `--dbmem` (`usemem=True`) branch were
removed. All annotation logic now lives in `eggnogmapper.annotator.e7`.
"""

import sqlite3
import array

from ..common import get_eggnogdb_file, existing_file


db = None


def get_eggnog_db():
    """Return the process-wide singleton AnnotDB."""
    global db
    if db is None or db.conn is None:
        db = AnnotDB()
    return db


def get_fresh_eggnog_db():
    """Return a brand-new AnnotDB (used by `get_db_version` and tests)."""
    return AnnotDB()


class AnnotDB:
    conn = None
    _int_mode = False  # True for v7+ integer-encoded databases

    def __init__(self):
        eggnogdb_file = get_eggnogdb_file()
        existing_file(eggnogdb_file)

        self.conn = sqlite3.connect(eggnogdb_file)

        curs = self.conn.cursor()
        curs.execute("PRAGMA synchronous=OFF;")
        curs.execute("PRAGMA journal_mode=OFF;")
        curs.execute("PRAGMA cache_size=2000;")

        # v7+ databases carry an `event_index` table; legacy v5 builds do not.
        # The v3 mapper requires v7 — `Annotator._annotate` raises if False.
        curs.execute(
            "SELECT name FROM sqlite_master "
            "WHERE type='table' AND name='event_index'"
        )
        self._int_mode = curs.fetchone() is not None

        if self._int_mode:
            from .codec import decode_intlist
            self._decode_intlist = decode_intlist
            curs.execute("SELECT MAX(id) FROM protein_names")
            max_id = curs.fetchone()[0] or 0
            self._taxids = array.array("i", [0]) * (max_id + 1)
            curs.execute("SELECT id, taxid FROM protein_names")
            for pid, taxid in curs:
                self._taxids[pid] = taxid
            print(
                f"  {max_id + 1} protein taxids cached "
                f"({len(self._taxids) * 4 // 1024 // 1024} MB)",
                flush=True,
            )
        else:
            self._taxids = None

        curs.close()

    def close(self):
        if self.conn is not None:
            self.conn.close()
            self.conn = None

    def get_db_version(self):
        cur = self.conn.cursor()
        cur.execute("SELECT * FROM version LIMIT 1")
        return cur.fetchone()[0]

    def decode_protein_ids(self, int_ids):
        """Bulk decode integer IDs to display names. Returns {id: name}.
        Used only by `output.py` when writing the .annotations TSV."""
        if not int_ids:
            return {}
        id_to_name = {}
        cur = self.conn.cursor()
        BATCH = 900
        ids = list(int_ids)
        for i in range(0, len(ids), BATCH):
            chunk = ids[i:i + BATCH]
            ph = ",".join(["?"] * len(chunk))
            cur.execute(
                f"SELECT id, name FROM protein_names WHERE id IN ({ph})",
                chunk,
            )
            id_to_name.update(cur.fetchall())
        return id_to_name

    def get_protein_name(self, protein_id):
        """Decode a single integer protein ID to its name (str)."""
        cur = self.conn.cursor()
        cur.execute(
            "SELECT name FROM protein_names WHERE id = ?",
            (int(protein_id),),
        )
        row = cur.fetchone()
        return row[0] if row else str(protein_id)
