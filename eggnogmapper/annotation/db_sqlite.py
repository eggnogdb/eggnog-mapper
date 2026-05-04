"""Thin sqlite3 connection holder for eggnog-mapper output and version display.

In v3 the only supported schema is v7 (integer-encoded). This module
provides:

- The shared singleton `AnnotDB` (open once, kept on a module global so
  output.py and the annotator share one connection).
- A 226 MB `_taxids` array preloaded once for O(1) species lookup; the
  v7 batch shim hands this array straight to
  `eggnogmapper.annotator.e7.EggnogDB` via `from_connection(...)` so
  memory is not duplicated.
- `decode_protein_ids` / `get_protein_name`: integer-id → display-name
  helpers used by `output.py` to write the .annotations TSV.
- `get_db_version`: read for the `--version` banner.

If the opened DB is not v7 (no `event_index` table), `__init__` raises
`EmapperException` immediately — there is no legacy v5 fallback path.
Use eggnog-mapper 2.x for v5 databases.
"""

import sqlite3
import array

from ..common import get_eggnogdb_file, existing_file
from ..emapperException import EmapperException


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

    def __init__(self):
        eggnogdb_file = get_eggnogdb_file()
        existing_file(eggnogdb_file)

        self.conn = sqlite3.connect(eggnogdb_file)

        curs = self.conn.cursor()
        curs.execute("PRAGMA synchronous=OFF;")
        curs.execute("PRAGMA journal_mode=OFF;")
        curs.execute("PRAGMA cache_size=2000;")

        # v7+ databases carry an `event_index` table. v3 requires v7;
        # raise loudly if a legacy v5 build is opened.
        curs.execute(
            "SELECT name FROM sqlite_master "
            "WHERE type='table' AND name='event_index'"
        )
        if curs.fetchone() is None:
            raise EmapperException(
                f"{eggnogdb_file}: not a v7 eggnog.db (missing event_index "
                "table). eggnog-mapper v3 requires a v7+ database. Rebuild "
                "with eggnog-builder build-emapper, or use eggnog-mapper 2.x "
                "for legacy v5 databases."
            )

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
