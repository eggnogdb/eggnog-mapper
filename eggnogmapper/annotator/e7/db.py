"""Database access for v7+ integer-encoded databases.

v7 databases use integer protein IDs with delta-varint encoded BLOBs.
This module provides bulk query methods for high-throughput annotation.
"""

import array
import hashlib
import logging
import os
import sqlite3
import tempfile
from typing import Dict, Iterable, List, Optional, Set, Tuple
from ..codec import decode_intlist

logger = logging.getLogger(__name__)


class EggnogDBError(Exception):
    """Raised when the eggnog database cannot be opened or queried."""


class EggnogDB:
    """Database access for v7+ eggNOG databases.

    Provides optimized bulk queries for batch annotation:
    - Pre-loaded taxid array for fast protein -> species lookup
    - Bulk event queries
    - Bulk protein annotation queries
    """

    def __init__(self, db_path: str, load_taxids: bool = True):
        """Initialize database connection.

        Args:
            db_path: Path to eggnog.db SQLite file
            load_taxids: If True, load protein_names.taxid into memory array
        """
        self.db_path = db_path
        try:
            self.conn = sqlite3.connect(
                f"file:{db_path}?mode=ro", uri=True, check_same_thread=False
            )
            self.conn.row_factory = sqlite3.Row
            self.conn.execute("PRAGMA mmap_size=2147483648")   # 2 GB mmap window
            self.conn.execute("PRAGMA cache_size=-131072")      # 128 MB page cache

            self._taxids = None
            self._name_cache: Dict[str, int] = {}
            if load_taxids:
                self._load_taxid_array()
        except (sqlite3.OperationalError, sqlite3.DatabaseError) as exc:
            raise EggnogDBError(
                f"Cannot open eggnog database at '{db_path}': {exc}"
            ) from exc

    def reopen_connection(self):
        """Reopen the SQLite connection — for fork-pool workers.

        sqlite3.Connection objects are not safe across `fork()`: shared
        file descriptors / locks can deadlock or corrupt. After fork,
        each worker calls this once to drop the inherited (now-stale)
        connection and open its own. The taxid_array, db_path, and any
        other in-memory state are preserved (inherited from the parent
        via copy-on-write — no reload of the 226 MB array).

        Safe to call from `from_connection`-instantiated DBs only when
        a `db_path` is known. Adopted connections without a path raise.
        """
        if not self.db_path:
            raise EggnogDBError(
                "reopen_connection() requires a db_path; this DB was "
                "created via from_connection() and the original path "
                "is not retained."
            )
        try:
            self.conn.close()
        except Exception:
            pass
        try:
            self.conn = sqlite3.connect(
                f"file:{self.db_path}?mode=ro", uri=True, check_same_thread=False
            )
            self.conn.row_factory = sqlite3.Row
            # Workers do random sub-batch lookups, not full scans.
            # 256 MB mmap + 32 MB page cache is sufficient; keeping these
            # small reduces per-worker virtual and physical memory
            # significantly compared to the parent's 2 GB / 128 MB settings.
            self.conn.execute("PRAGMA mmap_size=268435456")
            self.conn.execute("PRAGMA cache_size=-32768")
        except (sqlite3.OperationalError, sqlite3.DatabaseError) as exc:
            raise EggnogDBError(
                f"Cannot reopen eggnog database at '{self.db_path}': {exc}"
            ) from exc

    @classmethod
    def from_connection(cls, conn, taxid_array=None, db_path=None):
        """Adopt an existing sqlite3 connection (no reopen, no reload).

        Emapper already opens the DB and preloads the taxid array; sharing
        that state avoids duplicating the 226 MB in-memory array.

        Args:
            conn: Pre-opened sqlite3.Connection (row_factory will be set)
            taxid_array: Optional pre-loaded taxid list (index=protein_id).
                If None, it is lazily loaded on first access.
            db_path: Optional file path the connection was opened from.
                When omitted, the path is recovered from the connection
                via ``PRAGMA database_list``. Retaining the path keeps
                ``reopen_connection()`` available for fork-pool workers.
        """
        self = cls.__new__(cls)
        if db_path is None:
            try:
                row = conn.execute("PRAGMA database_list").fetchone()
                # row = (seq, name, file) for sqlite3.Row or tuple
                if row is not None:
                    file_field = row["file"] if hasattr(row, "keys") else row[2]
                    if file_field:
                        db_path = file_field
            except sqlite3.Error:
                pass
        self.db_path = db_path
        self.conn = conn
        self.conn.row_factory = sqlite3.Row
        self._taxids = taxid_array
        self._name_cache: Dict[str, int] = {}
        return self

    @staticmethod
    def _taxid_cache_path(db_path: str) -> str:
        """Return the preferred path for the binary taxid cache file.

        Tries a sibling file alongside the DB first. Falls back to /tmp
        when the DB directory is not writable (e.g. a read-only Docker mount).
        """
        preferred = db_path + ".taxids.bin"
        db_dir = os.path.dirname(db_path) or "."
        if os.access(db_dir, os.W_OK):
            return preferred
        md5 = hashlib.md5(db_path.encode()).hexdigest()[:12]
        return os.path.join(tempfile.gettempdir(), f"eggnog_taxids_{md5}.bin")

    def _load_taxid_array(self):
        """Load taxid array for fast protein ID -> species lookup.

        Creates an array.array('i') where index = protein_id, value = taxid.
        Uses a flat binary cache (~236 MB on disk) so warm restarts load in
        ~2 s instead of ~60 s for the full SQL scan over 59 M rows.
        """
        cursor = self.conn.execute("SELECT MAX(id) FROM protein_names")
        max_id = cursor.fetchone()[0] or 0
        n = max_id + 1

        cache_path = self._taxid_cache_path(self.db_path)
        expected_bytes = n * array.array('i').itemsize

        # Try loading from binary cache.
        if os.path.exists(cache_path):
            try:
                if os.path.getsize(cache_path) == expected_bytes:
                    a = array.array('i')
                    with open(cache_path, 'rb') as fh:
                        a.fromfile(fh, n)
                    self._taxids = a
                    logger.info(
                        "_load_taxid_array: loaded %d taxids from cache %s",
                        n, cache_path,
                    )
                    return
            except Exception as exc:
                logger.warning(
                    "_load_taxid_array: cache read failed (%s); rebuilding from DB", exc
                )

        # Build from SQL: allocate zero-initialised buffer (no Python list).
        logger.info(
            "_load_taxid_array: building taxid array from DB (n=%d) …", n
        )
        buf = bytearray(expected_bytes)
        a = array.array('i', buf)

        cursor = self.conn.execute("SELECT id, taxid FROM protein_names")
        for row in cursor:
            a[row["id"]] = row["taxid"]
        self._taxids = a

        # Write cache for future restarts.
        try:
            with open(cache_path, 'wb') as fh:
                a.tofile(fh)
            logger.info(
                "_load_taxid_array: cache written to %s", cache_path
            )
        except Exception as exc:
            logger.warning(
                "_load_taxid_array: could not write cache to %s: %s", cache_path, exc
            )

    @property
    def taxid_array(self):
        """Array mapping protein ID -> species taxid."""
        return self._taxids

    def get_protein_name(self, protein_id: int) -> Optional[str]:
        """Get protein name for an integer ID."""
        cursor = self.conn.execute(
            "SELECT name FROM protein_names WHERE id = ?",
            (protein_id,)
        )
        row = cursor.fetchone()
        return row["name"] if row else None

    def get_protein_id(self, name: str) -> Optional[int]:
        """Get integer ID for a protein name (cached)."""
        cached = self._name_cache.get(name)
        if cached is not None:
            return cached
        cursor = self.conn.execute(
            "SELECT id FROM protein_names WHERE name = ?",
            (name,)
        )
        row = cursor.fetchone()
        result = row["id"] if row else None
        # Only cache hits; bound to ~100 k entries (~10 MB) to avoid unlimited growth.
        if result is not None and len(self._name_cache) < 100_000:
            self._name_cache[name] = result
        return result

    def get_protein_names_bulk(self, protein_ids: List[int]) -> Dict[int, str]:
        """Get protein names for multiple IDs.

        Args:
            protein_ids: List of integer protein IDs

        Returns:
            Dict mapping protein_id -> name
        """
        if not protein_ids:
            return {}

        placeholders = ",".join("?" * len(protein_ids))
        cursor = self.conn.execute(
            f"SELECT id, name FROM protein_names WHERE id IN ({placeholders})",
            protein_ids
        )
        return {row["id"]: row["name"] for row in cursor}

    def get_events_for_protein(self, protein_id: int) -> List[int]:
        """Get event IDs for a protein from event_index.

        Args:
            protein_id: Integer protein ID

        Returns:
            List of event IDs (integers)
        """
        cursor = self.conn.execute(
            "SELECT events FROM event_index WHERE protein_id = ?",
            (protein_id,)
        )
        row = cursor.fetchone()
        if not row or not row["events"]:
            return []
        return decode_intlist(row["events"])

    def get_event_indices_bulk(self, protein_ids: List[int]) -> Dict[int, List[int]]:
        """Get event IDs for multiple proteins.

        Args:
            protein_ids: List of integer protein IDs

        Returns:
            Dict mapping protein_id -> list of event IDs
        """
        if not protein_ids:
            return {}

        result = {}
        BATCH = 900
        for i in range(0, len(protein_ids), BATCH):
            chunk = protein_ids[i:i + BATCH]
            placeholders = ",".join("?" * len(chunk))
            cursor = self.conn.execute(
                f"SELECT protein_id, events FROM event_index WHERE protein_id IN ({placeholders})",
                chunk
            )
            for row in cursor:
                result[row["protein_id"]] = decode_intlist(row["events"])
        return result

    def get_events_bulk(self, event_ids: List[int]) -> Dict[int, dict]:
        """Get event data for multiple event IDs.

        Args:
            event_ids: List of event IDs

        Returns:
            Dict mapping event_id -> event dict with keys:
            - name: cluster name
            - og: OG name
            - og_lca: OG LCA taxid
            - ev_lca: event LCA taxid
            - sp_overlap: species Jaccard overlap between the two sides
            - side1: set of protein IDs (for fast membership test; order not preserved)
            - side2: set of protein IDs (for fast membership test; order not preserved)
        """
        if not event_ids:
            return {}

        result = {}
        BATCH = 900
        ids = list(event_ids)
        for i in range(0, len(ids), BATCH):
            chunk = ids[i:i + BATCH]
            placeholders = ",".join("?" * len(chunk))
            cursor = self.conn.execute(
                f"SELECT i, name, og, og_lca, ev_lca, sp_overlap, side1, side2 "
                f"FROM sp_events WHERE i IN ({placeholders})",
                chunk
            )

            for row in cursor:
                result[row["i"]] = {
                    "name": row["name"],
                    "og": row["og"],
                    "og_lca": row["og_lca"],
                    "ev_lca": row["ev_lca"],
                    "sp_overlap": row["sp_overlap"],
                    "side1": set(decode_intlist(row["side1"])),
                    "side2": set(decode_intlist(row["side2"])),
                }
        return result

    def get_protein_annotations(self, protein_id: int) -> Optional[dict]:
        """Get functional annotations for a protein.

        Args:
            protein_id: Integer protein ID

        Returns:
            Dict with annotation columns, or None if not found
        """
        cursor = self.conn.execute(
            "SELECT * FROM prots WHERE id = ?",
            (protein_id,)
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def get_protein_annotations_bulk(self, protein_ids: List[int]) -> Dict[int, dict]:
        """Get functional annotations for multiple proteins.

        Args:
            protein_ids: List of integer protein IDs

        Returns:
            Dict mapping protein_id -> annotation dict
        """
        if not protein_ids:
            return {}

        result = {}
        BATCH = 900
        ids = list(protein_ids)
        for i in range(0, len(ids), BATCH):
            chunk = ids[i:i + BATCH]
            placeholders = ",".join("?" * len(chunk))
            cursor = self.conn.execute(
                f"SELECT * FROM prots WHERE id IN ({placeholders})",
                chunk
            )
            for row in cursor:
                result[row["id"]] = dict(row)
        return result

    def get_protein_ogs_bulk(self, protein_ids: List[int]) -> Dict[int, str]:
        """Get OG assignments for multiple proteins.

        Args:
            protein_ids: List of integer protein IDs

        Returns:
            Dict mapping protein_id -> ogs string (comma-separated OG names)
        """
        if not protein_ids:
            return {}

        result = {}
        BATCH = 900
        ids = list(protein_ids)
        for i in range(0, len(ids), BATCH):
            chunk = ids[i:i + BATCH]
            placeholders = ",".join("?" * len(chunk))
            cursor = self.conn.execute(
                f"SELECT id, ogs FROM prots WHERE id IN ({placeholders})",
                chunk
            )
            for row in cursor:
                if row["ogs"]:
                    result[row["id"]] = row["ogs"]
        return result

    def get_og_info(self, og_name: str) -> Optional[dict]:
        """Get OG metadata.

        Args:
            og_name: Full OG name (e.g., "cluster@taxid|clade")

        Returns:
            Dict with OG columns, or None if not found
        """
        cursor = self.conn.execute(
            "SELECT * FROM ogs WHERE name = ?",
            (og_name,)
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def get_og_info_bulk(self, og_pairs: Iterable[Tuple[str, str]]) -> List[Tuple]:
        """Get OG metadata for multiple (name, level) pairs.

        Uses SQLite IN (?) on name (leveraging the (name, level) composite
        index) and filters by level client-side. Chunked at 900 to stay
        below the SQLite parameter limit.

        Args:
            og_pairs: Iterable of (og_name, level) tuples

        Returns:
            Dict mapping (og_name, level) -> OG row dict
        """
        wanted = set(og_pairs)
        if not wanted:
            return {}

        names = sorted({n for n, _ in wanted})
        result = {}
        BATCH = 900
        for i in range(0, len(names), BATCH):
            chunk = names[i:i + BATCH]
            placeholders = ",".join("?" * len(chunk))
            cursor = self.conn.execute(
                f"SELECT * FROM ogs WHERE name IN ({placeholders})",
                chunk,
            )
            for row in cursor:
                key = (row["name"], str(row["level"]))
                if key in wanted:
                    result[key] = dict(row)
        return result

    def get_ref_term(self, term_type: str, term: str) -> Optional[dict]:
        """Get reference term description (KO, GO, etc.).

        Args:
            term_type: Term type (e.g., 'kegg_ko', 'go')
            term: Term ID

        Returns:
            Dict with term info, or None if not found
        """
        cursor = self.conn.execute(
            "SELECT * FROM ref_terms WHERE type = ? AND term = ?",
            (term_type, term)
        )
        row = cursor.fetchone()
        return dict(row) if row else None

    def get_version(self) -> Optional[str]:
        """Get database version string."""
        cursor = self.conn.execute("SELECT version FROM version LIMIT 1")
        row = cursor.fetchone()
        return row["version"] if row else None

    def close(self):
        """Close database connection."""
        if self.conn:
            self.conn.close()
            self.conn = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
