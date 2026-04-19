"""Lineage cache for taxonomy lookups.

Provides efficient access to species lineages from taxa.db or traverse.pkl.
Used by tax_scope filtering to determine if species belong to a taxonomic clade.
"""

import os
import pickle
import sqlite3


class LineageCache:
    """Cache for species lineages.

    Maps species taxid (str) -> set of ancestor taxids (str).
    Supports loading from:
    - taxa.db SQLite (species.track column)
    - traverse.pkl (pre-computed pickle)

    Loaded lazily on first access.
    """

    def __init__(self, taxa_db_path=None, traverse_pkl_path=None):
        """Initialize LineageCache.

        Args:
            taxa_db_path: Path to eggnog.taxa.db SQLite file
            traverse_pkl_path: Path to traverse.pkl file (faster if available)

        At least one path must be provided, or set via environment:
        - EGGNOG_DATA_DIR: Base directory containing taxa.db
        """
        self._cache = None
        self._taxa_db_path = taxa_db_path
        self._traverse_pkl_path = traverse_pkl_path

    def _load(self):
        """Load lineage data from file."""
        if self._cache is not None:
            return

        # Load from taxa.db (traverse.pkl is for prepostorder, not lineage)
        db_path = self._taxa_db_path
        if not db_path:
            # Try to find taxa.db from environment
            data_dir = os.environ.get("EGGNOG_DATA_DIR")
            if data_dir:
                db_path = os.path.join(data_dir, "eggnog.taxa.db")

        if db_path and os.path.exists(db_path):
            self._load_from_sqlite(db_path)
            return

        raise FileNotFoundError(
            "No taxa database found. Provide taxa_db_path, traverse_pkl_path, "
            "or set EGGNOG_DATA_DIR environment variable."
        )

    def _load_from_pickle(self, path):
        """Load lineage data from traverse.pkl."""
        with open(path, "rb") as f:
            data = pickle.load(f)

        self._cache = {}
        # traverse.pkl format: {taxid: [lineage_taxids]}
        for taxid, lineage in data.items():
            if lineage:
                self._cache[str(taxid)] = set(str(t) for t in lineage)

    def _load_from_sqlite(self, db_path):
        """Load lineage data from taxa.db SQLite."""
        self._cache = {}

        conn = sqlite3.connect(db_path)
        try:
            cursor = conn.execute(
                "SELECT taxid, track FROM species WHERE track IS NOT NULL"
            )
            for taxid, track in cursor:
                if track:
                    # track is comma-separated lineage: "1,131567,2759,..."
                    lineage_set = set(track.split(","))
                    self._cache[str(taxid)] = lineage_set
        finally:
            conn.close()

    def get(self, taxid, default=None):
        """Get lineage set for a species taxid.

        Args:
            taxid: Species taxid (int or str)
            default: Value to return if taxid not found

        Returns:
            set of ancestor taxid strings, or default if not found
        """
        self._load()
        return self._cache.get(str(taxid), default)

    def __contains__(self, taxid):
        """Check if taxid is in cache."""
        self._load()
        return str(taxid) in self._cache

    def __len__(self):
        """Return number of cached species."""
        self._load()
        return len(self._cache)

    def items(self):
        """Iterate (taxid_str, lineage_set) pairs for all cached species."""
        self._load()
        return self._cache.items()

    def is_descendant_of(self, species_taxid, clade_taxid):
        """Check if species is a descendant of a clade.

        Args:
            species_taxid: Species taxid to check
            clade_taxid: Clade taxid to check membership in

        Returns:
            True if species lineage contains clade_taxid
        """
        lineage = self.get(species_taxid)
        if lineage is None:
            return False
        return str(clade_taxid) in lineage

    def filter_by_clade(self, taxids, clade_taxid):
        """Filter taxids to those belonging to a clade.

        Args:
            taxids: Iterable of taxids to filter
            clade_taxid: Clade taxid to filter by

        Returns:
            List of taxids that are descendants of clade_taxid
        """
        clade_str = str(clade_taxid)
        result = []
        for taxid in taxids:
            lineage = self.get(taxid)
            if lineage and clade_str in lineage:
                result.append(taxid)
        return result
