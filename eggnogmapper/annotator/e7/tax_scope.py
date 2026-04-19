"""Tax scope filtering for v7+ integer-encoded databases.

v7 databases use NCBI taxids in OG names, not the v5-style taxonomic level names.
This module provides lineage-based filtering of orthologs:
- Orthologs are kept only if their species lineage includes the scope taxid(s)
- "auto" mode selects appropriate scope based on the seed ortholog's taxonomy
"""

import sqlite3
from typing import Optional

from ..lineage import LineageCache

# Domain taxids
BACTERIA = "2"
ARCHAEA = "2157"
EUKARYOTA = "2759"
CELLULAR_ORGANISMS = "131567"

# Major eukaryotic clades for auto mode
METAZOA = "33208"
FUNGI = "4751"
VIRIDIPLANTAE = "33090"

# Order matters - check specific before general
AUTO_SCOPE_CLADES = [
    (METAZOA, {METAZOA}),              # Metazoa
    (FUNGI, {FUNGI}),                  # Fungi
    (VIRIDIPLANTAE, {VIRIDIPLANTAE}),  # Viridiplantae
]


class LineageFilter:
    """Filters orthologs by taxonomic scope using lineage information.

    Wraps a LineageCache and provides ortholog filtering methods.
    """

    def __init__(self, lineage_cache, taxid_array=None):
        """Initialize LineageFilter.

        Args:
            lineage_cache: LineageCache instance for lineage lookups
            taxid_array: Optional array mapping protein ID -> species taxid.
                         If provided, enables filter_ortholog_ids().
        """
        self.lineage_cache = lineage_cache
        self.taxid_array = taxid_array
        self._scope_taxids = None
        self._is_auto = False

    def set_scope(self, tax_scope):
        """Set the taxonomic scope for filtering.

        Args:
            tax_scope: Scope specification:
                - None or "none": no filtering
                - "auto": automatic scope per query based on seed taxonomy
                - taxid string: single taxid
                - comma-separated taxids: multiple taxids
        """
        self._scope_taxids, self._is_auto = parse_tax_scope(tax_scope)

    @property
    def is_auto(self):
        """True if using auto mode (scope depends on each seed ortholog)."""
        return self._is_auto

    @property
    def scope_taxids(self):
        """Current scope taxids, or None for no filter."""
        return self._scope_taxids

    def get_auto_scope_taxids(self, seed_taxid):
        """Determine appropriate tax scope taxids for auto mode.

        Returns a set of taxids - orthologs must have at least one in their lineage.

        Strategy:
        - Prokaryotes (Bacteria/Archaea): Allow both domains due to common HGT
        - Metazoa: Restrict to Metazoa only
        - Viridiplantae: Restrict to Viridiplantae only
        - Fungi: Restrict to Fungi only
        - Other Eukaryotes: Allow all Eukaryota
        """
        lineage = self.lineage_cache.get(seed_taxid)
        if lineage is None:
            return None  # Unknown taxid, no filter

        # Prokaryotes: allow both Bacteria and Archaea (HGT is common)
        if BACTERIA in lineage or ARCHAEA in lineage:
            return {BACTERIA, ARCHAEA}

        # Eukaryotes: check major clades in order of specificity
        for clade_taxid, scope_set in AUTO_SCOPE_CLADES:
            if clade_taxid in lineage:
                return scope_set

        # Other Eukaryotes
        if EUKARYOTA in lineage:
            return {EUKARYOTA}

        # Unknown domain - no filter
        return None

    def get_effective_scope(self, seed_taxid=None):
        """Get effective scope taxids for filtering.

        In auto mode, computes scope based on seed_taxid.
        Otherwise returns the configured scope.

        Args:
            seed_taxid: Seed ortholog taxid (required for auto mode)

        Returns:
            Set of taxid strings, or None for no filter
        """
        if self._is_auto:
            if seed_taxid is None:
                return None  # Can't auto-scope without seed
            return self.get_auto_scope_taxids(seed_taxid)
        return self._scope_taxids

    def filter_ortholog_ids(self, ortholog_ids, seed_taxid=None):
        """Filter orthologs to those matching the taxonomic scope.

        Args:
            ortholog_ids: List of integer protein IDs
            seed_taxid: Seed ortholog taxid (required for auto mode)

        Returns:
            Filtered list of ortholog IDs

        Raises:
            ValueError: If taxid_array was not provided at init
        """
        if self.taxid_array is None:
            raise ValueError(
                "taxid_array required for filter_ortholog_ids(). "
                "Provide it when creating LineageFilter."
            )

        scope_taxids = self.get_effective_scope(seed_taxid)
        if not scope_taxids:
            return ortholog_ids

        filtered = []
        for oid in ortholog_ids:
            species_taxid = str(self.taxid_array[oid])
            lineage = self.lineage_cache.get(species_taxid)
            if lineage and (lineage & scope_taxids):  # intersection
                filtered.append(oid)
        return filtered

    def get_valid_species_ids(self, scope_taxids) -> Optional[frozenset]:
        """Return the frozenset of species taxids whose lineage intersects `scope_taxids`.

        Results are memoized per frozenset of scope taxids — useful when
        filtering many orthologs against a small number of distinct scopes.
        Returns a frozenset so callers cannot accidentally corrupt the cache
        via mutation.

        Args:
            scope_taxids: Iterable of taxid strings (e.g., {"33090"} for Viridiplantae)

        Returns:
            Frozenset of species taxid strings that qualify under the scope,
            or None if scope_taxids is empty.
        """
        if not scope_taxids:
            return None
        key = frozenset(scope_taxids)
        cache = getattr(self, "_valid_species_cache", None)
        if cache is None:
            cache = self._valid_species_cache = {}
        if key in cache:
            return cache[key]
        out = set()
        for taxid, lineage in self.lineage_cache.items():
            if lineage & key:
                out.add(taxid)
        result = frozenset(out)
        cache[key] = result
        return result

    def filter_by_taxids(self, species_taxids, seed_taxid=None):
        """Filter species taxids to those matching the taxonomic scope.

        Args:
            species_taxids: Iterable of species taxid strings
            seed_taxid: Seed ortholog taxid (required for auto mode)

        Returns:
            List of taxids that match the scope
        """
        scope_taxids = self.get_effective_scope(seed_taxid)
        if not scope_taxids:
            return list(species_taxids)

        filtered = []
        for taxid in species_taxids:
            lineage = self.lineage_cache.get(taxid)
            if lineage and (lineage & scope_taxids):
                filtered.append(taxid)
        return filtered


def parse_tax_scope(tax_scope):
    """Parse tax_scope argument for v7 databases.

    Returns:
        (scope_taxids, is_auto) where:
        - scope_taxids: set of taxid strings to filter by, or None for no filter
        - is_auto: True if auto mode (scope depends on each seed ortholog)

    Accepts:
        - None or "none": no filtering
        - "auto": automatic scope per query based on seed taxonomy
        - taxid or comma-separated taxids: filter to these clades
    """
    if tax_scope is None or (isinstance(tax_scope, str) and tax_scope.lower() == "none"):
        return None, False

    if isinstance(tax_scope, str) and tax_scope.lower() == "auto":
        return None, True  # Will be determined per-query

    # Parse as taxid(s)
    if isinstance(tax_scope, str):
        parts = [p.strip() for p in tax_scope.split(",")]
    else:
        parts = [str(tax_scope)]

    scope_taxids = set()
    for part in parts:
        if part.isdigit():
            scope_taxids.add(part)
        # Note: taxonomic name resolution requires taxa.db access
        # For now, only numeric taxids are supported

    return scope_taxids if scope_taxids else None, False


def resolve_name_to_taxid(name, taxa_db_path):
    """Resolve a taxonomic name to its taxid.

    Args:
        name: Taxonomic name (e.g., "Metazoa", "Homo sapiens")
        taxa_db_path: Path to eggnog.taxa.db

    Returns:
        Taxid string, or None if not found
    """
    try:
        conn = sqlite3.connect(taxa_db_path)
        cursor = conn.execute(
            "SELECT taxid FROM species WHERE spname = ? COLLATE NOCASE LIMIT 1",
            (name,)
        )
        row = cursor.fetchone()
        conn.close()
        return str(row[0]) if row else None
    except Exception:
        return None
