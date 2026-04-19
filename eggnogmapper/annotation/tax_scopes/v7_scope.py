"""Tax scope filtering for v7+ integer-encoded databases.

This module re-exports shared tax_scope logic from eggnog_annotator
and provides compatibility wrappers for the eggnog-mapper API.

v7 databases use NCBI taxids in OG names, not the v5-style taxonomic level names.
This module provides lineage-based filtering of orthologs:
- Orthologs are kept only if their species lineage includes the scope taxid(s)
- "auto" mode selects appropriate scope based on the seed ortholog's taxonomy
"""

# Re-export shared modules
from eggnog_annotator.lineage import LineageCache
from eggnog_annotator.e7.tax_scope import (
    parse_tax_scope,
    LineageFilter,
    resolve_name_to_taxid,
    # Constants
    BACTERIA,
    ARCHAEA,
    EUKARYOTA,
    CELLULAR_ORGANISMS,
    METAZOA,
    FUNGI,
    VIRIDIPLANTAE,
    AUTO_SCOPE_CLADES,
)

__all__ = [
    # Classes
    "LineageCache",
    "LineageFilter",
    # Functions
    "parse_v7_tax_scope",
    "get_auto_scope_taxids",
    "filter_orthologs_by_lineage",
    "resolve_name_to_taxid",
    # Constants
    "BACTERIA",
    "ARCHAEA",
    "EUKARYOTA",
    "CELLULAR_ORGANISMS",
    "METAZOA",
    "FUNGI",
    "VIRIDIPLANTAE",
    "AUTO_SCOPE_CLADES",
]


# Compatibility wrappers for eggnog-mapper's function-based API


def parse_v7_tax_scope(tax_scope, lineage_cache):
    """Parse tax_scope argument for v7 databases.

    This is a compatibility wrapper for eggnog_annotator.e7.tax_scope.parse_tax_scope.
    The eggnog-mapper API passes lineage_cache, but it's not used by parse_tax_scope.

    Returns:
        (scope_taxids, is_auto) where:
        - scope_taxids: set of taxid strings to filter by, or None for no filter
        - is_auto: True if auto mode (scope depends on each seed ortholog)
    """
    # lineage_cache parameter kept for API compatibility but not used
    return parse_tax_scope(tax_scope)


def get_auto_scope_taxids(seed_taxid, lineage_cache):
    """Determine appropriate tax scope taxids for auto mode.

    This is a compatibility wrapper that uses LineageFilter internally.

    Args:
        seed_taxid: Taxid of the seed ortholog
        lineage_cache: LineageCache instance

    Returns:
        Set of taxid strings, or None for no filter
    """
    filter_obj = LineageFilter(lineage_cache)
    return filter_obj.get_auto_scope_taxids(seed_taxid)


def filter_orthologs_by_lineage(ortholog_ids, scope_taxids, taxid_array, lineage_cache):
    """Filter orthologs to those whose species lineage includes at least one scope taxid.

    This is a compatibility wrapper that uses LineageFilter internally.

    Args:
        ortholog_ids: list of integer protein IDs
        scope_taxids: set of taxid strings to filter by (None = no filter)
        taxid_array: array mapping protein ID -> species taxid
        lineage_cache: LineageCache instance

    Returns:
        Filtered list of ortholog IDs
    """
    if not scope_taxids:
        return ortholog_ids

    filter_obj = LineageFilter(lineage_cache, taxid_array=taxid_array)
    # Set scope manually (not using set_scope since we have pre-parsed taxids)
    filter_obj._scope_taxids = scope_taxids
    filter_obj._is_auto = False
    return filter_obj.filter_ortholog_ids(ortholog_ids)
