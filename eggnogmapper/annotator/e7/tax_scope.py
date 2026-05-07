"""Tax scope filtering for v7+ integer-encoded databases.

v7 databases use NCBI taxids in OG names, not the v5-style taxonomic
level names.  This module provides a thin ``LineageFilter`` delegate that
forwards all ceiling resolution and species-set queries to a
:class:`~eggnogmapper.annotator.e7.ceiling.TaxScopeCeilingResolver`.

The old ``AUTO_SCOPE_CLADES``, ``parse_tax_scope``, ``set_scope``,
``is_auto``, and ``describe_scope_for_seed`` APIs have been removed.
Use ``TaxScopeCeilingResolver`` for ceiling resolution and the new
``LineageFilter`` for ortholog filtering.
"""

import logging
import sqlite3
from typing import TYPE_CHECKING, Optional

from ..lineage import LineageCache

if TYPE_CHECKING:
    from .ceiling import TaxScopeCeilingResolver

# ---------------------------------------------------------------------------
# Domain / major-clade taxid constants (kept for external consumers)
# ---------------------------------------------------------------------------
BACTERIA: str = "2"
ARCHAEA: str = "2157"
EUKARYOTA: str = "2759"
CELLULAR_ORGANISMS: str = "131567"

METAZOA: str = "33208"
FUNGI: str = "4751"
VIRIDIPLANTAE: str = "33090"

logger = logging.getLogger(__name__)


class LineageFilter:
    """Thin delegate that forwards ceiling / species-set queries to a resolver.

    All annotation logic that used to live here has been moved to
    :class:`~eggnogmapper.annotator.e7.ceiling.TaxScopeCeilingResolver`.
    ``LineageFilter`` now serves purely as a convenience façade used by
    ``AnnotationEngine`` when it needs to:

    1. Resolve the ceiling taxid for a given seed.
    2. Obtain the set of valid species for that ceiling (for the
       per-protein species filter in ``_collect_orthologs``).
    3. Obtain the set of OG-LCA descendants for that ceiling (for the
       strict-OG gate in the same method).

    Attributes:
        lineage_cache: ``LineageCache`` used by the resolver.
        taxid_array: Optional protein-ID → species-taxid mapping (same
            as what ``AnnotationEngine`` holds for fast lookups).
        _ceiling_resolver: Bound ``TaxScopeCeilingResolver`` instance.
    """

    def __init__(
        self,
        lineage_cache: LineageCache,
        ceiling_resolver: "TaxScopeCeilingResolver",
        taxid_array=None,
    ) -> None:
        """Initialise LineageFilter.

        Args:
            lineage_cache: ``LineageCache`` instance for lineage lookups.
            ceiling_resolver: Required
                :class:`~eggnogmapper.annotator.e7.ceiling.TaxScopeCeilingResolver`
                instance.  Must not be ``None``.
            taxid_array: Optional array mapping protein integer ID to
                species taxid integer.  When provided, enables
                :meth:`filter_ortholog_ids`.

        Raises:
            ValueError: When ``ceiling_resolver`` is ``None``.
        """
        if ceiling_resolver is None:
            raise ValueError(
                "LineageFilter requires a TaxScopeCeilingResolver. "
                "Build one via TaxScopeCeilingResolver.build() and pass it here."
            )
        self.lineage_cache = lineage_cache
        self.taxid_array = taxid_array
        self._ceiling_resolver = ceiling_resolver

    # ------------------------------------------------------------------
    # Ceiling delegation
    # ------------------------------------------------------------------

    def get_ceiling_for_seed(self, seed_taxid: str) -> str:
        """Return the resolved ceiling taxid for a seed's species.

        Delegates to :meth:`TaxScopeCeilingResolver.resolve_ceiling`.

        Args:
            seed_taxid: NCBI species taxid of the seed (string).

        Returns:
            Ceiling taxid string or ``PROKARYOTA_SYNTHETIC`` or ``"root"``.
        """
        return self._ceiling_resolver.resolve_ceiling(seed_taxid)

    def ceiling_name(self, ceiling_taxid: str) -> str:
        """Human-readable name for the ``tax_ceiling`` output column.

        Args:
            ceiling_taxid: Taxid string returned by :meth:`get_ceiling_for_seed`.

        Returns:
            Display name string (e.g. ``"Fungi"``, ``"Prokaryota"``).
        """
        return self._ceiling_resolver.ceiling_name(ceiling_taxid)

    # ------------------------------------------------------------------
    # Species / OG-descendant sets (memoized in resolver)
    # ------------------------------------------------------------------

    def get_valid_species_ids(self, ceiling_taxid: str) -> Optional[frozenset]:
        """Return the frozenset of species taxids at-or-below the ceiling.

        Delegates to :meth:`TaxScopeCeilingResolver.get_valid_species`.
        Results are memoized in the resolver.

        Args:
            ceiling_taxid: Taxid string from :meth:`get_ceiling_for_seed`.

        Returns:
            Frozenset of species taxid strings, or ``None`` for no filter.
        """
        return self._ceiling_resolver.get_valid_species(ceiling_taxid)

    def get_scope_og_descendants(self, ceiling_taxid: str) -> Optional[frozenset]:
        """Return the frozenset of taxids at-or-below the ceiling.

        Used as an ``og_lca`` whitelist. Delegates to
        :meth:`TaxScopeCeilingResolver.get_og_descendants`.

        Args:
            ceiling_taxid: Taxid string from :meth:`get_ceiling_for_seed`.

        Returns:
            Frozenset of taxid strings, or ``None`` for no filter.
        """
        return self._ceiling_resolver.get_og_descendants(ceiling_taxid)

    # ------------------------------------------------------------------
    # ev_lca gate (used directly by engine)
    # ------------------------------------------------------------------

    def ev_lca_passes_ceiling(self, ev_lca: str, ceiling_taxid: str) -> bool:
        """Return True iff the event's LCA is at or below the ceiling.

        Delegates to :meth:`TaxScopeCeilingResolver.ev_lca_passes_ceiling`.

        Args:
            ev_lca: Speciation event LCA taxid string.
            ceiling_taxid: Resolved ceiling taxid.

        Returns:
            ``True`` to keep the event; ``False`` to discard.
        """
        return self._ceiling_resolver.ev_lca_passes_ceiling(ev_lca, ceiling_taxid)


# ---------------------------------------------------------------------------
# Standalone utility — not called internally; exposed for external scripts.
# ---------------------------------------------------------------------------


def resolve_name_to_taxid(name: str, taxa_db_path: str) -> Optional[str]:
    """Resolve a taxonomic name to its NCBI taxid.

    Args:
        name: Taxonomic name (e.g. ``"Metazoa"``, ``"Homo sapiens"``).
        taxa_db_path: Path to ``eggnog.taxa.db``.

    Returns:
        Taxid string, or ``None`` if not found.
    """
    try:
        conn = sqlite3.connect(taxa_db_path)
        cursor = conn.execute(
            "SELECT taxid FROM species WHERE spname = ? COLLATE NOCASE LIMIT 1",
            (name,),
        )
        row = cursor.fetchone()
        conn.close()
        return str(row[0]) if row else None
    except sqlite3.OperationalError as exc:
        logger.error(
            "resolve_name_to_taxid failed for %r in %s: %s",
            name,
            taxa_db_path,
            exc,
        )
        return None
