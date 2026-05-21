"""Per-seed ev_lca ceiling resolver for v7+ tax-scope filtering.

Replaces the old ``--tax_scope inner_narrowest`` / ``AUTO_SCOPE_CLADES``
logic with an ``ev_lca ≤ ceiling`` gate applied per speciation event.

Two built-in auto modes are provided:

- ``auto`` (default): Eukaryota ceiling for all eukaryotes;
  Prokaryota (Bacteria ∪ Archaea) for all prokaryotes. No sub-kingdom
  tiers — avoids over-narrow ceilings for fungi, plants, metazoans.
- ``auto-broad``: broad → narrow, Metazoa-first for eukaryotes.

A named clade or numeric taxid may also be supplied as a fixed ceiling,
including ``131567`` (cellular organisms / LUCA) to disable all filtering.

Microsporidia (taxid 6029) are subtracted from the Fungi set in
auto-broad mode.  Prokaryota is a synthetic sentinel (``_prokaryota``)
covering all Bacteria ∪ Archaea − Eukaryota species.
"""

from __future__ import annotations

import logging
from typing import FrozenSet, Optional

from ...emapperException import EmapperException
from ..lineage import LineageCache
from .tax_scope import (
    BACTERIA,
    ARCHAEA,
    EUKARYOTA,
    FUNGI,
    METAZOA,
    VIRIDIPLANTAE,
    resolve_name_to_taxid,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Taxid constants added / extended beyond tax_scope.py
# ---------------------------------------------------------------------------
OPISTHOKONTA: str = "33154"
MICROSPORIDIA: str = "6029"
LUCA: str = "131567"  # cellular organisms — use as ceiling to disable all ev_lca filtering

# Synthetic sentinel used when no real taxid covers "all prokaryotes".
PROKARYOTA_SYNTHETIC: str = "_prokaryota"

# ---------------------------------------------------------------------------
# Priority lists
# ---------------------------------------------------------------------------
# Each list entry is a taxid string (or PROKARYOTA_SYNTHETIC).
# resolve_ceiling() walks the seed's lineage and returns the FIRST entry
# in this list that is present in that lineage.

# auto (default): domain-level only — Eukaryota or Prokaryota.
# No sub-kingdom tiers: a fungal seed gets Eukaryota ceiling, not Fungi;
# a mammalian seed gets Eukaryota, not Metazoa. This maximises ortholog
# recall at the cost of including more distant (but still same-domain)
# functional transfers.
AUTO_NARROW_PRIORITY: list[str] = [
    EUKARYOTA,             # 2759 — all eukaryotes
    PROKARYOTA_SYNTHETIC,  # all bacteria + archaea
]

# auto-broad: broadest eukaryotic match wins first.
AUTO_BROAD_PRIORITY: list[str] = [
    METAZOA,         # 33208
    VIRIDIPLANTAE,   # 33090
    FUNGI,           # 4751  — excl. Microsporidia (subtracted below)
    EUKARYOTA,       # 2759
    PROKARYOTA_SYNTHETIC,
]

# ---------------------------------------------------------------------------
# Human-readable names for the tax_ceiling output column
# ---------------------------------------------------------------------------
CEILING_NAMES: dict[str, str] = {
    BACTERIA:             "Bacteria",
    ARCHAEA:              "Archaea",
    EUKARYOTA:            "Eukaryota",
    METAZOA:              "Metazoa",
    FUNGI:                "Fungi",
    VIRIDIPLANTAE:        "Viridiplantae",
    OPISTHOKONTA:         "Opisthokonta",
    PROKARYOTA_SYNTHETIC: "Prokaryota",
    LUCA:                 "cellular_organisms",
}


class TaxScopeCeilingResolver:
    """Resolves a per-seed ev_lca ceiling for tax-scope filtering.

    The resolver is built once per mapper run via :meth:`build` and then
    called per-seed inside ``AnnotationEngine._annotate_batch_inproc``.

    Attributes:
        mode: ``"auto"``, ``"auto-broad"``, or a fixed taxid/name.
        _priority: Ordered list of ceiling candidates used in auto modes.
        _prokaryota_taxids: Frozenset of all prokaryotic species taxids.
        _fungi_minus_microsporidia: Frozenset of Fungi species taxids
            excluding Microsporidia descendants.
        _lineage_cache: LineageCache used for lineage lookups.
        _fixed_ceiling: Resolved taxid for non-auto modes (str or None).
        _ceiling_cache: Memoized seed-taxid → ceiling-taxid results.
        _valid_species_cache: Memoized ceiling-taxid → valid species set.
        _og_descendants_cache: Memoized ceiling-taxid → descendant set.
    """

    def __init__(
        self,
        lineage_cache: LineageCache,
        mode: str,
        prokaryota_taxids: FrozenSet[str],
        fungi_minus_microsporidia: FrozenSet[str],
        priority: list[str],
        fixed_ceiling: Optional[str],
    ) -> None:
        """Internal constructor — use :meth:`build` instead."""
        self._lineage_cache = lineage_cache
        self.mode = mode
        self._prokaryota_taxids = prokaryota_taxids
        self._fungi_minus_microsporidia = fungi_minus_microsporidia
        self._priority = priority
        self._fixed_ceiling = fixed_ceiling
        self._ceiling_cache: dict[str, str] = {}
        self._valid_species_cache: dict[str, Optional[FrozenSet[str]]] = {}
        self._og_descendants_cache: dict[str, Optional[FrozenSet[str]]] = {}

    # ------------------------------------------------------------------
    # Public properties
    # ------------------------------------------------------------------

    @property
    def lineage_cache(self) -> LineageCache:
        """The ``LineageCache`` instance used by this resolver."""
        return self._lineage_cache

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------

    @classmethod
    def build(
        cls,
        lineage_cache: LineageCache,
        mode: str,
        taxa_db_path: str,
    ) -> "TaxScopeCeilingResolver":
        """Build a fully initialised resolver.

        Pre-computes Microsporidia descendants, prokaryota taxids, and
        the Fungi-minus-Microsporidia set. Resolves named clades via
        the taxa DB. Hard-fails on an unresolvable clade name.

        Args:
            lineage_cache: Loaded ``LineageCache`` instance.
            mode: ``"auto"``, ``"auto-broad"``, or a clade name /
                numeric taxid string for a fixed ceiling.
            taxa_db_path: Path to ``eggnog.taxa.db`` — used only when
                ``mode`` is a named clade (not purely numeric).

        Returns:
            A ready-to-use ``TaxScopeCeilingResolver``.

        Raises:
            EmapperException: When ``mode`` is a named clade that cannot
                be resolved to a taxid in ``eggnog.taxa.db``.
        """
        logger.info("Building TaxScopeCeilingResolver (mode=%r)", mode)

        prokaryota_taxids = cls._compute_prokaryota(lineage_cache)
        fungi_minus_microsporidia = cls._compute_fungi_minus_microsporidia(lineage_cache)

        is_auto = mode in ("auto", "auto-broad")
        priority = (
            AUTO_NARROW_PRIORITY
            if mode == "auto"
            else AUTO_BROAD_PRIORITY
            if mode == "auto-broad"
            else []
        )
        fixed_ceiling: Optional[str] = None

        if not is_auto:
            # Named or numeric fixed ceiling.
            if mode.lstrip("-").isdigit():
                fixed_ceiling = mode
            else:
                resolved = resolve_name_to_taxid(mode, taxa_db_path)
                if resolved is None:
                    raise EmapperException(
                        f"--tax_scope: cannot resolve clade name {mode!r} "
                        f"to a taxid in {taxa_db_path}. "
                        "Use a valid NCBI taxonomic name or a numeric taxid."
                    )
                fixed_ceiling = resolved
            logger.info("Fixed tax_scope ceiling: %s (from %r)", fixed_ceiling, mode)

        return cls(
            lineage_cache=lineage_cache,
            mode=mode,
            prokaryota_taxids=prokaryota_taxids,
            fungi_minus_microsporidia=fungi_minus_microsporidia,
            priority=priority,
            fixed_ceiling=fixed_ceiling,
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def resolve_ceiling(self, seed_taxid: str) -> str:
        """Return the ceiling taxid for a seed's species taxid.

        Memoized per seed taxid. In auto modes, walks the seed's NCBI
        lineage against the priority list; the first match wins. In
        fixed-ceiling mode, returns the pre-resolved taxid immediately.

        Args:
            seed_taxid: NCBI species taxid of the seed protein (string).

        Returns:
            A taxid string or ``PROKARYOTA_SYNTHETIC``.  Returns
            ``PROKARYOTA_SYNTHETIC`` when the seed lineage matches none
            of the eukaryotic tiers (i.e. it is prokaryotic).
            Returns ``"root"`` when the lineage is entirely unknown, so
            callers treat it as "no ceiling" (all events pass).
        """
        if seed_taxid in self._ceiling_cache:
            return self._ceiling_cache[seed_taxid]

        if self._fixed_ceiling is not None:
            self._ceiling_cache[seed_taxid] = self._fixed_ceiling
            return self._fixed_ceiling

        lineage = self._lineage_cache.get(seed_taxid)
        if lineage is None:
            # Unknown seed — no ceiling applied (conservative).
            result = "root"
            self._ceiling_cache[seed_taxid] = result
            return result

        result = self._walk_priority(lineage, seed_taxid)
        self._ceiling_cache[seed_taxid] = result
        return result

    def ceiling_name(self, ceiling_taxid: str) -> str:
        """Return a human-readable name for the ``tax_ceiling`` output column.

        Args:
            ceiling_taxid: Taxid string returned by :meth:`resolve_ceiling`.

        Returns:
            Display name (e.g. ``"Fungi"``, ``"Prokaryota"``, ``"root"``).
        """
        return CEILING_NAMES.get(ceiling_taxid, ceiling_taxid)

    def ev_lca_passes_ceiling(self, ev_lca: str, ceiling_taxid: str) -> bool:
        """Return True iff the event's LCA is at or below the ceiling.

        For the ``PROKARYOTA_SYNTHETIC`` sentinel the test is membership
        in the pre-computed prokaryota set. For all other ceilings we
        check that ``ceiling_taxid`` is on the path from root to
        ``ev_lca`` — i.e. ``ceiling_taxid`` appears in the lineage of
        ``ev_lca``.  This means ``ev_lca == ceiling`` also passes
        (ceiling is ancestor of itself).

        Args:
            ev_lca: The speciation event's LCA taxid string.
            ceiling_taxid: Resolved ceiling from :meth:`resolve_ceiling`.

        Returns:
            ``True`` when the event should be kept; ``False`` to discard.
        """
        if not ev_lca:
            return True  # Missing ev_lca: keep (conservative)
        if ceiling_taxid == "root":
            return True  # Unknown seed: no ceiling
        if ceiling_taxid == PROKARYOTA_SYNTHETIC:
            return ev_lca in self._prokaryota_taxids
        # ev_lca passes if ceiling_taxid is an ancestor of (or equal to) ev_lca.
        lineage = self._lineage_cache.get(ev_lca)
        if lineage is None:
            # ev_lca not in lineage cache — internal node not in species table.
            # Be conservative: keep it.
            return True
        return ceiling_taxid in lineage

    def get_valid_species(self, ceiling_taxid: str) -> Optional[FrozenSet[str]]:
        """Return all species taxids at-or-below the ceiling, memoized.

        Args:
            ceiling_taxid: Taxid string or ``PROKARYOTA_SYNTHETIC``.

        Returns:
            Frozenset of species taxid strings, or ``None`` when
            ``ceiling_taxid`` is ``"root"`` (no filter).
        """
        if ceiling_taxid == "root":
            return None
        if ceiling_taxid in self._valid_species_cache:
            return self._valid_species_cache[ceiling_taxid]

        if ceiling_taxid == PROKARYOTA_SYNTHETIC:
            result: FrozenSet[str] = self._prokaryota_taxids
        elif ceiling_taxid == FUNGI:
            # Fungi without Microsporidia.
            result = self._fungi_minus_microsporidia
        else:
            out: set[str] = set()
            for taxid, lineage in self._lineage_cache.items():
                if ceiling_taxid in lineage:
                    out.add(taxid)
            result = frozenset(out)

        self._valid_species_cache[ceiling_taxid] = result
        return result

    def get_og_descendants(self, ceiling_taxid: str) -> Optional[FrozenSet[str]]:
        """Return all taxids at-or-below the ceiling, memoized.

        Used as an ``og_lca`` whitelist: events whose containing OG LCA
        is not in this set are dropped before orthologs are fetched.

        Implementation mirrors the legacy
        ``LineageFilter.get_scope_og_descendants``: walks each species'
        ordered track (root → leaf), finds the first occurrence of
        ``ceiling_taxid``, and adds every node from that position to the
        leaf.  The union over all in-scope species includes all internal
        clade nodes (Brassicaceae, Magnoliopsida, …) that appear as
        ``og_lca`` values but are not leaf keys in ``LineageCache``.

        Args:
            ceiling_taxid: Taxid string or ``PROKARYOTA_SYNTHETIC``.

        Returns:
            Frozenset of taxid strings (ceiling + all descendants), or
            ``None`` for ``"root"`` (no filter).
        """
        if ceiling_taxid == "root":
            return None
        if ceiling_taxid in self._og_descendants_cache:
            return self._og_descendants_cache[ceiling_taxid]

        if ceiling_taxid == PROKARYOTA_SYNTHETIC:
            # For prokaryotes the "OG ceiling" set is the full prokaryota
            # leaf set plus every internal node reachable from a prokaryotic
            # species' lineage track (these nodes appear as og_lca values).
            # We only walk tracks of prokaryotic species (those lacking
            # EUKARYOTA) so no eukaryotic internal nodes bleed in.
            out_set: set[str] = set(self._prokaryota_taxids)
            for taxid, track in self._lineage_cache.tracks():
                if taxid in self._prokaryota_taxids:
                    out_set.update(track)
            result: FrozenSet[str] = frozenset(out_set)
        else:
            key_set: set[str] = {ceiling_taxid}
            for _, track in self._lineage_cache.tracks():
                try:
                    idx = track.index(ceiling_taxid)
                except ValueError:
                    continue
                key_set.update(track[idx:])
            result = frozenset(key_set)

        self._og_descendants_cache[ceiling_taxid] = result
        return result

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _walk_priority(self, lineage: set[str], seed_taxid: str) -> str:
        """Walk the priority list and return the first matching ceiling.

        Special case: when FUNGI appears in the priority list, the check
        is against ``_fungi_minus_microsporidia`` (the real Fungi species
        set minus Microsporidia descendants), not the raw Fungi taxid
        presence in the lineage.  This avoids assigning a Microsporidia
        seed to the Fungi ceiling.

        Args:
            lineage: The seed's lineage set.
            seed_taxid: Seed species taxid (used for Fungi membership check).

        Returns:
            Matched ceiling taxid string, or ``PROKARYOTA_SYNTHETIC``.
        """
        for candidate in self._priority:
            if candidate == PROKARYOTA_SYNTHETIC:
                # Catch-all for prokaryotes (and unknowns that reach here).
                if seed_taxid in self._prokaryota_taxids:
                    return PROKARYOTA_SYNTHETIC
                continue
            if candidate == FUNGI:
                # Exclude Microsporidia from the Fungi ceiling.
                if seed_taxid in self._fungi_minus_microsporidia:
                    return FUNGI
                continue
            if candidate in lineage:
                return candidate

        # Seed matched none of the priority entries — no ceiling.
        return "root"

    # ------------------------------------------------------------------
    # Static helpers (called once in build())
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_prokaryota(lineage_cache: LineageCache) -> FrozenSet[str]:
        """Compute ``bacteria_descendants ∪ archaea_descendants − eukaryota_descendants``.

        This is the "synthetic Prokaryota" set used for the
        PROKARYOTA_SYNTHETIC sentinel.  Reads species table once.

        Args:
            lineage_cache: Loaded ``LineageCache``.

        Returns:
            Frozenset of prokaryotic species taxid strings.
        """
        result: set[str] = set()
        for taxid, lineage in lineage_cache.items():
            if (BACTERIA in lineage or ARCHAEA in lineage) and EUKARYOTA not in lineage:
                result.add(taxid)
        prok = frozenset(result)
        logger.debug("Prokaryota set: %d species", len(prok))
        return prok

    @staticmethod
    def _compute_fungi_minus_microsporidia(
        lineage_cache: LineageCache,
    ) -> FrozenSet[str]:
        """Compute Fungi species taxids excluding Microsporidia descendants.

        Args:
            lineage_cache: Loaded ``LineageCache``.

        Returns:
            Frozenset of non-Microsporidia fungal species taxid strings.
        """
        result: set[str] = set()
        for taxid, lineage in lineage_cache.items():
            if FUNGI in lineage and MICROSPORIDIA not in lineage:
                result.add(taxid)
        fmm = frozenset(result)
        logger.debug("Fungi-minus-Microsporidia set: %d species", len(fmm))
        return fmm
