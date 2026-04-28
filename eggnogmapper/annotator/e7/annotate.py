"""Batch annotation logic for v7+ databases.

High-throughput annotation using bulk queries and simplified ortholog collection.
Eliminates species grouping bottleneck for ~100 q/s performance.

Usage:
    from eggnog_annotator.e7.annotate import AnnotationEngine

    engine = AnnotationEngine(db_path)
    result = engine.annotate(seed_ortholog_id)
    # result = {
    #     "orthologs": [id1, id2, ...],
    #     "annotations": {"GOs": [...], "KEGG_ko": [...], ...},
    #     "og_info": {"name": ..., "description": ...},
    # }
"""

import logging
import os
import time
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Set, Tuple, Any

from .db import EggnogDB
from ..codec import decode_intlist

logger = logging.getLogger(__name__)

# Enable per-phase timing by setting EGGNOG_ANNOTATOR_PROFILE=1
_PROFILE = os.environ.get("EGGNOG_ANNOTATOR_PROFILE", "0") == "1"

# C++/Cython acceleration of the inner ortholog-meta loop. Falls back
# to the pure-Python path if the .so didn't compile (mirrors the codec
# fallback pattern). The flag is module-level so each `_collect_orthologs`
# call doesn't pay an import-time check.
try:
    from .._collect_inner import OrthologCollector as _OrthologCollector
    from .._collect_inner import pack_sort_key as _pack_sort_key
    _COLLECTOR_AVAILABLE = True
except ImportError:  # pragma: no cover — falls back transparently
    _OrthologCollector = None
    _pack_sort_key = None
    _COLLECTOR_AVAILABLE = False


class AnnotationEngine:
    """Unified annotation engine for v7+ eggNOG databases.

    Used by both eggnog-mapper and eggnog-website for consistent annotation logic.
    Optimized for high throughput with bulk queries and simplified ortholog collection.
    """

    # Annotation fields to aggregate from orthologs
    ANNOTATION_FIELDS = [
        "pname", "gos", "kegg_ko", "kegg_ec", "kegg_pathway",
        "kegg_module", "kegg_reaction", "kegg_rclass", "kegg_brite",
        "kegg_tc", "kegg_cazy", "bigg_reaction", "pfam"
    ]

    # Mapping from ortholog-relationship type to cascade tier. The cascade
    # walks tiers in order 0 → 2; within a tier, donors are sorted by
    # ev_lca proximity to the seed. Confidence label is derived from the
    # winning tier: 0 → "high", 1 → "medium", 2 → "low".
    TYPE_TIERS = {
        "one2one":  0,
        "one2many": 1,
        "many2one": 1,
        "many2many": 2,
    }
    TIER_CONFIDENCE = {0: "high", 1: "medium", 2: "low"}

    def __init__(self, db: EggnogDB, lineage_filter=None, lineage_cache=None):
        """Initialize annotation engine.

        Args:
            db: EggnogDB instance with loaded taxid array
            lineage_filter: Optional LineageFilter for tax-scope filtering.
                When set (typically to auto mode), orthologs are pruned
                per-seed to the seed's taxonomic domain before expensive
                annotation fetch — the single largest speedup for
                proteome-scale runs.
            lineage_cache: Optional LineageCache used by the cascade engine
                to rank events by ev_lca depth. If None and lineage_filter
                is set, the filter's own cache is reused. If both are
                absent, the cascade falls back to type-tier priority only
                (no ev_lca distance ordering).
        """
        self.db = db
        self.lineage_filter = lineage_filter
        self.lineage_cache = lineage_cache or (
            lineage_filter.lineage_cache if lineage_filter is not None else None
        )
        # Cross-batch cache for OG metadata. Most eukaryotic proteomes
        # reuse a small set of common OGs across thousands of proteins,
        # so caching them avoids repeated SQL + JSON parsing.
        self._og_cache: Dict[Tuple[str, str], Optional[dict]] = {}
        # Cache of valid-species sets keyed by the seed's species taxid.
        # Reused across seeds that share the same auto-scope outcome
        # (e.g. an all-Arabidopsis proteome => one scope computed once).
        self._valid_species_by_seed: Dict[str, Optional[set]] = {}
        # Cache of seed → set(lineage_taxid_str). Cuts repeated taxa.db
        # work for proteomes where many seeds share a species.
        self._seed_lineage_set_cache: Dict[int, Optional[frozenset]] = {}

    def annotate(
        self,
        seed_id: int,
        target_taxa: Optional[Set[int]] = None,
        excluded_taxa: Optional[Set[int]] = None,
    ) -> Dict[str, Any]:
        """Annotate a single seed ortholog.

        Args:
            seed_id: Integer protein ID of seed ortholog
            target_taxa: Optional set of taxids to include (filter orthologs)
            excluded_taxa: Optional set of taxids to exclude

        Returns:
            Dict with keys:
            - orthologs: list of ortholog protein IDs
            - annotations: dict of annotation field -> values
            - og_info: dict with OG name, description, category
        """
        # Get events for this protein
        event_ids = self.db.get_events_for_protein(seed_id)
        if not event_ids:
            return {
                "orthologs": [],
                "annotations": {},
                "annotations_confidence": {},
                "tax_scope_used": self._describe_seed_scope(seed_id),
                "og_info": None,
                "ortholog_types": {
                    "one2one": set(), "one2many": set(),
                    "many2one": set(), "many2many": set(), "all": set(),
                },
                "all_ogs": [],
            }

        # Fetch events
        events = self.db.get_events_bulk(event_ids)

        # Collect orthologs from events (simplified: no species grouping)
        orthologs, ortholog_types, ortholog_meta = self._collect_orthologs(
            seed_id, events, target_taxa, excluded_taxa
        )

        # Fetch and summarize annotations through the cascade
        annotations: Dict[str, Any] = {}
        annotations_confidence: Dict[str, str] = {}
        if orthologs:
            annot_data = self.db.get_protein_annotations_bulk(list(orthologs))
            annotations, annotations_confidence = self._summarize_annotations(
                annot_data, ortholog_meta
            )

        # Preferred_name: seed's own pname takes priority when it is an informative
        # gene name. Locus IDs and multi-alias entries fall through to consensus.
        seed_row = self.db.get_protein_annotations_bulk([seed_id]).get(seed_id)
        if seed_row:
            seed_pname = (seed_row.get("pname") or "").strip()
            if self._is_informative_pname(seed_pname):
                annotations["Preferred_name"] = [seed_pname]
                annotations_confidence["Preferred_name"] = "high"

        # Get OG info from events
        og_info = self._get_og_info(events)

        return {
            "orthologs": list(orthologs),
            "annotations": annotations,
            "annotations_confidence": annotations_confidence,
            "tax_scope_used": self._describe_seed_scope(seed_id),
            "og_info": og_info,
            "ortholog_types": ortholog_types,
            # Single-protein path has no prots.ogs lookup; batch path fills this.
            "all_ogs": [],
        }

    def annotate_batch(
        self,
        seed_ids: List[int],
        target_taxa: Optional[Set[int]] = None,
        excluded_taxa: Optional[Set[int]] = None,
        target_orthologs: str = "all",
        scope_strict_og: bool = False,
    ) -> Dict[int, Dict[str, Any]]:
        """Annotate multiple seed orthologs efficiently.

        Uses bulk queries to minimize database round-trips.

        Args:
            seed_ids: List of integer protein IDs
            target_taxa: Optional set of taxids to include
            excluded_taxa: Optional set of taxids to exclude
            target_orthologs: Cascade type floor; see ``TARGET_ORTHOLOGS_FLOORS``.
            scope_strict_og: When True, drop events whose containing OG
                (``sp_events.og_lca``) is broader than the seed's resolved
                tax-scope ceiling. Default False preserves the legacy
                behaviour where above-scope OG events stay in (and are
                pruned per-protein by the species filter only). The
                strict mode is materially faster on auto-scope plant /
                animal proteomes — most cross-kingdom paralog noise
                lives in cellular-organisms / Eukaryota OGs and the same
                in-scope orthologs already appear in lower OGs at higher
                cascade priority.

        Returns:
            Dict mapping seed_id -> annotation result
        """
        if not seed_ids:
            return {}

        def _t(label, t0):
            if _PROFILE:
                print(f"  [annotate_batch] {label}: {time.time() - t0:.2f}s",
                      flush=True)

        # Phase 1: Bulk fetch event indices
        t0 = time.time()
        event_index = self.db.get_event_indices_bulk(seed_ids)
        _t("p1 event_indices", t0)

        # Collect all event IDs
        all_event_ids = set()
        for eids in event_index.values():
            all_event_ids.update(eids)

        # Phase 2: Bulk fetch events
        t0 = time.time()
        events_cache = self.db.get_events_bulk(list(all_event_ids))
        _t(f"p2 events ({len(all_event_ids)})", t0)

        # Phase 3: Collect orthologs for each seed (with tax-scope prune)
        t0 = time.time()
        seed_orthologs: Dict[int, Set[int]] = {}
        seed_ortholog_types: Dict[int, Dict[str, Set[int]]] = {}
        seed_ortholog_meta: Dict[int, Dict[int, Dict[str, Any]]] = {}
        seed_tax_scope_used: Dict[int, str] = {}
        all_orthologs: Set[int] = set()
        for seed_id in seed_ids:
            eids = event_index.get(seed_id, [])
            # Filter events to those in our cache
            seed_events = {eid: events_cache[eid] for eid in eids if eid in events_cache}
            valid_species = self._resolve_valid_species(seed_id)
            allowed_og_lcas = (
                self._resolve_allowed_og_lcas(seed_id) if scope_strict_og else None
            )
            orthologs, ortholog_types, ortholog_meta = self._collect_orthologs(
                seed_id, seed_events, target_taxa, excluded_taxa, valid_species,
                allowed_og_lcas=allowed_og_lcas,
            )
            seed_orthologs[seed_id] = orthologs
            seed_ortholog_types[seed_id] = ortholog_types
            seed_ortholog_meta[seed_id] = ortholog_meta
            seed_tax_scope_used[seed_id] = self._describe_seed_scope(seed_id)
            all_orthologs.update(orthologs)
        _t(f"p3 collect_orthologs ({len(all_orthologs)})", t0)

        # Phase 4: Bulk fetch annotations for all orthologs and seeds.
        # Seeds are fetched so their pname can be used as the primary Preferred_name.
        # Pre-parse once per ortholog: the cascade walk is O(seeds × buckets ×
        # fields), so re-splitting comma-strings inside that loop dominated
        # phase 6 on plant proteomes.
        t0 = time.time()
        annot_cache = {}
        all_to_fetch = all_orthologs | set(seed_ids)
        if all_to_fetch:
            annot_cache = self.db.get_protein_annotations_bulk(list(all_to_fetch))
        parsed_cache = self._pre_parse_batch(annot_cache)
        _t("p4 annotations", t0)

        # Phase 5: Bulk fetch OG info from prots.ogs field
        t0 = time.time()
        ogs_map = self.db.get_protein_ogs_bulk(seed_ids)

        # Parse OG membership strings into (name, level) pairs
        seed_parsed_ogs = {}
        all_pairs = set()
        for seed_id in seed_ids:
            ogs_str = ogs_map.get(seed_id)
            if ogs_str:
                parsed = self._parse_ogs_string(ogs_str)
                seed_parsed_ogs[seed_id] = parsed
                all_pairs.update(parsed)

        # Resolve OG descriptions with cross-batch cache
        og_info_cache = self._lookup_og_info(all_pairs)
        _t(f"p5 ogs ({len(all_pairs)})", t0)

        # Phase 6: Build results
        t0 = time.time()
        results = {}
        for seed_id in seed_ids:
            orthologs = seed_orthologs.get(seed_id, set())

            # Filter annot_cache to this seed's orthologs (seed itself excluded from consensus)
            seed_annots = {oid: annot_cache[oid] for oid in orthologs if oid in annot_cache}
            seed_meta = seed_ortholog_meta.get(seed_id, {})
            annotations, annotations_confidence = (
                self._summarize_annotations(
                    seed_annots, seed_meta, target_orthologs,
                    parsed=parsed_cache,
                )
                if seed_annots else ({}, {})
            )

            # Preferred_name: use the seed's own pname when it is an informative
            # gene name. The direct DIAMOND hit is the most specific reference;
            # ortholog consensus can assign a wrong family-level name (e.g. SDC1
            # to a GAD protein). Locus IDs and multi-alias entries fall through to
            # the ortholog consensus. The seed-derived name is the highest-
            # confidence source by definition.
            seed_annot = annot_cache.get(seed_id)
            if seed_annot:
                seed_pname = (seed_annot.get("pname") or "").strip()
                if self._is_informative_pname(seed_pname):
                    annotations["Preferred_name"] = [seed_pname]
                    annotations_confidence["Preferred_name"] = "high"

            # Pick the most specific OG we have description for (deepest level
            # appears last in the parsed list)
            og_info = None
            for og_name, level in reversed(seed_parsed_ogs.get(seed_id, [])):
                og_data = og_info_cache.get((og_name, level))
                if og_data is None:
                    continue
                # Parse fprof_sum JSON if present
                fprof = None
                fprof_str = og_data.get("fprof_sum")
                if fprof_str:
                    import json
                    try:
                        fprof = json.loads(fprof_str)
                    except (json.JSONDecodeError, TypeError):
                        pass
                og_info = {
                    "name": og_name,
                    "description": og_data.get("description") or "-",
                    "cog_cat": og_data.get("cog_cat") or "-",
                    "level": level,
                    "fprof": fprof,
                }
                break

            # Preserve OG name insertion order from _parse_ogs_string
            all_ogs = [name for name, _level in seed_parsed_ogs.get(seed_id, [])]

            results[seed_id] = {
                "orthologs": list(orthologs),
                "annotations": annotations,
                "annotations_confidence": annotations_confidence,
                "tax_scope_used": seed_tax_scope_used.get(seed_id, "none"),
                "og_info": og_info,
                "ortholog_types": seed_ortholog_types.get(
                    seed_id,
                    {
                        "one2one":  set(),
                        "one2many": set(),
                        "many2one": set(),
                        "many2many": set(),
                        "all":      set(),
                    },
                ),
                "all_ogs": all_ogs,
            }

        _t("p6 build_results", t0)
        return results

    def _lookup_og_info(self, pairs):
        """Resolve OG descriptions with a persistent cache.

        Any (name, level) not in `self._og_cache` is fetched via a single
        bulk query; entries (including misses, stored as None) are cached
        so later batches in the same process reuse them.
        """
        missing = [p for p in pairs if p not in self._og_cache]
        if missing:
            fetched = self.db.get_og_info_bulk(missing)
            for p in missing:
                self._og_cache[p] = fetched.get(p)  # may be None (miss)
        return {p: self._og_cache[p] for p in pairs if self._og_cache[p] is not None}

    def _describe_seed_scope(self, seed_id: int) -> str:
        """Return the human-readable tax_scope decision for one seed
        (Phase 7.1b). ``"none"`` when no filter is active or the seed's
        species is unknown."""
        if self.lineage_filter is None:
            return "none"
        taxids = self.db.taxid_array
        if not taxids or seed_id >= len(taxids):
            return "none"
        seed_taxid = taxids[seed_id]
        if seed_taxid == 0:
            return "none"
        return self.lineage_filter.describe_scope_for_seed(seed_taxid)

    def _seed_lineage(self, seed_id: int) -> Optional[frozenset]:
        """Cached lookup of the seed's lineage as a frozenset of taxid strings.

        Used by the cascade to test whether an event's `ev_lca` is on the
        path between root and the seed. None when the lineage cache is
        absent or the seed's species is unknown."""
        if seed_id in self._seed_lineage_set_cache:
            return self._seed_lineage_set_cache[seed_id]
        result: Optional[frozenset] = None
        if self.lineage_cache is not None:
            taxids = self.db.taxid_array
            if taxids is not None and seed_id < len(taxids):
                seed_taxid = str(taxids[seed_id])
                lineage = self.lineage_cache.get(seed_taxid)
                if lineage is not None:
                    # Defensive copy as a frozenset so callers can't mutate
                    # the cache's internal state.
                    result = frozenset(lineage)
        self._seed_lineage_set_cache[seed_id] = result
        return result

    def _resolve_valid_species(self, seed_id):
        """Return the cached set of species taxids (as strings) allowed
        by the lineage_filter for this seed, or None for no filter.

        Result is memoized by the seed's species taxid so we only run
        the lineage intersection once per unique species in a batch.
        """
        if self.lineage_filter is None:
            return None
        taxids = self.db.taxid_array
        if not taxids or seed_id >= len(taxids):
            return None
        taxid = taxids[seed_id]
        if taxid == 0:
            logger.debug("seed_id %d has taxid=0 (possibly uninitialized slot)", seed_id)
        seed_taxid = str(taxid)
        if seed_taxid in self._valid_species_by_seed:
            return self._valid_species_by_seed[seed_taxid]
        scope = self.lineage_filter.get_effective_scope(seed_taxid)
        valid = (
            self.lineage_filter.get_valid_species_ids(scope)
            if scope else None
        )
        self._valid_species_by_seed[seed_taxid] = valid
        return valid

    def _resolve_allowed_og_lcas(self, seed_id):
        """Return the frozenset of taxids at-or-below the seed's
        tax-scope ceiling, or None if no filter applies.

        Used as an `og_lca` whitelist by the strict-OG mode. Memoized
        per seed-species (same scope ⇒ same allowed set).
        """
        if self.lineage_filter is None:
            return None
        taxids = self.db.taxid_array
        if not taxids or seed_id >= len(taxids):
            return None
        seed_taxid = str(taxids[seed_id])
        cache = getattr(self, "_allowed_og_lcas_by_seed", None)
        if cache is None:
            cache = self._allowed_og_lcas_by_seed = {}
        if seed_taxid in cache:
            return cache[seed_taxid]
        scope = self.lineage_filter.get_effective_scope(seed_taxid)
        allowed = (
            self.lineage_filter.get_scope_og_descendants(scope)
            if scope else None
        )
        cache[seed_taxid] = allowed
        return allowed

    def _collect_orthologs(
        self,
        seed_id: int,
        events: Dict[int, dict],
        target_taxa: Optional[Set[int]] = None,
        excluded_taxa: Optional[Set[int]] = None,
        valid_species: Optional[frozenset] = None,
        allowed_og_lcas: Optional[frozenset] = None,
    ) -> Tuple[Set[int], Dict[str, Set[int]], Dict[int, Dict[str, Any]]]:
        """Collect ortholog IDs from events, classified by orthology type.

        For each speciation event where the seed appears, the proteins on the
        opposite side are collected and classified by the cardinality of each
        side *before* species filtering. A protein can appear in multiple typed
        buckets if it participates in events with different cardinalities.

        Args:
            seed_id: Seed protein ID
            events: Dict of event_id -> event data
            target_taxa: Optional taxids to include
            excluded_taxa: Optional taxids to exclude
            valid_species: Optional frozenset of taxid strings allowed by
                the lineage scope filter. Produced by
                ``LineageFilter.get_valid_species_ids()``.
            allowed_og_lcas: Optional frozenset of taxid strings naming
                the OGs whose containing taxonomic level (``og_lca``) is
                at-or-below the seed's tax-scope ceiling. When provided,
                events from broader OGs are skipped entirely — see
                ``annotate_batch(scope_strict_og=True)``.

        Returns:
            ``(orthologs, ortholog_types, ortholog_meta)`` where:

            - ``orthologs``: filtered set of all ortholog protein IDs.
            - ``ortholog_types``: dict with keys ``"one2one"``,
              ``"one2many"``, ``"many2one"``, ``"many2many"``, ``"all"`` —
              each mapping to the subset of filtered orthologs with that
              relationship to the seed.
            - ``ortholog_meta``: dict mapping each ortholog id to its
              priority metadata as the cascade engine sees it. The metadata
              describes the *best* event the ortholog participated in
              (best = lowest sort key in the cascade order: in-seed-lineage
              first, deeper ev_lca next, lower type tier last). Keys:
              ``event_id``, ``ev_lca``, ``type``, ``type_tier``, ``depth``,
              ``in_seed_lineage``. Empty dict if no events.
        """
        typed: Dict[str, Set[int]] = {
            "one2one":  set(),
            "one2many": set(),
            "many2one": set(),
            "many2many": set(),
            "all":      set(),
        }

        # Candidate orthologs per event, keyed by relationship type.
        # We defer species filtering to a single pass after collecting all
        # candidates so that out-of-range ID warnings fire only once.
        candidates: Dict[str, Set[int]] = {
            "one2one":  set(),
            "one2many": set(),
            "many2one": set(),
            "many2many": set(),
        }

        # Per-ortholog priority metadata. The Cython path stores
        # `(packed_sort_key << 32) | event_idx` in a C++ unordered_map
        # and materialises Python payload dicts only for filtered
        # candidates at the end. The pure-Python fallback keeps the
        # original `(sort_key, payload)` tuple form. The two paths are
        # gated by `_COLLECTOR_AVAILABLE` set at import time.
        seed_lineage = self._seed_lineage(seed_id)
        if _COLLECTOR_AVAILABLE:
            collector = _OrthologCollector(seed_id)
            ortholog_meta_raw = None
        else:
            collector = None
            ortholog_meta_raw: Dict[int, Tuple[Tuple[int, int, int], Dict[str, Any]]] = {}

        for event_id, event in events.items():
            # Hard scope-OG filter: drop events whose containing OG is
            # broader than the seed's tax-scope ceiling. The same in-scope
            # orthologs already appear in lower OGs at higher cascade
            # priority, so dropping these events is biologically lossless
            # under auto-scope and removes ~99 % of fetched orthologs on
            # plant proteomes. Opt-in via `annotate_batch(scope_strict_og=True)`.
            if allowed_og_lcas is not None:
                og_lca = event.get("og_lca")
                if og_lca is None or og_lca not in allowed_og_lcas:
                    continue

            side1 = event.get("side1")
            side2 = event.get("side2")
            if side1 is None or side2 is None:
                continue

            if seed_id in side1:
                seed_side, other_side = side1, side2
            elif seed_id in side2:
                seed_side, other_side = side2, side1
            else:
                continue

            # Classify by UNFILTERED side sizes
            s, o = len(seed_side), len(other_side)
            if s == 1 and o == 1:
                rel = "one2one"
            elif s == 1 and o > 1:
                rel = "one2many"
            elif s > 1 and o == 1:
                rel = "many2one"
            else:
                rel = "many2many"

            candidates[rel].update(other_side)

            # Cascade priority metadata for this event. Computed once per
            # event, applied to every ortholog drawn from this event.
            ev_lca = str(event.get("ev_lca") or "")
            in_lineage = bool(seed_lineage) and ev_lca in seed_lineage
            depth = (
                self.lineage_cache.depth(ev_lca)
                if self.lineage_cache is not None and ev_lca
                else 0
            )
            type_tier = self.TYPE_TIERS[rel]
            payload = {
                "event_id": event_id,
                "ev_lca": ev_lca,
                "type": rel,
                "type_tier": type_tier,
                "depth": depth,
                "in_seed_lineage": in_lineage,
            }
            if collector is not None:
                packed = _pack_sort_key(in_lineage, depth, type_tier)
                collector.add_event(other_side, packed, payload)
            else:
                sort_key = (0 if in_lineage else 1, -depth, type_tier)
                for oid in other_side:
                    if oid == seed_id:
                        continue
                    prev = ortholog_meta_raw.get(oid)
                    if prev is None or sort_key < prev[0]:
                        ortholog_meta_raw[oid] = (sort_key, payload)

        # Remove the seed itself from every candidate bucket
        for bucket in candidates.values():
            bucket.discard(seed_id)

        all_candidates: Set[int] = set()
        for bucket in candidates.values():
            all_candidates.update(bucket)

        if not all_candidates:
            return typed["all"], typed, {}

        taxids = self.db.taxid_array

        # No taxid array → return all candidates unfiltered
        if taxids is None:
            for rel, bucket in candidates.items():
                typed[rel].update(bucket)
                typed["all"].update(bucket)
            if collector is not None:
                ortholog_meta = collector.export_meta(typed["all"])
            else:
                ortholog_meta = {
                    oid: payload
                    for oid, (_, payload) in ortholog_meta_raw.items()
                    if oid in typed["all"]
                }
            return typed["all"], typed, ortholog_meta

        # Combined species filter: auto/scope lineage + manual taxa filters.
        # Scope lineage uses taxid strings; manual filters typically use ints.
        want_scope = valid_species is not None
        want_manual = bool(target_taxa or excluded_taxa)

        _warned_oor = False
        _oor_count = 0
        n_taxids = len(taxids)

        # Build a set of candidates that pass all filters, then classify
        # into typed buckets from that filtered set.
        filtered: Set[int] = set()
        for oid in all_candidates:
            if oid >= n_taxids:
                if not _warned_oor:
                    logging.getLogger(__name__).warning(
                        "Protein ID %d >= taxid_array length %d; discarding "
                        "out-of-range orthologs (seed %d). DB may be "
                        "inconsistent.",
                        oid,
                        n_taxids,
                        seed_id,
                    )
                    _warned_oor = True
                _oor_count += 1
                continue
            if want_scope or want_manual:
                taxid_int = taxids[oid]
                if want_scope and str(taxid_int) not in valid_species:
                    continue
                if want_manual:
                    if excluded_taxa and taxid_int in excluded_taxa:
                        continue
                    if target_taxa and taxid_int not in target_taxa:
                        continue
            filtered.add(oid)

        # Distribute filtered proteins into their typed buckets
        for rel, bucket in candidates.items():
            typed_bucket = filtered & bucket
            typed[rel].update(typed_bucket)
        typed["all"] = filtered

        if collector is not None:
            ortholog_meta = collector.export_meta(filtered)
        else:
            ortholog_meta = {
                oid: payload
                for oid, (_, payload) in ortholog_meta_raw.items()
                if oid in filtered
            }
        return typed["all"], typed, ortholog_meta

    # `target_orthologs` floor → set of allowed event types accepted by
    # the cascade. `one2one` (strictest) accepts only 1:1 events; `all`
    # / `many2many` (loosest) accepts every type. `one2many` and
    # `many2one` accept 1:1 plus the named asymmetric tier — the user
    # can still select which side of medium-tier ambiguity is OK.
    TARGET_ORTHOLOGS_FLOORS = {
        "all":       frozenset({"one2one", "one2many", "many2one", "many2many"}),
        "many2many": frozenset({"one2one", "one2many", "many2one", "many2many"}),
        "one2many":  frozenset({"one2one", "one2many"}),
        "many2one":  frozenset({"one2one", "many2one"}),
        "one2one":   frozenset({"one2one"}),
    }

    def _pre_parse_batch(
        self,
        annot_data: Dict[int, dict],
    ) -> Dict[int, Dict[str, Tuple[str, ...]]]:
        """Pre-parse annotation comma-strings into tuples of clean values
        once per ortholog, ahead of the cascade walk.

        Without this cache the cascade re-splits the same comma-string
        every time an ortholog appears as a candidate donor for a field
        — on plant proteomes (~1300 orthologs/seed × 13 fields × ~1000
        seeds per batch) that's the single hottest python loop in the
        engine. Pre-parsing tightens phase 6 substantially while
        preserving cascade outputs byte-for-byte.

        ``pname`` is normalized to a 1-tuple ``(stripped,)`` (or ``()``).
        ``gos`` is parsed via :meth:`_parse_gos` (strips evidence codes,
        keeps only ``GO:`` prefixed terms). Every other field becomes a
        distinct-preserving tuple — same insertion order, deduplicated.
        """
        parsed: Dict[int, Dict[str, Tuple[str, ...]]] = {}
        for oid, annot in annot_data.items():
            if not annot:
                parsed[oid] = {}
                continue
            row: Dict[str, Tuple[str, ...]] = {}
            for field in self.ANNOTATION_FIELDS:
                raw = annot.get(field)
                if not raw:
                    continue
                if field == "pname":
                    s = raw.strip()
                    if s:
                        row[field] = (s,)
                elif field == "gos":
                    gos = tuple(self._parse_gos(raw))
                    if gos:
                        row[field] = gos
                else:
                    seen: List[str] = []
                    seen_set: Set[str] = set()
                    for v in str(raw).split(","):
                        v = v.strip()
                        if v and v not in seen_set:
                            seen.append(v)
                            seen_set.add(v)
                    if seen:
                        row[field] = tuple(seen)
            parsed[oid] = row
        return parsed

    def _summarize_annotations(
        self,
        annot_data: Dict[int, dict],
        ortholog_meta: Optional[Dict[int, dict]] = None,
        target_orthologs: str = "all",
        parsed: Optional[Dict[int, Dict[str, Tuple[str, ...]]]] = None,
    ) -> Tuple[Dict[str, Any], Dict[str, str]]:
        """Cascade summary: per-source closest-ev_lca + type-priority winner.

        For each functional source (KEGG_ko, GOs, Pfam, ...) independently,
        donors are walked from closest+best-typed first. The first bucket
        — defined by the cascade key ``(in_seed_lineage, -ev_lca_depth,
        type_tier)`` — that has any donor with a non-empty value for that
        source wins, and the consensus is taken across only the donors in
        that winning bucket. This preserves "if no 1:1 donors annotate
        source S, fall through to 1:many or many:many — but with caution"
        as a per-source rule, while keeping the strict upper-clade ceiling
        from ``LineageFilter`` (already applied during phase 3).

        ``target_orthologs`` is a *floor*: types not in
        ``TARGET_ORTHOLOGS_FLOORS[target_orthologs]`` are excluded from the
        cascade entirely. The post-aggregation filter that the v2-style
        mapper shim applied is no longer needed.

        Args:
            annot_data: ``{protein_id: annotation_dict}`` for every
                ortholog plus optionally the seed.
            ortholog_meta: ``{protein_id: meta_dict}`` produced by
                :meth:`_collect_orthologs`. Required for cascade mode.
                When absent or empty the flat aggregation is used (same
                logic as v3 phase 0) and confidence is empty — needed for
                back-compat with the single-protein ``annotate()`` path
                and the existing test_annotate.py expectations.
            target_orthologs: ``"all"``, ``"many2many"``, ``"one2many"``,
                ``"many2one"`` or ``"one2one"``. Anything else is treated
                as ``"all"``.

        Returns:
            ``(annotations, confidence)`` where ``annotations`` is the
            same shape returned by the legacy summarizer (output-named
            keys: ``Preferred_name``, ``GOs``, ``PFAMs``, ``KEGG_ko``,
            ...) and ``confidence`` is ``{output_field: "high" | "medium"
            | "low"}`` for every source actually emitted.
        """
        if not annot_data:
            return {}, {}

        # Back-compat path: no metadata → flat aggregation, no confidence.
        if not ortholog_meta:
            return self._summarize_annotations_flat(annot_data), {}

        # Use the pre-parsed cache if the caller provided one (the batch
        # path always does). Otherwise build a tiny one-shot cache here.
        if parsed is None:
            parsed = self._pre_parse_batch(annot_data)

        allowed_types = self.TARGET_ORTHOLOGS_FLOORS.get(
            target_orthologs, self.TARGET_ORTHOLOGS_FLOORS["all"]
        )

        # Bucket orthologs by cascade priority key. Keys are tuples so
        # they sort lexicographically: smallest first = best donors.
        buckets: Dict[Tuple[int, int, int], List[int]] = defaultdict(list)
        for oid, meta in ortholog_meta.items():
            if oid not in parsed:
                continue
            if meta["type"] not in allowed_types:
                continue
            key = (
                0 if meta["in_seed_lineage"] else 1,
                -meta["depth"],
                meta["type_tier"],
            )
            buckets[key].append(oid)

        if not buckets:
            return {}, {}

        priority_order = sorted(buckets.keys())

        annotations: Dict[str, Any] = {}
        confidence: Dict[str, str] = {}

        for field in self.ANNOTATION_FIELDS:
            output_field = self._cascade_output_field(field)
            for prio_key in priority_order:
                # Pre-parsed lookup is O(1); membership check just asks
                # whether the field exists in the per-ortholog dict.
                contributors = [
                    oid for oid in buckets[prio_key]
                    if field in parsed.get(oid, ())
                ]
                if not contributors:
                    continue
                values = self._aggregate_field(field, contributors, parsed)
                if values:
                    annotations[output_field] = values
                    confidence[output_field] = self.TIER_CONFIDENCE[prio_key[2]]
                # First bucket with non-empty contributors is the winner —
                # stop the cascade for this source even if `_aggregate_field`
                # filtered everything out (the bucket *had* signal; lower
                # buckets shouldn't be promoted on its behalf).
                break

        return annotations, confidence

    def _cascade_output_field(self, field: str) -> str:
        """Map an internal annotation field name to the output key the
        cascade emits."""
        if field == "pname":
            return "Preferred_name"
        if field == "gos":
            return "GOs"
        if field == "pfam":
            return "PFAMs"
        return self._field_to_output(field)

    def _aggregate_field(
        self,
        field: str,
        contributors: List[int],
        parsed: Dict[int, Dict[str, Tuple[str, ...]]],
    ) -> List[str]:
        """Aggregate one annotation field across the cascade winning bucket.

        Reads pre-parsed per-ortholog tuples instead of re-splitting comma
        strings on every cascade lookup (Phase 4 hot-path optimization).

        ``pname``: most-common across the bucket. No ≥2-occurrence rule —
        cascade selectivity already guarantees the bucket holds the
        highest-confidence donors. Counter is used because we need the
        winner.

        Other fields: distinct values via set union — order doesn't
        matter (downstream sorts at write time) and dedup is far faster
        than counting then taking ``.keys()``.
        """
        if field == "pname":
            counter: Counter = Counter()
            for oid in contributors:
                tup = parsed.get(oid, {}).get(field)
                if tup:
                    counter[tup[0]] += 1
            if not counter:
                return []
            return [counter.most_common(1)[0][0]]

        out: Set[str] = set()
        for oid in contributors:
            tup = parsed.get(oid, {}).get(field)
            if tup:
                out.update(tup)
        return list(out)

    def _summarize_annotations_flat(
        self,
        annot_data: Dict[int, dict],
    ) -> Dict[str, Any]:
        """Legacy flat aggregation used by the no-metadata path (single-
        protein ``annotate()`` and pre-cascade tests). Behavior matches
        v3-phase0-end exactly."""
        counters = defaultdict(Counter)
        for annot in annot_data.values():
            if not annot:
                continue
            for field in self.ANNOTATION_FIELDS:
                value = annot.get(field)
                if not value:
                    continue
                if field == "pname":
                    counters[field][value.strip()] += 1
                elif field == "gos":
                    for go in self._parse_gos(value):
                        counters[field][go] += 1
                else:
                    for v in str(value).split(","):
                        v = v.strip()
                        if v:
                            counters[field][v] += 1

        result = {}
        for field, counter in counters.items():
            if field == "pname":
                most_common = counter.most_common(1)
                if most_common and most_common[0][1] >= 2:
                    result["Preferred_name"] = [most_common[0][0]]
            elif field == "gos":
                result["GOs"] = list(counter.keys())
            elif field == "pfam":
                total = len(annot_data)
                result["PFAMs"] = [
                    k for k, c in counter.items()
                    if c > 1 and c / total > 0.05
                ]
            else:
                output_field = self._field_to_output(field)
                result[output_field] = list(counter.keys())
        return result

    def _parse_ogs_string(self, ogs_string: str) -> List[Tuple[str, str]]:
        """Parse OGs string from prots.ogs field.

        Format: "cluster@taxid|clade|taxon_name,cluster@taxid|clade|taxon_name,..."
        Returns list of (og_name, level) where og_name is "cluster@taxid|clade"
        and level is the taxid.
        """
        result = []
        for og in ogs_string.split(","):
            og = og.strip()
            if not og:
                continue
            # Split by | - format is cluster@taxid|clade|taxon_name
            parts = og.split("|")
            if len(parts) >= 2:
                # og_name is first two parts: cluster@taxid|clade
                og_name = "|".join(parts[:2])
                # level is the taxid part from cluster@taxid
                if "@" in parts[0]:
                    level = parts[0].split("@")[1]
                else:
                    level = "-"
                result.append((og_name, level))
        return result

    def _is_informative_pname(self, pname: str) -> bool:
        """Return True if pname is a real gene symbol worth transferring.

        Rejects:
        - empty string
        - multi-alias entries (contain a comma, e.g. "TP53,P53")
        - NCBI locus IDs  (LOC followed by digits, e.g. "LOC104234906")
        - purely numeric strings (e.g. "11423014")
        """
        if not pname:
            return False
        if "," in pname:
            return False
        if pname.startswith("LOC") and pname[3:].isdigit():
            return False
        if pname.isdigit():
            return False
        return True

    def _parse_gos(self, go_string: str) -> List[str]:
        """Parse GO terms from string, filtering by evidence if needed."""
        gos = []
        for term in go_string.split(","):
            term = term.strip()
            if term.startswith("GO:"):
                gos.append(term.split("|")[0])  # Strip evidence code if present
        return gos

    def _field_to_output(self, field: str) -> str:
        """Map internal field names to output format."""
        mapping = {
            "kegg_ko": "KEGG_ko",
            "kegg_ec": "EC",
            "kegg_pathway": "KEGG_Pathway",
            "kegg_module": "KEGG_Module",
            "kegg_reaction": "KEGG_Reaction",
            "kegg_rclass": "KEGG_rclass",
            "kegg_brite": "BRITE",
            "kegg_tc": "KEGG_TC",
            "kegg_cazy": "CAZy",
            "bigg_reaction": "BiGG_Reaction",
        }
        return mapping.get(field, field)

    def _get_og_info(self, events: Dict[int, dict]) -> Optional[dict]:
        """Get OG info from events."""
        for event in events.values():
            og_name = event.get("og")
            if og_name:
                og = self.db.get_og_info(og_name)
                if og:
                    return {
                        "name": og_name,
                        "description": og.get("description"),
                        "cog_cat": og.get("cog_cat"),
                        "level": og.get("level"),
                    }
        return None


# Convenience function for simple use cases
def annotate_protein(
    db_path: str,
    seed_id: int,
    target_taxa: Optional[Set[int]] = None,
    excluded_taxa: Optional[Set[int]] = None,
) -> Dict[str, Any]:
    """Annotate a single protein.

    Convenience wrapper that handles database connection.

    Args:
        db_path: Path to eggnog.db
        seed_id: Integer protein ID
        target_taxa: Optional taxids to include
        excluded_taxa: Optional taxids to exclude

    Returns:
        Annotation result dict
    """
    with EggnogDB(db_path) as db:
        engine = AnnotationEngine(db)
        return engine.annotate(seed_id, target_taxa, excluded_taxa)
