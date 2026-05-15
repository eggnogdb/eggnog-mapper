"""Batch annotation logic for v7+ databases.

High-throughput annotation using bulk queries and simplified ortholog collection.
Eliminates species grouping bottleneck for ~100 q/s performance.

Usage:
    from eggnogmapper.annotator.e7.annotate import AnnotationEngine

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

# GO namespace cascade — env var overrides the default OBO path. When the
# file is missing we fall back to the legacy combined-GO cascade (one
# winning tier across all namespaces) and emit a single warning.
_GO_OBO_ENV = "EGGNOG_GO_OBO"
_GO_OBO_DEFAULT = "/app/data/e7/full/source/reference/go-basic.obo"

# Internal per-namespace cascade keys. Each is parsed independently in the
# cascade (so cellular_component can fall through to a deeper tier when the
# winning MF/BP donors lack CC), but they all merge into the single output
# field "GOs".
_GO_NS_FIELDS = ("gos_mf", "gos_bp", "gos_cc")
_GO_NS_TO_FIELD = {
    "molecular_function": "gos_mf",
    "biological_process": "gos_bp",
    "cellular_component": "gos_cc",
}

# Module-level cache: { obo_path: { "GO:xxxxxxx": "molecular_function" | ... } }.
# Populated lazily on first parse; reused across AnnotationEngine instances.
_GO_NAMESPACE_CACHE: Dict[str, Dict[str, str]] = {}
_GO_OBO_MISSING_WARNED: Set[str] = set()


def _load_go_namespace_map(obo_path: str) -> Optional[Dict[str, str]]:
    """Parse a go-basic.obo file once and return ``{GO:xxxxxxx: namespace}``.

    Pure-Python OBO parser (no goatools dependency). Only ``[Term]``
    stanzas are read; ``[Typedef]`` and unknown stanzas are skipped.
    Obsolete terms are dropped. ``alt_id`` lines map to the same
    namespace as the parent ``id``.

    Returns ``None`` (and warns once) when the file cannot be opened
    so the engine can fall back to the legacy combined-GO cascade.
    """
    if obo_path in _GO_NAMESPACE_CACHE:
        return _GO_NAMESPACE_CACHE[obo_path]
    try:
        f = open(obo_path, "r", encoding="utf-8")
    except OSError as exc:
        if obo_path not in _GO_OBO_MISSING_WARNED:
            logger.warning(
                "GO namespace OBO not found at %s (%s) — falling back "
                "to legacy combined-GO cascade. Set %s to override.",
                obo_path, exc, _GO_OBO_ENV,
            )
            _GO_OBO_MISSING_WARNED.add(obo_path)
        return None

    mapping: Dict[str, str] = {}
    in_term = False
    cur_id: Optional[str] = None
    cur_alt: List[str] = []
    cur_ns: Optional[str] = None
    cur_obsolete = False

    def _commit():
        if cur_id and cur_ns and not cur_obsolete:
            mapping[cur_id] = cur_ns
            for alt in cur_alt:
                mapping[alt] = cur_ns

    with f:
        for line in f:
            line = line.rstrip()
            if not line:
                if in_term:
                    _commit()
                    in_term = False
                    cur_id = None
                    cur_alt = []
                    cur_ns = None
                    cur_obsolete = False
                continue
            if line.startswith("["):
                if in_term:
                    _commit()
                in_term = (line == "[Term]")
                cur_id = None
                cur_alt = []
                cur_ns = None
                cur_obsolete = False
                continue
            if not in_term:
                continue
            if line.startswith("id: "):
                cur_id = line[4:].strip()
            elif line.startswith("alt_id: "):
                cur_alt.append(line[8:].strip())
            elif line.startswith("namespace: "):
                ns = line[11:].strip()
                if ns in _GO_NS_TO_FIELD:
                    cur_ns = ns
            elif line.startswith("is_obsolete:"):
                if line.split(":", 1)[1].strip().lower() == "true":
                    cur_obsolete = True
        # Trailing term (no blank line at EOF)
        if in_term:
            _commit()

    logger.info(
        "Loaded GO namespace map from %s: %d terms across %d namespaces",
        obo_path, len(mapping), 3,
    )
    _GO_NAMESPACE_CACHE[obo_path] = mapping
    return mapping

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


# Process-pool worker globals.
#
# When `annotate_batch(pool=…)` is called with a pool created via the
# fork start method, workers inherit the parent's `AnnotationEngine`
# state — taxid_array (~226 MB numpy), lineage_cache, og_cache, etc. —
# via copy-on-write. They only need to (a) drop the inherited sqlite3
# connection and reopen one of their own (sqlite3 connections are not
# fork-safe) and (b) keep a module-level pointer to the engine so the
# pickle-by-name worker function can find it.
#
# The parent calls `_register_worker_engine(self)` *before* the Pool is
# forked; the post-fork initializer then reopens SQLite. The worker fn
# below dispatches each sub-batch to the inherited engine.
_WORKER_ENGINE: "AnnotationEngine | None" = None


def _register_worker_engine(engine: "AnnotationEngine") -> None:
    """Set the module-global engine reference visible to forked workers.

    Called by the parent process *before* a fork-context Pool is
    created so each worker inherits the same `AnnotationEngine`
    instance (with taxid_array, lineage_cache, og_cache already loaded).
    """
    global _WORKER_ENGINE
    _WORKER_ENGINE = engine


def _worker_init_after_fork() -> None:
    """Pool initializer — runs once per worker after fork().

    The inherited `_WORKER_ENGINE.db.conn` is a copy of the parent's
    sqlite3 connection. That object is *not* fork-safe (shared fds /
    locks can deadlock the next executable statement). Reopen it here.
    """
    if _WORKER_ENGINE is not None and _WORKER_ENGINE.db is not None:
        _WORKER_ENGINE.db.reopen_connection()


def _worker_annotate_subbatch(args):
    """Pool worker — annotate one sub-batch of seed IDs.

    Picklable by name. The engine state is reached via the module-global
    `_WORKER_ENGINE` set in the parent before fork. Returns the per-seed
    result dict so the parent can merge.
    """
    seed_ids, target_taxa, excluded_taxa, target_orthologs = args
    return _WORKER_ENGINE._annotate_batch_inproc(
        seed_ids,
        target_taxa=target_taxa,
        excluded_taxa=excluded_taxa,
        target_orthologs=target_orthologs,
    )


class AnnotationEngine:
    """Unified annotation engine for v7+ eggNOG databases.

    Used by both eggnog-mapper and eggnog-website for consistent annotation logic.
    Optimized for high throughput with bulk queries and simplified ortholog collection.
    """

    # Annotation fields to aggregate from orthologs.
    #
    # ``gos`` is split into three pseudo-sources — gos_mf, gos_bp,
    # gos_cc — so each GO namespace runs its own per-tier "first
    # non-empty wins" cascade. They all merge into the single output
    # key "GOs". When the OBO map is unavailable the engine falls back
    # to the legacy combined ``gos`` field at runtime (see
    # :meth:`_pre_parse_batch` and :meth:`_summarize_annotations`).
    ANNOTATION_FIELDS = [
        "pname", "gos_mf", "gos_bp", "gos_cc",
        "kegg_ko", "kegg_ec", "kegg_pathway",
        "kegg_module", "kegg_reaction", "kegg_rclass", "kegg_brite",
        "kegg_tc", "kegg_cazy", "bigg_reaction", "pfam"
    ]
    # Legacy combined-GO field name. Used by the flat fallback path and
    # by :meth:`_pre_parse_batch` when the OBO map is unavailable.
    LEGACY_GO_FIELD = "gos"

    # Per-worker cache size limits. Workers are long-lived (no
    # maxtasksperchild), so these bounds prevent unbounded RSS growth
    # over a 24h+ run. When either limit is hit, the oldest half of the
    # cache is evicted (dict insertion order, Python 3.7+).
    _OG_CACHE_MAX = 100_000
    _LINEAGE_CACHE_MAX = 100_000

    # Mapping from ortholog-relationship type to cascade tier. The cascade
    # walks tiers in order 0 → 2; within a tier, donors are sorted by
    # ev_lca proximity to the seed. Confidence label is derived from the
    # winning tier: 0 → "high", 1 → "medium", 2 → "low".
    TYPE_TIERS = {
        "self":     0,
        "one2one":  0,
        "one2many": 1,
        "many2one": 1,
        "many2many": 2,
    }
    TIER_CONFIDENCE = {0: "high", 1: "medium", 2: "low"}

    def __init__(
        self,
        db: EggnogDB,
        lineage_filter=None,
        lineage_cache=None,
        go_obo_path: Optional[str] = None,
        donor_pool: str = "closest",
    ):
        """Initialize annotation engine.

        Args:
            db: EggnogDB instance with loaded taxid array.
            lineage_filter: Optional LineageFilter for tax-scope filtering.
                When set, orthologs are pruned per-seed to the seed's
                taxonomic ceiling (``ev_lca ≤ ceiling``) before the
                expensive annotation-fetch phase — the single largest
                speedup for proteome-scale runs.
            lineage_cache: Optional LineageCache used by the cascade engine
                to rank events by ev_lca depth.  If None and
                lineage_filter is set, the filter's own cache is reused.
                If both are absent the cascade falls back to type-tier
                priority only (no ev_lca distance ordering).
            go_obo_path: Optional path to a ``go-basic.obo`` file.  When
                provided (or when the ``EGGNOG_GO_OBO`` env var is set),
                the cascade splits GO terms by namespace (MF / BP / CC)
                and runs an independent first-non-empty-tier walk per
                namespace.  All three subnamespaces merge into the single
                output key ``"GOs"``.  If the file is missing the engine
                logs a warning once and falls back to the legacy combined
                cascade.  Default: ``$EGGNOG_GO_OBO`` or the canonical path
                under ``data/e7/full/source/reference/go-basic.obo``.
            donor_pool: ``"closest"`` (default) or ``"union"``.

                - ``"closest"``: walk tiers in priority order; the first
                  non-empty tier for each annotation source wins (original
                  behaviour).
                - ``"union"``: walk *all* tiers, union values across every
                  tier; confidence is the best (smallest) tier seen for
                  that source.
        """
        self.db = db
        self.lineage_filter = lineage_filter
        self.donor_pool = donor_pool
        self.lineage_cache = lineage_cache or (
            lineage_filter.lineage_cache if lineage_filter is not None else None
        )
        # Lazy-load the GO namespace map. The mapping is module-cached so
        # repeated engine construction in tests / batch_annotate doesn't
        # re-parse the 32 MB OBO file.
        self.go_obo_path = (
            go_obo_path or os.environ.get(_GO_OBO_ENV) or _GO_OBO_DEFAULT
        )
        self._go_namespace_map: Optional[Dict[str, str]] = (
            _load_go_namespace_map(self.go_obo_path)
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
        # Resolve ceiling for this seed.
        ceiling_taxid = self._resolve_ceiling(seed_id)

        # Get events for this protein
        event_ids = self.db.get_events_for_protein(seed_id)
        if not event_ids:
            seed_taxid_str = self._get_seed_taxid_str(seed_id)
            lineage_str = self._lineage_str_for_taxid(seed_taxid_str)
            return {
                "orthologs": [],
                "annotations": {},
                "annotations_confidence": {},
                "tax_ceiling": self._ceiling_name(ceiling_taxid),
                "farthest_donor_taxid": seed_taxid_str or "-",
                "farthest_donor_lineage": lineage_str or "-",
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
            seed_id, events, target_taxa, excluded_taxa,
            ceiling_taxid=ceiling_taxid,
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

        farthest_taxid, farthest_lineage = self._compute_farthest_donor(
            ortholog_meta, seed_id
        )

        return {
            "orthologs": list(orthologs),
            "annotations": annotations,
            "annotations_confidence": annotations_confidence,
            "tax_ceiling": self._ceiling_name(ceiling_taxid),
            "farthest_donor_taxid": farthest_taxid,
            "farthest_donor_lineage": farthest_lineage,
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
        pool=None,
        sub_batch_size: int = 125,
    ) -> Dict[int, Dict[str, Any]]:
        """Annotate multiple seed orthologs efficiently.

        Uses bulk queries to minimize database round-trips.

        Args:
            seed_ids: List of integer protein IDs.
            target_taxa: Optional set of taxids to include.
            excluded_taxa: Optional set of taxids to exclude.
            target_orthologs: Cascade type floor; see
                ``TARGET_ORTHOLOGS_FLOORS``.
            pool: Optional ``multiprocessing.Pool`` (created with the
                fork start method) for parallel sub-batch dispatch.  The
                caller must have called
                :func:`_register_worker_engine` on this engine before
                forking the pool, and the pool must be created with
                ``initializer=_worker_init_after_fork`` so each worker
                reopens its own SQLite connection.  When ``None``
                (default), the batch is processed in-process.
            sub_batch_size: Sub-batch size for parallel dispatch.
                Smaller = better load balance, more dispatch overhead.
                The default of 125 fits comfortably under the cascade's
                event-cache reuse window.

        Returns:
            Dict mapping seed_id → annotation result.  Each result
            contains the keys ``"tax_ceiling"``,
            ``"farthest_donor_taxid"``, and
            ``"farthest_donor_lineage"`` in addition to the standard
            annotation fields.
        """
        if not seed_ids:
            return {}

        # Parallel path — slice and dispatch to the pool. Each worker
        # owns its own SQLite connection (fork-safe via the post-fork
        # initializer); the heavy in-memory state (taxid_array,
        # lineage_cache) is COW-shared with this process.
        if pool is not None and len(seed_ids) > sub_batch_size:
            sub_batches = [
                seed_ids[i:i + sub_batch_size]
                for i in range(0, len(seed_ids), sub_batch_size)
            ]
            args_list = [
                (sb, target_taxa, excluded_taxa, target_orthologs)
                for sb in sub_batches
            ]
            merged: Dict[int, Dict[str, Any]] = {}
            # apply_async + get(timeout) instead of imap_unordered so that
            # a single stuck worker (blocked on SQLite, mmap page fault, or
            # IPC pipe) never hangs the parent forever. imap_unordered waits
            # for ALL results with no escape; apply_async lets us set a
            # per-task deadline and let the caller's finally block clean up.
            import multiprocessing as _mp
            _TASK_TIMEOUT = 120  # seconds — generous for a 125-seed sub-batch
            pending = [
                pool.apply_async(_worker_annotate_subbatch, (a,))
                for a in args_list
            ]
            for ar in pending:
                try:
                    merged.update(ar.get(timeout=_TASK_TIMEOUT))
                except _mp.TimeoutError:
                    raise RuntimeError(
                        f"Annotation worker timed out after {_TASK_TIMEOUT}s "
                        "on a 125-seed sub-batch — pool will be terminated."
                    )
            return merged

        return self._annotate_batch_inproc(
            seed_ids,
            target_taxa=target_taxa,
            excluded_taxa=excluded_taxa,
            target_orthologs=target_orthologs,
        )

    def _annotate_batch_inproc(
        self,
        seed_ids: List[int],
        target_taxa: Optional[Set[int]] = None,
        excluded_taxa: Optional[Set[int]] = None,
        target_orthologs: str = "all",
    ) -> Dict[int, Dict[str, Any]]:
        """In-process implementation of ``annotate_batch``.

        The public method dispatches to either this or the parallel path;
        this also serves as the worker entrypoint after fork.
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

        # Phase 3: Collect orthologs for each seed (with ev_lca ceiling prune)
        t0 = time.time()
        seed_orthologs: Dict[int, Set[int]] = {}
        seed_ortholog_types: Dict[int, Dict[str, Set[int]]] = {}
        seed_ortholog_meta: Dict[int, Dict[int, Dict[str, Any]]] = {}
        seed_tax_ceiling: Dict[int, str] = {}
        all_orthologs: Set[int] = set()
        for seed_id in seed_ids:
            eids = event_index.get(seed_id, [])
            # Filter events to those in our cache
            seed_events = {eid: events_cache[eid] for eid in eids if eid in events_cache}
            ceiling_taxid = self._resolve_ceiling(seed_id)
            valid_species = self._resolve_valid_species_for_ceiling(
                seed_id, ceiling_taxid
            )
            allowed_og_lcas = self._resolve_allowed_og_lcas_for_ceiling(ceiling_taxid)
            orthologs, ortholog_types, ortholog_meta = self._collect_orthologs(
                seed_id,
                seed_events,
                target_taxa,
                excluded_taxa,
                valid_species,
                ceiling_taxid=ceiling_taxid,
                allowed_og_lcas=allowed_og_lcas,
            )
            seed_orthologs[seed_id] = orthologs
            seed_ortholog_types[seed_id] = ortholog_types
            seed_ortholog_meta[seed_id] = ortholog_meta
            seed_tax_ceiling[seed_id] = self._ceiling_name(ceiling_taxid)
            all_orthologs.update(orthologs)
        _t(f"p3 collect_orthologs ({len(all_orthologs)})", t0)

        # Phase 4: Bulk fetch annotations for all orthologs and seeds.
        # Seeds are fetched so they can act as the cascade's tier-0 self-donor
        # for every functional source (and so their pname can be promoted to
        # the primary Preferred_name). Pre-parse once per ortholog: the
        # cascade walk is O(seeds × buckets × fields), so re-splitting
        # comma-strings inside that loop dominated phase 6 on plant proteomes.
        t0 = time.time()
        annot_cache = {}
        all_to_fetch = all_orthologs | set(seed_ids)
        if all_to_fetch:
            annot_cache = self.db.get_protein_annotations_bulk(list(all_to_fetch))
        parsed_cache = self._pre_parse_batch(annot_cache)
        # The seed contributes every functional source through the cascade
        # except pname — pname stays on the post-cascade promotion path
        # below so an uninformative seed pname (locus IDs, multi-aliases,
        # numeric strings) still falls through to ortholog consensus.
        for seed_id in seed_ids:
            seed_parsed = parsed_cache.get(seed_id)
            if seed_parsed:
                seed_parsed.pop("pname", None)
        _t("p4 annotations", t0)

        # Inject seed-as-self-donor into each seed's ortholog_meta. The
        # diamond hit is the closest possible reference and, for SwissProt-
        # curated seeds, carries experimentally-validated annotations the
        # cascade should prefer over inferred OG-paralog consensus. The
        # synthetic depth (seed_taxid lineage depth + 1) is strictly deeper
        # than any real ev_lca on the seed's lineage, so the seed's bucket
        # wins the cascade for every source it provides; for sources where
        # the seed is empty (e.g. KEGG_ko on a Pfam-only seed), the bucket
        # has no contributors and the cascade falls through to orthologs.
        taxids = self.db.taxid_array
        for seed_id in seed_ids:
            if seed_id not in annot_cache:
                continue
            seed_depth = 0
            if (
                self.lineage_cache is not None
                and taxids is not None
                and seed_id < len(taxids)
            ):
                seed_taxid = str(taxids[seed_id])
                if seed_taxid:
                    seed_depth = self.lineage_cache.depth(seed_taxid)
            seed_ortholog_meta.setdefault(seed_id, {})[seed_id] = {
                "event_id": -1,
                "ev_lca": "",
                # Empty og_lca on the seed-self-donor is the sentinel that
                # bypasses the GO scope filter — the seed is always in scope
                # for its own annotation.
                "og_lca": "",
                "type": "self",
                "type_tier": 0,
                "depth": seed_depth + 1,
                "in_seed_lineage": True,
            }

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

            # Filter annot_cache to this seed's orthologs PLUS the seed itself.
            # The seed enters the cascade as a synthetic tier-0 self-donor (see
            # phase 4 above); for every source the seed has, its bucket wins
            # the cascade with confidence "high" and overrides consensus that
            # would otherwise be drawn from OG-paralog inference.
            cascade_ids = orthologs | {seed_id} if seed_id in annot_cache else orthologs
            seed_annots = {oid: annot_cache[oid] for oid in cascade_ids if oid in annot_cache}
            seed_meta = seed_ortholog_meta.get(seed_id, {})
            annotations, annotations_confidence = (
                self._summarize_annotations(
                    seed_annots,
                    seed_meta,
                    target_orthologs,
                    parsed=parsed_cache,
                    donor_pool=self.donor_pool,
                )
                if seed_annots
                else ({}, {})
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

            farthest_taxid, farthest_lineage = self._compute_farthest_donor(
                seed_meta, seed_id
            )

            results[seed_id] = {
                "orthologs": list(orthologs),
                "annotations": annotations,
                "annotations_confidence": annotations_confidence,
                "tax_ceiling": seed_tax_ceiling.get(seed_id, "-"),
                "farthest_donor_taxid": farthest_taxid,
                "farthest_donor_lineage": farthest_lineage,
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
            if len(self._og_cache) > self._OG_CACHE_MAX:
                evict_n = len(self._og_cache) // 2
                for k in list(self._og_cache)[:evict_n]:
                    del self._og_cache[k]
        return {p: self._og_cache[p] for p in pairs if self._og_cache[p] is not None}

    def _get_seed_taxid_str(self, seed_id: int) -> str:
        """Return the species taxid string for ``seed_id``, or ``""``."""
        taxids = self.db.taxid_array
        if taxids is None or seed_id >= len(taxids):
            return ""
        taxid_int = taxids[seed_id]
        return str(taxid_int) if taxid_int else ""

    def _lineage_str_for_taxid(self, taxid_str: str) -> str:
        """Return the semicolon-separated lineage string for a taxid.

        Used to populate ``farthest_donor_lineage``.  Returns ``"-"``
        when the taxid is empty or not in the lineage cache.
        """
        if not taxid_str or self.lineage_cache is None:
            return "-"
        track = getattr(self.lineage_cache, "_tracks", {})
        if not track:
            # Fallback: use the set form sorted arbitrarily.
            lineage = self.lineage_cache.get(taxid_str)
            if lineage:
                return ";".join(sorted(lineage))
            return "-"
        ordered = track.get(taxid_str)
        if ordered:
            return ";".join(ordered)
        return "-"

    def _ceiling_name(self, ceiling_taxid: str) -> str:
        """Return the human-readable ceiling name for output.

        Delegates to ``lineage_filter.ceiling_name`` when a filter is
        present; otherwise returns the raw taxid (or ``"-"``).
        """
        if self.lineage_filter is not None:
            return self.lineage_filter.ceiling_name(ceiling_taxid)
        from .ceiling import CEILING_NAMES
        return CEILING_NAMES.get(ceiling_taxid, ceiling_taxid or "-")

    def _resolve_ceiling(self, seed_id: int) -> str:
        """Resolve the ev_lca ceiling taxid for a seed protein.

        Returns ``"root"`` (no filter) when no ``lineage_filter`` is
        configured or the seed taxid is unknown.

        Args:
            seed_id: Integer protein ID.

        Returns:
            Ceiling taxid string (or ``PROKARYOTA_SYNTHETIC`` / ``"root"``).
        """
        if self.lineage_filter is None:
            return "root"
        seed_taxid = self._get_seed_taxid_str(seed_id)
        if not seed_taxid or seed_taxid == "0":
            return "root"
        return self.lineage_filter.get_ceiling_for_seed(seed_taxid)

    def _compute_farthest_donor(
        self,
        ortholog_meta: Dict[int, Dict[str, Any]],
        seed_id: int,
    ) -> Tuple[str, str]:
        """Return (taxid_str, lineage_str) of the farthest used donor.

        "Farthest" = the donor whose ev_lca has the smallest depth
        (shallowest in the tree = broadest LCA = most evolutionarily
        distant from the seed).  Ties are broken arbitrarily.

        Falls back to the seed's own species values when no non-self
        donor orthologs are present in ``ortholog_meta``.

        Args:
            ortholog_meta: Mapping from ortholog protein ID to cascade
                metadata dict (as produced by ``_collect_orthologs``).
            seed_id: Integer protein ID of the seed.

        Returns:
            A ``(taxid_str, lineage_str)`` tuple where ``lineage_str``
            is semicolon-separated from root to leaf.
        """
        taxids = self.db.taxid_array
        best_depth: Optional[int] = None
        best_taxid_str: str = ""

        for oid, meta in ortholog_meta.items():
            if oid == seed_id:
                continue
            depth = meta.get("depth", 0)
            if best_depth is None or depth < best_depth:
                best_depth = depth
                if taxids is not None and oid < len(taxids):
                    best_taxid_str = str(taxids[oid])

        if best_taxid_str and best_taxid_str != "0":
            return best_taxid_str, self._lineage_str_for_taxid(best_taxid_str)

        # Fallback to seed's own species.
        seed_taxid_str = self._get_seed_taxid_str(seed_id)
        if seed_taxid_str and seed_taxid_str != "0":
            return seed_taxid_str, self._lineage_str_for_taxid(seed_taxid_str)
        return "-", "-"

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
        if len(self._seed_lineage_set_cache) > self._LINEAGE_CACHE_MAX:
            evict_n = len(self._seed_lineage_set_cache) // 2
            for k in list(self._seed_lineage_set_cache)[:evict_n]:
                del self._seed_lineage_set_cache[k]
        return result

    def _resolve_valid_species_for_ceiling(
        self,
        seed_id: int,
        ceiling_taxid: str,
    ) -> Optional[frozenset]:
        """Return the valid species frozenset for a pre-resolved ceiling.

        Memoized by ceiling taxid (multiple seeds may share the same
        ceiling on a single-species proteome run).

        Args:
            seed_id: Integer protein ID (used only for debug logging).
            ceiling_taxid: Pre-resolved ceiling from :meth:`_resolve_ceiling`.

        Returns:
            Frozenset of species taxid strings, or ``None`` for no filter.
        """
        if self.lineage_filter is None:
            return None
        if ceiling_taxid == "root":
            return None
        if ceiling_taxid in self._valid_species_by_seed:
            return self._valid_species_by_seed[ceiling_taxid]
        valid = self.lineage_filter.get_valid_species_ids(ceiling_taxid)
        self._valid_species_by_seed[ceiling_taxid] = valid
        return valid

    def _resolve_allowed_og_lcas_for_ceiling(
        self, ceiling_taxid: str
    ) -> Optional[frozenset]:
        """Return the ``og_lca`` whitelist for a pre-resolved ceiling.

        Used as the strict-OG gate in ``_collect_orthologs``.  Memoized
        per ceiling taxid.

        Args:
            ceiling_taxid: Pre-resolved ceiling from :meth:`_resolve_ceiling`.

        Returns:
            Frozenset of taxid strings (ceiling + all descendants), or
            ``None`` for no filter.
        """
        if self.lineage_filter is None:
            return None
        if ceiling_taxid == "root":
            return None
        cache = getattr(self, "_allowed_og_lcas_by_ceiling", None)
        if cache is None:
            cache = self._allowed_og_lcas_by_ceiling = {}
        if ceiling_taxid in cache:
            return cache[ceiling_taxid]
        allowed = self.lineage_filter.get_scope_og_descendants(ceiling_taxid)
        cache[ceiling_taxid] = allowed
        return allowed

    def _collect_orthologs(
        self,
        seed_id: int,
        events: Dict[int, dict],
        target_taxa: Optional[Set[int]] = None,
        excluded_taxa: Optional[Set[int]] = None,
        valid_species: Optional[frozenset] = None,
        ceiling_taxid: Optional[str] = None,
        allowed_og_lcas: Optional[frozenset] = None,
    ) -> Tuple[Set[int], Dict[str, Set[int]], Dict[int, Dict[str, Any]]]:
        """Collect ortholog IDs from events, classified by orthology type.

        For each speciation event where the seed appears, the proteins on the
        opposite side are collected and classified by the cardinality of each
        side *before* species filtering. A protein can appear in multiple typed
        buckets if it participates in events with different cardinalities.

        Args:
            seed_id: Seed protein ID.
            events: Dict of event_id → event data.
            target_taxa: Optional taxids to include.
            excluded_taxa: Optional taxids to exclude.
            valid_species: Optional frozenset of taxid strings allowed by
                the ceiling resolver.  Produced by
                :meth:`_resolve_valid_species_for_ceiling`.
            ceiling_taxid: Pre-resolved ceiling taxid string (or
                ``"root"``).  When provided, events whose ``ev_lca`` does
                not pass the ceiling are discarded before orthologs are
                collected.  ``None`` behaves as ``"root"`` (no ev_lca
                filter).
            allowed_og_lcas: Optional frozenset of taxid strings naming
                the OGs whose containing taxonomic level (``og_lca``) is
                at-or-below the ceiling.  When provided, events from
                broader OGs are skipped entirely — the strict-OG gate.

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

        # Resolve the ev_lca ceiling check once per call (not per event).
        effective_ceiling = ceiling_taxid or "root"
        do_ev_lca_filter = (
            effective_ceiling != "root"
            and self.lineage_filter is not None
        )

        for event_id, event in events.items():
            # Hard scope-OG filter: drop events whose containing OG is
            # broader than the seed's tax-scope ceiling. The same in-scope
            # orthologs already appear in lower OGs at higher cascade
            # priority, so dropping these events is biologically lossless
            # under auto-scope and removes ~99 % of fetched orthologs on
            # plant proteomes.
            #
            # Empty / missing `og_lca` (rare — events whose source tree
            # node had no LCA prop): keep the event through this gate
            # rather than silently drop. Out-of-scope orthologs from
            # such events are still excluded by the per-protein species
            # filter (`valid_species`) below. Dropping them would mean
            # silent biological data loss.
            if allowed_og_lcas is not None:
                og_lca = event.get("og_lca")
                if og_lca and og_lca not in allowed_og_lcas:
                    continue

            # ev_lca ceiling gate: drop events whose ev_lca is above the
            # per-seed ceiling (i.e. the LCA taxid is not at-or-below the
            # ceiling in the NCBI tree).  This is the primary new filter
            # replacing the old `allowed_og_lcas` strict-OG approach.
            if do_ev_lca_filter:
                ev_lca_for_check = str(event.get("ev_lca") or "")
                if ev_lca_for_check and not self.lineage_filter.ev_lca_passes_ceiling(
                    ev_lca_for_check, effective_ceiling
                ):
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
            # og_lca: the LCA taxid of the OG that produced this event.
            # Threaded through to the cascade so the GO sub-namespace cross-
            # tier union can scope-cap contributions to the seed's auto-scope
            # ceiling regardless of the global `scope_strict_og` flag.
            og_lca_str = str(event.get("og_lca") or "")
            payload = {
                "event_id": event_id,
                "ev_lca": ev_lca,
                "og_lca": og_lca_str,
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
        "all":       frozenset({"self", "one2one", "one2many", "many2one", "many2many"}),
        "many2many": frozenset({"self", "one2one", "one2many", "many2one", "many2many"}),
        "one2many":  frozenset({"self", "one2one", "one2many"}),
        "many2one":  frozenset({"self", "one2one", "many2one"}),
        "one2one":   frozenset({"self", "one2one"}),
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
        keeps only ``GO:`` prefixed terms) and then split by GO namespace
        into ``gos_mf`` / ``gos_bp`` / ``gos_cc`` so each subnamespace
        runs an independent cascade. Terms whose namespace is not in
        the OBO map are dropped (a single debug-level summary is logged
        per batch — see :meth:`_split_gos_by_namespace`). When the OBO
        map is unavailable the legacy combined ``gos`` field is emitted
        instead and the cascade falls back to the historical behaviour.
        Every other field becomes a distinct-preserving tuple — same
        insertion order, deduplicated.
        """
        parsed: Dict[int, Dict[str, Tuple[str, ...]]] = {}
        ns_map = self._go_namespace_map
        unmapped_count = 0
        for oid, annot in annot_data.items():
            if not annot:
                parsed[oid] = {}
                continue
            row: Dict[str, Tuple[str, ...]] = {}
            for field in self.ANNOTATION_FIELDS:
                # gos_mf / gos_bp / gos_cc come from the single source
                # field "gos" — handled outside the per-field loop below.
                if field in _GO_NS_FIELDS:
                    continue
                raw = annot.get(field)
                if not raw:
                    continue
                if field == "pname":
                    s = raw.strip()
                    if s:
                        row[field] = (s,)
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
            # GO terms — single parse, then namespace split if possible.
            raw_gos = annot.get(self.LEGACY_GO_FIELD)
            if raw_gos:
                gos = tuple(self._parse_gos(raw_gos))
                if gos:
                    if ns_map is None:
                        # Legacy combined-GO fallback: keep the single
                        # "gos" field; cascade walk uses it directly.
                        row[self.LEGACY_GO_FIELD] = gos
                    else:
                        mf: List[str] = []
                        bp: List[str] = []
                        cc: List[str] = []
                        for go in gos:
                            ns = ns_map.get(go)
                            if ns == "molecular_function":
                                mf.append(go)
                            elif ns == "biological_process":
                                bp.append(go)
                            elif ns == "cellular_component":
                                cc.append(go)
                            else:
                                # Unmapped GO (obsolete, alt namespace, or
                                # newer than the OBO snapshot). We drop
                                # rather than mis-bucket — over-counting
                                # any namespace would skew the per-ns
                                # short-circuit. Rare in practice (<<1 %
                                # on a current OBO release).
                                unmapped_count += 1
                        if mf:
                            row["gos_mf"] = tuple(mf)
                        if bp:
                            row["gos_bp"] = tuple(bp)
                        if cc:
                            row["gos_cc"] = tuple(cc)
            parsed[oid] = row
        if unmapped_count and ns_map is not None:
            logger.debug(
                "_pre_parse_batch: dropped %d GO terms with unknown namespace",
                unmapped_count,
            )
        return parsed

    def _summarize_annotations(
        self,
        annot_data: Dict[int, dict],
        ortholog_meta: Optional[Dict[int, dict]] = None,
        target_orthologs: str = "all",
        parsed: Optional[Dict[int, Dict[str, Tuple[str, ...]]]] = None,
        scope_og_lcas: Optional[frozenset] = None,
        donor_pool: str = "closest",
    ) -> Tuple[Dict[str, Any], Dict[str, str]]:
        """Cascade summary: per-source closest-ev_lca + type-priority winner.

        For each functional source (KEGG_ko, GOs, Pfam, …) independently,
        donors are walked from closest + best-typed first.

        ``donor_pool`` controls the walk strategy:

        - ``"closest"`` (default): The first bucket — defined by the
          cascade key ``(in_seed_lineage, -ev_lca_depth, type_tier)`` —
          that has any donor with a non-empty value for that source wins,
          and the consensus is taken across *only* the donors in that
          winning bucket.  This preserves the original semantics.
        - ``"union"``: Walk *all* tiers; union the values across every
          tier that contributes any donor for each source.  The confidence
          reported is the best (smallest ``type_tier``) seen across the
          contributing tiers.

        ``target_orthologs`` is a *floor*: types not in
        ``TARGET_ORTHOLOGS_FLOORS[target_orthologs]`` are excluded from
        the cascade entirely.

        Args:
            annot_data: ``{protein_id: annotation_dict}`` for every
                ortholog plus optionally the seed.
            ortholog_meta: ``{protein_id: meta_dict}`` produced by
                :meth:`_collect_orthologs`.  Required for cascade mode.
                When absent or empty the flat aggregation is used (same
                logic as v3 phase 0) and confidence is empty.
            target_orthologs: ``"all"``, ``"many2many"``, ``"one2many"``,
                ``"many2one"`` or ``"one2one"``.  Anything else is treated
                as ``"all"``.
            parsed: Optional pre-parsed annotation cache from
                :meth:`_pre_parse_batch`.
            scope_og_lcas: Unused; kept for call-site compatibility only.
            donor_pool: ``"closest"`` or ``"union"``; see above.

        Returns:
            ``(annotations, confidence)`` where ``annotations`` is the
            same shape returned by the legacy summariser (output-named
            keys: ``Preferred_name``, ``GOs``, ``PFAMs``, ``KEGG_ko``,
            …) and ``confidence`` is ``{output_field: "high" | "medium"
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
        # Per-output-field accumulator for sources that share an output
        # key (currently only the three GO sub-namespaces, which all
        # merge into "GOs"). We accumulate the union of values and keep
        # the *best* (lowest-tier-int) confidence among contributing
        # sub-namespaces, so the GO row ends up with the strongest
        # confidence of the three independent cascade winners.
        merged_values: Dict[str, Set[str]] = {}
        merged_tier: Dict[str, int] = {}
        # Output fields that aggregate from more than one source. Used to
        # branch the assignment below from "overwrite" to "union+merge".
        multi_source_outputs = {self._cascade_output_field(f) for f in _GO_NS_FIELDS}

        use_union = donor_pool == "union"

        for field in self.ANNOTATION_FIELDS:
            # Skip the per-namespace GO sub-fields entirely when the OBO
            # map was unavailable: in that case `_pre_parse_batch` writes
            # the legacy combined "gos" field instead, which is appended
            # to the walk list below.
            if field in _GO_NS_FIELDS and self._go_namespace_map is None:
                continue
            output_field = self._cascade_output_field(field)
            # Both modes walk tiers in priority order. "closest" stops at
            # the first non-empty tier (the original behaviour). "union"
            # continues through all tiers and accumulates values.
            field_union_values: Set[str] = set()
            field_best_tier: Optional[int] = None

            for prio_key in priority_order:
                contributors = [
                    oid for oid in buckets[prio_key]
                    if field in parsed.get(oid, ())
                ]
                if not contributors:
                    continue
                values = self._aggregate_field(field, contributors, parsed)
                if not values:
                    continue

                if output_field in multi_source_outputs:
                    # GO sub-namespaces always union into "GOs" regardless
                    # of donor_pool, to preserve per-namespace completeness.
                    bucket = merged_values.setdefault(output_field, set())
                    bucket.update(values)
                    cur_tier = merged_tier.get(output_field)
                    if cur_tier is None or prio_key[2] < cur_tier:
                        merged_tier[output_field] = prio_key[2]
                    if not use_union:
                        break
                elif use_union:
                    # Accumulate all tiers; record the best (smallest) tier.
                    field_union_values.update(values)
                    if field_best_tier is None or prio_key[2] < field_best_tier:
                        field_best_tier = prio_key[2]
                else:
                    # "closest": first non-empty bucket wins.
                    annotations[output_field] = values
                    confidence[output_field] = self.TIER_CONFIDENCE[prio_key[2]]
                    break

            # Materialise union-mode non-GO fields.
            if use_union and field_union_values and output_field not in multi_source_outputs:
                annotations[output_field] = list(field_union_values)
                confidence[output_field] = self.TIER_CONFIDENCE[
                    field_best_tier  # type: ignore[index]
                ]

        # Legacy combined-GO walk when the OBO map was unavailable. The
        # three sub-fields above were no-ops; "gos" runs the historical
        # single-tier short-circuit and emits one "GOs" output.
        if self._go_namespace_map is None:
            field = self.LEGACY_GO_FIELD
            output_field = self._cascade_output_field(field)
            gos_union: Set[str] = set()
            gos_best_tier: Optional[int] = None
            for prio_key in priority_order:
                contributors = [
                    oid for oid in buckets[prio_key]
                    if field in parsed.get(oid, ())
                ]
                if not contributors:
                    continue
                values = self._aggregate_field(field, contributors, parsed)
                if values:
                    if use_union:
                        gos_union.update(values)
                        if gos_best_tier is None or prio_key[2] < gos_best_tier:
                            gos_best_tier = prio_key[2]
                    else:
                        annotations[output_field] = values
                        confidence[output_field] = self.TIER_CONFIDENCE[prio_key[2]]
                        break
            if use_union and gos_union:
                annotations[output_field] = list(gos_union)
                confidence[output_field] = self.TIER_CONFIDENCE[
                    gos_best_tier  # type: ignore[index]
                ]

        # Materialise multi-source merged outputs (currently just GOs).
        for output_field, values_set in merged_values.items():
            if values_set:
                annotations[output_field] = list(values_set)
                confidence[output_field] = self.TIER_CONFIDENCE[
                    merged_tier[output_field]
                ]

        return annotations, confidence

    def _cascade_output_field(self, field: str) -> str:
        """Map an internal annotation field name to the output key the
        cascade emits.

        Note: ``gos_mf``, ``gos_bp``, ``gos_cc`` and the legacy combined
        ``gos`` all map to the same output key ``"GOs"`` — the cascade
        engine merges them via :meth:`_summarize_annotations`.
        """
        if field == "pname":
            return "Preferred_name"
        if field == "gos" or field in _GO_NS_FIELDS:
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
        v3-phase0-end exactly.

        Operates on raw DB rows whose GO terms are stored under the
        single key ``"gos"`` — the per-namespace split only happens in
        the cascade path. The flat path therefore iterates over the
        legacy field set rather than ``ANNOTATION_FIELDS`` (which holds
        the per-namespace cascade keys).
        """
        # Replace the three per-namespace sub-keys with the legacy
        # combined "gos" field for raw-row reads.
        legacy_fields = [
            self.LEGACY_GO_FIELD if f == _GO_NS_FIELDS[0] else f
            for f in self.ANNOTATION_FIELDS
            if f == _GO_NS_FIELDS[0] or f not in _GO_NS_FIELDS
        ]

        counters = defaultdict(Counter)
        for annot in annot_data.values():
            if not annot:
                continue
            for field in legacy_fields:
                value = annot.get(field)
                if not value:
                    continue
                if field == "pname":
                    counters[field][value.strip()] += 1
                elif field == self.LEGACY_GO_FIELD:
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
            elif field == self.LEGACY_GO_FIELD:
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
