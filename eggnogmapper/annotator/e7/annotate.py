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

    def __init__(self, db: EggnogDB, lineage_filter=None):
        """Initialize annotation engine.

        Args:
            db: EggnogDB instance with loaded taxid array
            lineage_filter: Optional LineageFilter for tax-scope filtering.
                When set (typically to auto mode), orthologs are pruned
                per-seed to the seed's taxonomic domain before expensive
                annotation fetch — the single largest speedup for
                proteome-scale runs.
        """
        self.db = db
        self.lineage_filter = lineage_filter
        # Cross-batch cache for OG metadata. Most eukaryotic proteomes
        # reuse a small set of common OGs across thousands of proteins,
        # so caching them avoids repeated SQL + JSON parsing.
        self._og_cache: Dict[Tuple[str, str], Optional[dict]] = {}
        # Cache of valid-species sets keyed by the seed's species taxid.
        # Reused across seeds that share the same auto-scope outcome
        # (e.g. an all-Arabidopsis proteome => one scope computed once).
        self._valid_species_by_seed: Dict[str, Optional[set]] = {}

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
        orthologs, ortholog_types = self._collect_orthologs(
            seed_id, events, target_taxa, excluded_taxa
        )

        # Fetch and summarize annotations
        annotations = {}
        if orthologs:
            annot_data = self.db.get_protein_annotations_bulk(list(orthologs))
            annotations = self._summarize_annotations(annot_data)

        # Get OG info from events
        og_info = self._get_og_info(events)

        return {
            "orthologs": list(orthologs),
            "annotations": annotations,
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
    ) -> Dict[int, Dict[str, Any]]:
        """Annotate multiple seed orthologs efficiently.

        Uses bulk queries to minimize database round-trips.

        Args:
            seed_ids: List of integer protein IDs
            target_taxa: Optional set of taxids to include
            excluded_taxa: Optional set of taxids to exclude

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
        all_orthologs: Set[int] = set()
        for seed_id in seed_ids:
            eids = event_index.get(seed_id, [])
            # Filter events to those in our cache
            seed_events = {eid: events_cache[eid] for eid in eids if eid in events_cache}
            valid_species = self._resolve_valid_species(seed_id)
            orthologs, ortholog_types = self._collect_orthologs(
                seed_id, seed_events, target_taxa, excluded_taxa, valid_species,
            )
            seed_orthologs[seed_id] = orthologs
            seed_ortholog_types[seed_id] = ortholog_types
            all_orthologs.update(orthologs)
        _t(f"p3 collect_orthologs ({len(all_orthologs)})", t0)

        # Phase 4: Bulk fetch annotations for all orthologs
        t0 = time.time()
        annot_cache = {}
        if all_orthologs:
            annot_cache = self.db.get_protein_annotations_bulk(list(all_orthologs))
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

            # Filter annot_cache to this seed's orthologs
            seed_annots = {oid: annot_cache[oid] for oid in orthologs if oid in annot_cache}
            annotations = self._summarize_annotations(seed_annots) if seed_annots else {}

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

    def _collect_orthologs(
        self,
        seed_id: int,
        events: Dict[int, dict],
        target_taxa: Optional[Set[int]] = None,
        excluded_taxa: Optional[Set[int]] = None,
        valid_species: Optional[frozenset] = None,
    ) -> Tuple[Set[int], Dict[str, Set[int]]]:
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

        Returns:
            (orthologs, ortholog_types) where:
            - orthologs: filtered set of all ortholog protein IDs
            - ortholog_types: dict with keys "one2one", "one2many",
              "many2one", "many2many", "all" — each mapping to the subset
              of filtered orthologs with that relationship to the seed.
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

        for event in events.values():
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

        # Remove the seed itself from every candidate bucket
        for bucket in candidates.values():
            bucket.discard(seed_id)

        all_candidates: Set[int] = set()
        for bucket in candidates.values():
            all_candidates.update(bucket)

        if not all_candidates:
            return typed["all"], typed

        taxids = self.db.taxid_array

        # No taxid array → return all candidates unfiltered
        if taxids is None:
            for rel, bucket in candidates.items():
                typed[rel].update(bucket)
                typed["all"].update(bucket)
            return typed["all"], typed

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

        return typed["all"], typed

    def _summarize_annotations(
        self,
        annot_data: Dict[int, dict],
    ) -> Dict[str, Any]:
        """Summarize annotations from multiple orthologs.

        Aggregates annotations across orthologs using frequency counting.

        Args:
            annot_data: Dict of protein_id -> annotation dict

        Returns:
            Dict of annotation field -> aggregated values
        """
        counters = defaultdict(Counter)

        for annot in annot_data.values():
            if not annot:
                continue

            for field in self.ANNOTATION_FIELDS:
                value = annot.get(field)
                if not value:
                    continue

                if field == "pname":
                    # Preferred name: count occurrences
                    counters[field][value.strip()] += 1
                elif field == "gos":
                    # GO terms: parse and count
                    for go in self._parse_gos(value):
                        counters[field][go] += 1
                else:
                    # Other fields: split by comma and count
                    for v in str(value).split(","):
                        v = v.strip()
                        if v:
                            counters[field][v] += 1

        # Build result
        result = {}
        for field, counter in counters.items():
            if field == "pname":
                # Most common preferred name (require at least 2 occurrences).
                # Wrapped in a list so downstream ",".join(sorted(list(...)))
                # works the same as for the other multi-value fields.
                most_common = counter.most_common(1)
                if most_common and most_common[0][1] >= 2:
                    result["Preferred_name"] = [most_common[0][0]]
            elif field == "gos":
                result["GOs"] = list(counter.keys())
            elif field == "pfam":
                # Filter PFAMs by frequency
                total = len(annot_data)
                result["PFAMs"] = [k for k, c in counter.items()
                                   if c > 1 and c / total > 0.05]
            else:
                # Map field names to output format
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
