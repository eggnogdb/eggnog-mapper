"""Batch annotation for v7+ integer-encoded eggnog.db.

Thin adapter over eggnog_annotator.e7.AnnotationEngine: emapper parses hits
from DIAMOND/MMseqs, calls the engine for all annotation logic (orthologs,
annotations, OG resolution, taxonomic filters, caches), and re-shapes the
result into the tuple that emapper's output formatter consumes.

All actual annotation logic lives in eggnog-annotator. This module owns
only the hit parsing and the emapper-specific tuple packaging.
"""

import logging

from eggnog_annotator.e7 import AnnotationEngine, EggnogDB, LineageFilter

from ..emapperException import EmapperException

logger = logging.getLogger(__name__)

from . import output as output_mod


def filter_out(hit_name, hit_evalue, hit_score, threshold_evalue, threshold_score):
    """Drop a hit before annotation if the seed name is sentinel or it
    fails score/e-value thresholds. Inlined from the now-deleted
    annotator_worker module (Phase 2 commit 3)."""
    if hit_name == "-" or hit_name == "ERROR":
        return True
    if threshold_evalue is not None and hit_evalue is not None and hit_evalue > threshold_evalue:
        return True
    if threshold_score is not None and hit_score is not None and hit_score < threshold_score:
        return True
    return False

ANNOTATIONS_HEADER = output_mod.ANNOTATIONS_HEADER

# Module-level lineage cache (loaded once, shared across batches)
_lineage_cache = None

# Module-level AnnotationEngine singleton, keyed by db connection id.
# Reusing the engine across batches preserves its OG-description cache.
_engine = None
_engine_conn_id = None
_engine_scope_key = None


def get_lineage_cache():
    """Get or create the module-level lineage cache."""
    global _lineage_cache
    if _lineage_cache is None:
        from ..common import get_ncbitaxadb_file
        from .tax_scopes.v7_scope import LineageCache
        taxa_db = get_ncbitaxadb_file()
        try:
            _lineage_cache = LineageCache(taxa_db_path=taxa_db)
        except Exception as exc:
            raise EmapperException(
                f"Failed to load lineage cache from {taxa_db}: {exc}. "
                "Ensure eggnog.taxa.db is present in --data_dir."
            ) from exc
    return _lineage_cache


def _build_lineage_filter(v7_tax_scope, v7_tax_scope_auto, taxid_array):
    """Build a LineageFilter for scope-based ortholog pruning.

    Returns None when no filtering was requested. Without this filter
    the engine degrades to collecting every ortholog across all events,
    which is orders of magnitude more work than needed for a single
    taxonomic domain.
    """
    if not v7_tax_scope_auto and not v7_tax_scope:
        return None
    lf = LineageFilter(get_lineage_cache(), taxid_array=taxid_array)
    if v7_tax_scope_auto:
        lf.set_scope("auto")
    else:
        lf.set_scope(",".join(sorted(v7_tax_scope)))
    return lf


def _get_engine(annot_db, v7_tax_scope, v7_tax_scope_auto):
    """Return the singleton AnnotationEngine bound to this DB connection.

    Reuses emapper's already-open sqlite3 connection and its pre-loaded
    taxid array, so no memory is duplicated. Recreated if scope changes.
    """
    global _engine, _engine_conn_id, _engine_scope_key
    # scope_key distinguishes: no scope (None), auto scope ("auto"),
    # and explicit scope (sorted tuple of taxid strings).
    # Empty v7_tax_scope list maps to None (same as no scope) — correct,
    # since an empty explicit scope list is effectively unconstrained.
    scope_key = (
        "auto" if v7_tax_scope_auto
        else tuple(sorted(v7_tax_scope)) if v7_tax_scope
        else None
    )
    if (_engine is None
            or _engine_conn_id != id(annot_db.conn)
            or _engine_scope_key != scope_key):
        db = EggnogDB.from_connection(annot_db.conn, taxid_array=annot_db._taxids)
        lf = _build_lineage_filter(v7_tax_scope, v7_tax_scope_auto, annot_db._taxids)
        _engine = AnnotationEngine(db, lineage_filter=lf)
        _engine_conn_id = id(annot_db.conn)
        _engine_scope_key = scope_key
    return _engine


def annotate_batch(batch, eggnog_db, annot, target_orthologs,
                   target_taxa, excluded_taxa, tax_scope_mode,
                   tax_scope_ids, go_evidence, go_excluded,
                   seed_ortholog_score, seed_ortholog_evalue,
                   pool=None, v7_tax_scope=None, v7_tax_scope_auto=False):
    """Annotate a batch of hits using eggnog-annotator.

    batch: list of (hit, ...) argument tuples from iter_hit_lines
    pool: ignored (kept for signature compatibility)
    Yields ((hit, annotation), False) tuples, same interface as the
    legacy per-hit annotator.
    """
    # Parse hits, drop those that fail score/evalue thresholds
    valid_hits = []
    for args in batch:
        hit = args[0]
        best_hit_name = hit[1]
        best_hit_evalue = float(hit[2])
        best_hit_score = float(hit[3])

        if filter_out(best_hit_name, best_hit_evalue, best_hit_score,
                       seed_ortholog_evalue, seed_ortholog_score):
            yield ((hit, None), False)
            continue

        # v7 int_mode: DIAMOND/MMseqs emit integer IDs as strings
        try:
            seed_id = int(best_hit_name)
        except ValueError:
            logger.warning(
                "Cannot parse integer seed ID from '%s'; skipping hit"
                " (DB version mismatch?)",
                best_hit_name,
            )
            yield ((hit, None), False)
            continue

        valid_hits.append((hit, hit[0], seed_id, best_hit_evalue, best_hit_score))

    if not valid_hits:
        return

    engine = _get_engine(eggnog_db, v7_tax_scope, v7_tax_scope_auto)
    seed_ids = [s for _, _, s, _, _ in valid_hits]
    # `target_orthologs` is now a *floor* on the cascade — the engine
    # restricts which donor types contribute annotations. The orthologs
    # list shown in the .emapper.orthologs file is still derived from the
    # full classified `ortholog_types` dict and post-filtered there, so
    # the user's requested target is respected end-to-end.
    results = engine.annotate_batch(
        seed_ids,
        target_taxa=target_taxa,
        excluded_taxa=excluded_taxa,
        target_orthologs=target_orthologs,
    )

    # Re-shape engine results into the tuple consumed by output.py
    for hit, query_name, seed_id, evalue, score in valid_hits:
        r = results.get(seed_id)
        if not r or (not r.get("orthologs") and not r.get("annotations")):
            yield ((hit, None), False)
            continue

        annotations = r.get("annotations") or {}
        annotations_confidence = r.get("annotations_confidence") or {}
        tax_scope_used = r.get("tax_scope_used") or "none"
        og_info = r.get("og_info") or {}
        orthologs = r.get("orthologs") or []

        og_name = og_info.get("name")
        og_cat = og_info.get("cog_cat", "-")
        og_desc = og_info.get("description", "-")
        max_annot_lvl = og_info.get("level", "-")
        match_nog_names = r.get("all_ogs") or ([og_name] if og_name else [])
        all_orthologies = r.get("ortholog_types", {"all": set(orthologs)})
        annot_orthologs = sorted(
            all_orthologies.get(target_orthologs, all_orthologies.get("all", []))
        )

        annotation = (
            query_name,
            str(seed_id),
            evalue,
            score,
            annotations,
            (og_name, og_cat, og_desc),
            max_annot_lvl,
            match_nog_names,
            all_orthologies,
            annot_orthologs,
            annotations_confidence,
            tax_scope_used,
        )

        yield ((hit, annotation), False)
