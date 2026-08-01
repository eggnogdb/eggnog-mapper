"""Batch annotation for v7+ integer-encoded eggnog.db.

Thin adapter over eggnogmapper.annotator.e7.AnnotationEngine: emapper parses
hits from DIAMOND/MMseqs, calls the engine for all annotation logic
(orthologs, annotations, OG resolution, taxonomic filters, caches), and
re-shapes the result into the tuple that emapper's output formatter consumes.

All actual annotation logic lives in eggnog-annotator. This module owns
only the hit parsing and the emapper-specific tuple packaging.
"""

import logging
import os
import sys

from eggnogmapper.annotator.e7 import AnnotationEngine, EggnogDB, LineageFilter
from eggnogmapper.annotator.e7.ceiling import TaxScopeCeilingResolver

from ..emapperException import EmapperException
from ..common import get_data_path
from ..utils import colorify

logger = logging.getLogger(__name__)

from . import output as output_mod


# Relative locations (under the resolved --data_dir) where go-basic.obo is
# looked up, in order. The flat location matches "shipped next to the DB
# files"; the source-subtree matches the builder/dev layout.
_GO_OBO_REL_CANDIDATES = (
    "go-basic.obo",
    os.path.join("e7", "full", "source", "reference", "go-basic.obo"),
)


def resolve_go_obo_path():
    """Resolve the go-basic.obo path for the current run.

    Precedence: the ``EGGNOG_GO_OBO`` env var (explicit override) wins;
    otherwise look under the resolved data dir (``--data_dir``) at the
    candidate locations and return the first that exists. When none exists,
    return the primary data-dir location anyway so the engine's "not found"
    warning points the user at where to drop the file.

    Returns:
        An absolute path string, or ``None`` when the data dir is unknown
        (the engine then falls back to its own default / legacy GO cascade).
    """
    env = os.environ.get("EGGNOG_GO_OBO")
    if env:
        return env
    try:
        data_dir = get_data_path()
    except Exception:
        return None
    if not data_dir:
        return None
    candidates = [os.path.join(data_dir, rel) for rel in _GO_OBO_REL_CANDIDATES]
    for cand in candidates:
        if os.path.isfile(cand):
            return cand
    return candidates[0]


def filter_out(
    hit_name: str,
    hit_evalue: float,
    hit_score: float,
    threshold_evalue: float,
    threshold_score: float,
) -> None | str:
    """Decide whether a hit should be dropped before annotation.

    Returns ``None`` when the hit survives, otherwise a short reason
    string used by the optional drop log (Phase 7.1c).  The reasons are
    stable identifiers so users can grep/filter the .dropped file
    programmatically.

    Args:
        hit_name: Seed ortholog name string (``"-"`` or ``"ERROR"`` = drop).
        hit_evalue: E-value of the DIAMOND/MMseqs hit.
        hit_score: Bit-score of the DIAMOND/MMseqs hit.
        threshold_evalue: Maximum allowed E-value (``None`` = no threshold).
        threshold_score: Minimum required bit-score (``None`` = no threshold).

    Returns:
        ``None`` if the hit passes, or a reason string if it should be dropped.
    """
    if hit_name == "-" or hit_name == "ERROR":
        return "error_seed"
    if (
        threshold_evalue is not None
        and hit_evalue is not None
        and hit_evalue > threshold_evalue
    ):
        return "evalue_above_threshold"
    if (
        threshold_score is not None
        and hit_score is not None
        and hit_score < threshold_score
    ):
        return "score_below_threshold"
    return None


ANNOTATIONS_HEADER = output_mod.ANNOTATIONS_HEADER

# Module-level lineage cache (loaded once, shared across batches)
_lineage_cache = None

# Module-level AnnotationEngine singleton, keyed by db connection id.
# Reusing the engine across batches preserves its OG-description cache.
_engine = None
_engine_conn_id = None
_engine_scope_key = None


def get_lineage_cache():
    """Get or create the module-level lineage cache.

    Returns:
        Loaded ``LineageCache`` instance backed by the project taxa DB.

    Raises:
        EmapperException: If the taxa DB cannot be found or loaded.
    """
    global _lineage_cache
    if _lineage_cache is None:
        from ..common import get_ncbitaxadb_file
        from eggnogmapper.annotator.lineage import LineageCache

        taxa_db = get_ncbitaxadb_file()
        try:
            _lineage_cache = LineageCache(taxa_db_path=taxa_db)
        except Exception as exc:
            raise EmapperException(
                f"Failed to load lineage cache from {taxa_db}: {exc}. "
                "Ensure eggnog.taxa.db is present in --data_dir."
            ) from exc
    return _lineage_cache


def _build_lineage_filter(
    ceiling_resolver: TaxScopeCeilingResolver,
    taxid_array,
) -> LineageFilter:
    """Wrap a ceiling resolver in a LineageFilter.

    Args:
        ceiling_resolver: Pre-built ``TaxScopeCeilingResolver`` instance.
        taxid_array: Protein-ID → species-taxid mapping from the DB.

    Returns:
        A ``LineageFilter`` bound to the same lineage cache and resolver.
    """
    return LineageFilter(
        lineage_cache=ceiling_resolver.lineage_cache,
        ceiling_resolver=ceiling_resolver,
        taxid_array=taxid_array,
    )


def _get_engine(
    annot_db,
    ceiling_resolver: TaxScopeCeilingResolver,
    donor_pool: str,
    lazy_cascade: bool = False,
) -> AnnotationEngine:
    """Return the singleton AnnotationEngine bound to this DB connection.

    Reuses emapper's already-open sqlite3 connection and its pre-loaded
    taxid array, so no memory is duplicated.  Recreated when the
    ceiling resolver mode, donor_pool or lazy_cascade toggle changes.

    Args:
        annot_db: Open annotation DB object with ``.conn`` and
            ``._taxids`` attributes.
        ceiling_resolver: Pre-built ``TaxScopeCeilingResolver``.
        donor_pool: ``"closest"`` or ``"union"``.
        lazy_cascade: Enable the tier-staged lazy cascade (closest only).

    Returns:
        Singleton ``AnnotationEngine`` instance.
    """
    global _engine, _engine_conn_id, _engine_scope_key
    scope_key = (ceiling_resolver.mode, donor_pool, bool(lazy_cascade))
    if (
        _engine is None
        or _engine_conn_id != id(annot_db.conn)
        or _engine_scope_key != scope_key
    ):
        db = EggnogDB.from_connection(annot_db.conn, taxid_array=annot_db._taxids)
        lf = _build_lineage_filter(ceiling_resolver, annot_db._taxids)
        _engine = AnnotationEngine(
            db,
            lineage_filter=lf,
            donor_pool=donor_pool,
            lazy_cascade=lazy_cascade,
            go_obo_path=resolve_go_obo_path(),
        )
        # Build the per-protein field-presence mask once, up front, so it is
        # COW-shared across fork workers (like the taxid array). On failure the
        # engine silently uses the byte-identical tier-staged fallback.
        if lazy_cascade and donor_pool == "closest":
            _engine.load_field_presence()
        _engine_conn_id = id(annot_db.conn)
        _engine_scope_key = scope_key
    return _engine


def annotate_batch(
    batch,
    eggnog_db,
    annot: bool,
    target_orthologs: str,
    target_taxa,
    excluded_taxa,
    seed_ortholog_score,
    seed_ortholog_evalue,
    ceiling_resolver: TaxScopeCeilingResolver,
    donor_pool: str = "closest",
    lazy_cascade: bool = False,
    pool=None,
    dropped_writer=None,
    sub_batch_size: int = 125,
):
    """Annotate a batch of hits using eggnog-annotator.

    Args:
        batch: List of ``(hit, ...)`` argument tuples from
            ``iter_hit_lines``.
        eggnog_db: Open annotation DB object.
        annot: Whether annotation is requested.
        target_orthologs: Cascade type floor (``"all"``, ``"one2one"``,
            etc.).
        target_taxa: Optional set of taxids to include.
        excluded_taxa: Optional set of taxids to exclude.
        seed_ortholog_score: Minimum bit-score threshold.
        seed_ortholog_evalue: Maximum E-value threshold.
        ceiling_resolver: Pre-built ``TaxScopeCeilingResolver`` for
            per-seed ev_lca ceiling resolution.
        donor_pool: ``"closest"`` (default) or ``"union"``.
        lazy_cascade: Enable the tier-staged lazy closest cascade
            (byte-identical output, fewer ortholog fetches). Ignored for
            ``donor_pool="union"``. Default ``False``.
        pool: Optional ``multiprocessing.Pool`` (fork start method).
            When set, the engine slices each batch into sub-batches and
            dispatches them to pool workers.  Caller-managed lifecycle.
        dropped_writer: Optional callable
            ``dropped_writer(query, reason, seed, evalue, score)``
            invoked for every hit dropped before or after the engine
            runs.  ``None`` disables the drop log.

    Yields:
        ``((hit, annotation), False)`` tuples, same interface as the
        legacy per-hit annotator.
    """
    # Parse hits, drop those that fail score/evalue thresholds.
    valid_hits = []
    for args in batch:
        hit = args[0]
        query_name = hit[0]
        best_hit_name = hit[1]
        best_hit_evalue = float(hit[2])
        best_hit_score = float(hit[3])

        reason = filter_out(
            best_hit_name,
            best_hit_evalue,
            best_hit_score,
            seed_ortholog_evalue,
            seed_ortholog_score,
        )
        if reason is not None:
            if dropped_writer is not None:
                dropped_writer(
                    query_name, reason, best_hit_name, best_hit_evalue, best_hit_score
                )
            yield ((hit, None), False)
            continue

        # v7: DIAMOND/MMseqs emit integer protein IDs as strings.
        try:
            seed_id = int(best_hit_name)
        except ValueError:
            logger.warning(
                "Cannot parse integer seed ID from '%s'; skipping hit"
                " (DB version mismatch?)",
                best_hit_name,
            )
            if dropped_writer is not None:
                dropped_writer(
                    query_name,
                    "non_integer_seed_id",
                    best_hit_name,
                    best_hit_evalue,
                    best_hit_score,
                )
            yield ((hit, None), False)
            continue

        valid_hits.append((hit, query_name, seed_id, best_hit_evalue, best_hit_score))

    if not valid_hits:
        return

    engine = _get_engine(eggnog_db, ceiling_resolver, donor_pool, lazy_cascade)
    seed_ids = [s for _, _, s, _, _ in valid_hits]
    # Report the seed dedup ratio for this batch: how many query lines map to
    # how many unique seed orthologs (the engine computes each unique seed
    # once and fans the result back to its duplicates). The ratio is ~1x
    # without --sort_entries (scattered duplicates rarely share a batch) and
    # approaches the global redundancy with it (adjacent duplicates collapse).
    _n_queries = len(seed_ids)
    _n_unique = len(set(seed_ids))
    if _n_unique:
        print(colorify(
            f"  batch: {_n_queries} queries → {_n_unique} unique seeds "
            f"({_n_queries / _n_unique:.1f}x dedup)", "blue"), file=sys.stderr)
    # ``target_orthologs`` is a *floor* on the cascade — the engine
    # restricts which donor types contribute annotations.
    results = engine.annotate_batch(
        seed_ids,
        target_taxa=target_taxa,
        excluded_taxa=excluded_taxa,
        target_orthologs=target_orthologs,
        pool=pool,
        sub_batch_size=sub_batch_size,
    )

    # Re-shape engine results into the tuple consumed by output.py
    for hit, query_name, seed_id, evalue, score in valid_hits:
        r = results.get(seed_id)
        if not r or (not r.get("orthologs") and not r.get("annotations")):
            if dropped_writer is not None:
                # Engine ran but found nothing useful — could be a tax_scope
                # that excluded every donor, target_orthologs floor that
                # filtered the cascade empty, or a seed with no events.
                dropped_writer(
                    query_name,
                    "no_donor_orthologs",
                    str(seed_id),
                    evalue,
                    score,
                )
            yield ((hit, None), False)
            continue

        annotations = r.get("annotations") or {}
        annotations_confidence = r.get("annotations_confidence") or {}
        tax_ceiling = r.get("tax_ceiling") or "-"
        farthest_donor_taxid = r.get("farthest_donor_taxid") or "-"
        farthest_donor_lineage = r.get("farthest_donor_lineage") or "-"
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
            tax_ceiling,
            farthest_donor_taxid,
            farthest_donor_lineage,
        )

        yield ((hit, annotation), False)
