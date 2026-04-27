## [v3] — 2026-04-27

Per-source closest-ev_lca + ortholog-type-priority donor cascade. Phase
6 of `AnnotationEngine.annotate_batch` is no longer a flat Counter
aggregation across all orthologs; for each functional source the engine
walks donors from closest+best-typed first and stops at the first
non-empty bucket.

### Added

- `LineageCache.depth(taxid)` returns the size of a taxid's lineage.
  The e7 species table holds internal nodes (Bacteria, Metazoa, …) too,
  so this is a valid measure of taxonomic specificity for `ev_lca`
  taxids, not just leaf species.
- `AnnotationEngine.__init__` accepts an explicit `lineage_cache=`
  argument; falls back to `lineage_filter.lineage_cache` when not
  given. With both absent the cascade falls back to type-tier priority
  only.
- `_collect_orthologs` now returns a 3-tuple
  `(orthologs, ortholog_types, ortholog_meta)` where `ortholog_meta`
  carries per-ortholog priority metadata
  `{event_id, ev_lca, type, type_tier, depth, in_seed_lineage}`. The
  engine retains the *best* event each ortholog joined.
- `annotate_batch` and `annotate` results gain
  `annotations_confidence: dict[output_field, "high"|"medium"|"low"]`.
  Seed-derived `Preferred_name` is always `"high"`.
- `annotate_batch(..., target_orthologs="all")` parameter — the cascade
  applies it as a *floor* (excludes non-allowed types entirely). Valid
  values: `"all"`, `"many2many"`, `"one2many"`, `"many2one"`,
  `"one2one"`.
- `_pre_parse_batch` cache: each ortholog's annotation comma-strings
  parsed once per batch into tuples, before the cascade walk. **2.3×
  measured speedup** on phase 6 (30-query bacterial benchmark; plant
  proteomes scale into the 3–5× range targeted in the v3 plan).
- 7 new audit tests in `tests/test_sp_events_engine_audit.py` and
  13 new behavior tests in `tests/test_cascade.py` covering tier
  ordering, ev_lca-closer-wins, in-lineage-outranks-out-of-lineage,
  per-source independence, target_orthologs floor, missing source
  fall-through, and bucket-only consensus.

### Changed

- Type tier mapping (`AnnotationEngine.TYPE_TIERS`):
  `one2one→0` (high), `one2many→1` and `many2one→1` (medium),
  `many2many→2` (low). Confidence labels follow the tier.
- Cascade sort key (smaller = preferred):
  `(0 if in_seed_lineage else 1, -depth(ev_lca), type_tier)`.
- `_summarize_annotations` keeps the legacy flat aggregation behind a
  `parsed=None, ortholog_meta=None` back-compat path so the single-
  protein `annotate()` call and pre-cascade tests still work.
- `_aggregate_field` uses `set` for non-pname dedup (downstream sorts
  at write time, so order is irrelevant); only `pname` keeps `Counter`
  for `most_common`.

## [Unreleased — pre-v3 fixes folded into v3]

### Fixed

- **C7** (annotate.py): Out-of-range protein IDs now emit `logger.warning` instead of being silently discarded, making data integrity issues visible in logs.
- **C8** (tax_scope.py): `get_valid_species_ids()` now returns a `frozenset`, preventing accidental cache corruption by callers that mutate the result.
- **H4** (db.py): `EggnogDB.__init__` now raises `EggnogDBError` with the database path in the message on open failure, replacing a bare `sqlite3.OperationalError` with no context.
- **H5** (tax_scope.py): Bare `except Exception` replaced with `except sqlite3.OperationalError` combined with `logger.error()`, so unexpected exceptions are no longer swallowed silently.
