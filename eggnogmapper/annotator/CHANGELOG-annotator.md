## [v3.4] — 2026-04-28

Annotation phase performance + biology pass. End-to-end on the test
proteomes: araport (27,596 hits) and itag4 (32,692 hits) each finish
in ~10 min wall at `--cpu 10`, ~15× faster than the v3.3 single-thread
permissive baseline. Output is biologically self-consistent with the
seed ortholog's source annotations (mean Jaccard 0.94-0.97 on PFAMs,
0.7-0.8 on EC/KEGG_ko).

### Changed

- **`AnnotationEngine.annotate_batch(scope_strict_og=…)`** is now
  `True` by default. Events whose containing OG (`sp_events.og_lca`)
  is broader than the seed's resolved tax-scope ceiling are dropped
  before any orthologs are fetched. On auto-scope plant proteomes
  this discards 56 % of events that contributed 99 % of orthologs —
  the cross-kingdom paralog noise that was leaking through the
  per-protein species filter. Set `False` for the legacy permissive
  behaviour.
- **`AnnotationEngine.annotate_batch(pool=…)`** accepts a fork-context
  `multiprocessing.Pool`. When set, the batch is sliced into
  sub-batches (default 125 seeds each) and dispatched via
  `imap_unordered`. Workers inherit the parent's loaded engine state
  (taxid_array, lineage_cache, og_cache) via fork copy-on-write and
  reopen their own SQLite connection in the post-fork initializer.

### Added

- **Cython `_collect_inner.pyx`** — `OrthologCollector` cdef class
  backed by a C++ `unordered_map<int64, OrthEntry>`. Replaces the
  6.6 M-iteration Python inner loop in `_collect_orthologs` that built
  `ortholog_meta_raw[oid] = (sort_key_tuple, payload_dict)`. On full
  e7 plant batches this phase drops from ~185 s to ~9 s
  (single-thread). Falls back to the pure-Python loop when the .so
  didn't compile.
- **`LineageCache._tracks`** — ordered root → leaf species lineages,
  populated alongside the existing set form. Required to compute the
  descendant set of an internal clade (where set-form lineages lose
  the ordering needed to identify "below scope").
- **`LineageFilter.get_scope_og_descendants(scope_taxids)`** — returns
  the frozenset of taxids at-or-below any scope taxid; cached. Used as
  the `og_lca` whitelist for strict-OG mode.
- **`EggnogDB.reopen_connection()`** — drops the (fork-stale) inherited
  connection and opens a fresh one against the retained `db_path`.
  `from_connection()` now auto-recovers the path via
  `PRAGMA database_list` so adopted connections can be reopened
  post-fork without the caller passing a path.
- **Module-level pool plumbing** in `e7/annotate.py`: `_WORKER_ENGINE`
  global, `_register_worker_engine(engine)`, `_worker_init_after_fork()`,
  `_worker_annotate_subbatch(args)`. The mapper-side caller registers
  its engine in the parent before forking the pool.
- **`_annotate_batch_inproc()`** — refactored from the body of
  `annotate_batch`; called both directly (serial path) and from the
  pool worker function.

### Performance

1000-seed araport bench, single-thread Python baseline = 324 s:

| Step | Wall | Cumulative |
|---|---:|---:|
| Cython `_collect_orthologs` only | 143 s | 2.27× |
| + `scope_strict_og` default ON | 72 s | 4.5× |
| + 8-worker fork pool | **35 s** | **9.3×** |
| Annotation phase only (excl. 30 s init) | 9 s | **~36×** |

### Notes

- The strict-OG mode produces output that differs from permissive on
  ~35 % of common queries — and 7.5 % of seeds get no annotation at
  all if every event was above scope. All differences are in the
  expected direction: queries either lose data borrowed from
  cross-kingdom paralogs (honest auto-scope answer) or gain data the
  seed itself was missing (cascade discovers KEGG/EC/etc.).
- High-confidence disagreement examples (zero KEGG_ko or EC overlap
  between query and its seed) all localise to **close-paralog
  substitutions** within the same enzyme family or protein complex
  (PsbB/PsbD, PsaB/PsaC, related phosphatases, terpene synthases,
  methyltransferase variants). These would carry a `low` confidence
  label.

## [v3.2] — 2026-04-27

### Added

- **Cython codec** (`eggnog_annotator/_codec.pyx`). The hot delta-varint
  encode/decode loops compile to a native extension at install time;
  the public `eggnog_annotator.codec` module is now a thin import shim
  that prefers the compiled `_codec.so` and falls back transparently
  to `eggnog_annotator.codec_py` (byte-identical pure-Python) when the
  C compile is unavailable. Build-system requires Cython>=3.0; runtime
  requires nothing new. `_BACKEND` constant exposes the resolved
  backend (`"cython"` or `"python"`) for tests and debugging.
- 12 backend-selection tests in `tests/test_codec_backend.py` covering
  per-backend round-trip, edge cases (empty / single value / varint
  boundaries), and a cross-backend byte-identity check that runs
  whenever both backends are available.

### Changed

- `eggnog_annotator/codec.py` is now an import shim (was the canonical
  pure-Python implementation). The pure-Python implementation moved to
  `codec_py.py`; same byte format, same API.

Why: at full e7 the codec encodes 67 M side1+side2 BLOBs during
`build-web --step eggnogdb` Phase A. Cython drops that step from
~25 min to ~30 sec on the same hardware. The maintenance step-up is
contained — the codec spec is frozen and the file rarely changes.

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
