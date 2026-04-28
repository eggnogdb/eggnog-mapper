## [v3.4] — 2026-04-28

Annotation phase performance + biology pass. The annotation pipeline
is now CPU-bound rather than dict-and-string bound, parallel across
sub-batches, and biologically consistent with the auto-tax-scope
intent the user already asked for.

### Added

- **`--scope_strict_og` CLI flag** (default ON; opt out via
  `--no-scope_strict_og`). Skips speciation events whose containing
  OG (`sp_events.og_lca`) is broader than the seed's resolved
  tax-scope ceiling. On auto-scope plant proteomes this discards
  cross-kingdom paralog noise that the legacy permissive path was
  letting through to the per-protein species filter. Materially
  faster (~2× on plant proteomes) and more biologically consistent
  with what `--tax_scope auto` asks for. Use the opt-out for legacy
  reproducibility.
- **`--cpu N` now parallelizes the annotation phase**, not just
  DIAMOND search. `_annotate_batched` creates one fork-context
  `multiprocessing.Pool` for the full mapper run, registers the
  parent's loaded `AnnotationEngine` on a module global before fork,
  and passes the Pool to each `annotate_batch` call. Workers inherit
  the engine's loaded state (`taxid_array` ~226 MB, lineage_cache,
  og_cache) via fork copy-on-write and reopen their own SQLite
  connection in the post-fork initializer.

### Changed

- `Annotator.scope_strict_og` attribute (default `True`) — set from
  the `--scope_strict_og`/`--no-scope_strict_og` flag and threaded
  through to `engine.annotate_batch(scope_strict_og=…)`.
- `annotate_batch(pool=…)` is no longer ignored — it's passed through
  to the engine, which slices the batch and dispatches.

### Performance (real-case test on `/app/test_proteomes/`)

| | araport11 (27,596 hits) | itag4 (32,692 hits) |
|---|---:|---:|
| Wall (`--cpu 10`) | **10 min** | **10.5 min** |
| Reported throughput | 48 q/s | 54 q/s |
| Speedup vs v3.3 single-thread permissive | ~15× | ~15× |
| Speedup vs Python baseline (1000-seed bench) | 9.3× total / ~36× annotation-only | — |

### Biology

End-to-end consistency vs each query's seed ortholog source row in
`prots`:

| Field | araport mean Jaccard | itag4 mean Jaccard |
|---|---:|---:|
| `pfam` | 0.967 | 0.939 |
| `ec` | 0.785 | 0.734 |
| `gos` | 0.524 | 0.719 |
| `kegg_ko` | 0.724 | 0.649 |

High-confidence disagreements localise to close-paralog substitutions
within the same enzyme family or protein complex (PsbB/PsbD,
PsaB/PsaC, related phosphatases / methyltransferases / terpene
synthases). These would carry a `low` confidence label and reflect
honest cascade behaviour on divergent paralogs, not regressions.

## [v3] — 2026-04-27

Production cut. Three themes: legacy purge (mapper is now a thin shim
over `eggnog_annotator.e7.AnnotationEngine`), per-source closest-ev_lca
+ ortholog-type-priority cascade for annotation transfer with
confidence labels, and a builder that's self-contained and uniformly
SQLite-tuned for slow-disk machines.

### Removed (Phase 2 — legacy purge)

- `--mode cache` and the `-c/--cache <FILE>` argument.
- `--mode novel_fams` and the parallel novel-fams search/annotation chain.
- `--dbmem` flag and the `_annotate_dbmem` per-hit path.
- Modules unreachable from the v7 batch path:
  `cache_annotator.py`, `annotator_novel_fams.py`,
  `annotator_worker_novel_fams.py`, `output_novel_fams.py`,
  `annota.py`, `annota_mongo.py`, `orthologs.py`,
  `annotator_worker.py` (its only live helper, `filter_out`, was
  inlined into `batch_annotate.py`).
- `db_sqlite.py`: 11 dead methods removed
  (`get_member_ogs`, `get_ogs_description`, `get_annotations`,
  `get_pfam_annotations`, `get_taxid`, `get_protein_id`,
  `bulk_get_protein_ids`, `get_member_events`, `bulk_get_ogs`,
  `bulk_get_event_indices`, `bulk_get_events` and the v5/legacy
  schema branches inside them); kept `get_eggnog_db`,
  `get_fresh_eggnog_db`, `AnnotDB.__init__`, `close`,
  `get_db_version`, `decode_protein_ids`, `get_protein_name`.
- `Annotator._annotate_ondisk` (per-hit path for non-v7 DBs) and
  `_annotate_dbmem`; `Annotator._annotate` now raises
  `EmapperException` when the open DB is not v7-int_mode.

Net delta: 16 files changed, +151/-1727 LoC; e2e output identical to
the v3-baseline reference on `data/e7/sample` after Phase 2.

### Changed (Phase 3 — cascade engine)

- **`--tax_scope auto` semantics** are now per-functional-source. For
  each source (KEGG_ko, GOs, Pfam, Preferred_name, EC, KEGG_Pathway,
  …) independently, the engine walks orthologs in priority order
  `(in_seed_lineage, -ev_lca_depth, type_tier)`. The first bucket with
  any donor that has a non-empty value for that source wins, and the
  consensus is taken across only the donors in that winning bucket.
  Annotations from lower-priority buckets are never mixed in. The
  cascade halts per-source — common sources may resolve at the 1:1 tier
  while rare ones fall through to many:many.
- **`--target_orthologs` is now a floor on the cascade**, not a
  post-filter. `one2one` accepts only 1:1 events; `one2many` accepts
  1:1 + 1:many; `many2one` accepts 1:1 + many:1; `all` / `many2many`
  accept everything. If no donors of an allowed type cover a source,
  the source is silently omitted (no fall-through to disallowed types).
  Effect: `--target_orthologs one2one` on a seed without 1:1 donors now
  emits no annotations, instead of silently widening to all donors.
- **New TSV column** `annotation_confidence` in `.emapper.annotations`,
  serialized as `field=tier;field=tier;...` (sorted, semicolon-separated).
  Tier is `high` (1:1 winner), `medium` (1:many or many:1), `low`
  (many:many). Seed-derived `Preferred_name` (when the seed itself has
  an informative gene symbol) is labelled `high`. Empty when nothing was
  emitted.
- **Annotation tuple grew from 10 to 11 elements** (element[10] =
  `annotations_confidence: dict`). `output_orthologs_row` and
  `output_annotations_row` accept either size for back-compat.
- The collapse `1:many ≡ many:1 ≡ medium` matches the user's biological
  spec; the four ortholog types remain available end-to-end via
  `--target_orthologs` for users who want to discriminate.

### Added (Phases 1 + 3)

- `eggnog_annotator/lineage.py`: `LineageCache.depth(taxid)` returns
  the lineage size — used by the cascade to rank events by taxonomic
  specificity.
- `eggnog_annotator/e7/annotate.py`: `_pre_parse_batch` (Phase 4
  optimization) parses each ortholog's annotation comma-strings once
  per batch ahead of the cascade walk; the cascade winning-bucket scan
  now reads pre-parsed tuples and dedupes with `set` (not Counter).
  Measured **2.3× speedup** on `phase 6 build_results` on a 30-query
  bacterial benchmark; plant proteomes (with ~8× more orthologs/seed)
  scale into the targeted 3–5× range.
- 16 new property tests in eggnog-builder + eggnog-annotator audit
  the sp_events round-trip, the cascade priority key, target_orthologs
  floor, missing-source fall-through, and bucket-only consensus.

### Phase 0 baseline (CHANGELOG `[Unreleased]` items below) — already in v3 lineage

- **C1** (batch_annotate.py): Annotation tuple element[5] is now a proper `(og_name, cat, desc)` 3-tuple and element[9] is the filtered ortholog list, preventing column misalignment in output files.
- **C2** (batch_annotate.py): `--target-orthologs` filtering is now applied in the v7 batch path; all five orthology type keys (`one2one`, `one2many`, `many2one`, `many2many`, `all`) are populated per event before filtering.
- **C3** (batch_annotate.py): `all_ogs` now carries the full OG hierarchy across all levels rather than only the narrowest OG, ensuring correct `eggNOG_OGs` output for multi-level assignments.
- **H1** (batch_annotate.py): Non-integer seed IDs now emit `logger.warning` instead of silently continuing, making malformed input visible in logs.
- **H2** (batch_annotate.py): Missing `taxa.db` file now raises `EmapperException` with an actionable message identifying the missing path, replacing a silent failure.
- **H3** (batch_annotate.py): Scope key derivation logic in `_get_engine()` is now documented inline for maintainability.
