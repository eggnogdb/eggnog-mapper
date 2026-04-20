## [Unreleased]

### Fixed

- **C1** (batch_annotate.py): Annotation tuple element[5] is now a proper `(og_name, cat, desc)` 3-tuple and element[9] is the filtered ortholog list, preventing column misalignment in output files.
- **C2** (batch_annotate.py): `--target-orthologs` filtering is now applied in the v7 batch path; all five orthology type keys (`one2one`, `one2many`, `many2one`, `many2many`, `all`) are populated per event before filtering.
- **C3** (batch_annotate.py): `all_ogs` now carries the full OG hierarchy across all levels rather than only the narrowest OG, ensuring correct `eggNOG_OGs` output for multi-level assignments.
- **H1** (batch_annotate.py): Non-integer seed IDs now emit `logger.warning` instead of silently continuing, making malformed input visible in logs.
- **H2** (batch_annotate.py): Missing `taxa.db` file now raises `EmapperException` with an actionable message identifying the missing path, replacing a silent failure.
- **H3** (batch_annotate.py): Scope key derivation logic in `_get_engine()` is now documented inline for maintainability.
