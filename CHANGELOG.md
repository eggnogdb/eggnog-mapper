# Changelog

All notable changes to eggNOG-mapper are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Changed

- **`--tax_scope` now accepts `auto` (default), `auto-broad`, or a fixed
  clade name/taxid** (e.g. `Metazoa`, `33208`). The old predefined scope
  list files (`bacteria`, `eukaryota`, `all_narrow`, etc.) are removed.
  Each seed gets an individual `ev_lca`-based ceiling: only speciation
  events whose `ev_lca` is at or below the ceiling are used.
  Dropped: `--tax_scope_mode`, `--scope_strict_og`, file-path scopes.
- **`--donor_pool {closest,union}`** (default `closest`). `union` takes
  consensus across all priority tiers rather than stopping at the first
  non-empty tier.
- **Output column `tax_scope_used` renamed to `tax_ceiling`** — reflects
  the resolved ceiling clade per seed (e.g. `Viridiplantae`).
- **New output columns `farthest_donor_taxid` and `farthest_donor_lineage`**
  — taxid and full lineage of the most distant donor used in annotation
  transfer.

---

## [3.0.0-beta2] — 2026-05-04

### Added

- **Compressed FASTA input**: gzip (`.gz`) and bzip2 (`.bz2`) query files
  are autodetected by magic bytes — no flag needed.
- **`--resume` improvements**: all skipped search steps now emit `[resume]`
  log lines; skipped hits are excluded from the reported query/seed rates.
- **Apptainer/Singularity image**: `apptainer/build.sh` builds a
  self-contained HPC image with DIAMOND, MMseqs2, HMMER, and Prodigal
  bundled. See `apptainer/` for usage.

### Changed

- **CLI simplified**:
  - `--go_evidence` / `--go_excluded` removed (dead code since v3.0).
  - `--outfmt_short` removed.
  - `--tax_scope auto-narrow` renamed to `--tax_scope auto` (default).
- **`--resume` safety**: existing output files are now validated before
  reuse; corrupted or incomplete hits files are rejected with a clear error.

---

## [3.0.0-beta1] — 2026-04-29

First public beta of the v3 lineage. Requires **eggNOG v7** database.
Not compatible with eggNOG v5 — use eggNOG-mapper 2.x for v5 databases.

### Added

- **eggNOG v7 database support** with integer-encoded orthology and
  phylogeny-based speciation events across ~12 M proteins and ~10 k taxa.
  v5 databases are rejected at startup with an actionable error.
- **Curated-only functional donors**: only manually curated terms (SwissProt
  and equivalent curated sources) are used as annotation donors. This avoids
  propagating misannotations that were inherited from automated pipelines in
  earlier database versions. v3 achieves better annotation coverage than v2
  despite the stricter source filtering.
- **Per-seed taxonomic ceiling** (`--tax_scope auto`): each query seed
  gets its own `ev_lca`-derived ceiling automatically narrowed to the most
  informative phylogenetic level. Plant seeds resolve to `Viridiplantae`,
  bacterial seeds to their phylum, etc. Fixed clades still accepted.
- **Cascade annotation engine**: for each functional source (GO, KEGG_ko,
  Pfam, EC, …) independently, donors are walked from closest and
  best-typed first. The seed's own curated annotation is the strongest
  tier-0 donor, preventing paralog inferences from overwriting
  SwissProt-curated values.
- **`annotation_confidence` output column**: per-field confidence tier
  (`high` = 1:1 donor, `medium` = 1:many or many:1, `low` = many:many).
- **`tax_ceiling` output column**: resolved ceiling clade name per seed.
- **`--cpu N` parallelises the annotation phase** via fork-context
  multiprocessing. Workers share the loaded engine state via COW and
  reopen their own SQLite connection post-fork.
- **Cython-accelerated inner loops**: `_codec` (delta-varint encode/decode)
  and `_collect_inner` (`OrthologCollector`). Annotation phase: ~25 min →
  ~30 s on a full eggNOG v7 run; `_collect_orthologs`: ~185 s → ~9 s on
  plant batches.
- **`--dmnd_top {1,3}`** (default `1`): uses `--max-target-seqs 1` for a
  free ~20–30 % faster DIAMOND pass with no biological impact.
- **DIAMOND auto-tune**: `block_size`, `index_chunks`, and `--algo`
  selected automatically from host RAM and query size. User overrides take
  precedence.
- **`missing taxa.db` error**: raises `EmapperException` with the expected
  path when `eggnog.taxa.db` is absent, replacing a silent failure.
- **`require_tool()` helper**: raises `EmapperException` with a
  `conda install` hint when a required external binary is not on PATH.

### Changed

- **`--target_orthologs` is now a cascade floor**, not a post-filter:
  `one2one` accepts only 1:1 events; seeds without matching donor types
  emit no annotations rather than silently widening to all types.
- **`--sensmode`** is the *iterate-ceiling* sensitivity (default
  `sensitive`); divergent queries escalate, easy ones do not.
- **`python_requires >= 3.9`** enforced in `setup.cfg`.

### Removed

- **Bundled binaries** (~150 MB: DIAMOND, HMMER suite, MMseqs2, Prodigal)
  removed from the package. Install via conda or your system package manager.
  Wheel size: ~150 MB → ~5 MB. macOS/Windows installs now work.
- **`eggnog-annotator` merged** into `eggnogmapper.annotator` (imports
  rewritten from `eggnog_annotator.X` → `eggnogmapper.annotator.X`).
  The `eggnog-annotator` repo is archived.
- `--mode cache` / `-c` / `--cache`
- `--mode novel_fams` and the novel-fams search/annotation chain
- `--dbmem`
- `--go_evidence`, `--go_excluded`
- Predefined tax-scope files (`bacteria`, `eukaryota`, `all_narrow`, etc.)
- `eggnogmapper/annotation/tax_scopes/` directory

### Fixed

- `scope_strict_og` was silently dropped on every full 1000-hit batch,
  causing proteomes >1000 hits to run in hybrid strict/permissive mode.
  Fix: `--scope_strict_og` now defaults to `True` and is threaded
  consistently through all batch calls.
- Seed IDs that are not integers now emit a warning instead of silently
  continuing.

### Migration from v2

```bash
# 1. Install external tools (no longer bundled)
conda install -c bioconda diamond hmmer mmseqs2 prodigal

# 2. Upgrade eggNOG-mapper
pip install --upgrade eggnog-mapper   # or: pip install eggnog-mapper==3.0.0b2

# 3. Download eggNOG v7 database  (v5 databases are not compatible)
download_eggnog_data.py --data_dir /path/to/eggnog-data

# 4. Run
emapper.py --version  # should print 3.0.0-beta2
emapper.py -m diamond -i proteins.fa --itype proteins \
    --data_dir /path/to/eggnog-data -o out --output_dir results/ --cpu 20
```

**CLI changes from v2:**

| v2 flag | v3 equivalent |
|---------|--------------|
| `--tax_scope bacteria` | `--tax_scope Bacteria` (taxid or clade name) |
| `--tax_scope auto` + `--tax_scope_mode inner_narrowest` | `--tax_scope auto` (default) |
| `--tax_scope auto` + `--tax_scope_mode broadest` | `--tax_scope auto-broad` |
| `--mode cache` | removed |
| `--mode novel_fams` | removed |
| `--dbmem` | removed |
| `--go_evidence` | removed |

---

## [2.1.15] — 2026-05-22

Final v2 release. Synced with all patches from the `master` branch.
Use this version if you need eggNOG v5 database compatibility.

- Fixed download URL for eggNOG database files
  ([#591](https://github.com/eggnogdb/eggnog-mapper/pull/591))

---

## [2.1.13] — 2024-xx-xx

See the [`master` branch](https://github.com/eggnogdb/eggnog-mapper/tree/master)
for the full v2 changelog.
