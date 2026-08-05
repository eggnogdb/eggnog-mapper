# Changelog

All notable changes to eggNOG-mapper are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Removed / renamed (pre-stable CLI cleanup)

- **`--db BACKEND` removed.** Data location is set exclusively via `--data_dir`
  or the `EGGNOG_DATA_DIR` environment variable (both `emapper.py` and
  `download_eggnog_data.py`). The internal default backend still resolves when
  neither is set.
- **`--translate` removed; translation is now automatic.** `--itype cds`
  always translates before the protein search (previously required
  `--translate`); `--itype genome`/`metagenome` gene prediction is unchanged
  (predicted ORFs remain in their existing default representation).
- **`--mp_start_method` removed.** The annotation pool hard-codes the fork
  start method (already the effective default for eggNOG-mapper). Standalone
  `hmm_server.py` / `hmm_worker.py` still expose the flag.
- **`--unsorted_seeds` removed; file-based seed inputs are always sorted** by
  seed ortholog id before annotation. Byte-identical to the previous default;
  the sort-vs-unsorted resume mode guard was removed. The in-memory
  (genome/metagenome) hits path is unchanged. The `go_namespace_split` output
  header marker stays.
- **HMMER as a search mode removed** (`-m hmmer` no longer accepted). The
  associated HMMER-search-only options were removed: `-d/--database`,
  `--servers_list`, `--qtype`, `--dbtype`, `--usemem`, `--Z`, `--cut_ga`,
  `--clean_overlaps`, `--hmm_maxhits`, `--report_no_hits`, `--hmm_maxseqlen`,
  plus the `hmm_mapper.py` entrypoint. `--pfam_realign` (realign / denovo) is
  unaffected — it drives its own hmmpgmd server internally, and the
  `--port`, `--end_port`, `--num_servers`, `--num_workers`,
  `--timeout_load_server` options were moved into the Annotation group where
  they now serve `--pfam_realign` only. The standalone `hmm_server.py` /
  `hmm_worker.py` scripts (for pre-starting an hmmpgmd server) are kept.
- **Diamond flags renamed to a consistent `--dmnd_*` prefix** (flag == dest):
  `--sensmode`→`--dmnd_sensmode`, `--matrix`→`--dmnd_matrix`,
  `--gapopen`→`--dmnd_gapopen`, `--gapextend`→`--dmnd_gapextend`,
  `--block_size`→`--dmnd_block_size`, `--index_chunks`→`--dmnd_index_chunks`.
  `--dmnd_top` is unchanged.
- **Annotation-stage thresholds renamed for clarity**:
  `--seed_ortholog_evalue`→`--annot_evalue`,
  `--seed_ortholog_score`→`--annot_score`. The search-stage `--evalue` /
  `--score` (shared by diamond + mmseqs) are unchanged.

### Added

- **Seed sorting is the default** (and now the only mode for file-based seed
  inputs): the `.seed_orthologs` entries are sorted by seed ortholog id before
  annotation, so queries sharing a seed form contiguous blocks and each unique
  seed is annotated once (its result fanned out to the duplicates). On full
  UniProt this removes ~3× redundancy. The sort is atomic and reused on
  `--resume`.
- **The lazy closest-cascade is now on by default** (disable with
  **`--no-lazy_cascade`**): a byte-identical tier-lazy annotation cascade for
  `donor_pool=closest` (the default; automatically a no-op under `union`). A
  per-protein field-presence mask (built once from `prots`, cached, and shared
  across workers) lets the cascade select each field's winning donor bucket and
  skip absent/sparse fields with zero ortholog fetches, then fetch only the
  winning donors. Verified byte-identical by a fuzz + real-data parity test.
- **`--annot_batch_size N`**: tune the number of seed orthologs annotated per
  worker task (bulk-queried together); with sorted seeds (the default) the outer
  batch is sized by distinct seeds so all workers stay fed on redundant inputs.
- **`go-basic.obo` is now shipped by `download_eggnog_data.py`** and the active
  GO mode is recorded in the output header as `go_namespace_split` — so a run
  using the per-namespace (MF/BP/CC) cascade vs the legacy combined fallback is
  reproducibly distinguishable from the file alone, not just a log line.
- **The example-citation printed at the end of a run now includes the eggNOG DB
  version used** (e.g. `eggNOG-mapper (version X; eggNOG DB version 7.0.0)`).
- **Offline self-test dataset** (`test_datasets/`): a self-consistent subsample
  of the final v7 mapper DB — nif-operon query seeds spanning Archaea and
  Bacteria together with their full annotation-transfer closure
  (`event_index` → `sp_events` → ortholog donors) — plus per-`--itype` query
  fixtures (proteins, CDS, genome, metagenome) and golden outputs. A `pytest`
  runner (`test_selftest.py`) exercises the whole option surface offline against
  the mini DB, including auto-translate (CDS) and prodigal gene prediction
  (genome/metagenome), and checks the nifH query recovers the `Fer4_NifH` OG.
  The query input fixtures are committed; the subsampled DB and goldens ship as
  a downloadable bundle (paired with the reference full DB, which can change),
  rebuilt reproducibly by `make_test_dataset.py` + `gen_golden.py`.

### Changed

- **Annotation throughput** on proteome/UniProt-scale runs is substantially
  higher (≈90–300 → ≈800+ q/s on full UniProt) via seed-sorted global dedup,
  distinct-seed batching, and the lazy cascade (all on by default) — all
  producing byte-identical output.
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
- **`go-basic.obo` is resolved relative to `--data_dir`** (`$EGGNOG_GO_OBO` →
  `<data_dir>/go-basic.obo` → `<data_dir>/e7/full/source/reference/go-basic.obo`)
  instead of a hardcoded path, so shipping the OBO with the database works for
  any data dir.
- **Progress `q/s` is measured from the first produced result**, excluding the
  one-time setup (DB load, worker fork, sort), so it reflects real annotation
  throughput instead of being diluted by startup.
- **CI migrated from Travis to GitHub Actions** (Python 3.9–3.11); the
  `biopython` floor in `requirements.txt` was raised to `>=1.79` to match
  `setup.cfg`. The README no longer directs production users to install v2.

### Fixed

- **Malformed `--annotate_hits_table` lines** now raise a clear error naming
  the file/line/content, instead of silently re-emitting the previous line's
  hit (which misattributed its annotation) or raising `UnboundLocalError`.
- **`--resume` in report-orthologs mode** matched on the first *character* of
  each line rather than the query column, so the resume lockstep never matched:
  every query was re-annotated and every ortholog row was duplicated (append
  mode). It now matches on the query column.
- **Stale binary caches after an in-place DB rebuild**: the `taxid` array and
  the lazy-cascade field-presence caches were validated by protein count only,
  so a rebuild that kept the same count could silently reuse a stale cache
  (wrong taxonomy / dropped or misplaced annotation fields). Both are now keyed
  on a DB size+mtime fingerprint; an archived older DB keeps its own cache.

---

## [3.0.0-beta4] — 2026-07-31

### Added

- **`--resume` skips seed_orthologs regeneration when the file is complete**:
  on resume, if `.seed_orthologs` is already fully written (detected via the
  `## N queries scanned` footer sentinel), both the DIAMOND/MMseqs2 search and
  the seed-writing/`.hits` re-parse are skipped and annotation reads the
  existing file directly. Incomplete or unmarked files (e.g. written with
  `--no_file_comments`) fall back to full regeneration, so the check can never
  reuse a truncated file. Proteins/CDS only — genome/meta still regenerate
  because the seed pass also yields the contig-relative hits used downstream
  for GFF/FASTA.

### Fixed

- **Search/annotation pipeline decoupling**: the `.seed_orthologs` file is now
  written completely before the annotation phase begins. Previously,
  `output_seeds` was a lazy generator that wrote to the file as the annotation
  workers consumed hits, leaving the file incomplete until the full run
  finished. The file is now generated in a single sequential pass immediately
  after the search tool (DIAMOND or MMseqs2) exits.
- **Annotation crash on large runs (`KeyError` from OG-info cache)**: on runs
  with more than ~100k unique orthologous groups, the OG-description cache
  evicted its oldest half *before* the current lookup had read its results,
  so a group carried over from an earlier batch could be dropped while still
  needed — aborting the whole run mid-way (observed at ~5.5M queries). The
  result is now read before eviction. This also closes a latent correctness
  gap where an evicted-but-needed group would have silently lost its
  description and COG category.

### Changed

- **Compressed input is now streamed, not pre-decompressed**: gzipped query
  files are passed straight to DIAMOND/MMseqs2 (both read `.gz` natively) and
  to the FASTA iterators, instead of first inflating the whole file to a temp
  copy on disk. This removes the full on-disk decompression step for large
  gzipped inputs (e.g. UniProt). bzip2 is still decompressed (no search tool
  streams it), and Prodigal now detects compression by magic bytes rather than
  the `.gz` extension. Annotations are identical to plain-input runs.

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
emapper.py --version  # should print 3.0.0-beta3
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
