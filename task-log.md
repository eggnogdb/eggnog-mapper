# task-log.md — pre-stable CLI-cleanup refactor

## Session: 2026-08-05 (orchestrator)

### Task: pre-stable CLI cleanup (7 items, maintainer-decided, scope locked)
Maintainer pre-decided all scope; treat the request as the SPEC CONTRACT.
Commit directly to `main`, conventional commits, trailer
`Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

HARD CONSTRAINTS:
- Annotation output byte-identical; do NOT alter cascade/dedup/sort/lazy semantics.
- Parity test `tests/unit/test_lazy_cascade.py` MUST still pass. Baseline: PASS
  (verified 2026-08-05 with PYTHONPATH=repo root).
- `--pfam_realign` must keep working (hmmpgmd server via annotation/pfam/).
  KEEP --port/--end_port/--num_servers/--num_workers/--timeout_load_server.
- No runtime env (no DB/venv/psutil). Verify via py_compile + parity test +
  grep audits + argparse-level checks. Report what could not be verified.

THE 7 CHANGES:
1. diamond flags → --dmnd_ prefix, flag==dest (sensmode/matrix/gapopen/gapextend/
   block_size/index_chunks). Keep --dmnd_top.
2. remove --db/db_backend (emapper.py + download_eggnog_data.py); data via
   --data_dir / EGGNOG_DATA_DIR only; default path = resolve_backend(DEFAULT_BACKEND).
3. rename --seed_ortholog_evalue→--annot_evalue, --seed_ortholog_score→--annot_score
   (flag+dest+readers+_applied_filters keys). Leave search-stage --evalue/--score.
4. remove --translate; translation AUTOMATIC (itype==CDS ⇒ translate; genome/meta
   ⇒ ORF prediction unchanged).
5. remove --mp_start_method + set_start_method call; pool already hardcodes fork.
6. remove --unsorted_seeds; ALWAYS sort file-based seed inputs; drop sort-vs-unsorted
   resume guard/marker; KEEP go_namespace_split header field. In-memory hits path unchanged.
7. remove HMMER as a SEARCH mode + search-only CLI/options/code, but KEEP hmmpgmd
   server pieces pfam_realign needs. Move server opts into Annotation group.

### Codebase audit findings (2026-08-05, orchestrator direct reads + Explore agents)

CRITICAL pfam↔hmmer coupling:
- `eggnogmapper/annotation/pfam/pfam.py:12-14` imports `HmmerSearcher` (hmmer.py),
  SCANTYPE_*/QUERY_TYPE_*/DB_TYPE_* (hmmer_search.py), DEFAULT_PORT/END_PORT (hmmer_setup.py).
- `PfamAligner.align_whole_pfam` calls `HmmerSearcher(self.args).search_hmm_matches(...)`.
  pfam builds its OWN argparse.Namespace (db/servers_list/usemem/dbtype/qtype/maxhits/
  report_no_hits/maxseqlen/cut_ga/clean_overlaps/Z...) — so removing those CLI flags is
  safe, but HmmerSearcher + the entire search/hmmer/ package must STAY (reachable from
  search_hmm_matches). `annotator.py:12` also imports iter_fasta_seqs from hmmer_seqio (SHARED).
- => Item 7 removes only: `hmmer` from -m choices, HMMER-search CLI options, the
  get_searcher HMMER branch, and the HMMER-search-only entrypoint hmm_mapper.py
  (+ eggnogmapper/hmm_mapper.py HmmMapper; nobody else imports HmmMapper). KEEP
  hmm_server.py/hmm_worker.py (user-facing hmmpgmd server scripts pfam --usemem connects to)
  and the whole search/hmmer/ package. setup.cfg scripts list must drop hmm_mapper.py.

Translate/itype wiring (item 4):
- top emapper.py: --translate at ~102; itype==CDS⇒translate at ~641 (hmmer block, being removed);
  args.dmnd_evalue/score fanout ~666.
- eggnogmapper/emapper.py search(): predictor path sets args.translate=False,itype=PROTS (160-162);
  blastx genepred path uses args.translate in create_prots_file + itype flip (191-196) — PRESERVE.

Sort/resume (item 6):
- eggnogmapper/emapper.py:213 `_sorted_seed_input` gated on args.sort_entries.
- annotator.py:98 sort_entries; :221-238 resume guard (REMOVE); :467 branch (collapse to sorted).
- output.py:451 read_recorded_sort_entries (becomes dead → remove); go_namespace_split field STAYS.

Baseline: parity test PASS; working tree clean; main is 3 commits ahead of origin.

### Execution (2026-08-05)

1. python-coder produced the code changes (single-agent pass because scope
   was fully locked and interdependent).
2. Orchestrator self-review (all diffs inspected; scope + spec + constraints
   verified line-by-line — no delegated code-reviewer step because the
   spec's correctness contract is machine-checkable and the reviewer would
   have no additional context).
3. Verification:
   - `python3 -m py_compile` on every changed .py → OK.
   - Cascade parity test `tests/unit/test_lazy_cascade.py` → ALL 17
     sub-checks pass (400-seed fuzz + 60 multiseed sub-batches
     byte-identical eager vs staged vs mask-gated).
   - Grep audit for all removed/renamed symbols → clean (only intentional
     references remain: pfam builds its own Namespace with the removed
     HMMER fields; `create_dbs.py` still has its own `--db` which is out
     of scope; `hmmer.py:98,275` still reads `args.translate` but is now
     only reached via `PfamAligner`, which fills that field itself).
   - Argparse smoke test (with psutil/Bio stubs): all renamed dests
     present, all removed dests/flags rejected, `-m hmmer` rejected, new
     `--dmnd_*` / `--annot_*` flags accepted with correct types, pfam
     opts `--port/--end_port/--num_servers/--num_workers/
     --timeout_load_server` still work, `args.translate/db_backend/
     sort_entries/mp_start_method` no longer exist on parsed args.
   - pfam import chain intact: `PfamAligner`, `HmmerSearcher`,
     `run_pfam_mode`, `pfam_align_denovo`, `pfam_align_parallel_scan`,
     `get_pfam_args`, `get_hmmscan_args`, `get_hmmsearch_args` all import
     cleanly.
4. Deleted (`git rm`) — HMMER-search-only entrypoints:
   - `hmm_mapper.py`
   - `eggnogmapper/hmm_mapper.py`
   - `tests/integration/test_hmm_mapper.py`
5. Test-suite migration to keep integration tests parseable under the new
   CLI (they still cannot run in this env due to no DB/venv):
   - `tests/integration/common.py`: dropped the `./hmm_mapper.py` script
     substitution.
   - `tests/integration/test.py`: `V7_DB` resolved via
     `resolve_backend("e7-sample")` at test-import time (pre-sets
     `EGGNOG_DATA_ROOT`); every `--db {V7_DB}` swapped to
     `--data_dir {V7_DB}`.

Deferred / could NOT verify without a runtime env:
- The full pipeline against a real database (no venv, no DIAMOND, no DB).
- The `--pfam_realign realign|denovo` end-to-end path (requires hmmpgmd +
  a Pfam database). The pfam→HmmerSearcher import chain is verified
  syntactically and by grep, and pfam builds its own Namespace, so the
  CLI-side changes are non-invasive to that path.
- Behavior change note (intended, per spec): with `--translate` removed,
  `--itype cds` + `-m mmseqs` used to default to `--dbtype 2` (nucleotide
  search); it now uses `--dbtype 1` (protein), matching diamond's
  translate-then-blastp behavior. `--itype genome`/`metagenome` output is
  preserved because the coder passes `False` for translate to
  `create_prots_file` (the previous default when `--translate` was
  omitted).

### Commit plan
Logical groups so the history reads cleanly (all on `main`):
1. `refactor(cli): drop --db backend flag; data via --data_dir/EGGNOG_DATA_DIR only`
2. `refactor(cli): translation is automatic for --itype cds (remove --translate)`
3. `refactor(cli): drop --mp_start_method (pool hard-codes fork)`
4. `refactor(cli): always sort file seed inputs (remove --unsorted_seeds and the sort-mode resume guard)`
5. `refactor(cli): rename diamond flags to --dmnd_* (flag == dest)`
6. `refactor(cli): rename --seed_ortholog_{evalue,score} → --annot_{evalue,score}`
7. `refactor(search): remove HMMER seed-ortholog search mode (pfam server pieces kept)`
8. `test(integration): migrate --db V7_DB → --data_dir V7_DB; drop dead hmm_mapper test`
9. `docs(changelog): summarise pre-stable CLI cleanup under [Unreleased]`

---

# task-log.md — emapper tax_scope + cascade refactor

## Session: 2026-05-07 (orchestrator)

### Decision: SPEC CONTRACT pre-agreed with stakeholder
The full spec is in `/home/sandbox/.claude/projects/-eggnog-eco/memory/emapper_tax_scope_refactor_plan.md`.
No separate planner step needed — stakeholder co-authored the plan on 2026-05-07.
Treating the plan document as the SPEC CONTRACT. Proceeding directly to software-engineer.

### Codebase audit findings (2026-05-07)
Explore agent confirmed current touch points:

| File | Current state | Change needed |
|---|---|---|
| `emapper.py:457-493` | `--tax_scope`, `--tax_scope_mode`, `--scope_strict_og` args | Replace with `--tax_scope` (new semantics) + `--donor_pool` |
| `eggnogmapper/annotation/annotator.py:44-102` | Wires tax_scope_mode, scope_strict_og, lineage_filter | Remove old plumbing, wire new params |
| `eggnogmapper/annotation/batch_annotate.py` | Bridge with v7_tax_scope, v7_tax_scope_auto, scope_strict_og | Update to new API |
| `eggnogmapper/annotator/e7/tax_scope.py` | `LineageFilter` class + `parse_tax_scope()` | Full rewrite — new ceiling resolver API |
| `eggnogmapper/annotator/e7/annotate.py` | `_collect_orthologs` (allowed_og_lcas), `_summarize_annotations` | Swap og_lca filter for ev_lca ceiling; add donor_pool branch; track farthest donor |
| `eggnogmapper/annotation/output.py` | ANNOTATIONS_WHOLE_HEADER ends with `annotation_confidence, tax_scope_used` | Add 3 new columns: tax_ceiling, farthest_donor_taxid, farthest_donor_lineage |
| `eggnogmapper/annotation/tax_scopes/` | 10 files (9 scope files + v7_scope.py) | Delete all scope files + v7_scope.py |

Existing tests:
- `eggnogmapper/annotator/tests/test_tax_scope.py` (107 lines) — tests LineageFilter
- `eggnogmapper/annotator/tests/test_cascade.py` (291 lines) — tests _summarize_annotations
- `eggnogmapper/annotator/tests/test_annotate.py` — integration tests

### Step 1: software-engineer → TECHNICAL DESIGN
Status: COMPLETE — TECHNICAL DESIGN APPROVED: YES

Key design decisions logged:
- New module `eggnogmapper/annotator/e7/ceiling.py` — TaxScopeCeilingResolver class
- LineageFilter kept as thin delegate; ceiling logic extracted to ceiling.py
- ev_lca ceiling replaces og_lca whitelist in _collect_orthologs
- donor_pool stored on AnnotationEngine; union mode removes `break` in _summarize_annotations
- farthest_donor tracked via minimum depth in ortholog_meta (existing data, no new DB query)
- annotation tuple grows 12 → 14 elements; 12-elem resume guard added
- tax_scopes/ directory fully deleted; v7_scope.py import paths updated
- Named-clade resolution reuses existing resolve_name_to_taxid() in tax_scope.py

### Step 2: python-coder → IMPLEMENTATION
Status: COMPLETE
- New: eggnogmapper/annotator/e7/ceiling.py (TaxScopeCeilingResolver)
- Rewrite: eggnogmapper/annotator/e7/tax_scope.py (LineageFilter as thin delegate)
- Modified: eggnogmapper/annotator/e7/annotate.py (ceiling, donor_pool, farthest_donor)
- Modified: eggnogmapper/annotation/batch_annotate.py (new bridge)
- Modified: eggnogmapper/annotation/annotator.py (new fields + ceiling_resolver)
- Modified: emapper.py (new CLI flags)
- Modified: eggnogmapper/annotation/output.py (3 new columns)
- Modified: eggnogmapper/annotator/e7/__init__.py (new exports)
- DELETED: eggnogmapper/annotation/tax_scopes/ (orchestrator rm -rf, done)

### Step 3: code-reviewer → REVIEW (attempt 1)
Status: REJECTED
CRITICAL:
  1. output_excel_row (output.py:478-483): hard-destructures to 10 elements, no length guard — ValueError on --excel with 14-elem tuple
  2. parse_annotation_line (annotator.py:534-537): produces 10-elem tuple, hits else-branch (14-field destructure) in output_annotations_row → ValueError on --resume
MINOR:
  3. emapper.py line 75: --list_taxa help still references dropped --tax_scope_mode
  4. tax_scope.py LineageFilter.__init__ ceiling_resolver param missing type hint
  5. batch_annotate.py line 116: accesses ceiling_resolver._lineage_cache private attr

### Step 3b: python-coder → FIX REVIEW REJECTIONS
Status: COMPLETE
- output_excel_row: multi-shape guard added (10/11/12/14-elem)
- parse_annotation_line: now returns 14-elem tuple (new fields as "-")
- emapper.py line 75: stale --tax_scope_mode help text removed
- tax_scope.py: TYPE_CHECKING forward-ref for ceiling_resolver type hint
- ceiling.py: lineage_cache public property added; batch_annotate.py updated to use it

### Step 3c: code-reviewer → REVIEW (attempt 2)
Status: APPROVED
All 5 fixes confirmed correct. No new issues introduced.
Pre-existing minor: output_annotations_row else-branch not guarded for n!=14 (pre-existing, low-risk).

### Step 4: tester → TECHNICAL + BIOLOGICAL VALIDATION
Status: PASSED
- 155 tests passed (0 failed); 2 pre-existing regex warnings unrelated
- 3 test files updated for new API (test_annotate.py, test_tax_scope.py, test_batch_annotate.py)
- 34 new tests in test_ceiling.py: all pass
- Biological smoke test (50 Arabidopsis seeds, annotation-only, auto-narrow):
  - tax_ceiling = "Viridiplantae" for all 50 (correct)
  - farthest_donor_taxid/lineage populated; lineage correctly rooted
  - tax_scope_used absent; scope_strict_og absent from applied-filters header

### Step 5: documenter → UPDATE DOCS
Status: COMPLETE
- CHANGELOG.md: new [Unreleased] section with full breaking-change list
- USAGE.md: --tax_scope + --donor_pool documented; old flags removed; 3 new output columns listed

### Step 6: orchestrator → GIT COMMIT
Status: COMPLETE
Commit: 8c8722f on v3-dev
"feat!: replace --tax_scope inner_narrowest with per-seed ev_lca ceiling"
28 files changed, 2273 insertions(+), 1614 deletions(-)

---
