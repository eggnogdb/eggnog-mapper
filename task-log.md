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
Status: IN PROGRESS

---
