# eggnog-mapper Code Audit

Audit performed on v3 branch, 2026-03-24.

## Summary

| Category | Total | Fixed | Remaining |
|----------|-------|-------|-----------|
| Critical | 3 | 3 | 0 |
| High | 1 | 0 | 1 |
| Medium | 11 | 11 | 0 |
| Low / Hardening | 8 | 0 | 8 |

---

## Critical Issues

### C1. SQL Injection in `get_pfam_annotations()` — FIXED

- **File:** `eggnogmapper/annotation/db_sqlite.py:90-94`
- **Problem:** `seq_names` was interpolated directly into SQL via `%s` string
  formatting, while every other query in the file used parameterized `?`
  placeholders.
- **Impact:** A crafted sequence name could execute arbitrary SQL against the
  annotation database.
- **Fix:** Replaced with parameterized query using `?` placeholders and
  `WHERE name IN (...)`.

### C2. Undefined Variable Crashes Runtime — FIXED

- **File:** `eggnogmapper/annotation/cache_annotator.py:84`
- **Problem:** Referenced `cached_annot_fn` which does not exist; the parameter
  is named `cache_file`. Any run with `-m cache` where the cache file is missing
  would crash with `NameError`.
- **Fix:** Changed to `cache_file`.

### C3. N+1 Query Pattern in Hot Path — FIXED

- **File:** `eggnogmapper/annotation/db_sqlite.py:100-122`
- **Problem:** `get_member_events()` executed one SQL query per
  `(event_index, target_level)` pair. For a protein with N orthoindex entries
  and M target levels, this was N*M individual queries. This function is called
  for every hit during annotation — the single hottest database path.
- **Impact:** Severe performance degradation on large annotation runs.
- **Fix:** Replaced with a single batched query using
  `WHERE i IN (...) AND level IN (...)` with parameterized placeholders.

---

## High-Priority Issues

### H1. Shell Command Injection Surface — REMAINING

- **Files:**
  - `download_eggnog_data.py:28-31` — `os.system(cmd)` with f-string commands
  - `create_dbs.py:26-29` — same pattern
  - `eggnogmapper/search/diamond/diamond.py:262` — `subprocess.run(cmd, shell=True)`
  - `eggnogmapper/search/mmseqs/mmseqs.py:28,41,219,237,252,272` — same
  - `eggnogmapper/search/hmmer/hmmer_search.py:202,343` — `subprocess.call(cmd, shell=True)` and `os.system(cmd)`
  - `eggnogmapper/genepred/prodigal.py:96,142` — `subprocess.run(cmd, shell=True)`
  - `eggnogmapper/annotation/pfam/pfam_scan.py:72` — same
  - `eggnogmapper/annotation/pfam/pfam_common.py:55` — same
  - `eggnogmapper/annotation/pfam/pfam.py:228` — same
- **Problem:** All external tool invocations use `shell=True` with commands
  built via f-string interpolation of file paths. If a path contains shell
  metacharacters (`$()`, backticks, semicolons, quotes), arbitrary commands
  could be executed. While not exploitable in typical use (paths come from CLI
  args validated by argparse), it is fragile and violates defense-in-depth.
- **Recommended fix:** Replace `os.system()` and `subprocess.run(..., shell=True)`
  with `subprocess.run(cmd_list)` using list arguments. This is a large
  refactor touching every external tool call and should be done incrementally,
  one tool at a time, with integration tests after each change.

---

## Medium Issues — ALL FIXED

### M1. Non-standard SQL Operator `==`

- **File:** `eggnogmapper/annotation/db_sqlite.py:65,71,81,102,119`
- **Problem:** Used `==` instead of standard SQL `=` for equality comparisons.
  Works in SQLite but is non-standard and could break with other SQL engines.
- **Fix:** Replaced all `==` with `=` in SQL queries.

### M2. Bare `except:` Clause

- **File:** `eggnogmapper/annotation/annota.py:39`
- **Problem:** `except:` without an exception type catches `SystemExit`,
  `KeyboardInterrupt`, and other signals that should propagate.
- **Fix:** Changed to `except Exception:`.

### M3. Unclosed File Handle in `pfam.py`

- **File:** `eggnogmapper/annotation/pfam/pfam.py:106`
- **Problem:** `open(fasta_file)` used in a list comprehension without closing.
- **Fix:** Wrapped with `with open() as _f:` context manager.

### M4. Unclosed File Handle in `pfam_scan.py`

- **File:** `eggnogmapper/annotation/pfam/pfam_scan.py:92`
- **Problem:** `open(P.name, 'r').read()` inside a `with` block for the
  destination file, but the source file was never closed.
- **Fix:** Added nested `with open(P.name, 'r') as src:`.

### M5. Unclosed File Handle in `pfam_search.py`

- **File:** `eggnogmapper/annotation/pfam/pfam_search.py:98`
- **Problem:** Same pattern as M4.
- **Fix:** Same fix as M4.

### M6. Unclosed File Handle in `hmmer_search.py`

- **File:** `eggnogmapper/search/hmmer/hmmer_search.py:347`
- **Problem:** `for line in open(tempout):` — file never explicitly closed.
- **Fix:** Wrapped with `with open(tempout) as _f:`.

### M7. Missing Return Code Check for `hmmpress`

- **File:** `eggnogmapper/annotation/pfam/pfam_scan.py:72`
- **Problem:** `subprocess.run()` return code was not checked. If `hmmpress`
  failed, subsequent operations would fail with cryptic errors.
- **Fix:** Added `if cp.returncode != 0: raise EmapperException(...)`.

### M8. Missing Return Code Check for `hmmfetch`

- **File:** `eggnogmapper/annotation/pfam/pfam_common.py:55`
- **Problem:** Return code not checked; only file size was checked.
- **Fix:** Added `cp.returncode != 0` to the existing size check condition.

### M9. Missing Return Code Check for `esl-reformat`

- **File:** `eggnogmapper/annotation/pfam/pfam.py:228`
- **Problem:** Return code not checked.
- **Fix:** Added `if cp.returncode != 0: raise EmapperException(...)`.

### M10. Improper Multiprocessing Pool Cleanup in `annotator.py`

- **File:** `eggnogmapper/annotation/annotator.py:262-264`
- **Problem:** Called `pool.close()` then immediately `pool.terminate()`, which
  abruptly kills workers. Should wait for workers to finish.
- **Fix:** Changed to `pool.close(); pool.join()`.

### M11. Improper Multiprocessing Pool Cleanup in Pfam Modules

- **Files:** `eggnogmapper/annotation/pfam/pfam_scan.py:49-50`,
  `eggnogmapper/annotation/pfam/pfam_search.py:65-66`
- **Problem:** Called only `pool.terminate()` without `pool.close()` or
  `pool.join()`. Can orphan worker processes.
- **Fix:** Changed to `pool.close(); pool.join()`.

---

## Low Priority / Hardening — REMAINING

These are not bugs but improvements that would make the code more robust.

### L1. File Handle Leak in `cache_annotator.py` Generator

- **File:** `eggnogmapper/annotation/cache_annotator.py:59`
- **Problem:** `OUT = open(annot_file, "w")` is not in a context manager. Since
  the function is a generator (uses `yield`), the file must stay open across
  yields, making a `with` block difficult. If the generator is not fully
  consumed, `OUT.close()` on line 82 is never reached.
- **Recommendation:** Refactor to separate the file-writing concern from the
  generator, or use `try/finally` to ensure cleanup.

### L2. `get_annotations()` Iterates One Query Per Sequence

- **File:** `eggnogmapper/annotation/db_sqlite.py:75-88`
- **Problem:** For a comma-separated list of sequence names, executes one
  `SELECT` per name. Could be batched with `WHERE name IN (...)`.
- **Impact:** Minor — the list is typically co-orthologs for a single OG, so N
  is usually small (tens, not thousands).

### L3. SQLite Cache Size Could Be Larger

- **File:** `eggnogmapper/annotation/db_sqlite.py:41`
- **Problem:** `PRAGMA cache_size=2000` allocates ~8 MB of page cache. For a
  41 GB database with random access patterns, this causes excessive disk I/O.
- **Recommendation:** Consider `PRAGMA cache_size=-500000` (~500 MB) for
  on-disk mode, or make it configurable. The `--dbmem` option (load entire DB
  to RAM) is the alternative but requires 45 GB.

### L4. Inefficient HMM File Parsing

- **File:** `eggnogmapper/search/hmmer/hmmer_search.py:177-183`
- **Problem:** Parses HMM file line by line looking for NAME/LENG fields
  without early termination. No caching across calls.
- **Impact:** Minor for typical use.

### L5. `hmmer_search.py` Silently Continues on Subprocess Failure

- **File:** `eggnogmapper/search/hmmer/hmmer_search.py:202`
- **Problem:** `subprocess.call(cmd, shell=True)` return code is checked later
  but yields empty results on failure rather than raising an error.
- **Recommendation:** Raise or warn on non-zero return code.

### L6. Missing `hmmer_search.py` phmmer Return Code Handling

- **File:** `eggnogmapper/search/hmmer/hmmer_search.py:343`
- **Problem:** Uses `os.system(cmd)` and checks `status == 0` for the success
  path, but the `tempout` file may not exist if the command failed for reasons
  other than a non-zero exit code (e.g., signal).
- **Recommendation:** Replace with `subprocess.run()` and check both return
  code and file existence.

### L7. Temporary File Cleanup Could Use `atexit`

- **Files:** Various Pfam and HMMER modules
- **Problem:** Temporary files are cleaned up manually with `os.remove()`. If
  the process is interrupted, temporary files remain.
- **Recommendation:** Use `tempfile.NamedTemporaryFile(delete=True)` where
  possible, or register cleanup with `atexit`.

### L8. `download_eggnog_data.py` and `create_dbs.py` Use Global `args`

- **Files:** `download_eggnog_data.py:28-31`, `create_dbs.py:26-29`
- **Problem:** The `run()` function references `args.simulate` as a global
  variable. The `gunzip_flag()` function references `args.force` as a global.
  This makes the code harder to test and refactor.
- **Recommendation:** Pass `args` explicitly or refactor into a class.

---

## Change Log

| Commit | Issues Fixed |
|--------|-------------|
| `8de6689` (v3) | C1, C2, C3, M1–M11 |
