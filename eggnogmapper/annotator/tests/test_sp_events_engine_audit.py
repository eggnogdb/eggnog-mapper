"""Phase 1 correctness audit on the engine side.

Verifies that AnnotationEngine._collect_orthologs correctly transforms
sp_events into the (orthologs, ortholog_types) tuple downstream code relies
on, by recomputing the expected result directly from the SQLite tables and
comparing.

Properties under test
---------------------
1. For any seed, the union of `_collect_orthologs` typed buckets equals the
   union of (other_side ∖ {seed}) over every event the seed appears in,
   filtered by the active lineage scope. No filter is set in this test, so
   filtering is the identity.
2. Type classification depends only on side cardinality. We exhaustively
   verify every event the seed touches, against a brute-force recomputation
   from the same event dicts the engine sees.
3. Type assignment is invariant under set semantics: shuffling the
   construction order of side1 / side2 (using a Python set, which has
   non-deterministic iteration) does not change the typed buckets.
4. `get_events_bulk` returns side1 / side2 as `set` (not `list`), making
   ordering not part of the contract.
"""

from __future__ import annotations

import os
import random
import sqlite3
from collections import defaultdict
from pathlib import Path

import pytest

from eggnogmapper.annotator.codec import decode_intlist
from eggnogmapper.annotator.e7 import AnnotationEngine, EggnogDB

# Use the web-side DB rather than the mapper-side because it has the
# ogs/seqinfo tables — but for collect_orthologs only sp_events / event_index
# are needed and either DB works. We pick web for full data and the existing
# mapper-side tests cover the stripped DB path.
import os as _os
WEB_DB_PATH = Path(
    _os.environ.get("EGGNOG_SAMPLE_WEB_DB")
    or Path(__file__).resolve().parents[4] / "data" / "e7" / "sample" / "final" / "web" / "eggnog.db"
)


@pytest.fixture(scope="module")
def raw_conn():
    if not WEB_DB_PATH.exists():
        pytest.skip(f"Sample web DB not built: {WEB_DB_PATH}")
    conn = sqlite3.connect(str(WEB_DB_PATH))
    conn.row_factory = sqlite3.Row
    yield conn
    conn.close()


@pytest.fixture(scope="module")
def engine_no_filter(raw_conn):
    """Engine with NO lineage filter — pure collection from sp_events."""
    db = EggnogDB.from_connection(raw_conn, taxid_array=None)
    return AnnotationEngine(db, lineage_filter=None)


# ---------------------------------------------------------------------------
# Sample-DB helpers
# ---------------------------------------------------------------------------

def _seeds_with_events(conn, n: int = 25, seed: int = 1729) -> list[int]:
    """Return up to N protein IDs that appear in at least one event.

    Sampling is deterministic so failures are reproducible; we exclude
    proteins that touch suspiciously many events to keep the brute-force
    recomputation fast.
    """
    rows = conn.execute(
        "SELECT protein_id, length(events) AS n FROM event_index "
        "ORDER BY protein_id"
    ).fetchall()
    pids = [r["protein_id"] for r in rows if r["n"] < 200]
    rng = random.Random(seed)
    rng.shuffle(pids)
    return pids[:n]


def _expected_for_seed(conn, seed_id: int) -> tuple[set[int], dict[str, set[int]]]:
    """Brute-force recompute (orthologs, typed-buckets) for one seed
    directly from sp_events rows. The engine must agree exactly when run
    with no taxid array and no lineage filter.
    """
    ev_ids_blob = conn.execute(
        "SELECT events FROM event_index WHERE protein_id = ?",
        (seed_id,),
    ).fetchone()
    if ev_ids_blob is None:
        return set(), {k: set() for k in (
            "one2one", "one2many", "many2one", "many2many", "all"
        )}

    ev_ids = decode_intlist(ev_ids_blob["events"])
    if not ev_ids:
        return set(), {k: set() for k in (
            "one2one", "one2many", "many2one", "many2many", "all"
        )}

    placeholders = ",".join("?" * len(ev_ids))
    rows = conn.execute(
        f"SELECT i, side1, side2 FROM sp_events WHERE i IN ({placeholders})",
        ev_ids,
    ).fetchall()

    candidates: dict[str, set[int]] = {
        "one2one": set(), "one2many": set(),
        "many2one": set(), "many2many": set(),
    }
    for row in rows:
        s1 = set(decode_intlist(row["side1"]))
        s2 = set(decode_intlist(row["side2"]))
        if seed_id in s1:
            seed_side, other = s1, s2
        elif seed_id in s2:
            seed_side, other = s2, s1
        else:
            continue
        s, o = len(seed_side), len(other)
        if s == 1 and o == 1:
            rel = "one2one"
        elif s == 1 and o > 1:
            rel = "one2many"
        elif s > 1 and o == 1:
            rel = "many2one"
        else:
            rel = "many2many"
        candidates[rel].update(other)

    for bucket in candidates.values():
        bucket.discard(seed_id)
    all_ = set().union(*candidates.values())
    typed = dict(candidates)
    typed["all"] = all_
    return all_, typed


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_get_events_bulk_returns_sets(engine_no_filter, raw_conn):
    """The engine relies on set membership for `seed_id in side1`.
    Confirm the contract is upheld."""
    ev_ids = [r["i"] for r in raw_conn.execute(
        "SELECT i FROM sp_events ORDER BY i LIMIT 5"
    )]
    events = engine_no_filter.db.get_events_bulk(ev_ids)
    assert events
    for eid, ev in events.items():
        assert isinstance(ev["side1"], set), f"event {eid} side1 not set"
        assert isinstance(ev["side2"], set), f"event {eid} side2 not set"


def test_collect_orthologs_matches_brute_force_recomputation(
    engine_no_filter, raw_conn
):
    """For each sampled seed, engine collection must equal the brute-force
    expectation byte-for-byte. Run with NO lineage filter and no taxid
    array (engine bypasses species filtering in that case).
    """
    seeds = _seeds_with_events(raw_conn, n=25)
    assert seeds, "no events in sample DB"

    for seed in seeds:
        ev_ids = engine_no_filter.db.get_events_for_protein(seed)
        events = engine_no_filter.db.get_events_bulk(ev_ids)
        actual_all, actual_typed, _meta = engine_no_filter._collect_orthologs(
            seed, events
        )
        expected_all, expected_typed = _expected_for_seed(raw_conn, seed)

        assert actual_all == expected_all, (
            f"seed {seed}: ortholog set mismatch "
            f"(missing={sorted(expected_all-actual_all)[:5]}, "
            f"extra={sorted(actual_all-expected_all)[:5]})"
        )
        for rel in ("one2one", "one2many", "many2one", "many2many", "all"):
            assert actual_typed[rel] == expected_typed[rel], (
                f"seed {seed} bucket {rel}: mismatch "
                f"(missing={sorted(expected_typed[rel]-actual_typed[rel])[:5]}, "
                f"extra={sorted(actual_typed[rel]-expected_typed[rel])[:5]})"
            )


def test_seed_never_appears_in_its_own_ortholog_set(
    engine_no_filter, raw_conn
):
    """A seed must not appear in any of its own ortholog buckets."""
    for seed in _seeds_with_events(raw_conn, n=15):
        ev_ids = engine_no_filter.db.get_events_for_protein(seed)
        events = engine_no_filter.db.get_events_bulk(ev_ids)
        all_orth, typed, _ = engine_no_filter._collect_orthologs(seed, events)
        assert seed not in all_orth, (
            f"seed {seed} present in own ortholog set"
        )
        for rel in ("one2one", "one2many", "many2one", "many2many"):
            assert seed not in typed[rel]


def test_typed_buckets_are_subsets_of_all(engine_no_filter, raw_conn):
    """For every seed, typed[rel] ⊆ typed['all']. The engine's contract is
    that 'all' is the union of all typed buckets."""
    for seed in _seeds_with_events(raw_conn, n=15):
        ev_ids = engine_no_filter.db.get_events_for_protein(seed)
        events = engine_no_filter.db.get_events_bulk(ev_ids)
        all_orth, typed, _ = engine_no_filter._collect_orthologs(seed, events)
        for rel in ("one2one", "one2many", "many2one", "many2many"):
            assert typed[rel] <= typed["all"], (
                f"seed {seed}: typed[{rel}] not subset of typed['all']"
            )
        union = set()
        for rel in ("one2one", "one2many", "many2one", "many2many"):
            union |= typed[rel]
        assert union == typed["all"], (
            f"seed {seed}: typed['all'] not equal to union of typed buckets"
        )


def test_type_classification_invariant_to_side_construction_order(
    engine_no_filter,
):
    """Construct a synthetic event dict twice with sides that have the
    same membership but different *insertion* order (set ordering is
    unspecified). The engine must produce identical typed buckets.

    This pins down that type classification depends solely on cardinality
    after the codec's sorting, not on any incidental order.
    """
    # 1:many event: seed alone vs three others
    seed = 100
    other = [200, 201, 202]
    ev_a = {1: {"side1": {seed}, "side2": set(other)}}
    ev_b = {1: {"side1": set(reversed(other)) | {0xDEADBEEF},  # noise
                "side2": {seed}}}  # swapped sides + noise — should not affect type if seed is alone

    # We construct the canonical case (seed in side1, alone) with two
    # different *materializations* of side2 to ensure result is stable.
    other_set_a = set(other)
    other_set_b = set(reversed(other))
    assert other_set_a == other_set_b  # sanity
    ev_c = {1: {"side1": {seed}, "side2": other_set_b}}

    all_a, typed_a, _ = engine_no_filter._collect_orthologs(seed, ev_a)
    all_c, typed_c, _ = engine_no_filter._collect_orthologs(seed, ev_c)
    assert all_a == all_c
    assert typed_a == typed_c
    assert all_a == set(other)
    assert typed_a["one2many"] == set(other)
    for rel in ("one2one", "many2one", "many2many"):
        assert typed_a[rel] == set()


def test_type_classification_exhaustive_table(engine_no_filter):
    """Synthetic table covering all four cardinality classes."""
    cases = [
        # seed_side, other_side, expected_rel
        ({1}, {2}, "one2one"),
        ({1}, {2, 3, 4}, "one2many"),
        ({1, 5}, {2}, "many2one"),
        ({1, 5, 9}, {2, 3, 4, 11}, "many2many"),
    ]
    seed = 1
    for seed_side, other_side, expected_rel in cases:
        events = {0: {"side1": set(seed_side), "side2": set(other_side)}}
        all_orth, typed, _ = engine_no_filter._collect_orthologs(seed, events)
        assert all_orth == other_side, (
            f"seed_side={seed_side}, other={other_side}: orth mismatch"
        )
        assert typed[expected_rel] == other_side, (
            f"seed_side={seed_side}, other={other_side}: "
            f"expected typed[{expected_rel}]={other_side}, got {typed}"
        )
        for rel in ("one2one", "one2many", "many2one", "many2many"):
            if rel != expected_rel:
                assert typed[rel] == set(), (
                    f"seed_side={seed_side}, other={other_side}: "
                    f"typed[{rel}] should be empty, got {typed[rel]}"
                )


def test_ortholog_meta_records_best_event_per_ortholog(engine_no_filter):
    """When the same ortholog appears in multiple events the engine sees
    while collecting for one seed, `ortholog_meta` records the BEST
    event — lowest cascade sort key (in_seed_lineage first, then deeper
    ev_lca, then lower type tier). Without a lineage cache the in_lineage
    flag is False everywhere and depth is 0, so the tie-break collapses
    to type tier alone."""
    seed = 1
    other = 100
    # Three events, all covering (seed, other), with varying types and ev_lcas:
    #   - event 10: 1:1 (best type)
    #   - event 11: many:many (worst type), ev_lca closer
    #   - event 12: 1:many
    events = {
        10: {"side1": {seed}, "side2": {other},          "ev_lca": "9606"},
        11: {"side1": {seed, 5, 6}, "side2": {other, 7}, "ev_lca": "33208"},
        12: {"side1": {seed}, "side2": {other, 8, 9},    "ev_lca": "131567"},
    }
    all_orth, typed, meta = engine_no_filter._collect_orthologs(seed, events)
    assert other in meta
    # No lineage cache on engine_no_filter -> depth=0 for all, in_lineage=False
    # all events. Tie-break is purely on type tier: 1:1 wins.
    assert meta[other]["type"] == "one2one"
    assert meta[other]["type_tier"] == 0
    assert meta[other]["event_id"] == 10


def test_ortholog_meta_keys_match_filtered_orthologs(engine_no_filter):
    """`ortholog_meta` is keyed on exactly the same set as `typed['all']`.
    This is the contract phase 6 (cascade) relies on: every ortholog in
    `all` has metadata; nothing in `meta` is absent from `all`."""
    seed = 1
    events = {
        20: {"side1": {seed}, "side2": {200, 201, 202}, "ev_lca": "9606"},
        21: {"side1": {seed, 5}, "side2": {300}, "ev_lca": "131567"},
    }
    all_orth, typed, meta = engine_no_filter._collect_orthologs(seed, events)
    assert set(meta.keys()) == typed["all"]


def test_seed_in_both_sides_is_handled_deterministically(engine_no_filter):
    """Defensive test: the loader emits disjoint sides (verified by
    `test_side1_and_side2_are_disjoint` on the builder side). If a future
    bug let a seed appear on both sides of an event, _collect_orthologs
    picks side1 first (line 386-389 of annotate.py). Pin the behavior so
    a regression here changes the test, not silently the data."""
    seed = 7
    events = {0: {"side1": {seed, 100}, "side2": {seed, 200, 300}}}
    # If side1 wins, seed_side=={7,100}, other={7,200,300} -> after
    # discarding seed: {200, 300}, type=many2many (s=2, o=3).
    all_orth, typed, _ = engine_no_filter._collect_orthologs(seed, events)
    assert all_orth == {200, 300}
    assert typed["many2many"] == {200, 300}
    for rel in ("one2one", "one2many", "many2one"):
        assert typed[rel] == set()
