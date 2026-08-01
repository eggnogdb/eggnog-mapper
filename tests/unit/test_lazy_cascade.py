"""Byte-identity verification for the lazy ``closest`` cascade.

Proves that :meth:`AnnotationEngine._lazy_cascade_summarize_batch` returns the
*exact* ``(annotations, confidence)`` that the eager
:meth:`AnnotationEngine._summarize_annotations` produces for
``donor_pool="closest"`` — including dict key sets, value lists, and value
*order* — while fetching only the orthologs the cascade actually consumes.

``_summarize_annotations`` is pure (operates on in-memory dicts), so both the
eager and lazy paths are driven here over synthetic ``prots`` rows and
synthetic ``ortholog_meta`` with a mock, fetch-counting DB. No real eggnog.db,
psutil, or DIAMOND is required.

Runnable two ways::

    python3 tests/unit/test_lazy_cascade.py     # standalone (no pytest dep)
    pytest tests/unit/test_lazy_cascade.py       # if pytest is installed
"""

from __future__ import annotations

import array
import random
from typing import Any, Dict, List, Optional, Set, Tuple

from eggnogmapper.annotator.e7.annotate import (
    AnnotationEngine,
    _GO_NS_FIELDS,
    _LAZY_MAX_ROUNDS,
    _presence_mask_from_parsed,
)


# --------------------------------------------------------------------------- #
# Mock DB + engine harness
# --------------------------------------------------------------------------- #
class CountingDB:
    """Mock ``prots`` source that counts distinct ortholog fetches."""

    def __init__(self, prots: Dict[int, dict]) -> None:
        self.prots = prots
        self.distinct_fetched: Set[int] = set()
        self.n_queries = 0

    def get_protein_annotations_bulk(self, ids: List[int]) -> Dict[int, dict]:
        self.n_queries += 1
        self.distinct_fetched.update(ids)
        return {oid: dict(self.prots[oid]) for oid in ids if oid in self.prots}


def make_engine(go_ns_map: Optional[Dict[str, str]], db: CountingDB) -> AnnotationEngine:
    """Build an AnnotationEngine bypassing __init__ (no real DB needed)."""
    eng = AnnotationEngine.__new__(AnnotationEngine)
    eng.db = db
    eng.donor_pool = "closest"
    eng.lazy_cascade = True
    eng._go_namespace_map = go_ns_map
    return eng


# --------------------------------------------------------------------------- #
# Eager reference + lazy runner — mirror the engine's phase 4/6 orchestration
# --------------------------------------------------------------------------- #
def _pop_seed_pnames(eng: AnnotationEngine, parsed: dict, seeds: List[int]) -> None:
    for s in seeds:
        row = parsed.get(s)
        if row:
            row.pop("pname", None)


def run_eager(
    eng: AnnotationEngine,
    db: CountingDB,
    seeds: List[int],
    meta_map: Dict[int, Dict[int, dict]],
    target: str,
) -> Dict[int, Tuple[dict, dict]]:
    """Reproduce the eager phase-4/6 path exactly (fetch everything up front)."""
    all_orthologs: Set[int] = set()
    for s in seeds:
        all_orthologs.update(meta_map[s].keys())
    annot = db.get_protein_annotations_bulk(list(all_orthologs | set(seeds)))
    parsed = eng._pre_parse_batch(annot)
    _pop_seed_pnames(eng, parsed, seeds)

    out: Dict[int, Tuple[dict, dict]] = {}
    for s in seeds:
        meta = meta_map[s]
        seed_annots = {oid: annot[oid] for oid in meta.keys() if oid in annot}
        if seed_annots:
            out[s] = eng._summarize_annotations(
                seed_annots, meta, target, parsed=parsed, donor_pool="closest"
            )
        else:
            out[s] = ({}, {})
    return out


def run_lazy(
    eng: AnnotationEngine,
    db: CountingDB,
    seeds: List[int],
    meta_map: Dict[int, Dict[int, dict]],
    target: str,
) -> Dict[int, Tuple[dict, dict]]:
    """Reproduce the lazy phase-4/4b path (seed rows now, orthologs on demand)."""
    annot = db.get_protein_annotations_bulk(list(set(seeds)))
    parsed = eng._pre_parse_batch(annot)
    _pop_seed_pnames(eng, parsed, seeds)
    fetched_ids: Set[int] = set(seeds)
    return eng._lazy_cascade_summarize_batch(
        seeds, meta_map, target, annot, parsed, fetched_ids
    )


def build_masks(eng: AnnotationEngine, prots: Dict[int, dict]) -> "array.array":
    """Build the presence-mask array from the synthetic prots table.

    Uses the engine's own build path (which calls ``_pre_parse_batch`` with the
    engine's GO map), exactly as the load-time DB scan would. Mask building does
    NOT touch the counting DB, mirroring that it happens once at load time and
    is not part of the per-query fetch cost.
    """
    n = (max(prots) + 1) if prots else 0
    return eng.build_field_presence((dict(r) for r in prots.values()), n)


def run_masked(
    eng: AnnotationEngine,
    db: CountingDB,
    seeds: List[int],
    meta_map: Dict[int, Dict[int, dict]],
    target: str,
) -> Dict[int, Tuple[dict, dict]]:
    """Reproduce the mask-gated phase-4/4b path (single bulk fetch of winners)."""
    annot = db.get_protein_annotations_bulk(list(set(seeds)))
    parsed = eng._pre_parse_batch(annot)
    _pop_seed_pnames(eng, parsed, seeds)
    fetched_ids: Set[int] = set(seeds)
    return eng._masked_cascade_summarize_batch(
        seeds, meta_map, target, annot, parsed, fetched_ids
    )


def assert_identical(eager: dict, lazy: dict, label: str) -> None:
    """Assert byte-identity of the per-seed (annotations, confidence)."""
    assert set(eager) == set(lazy), f"[{label}] seed set mismatch"
    for s in eager:
        ea, ec = eager[s]
        la, lc = lazy[s]
        assert set(ea) == set(la), (
            f"[{label}] seed {s} annotation KEY set mismatch: "
            f"{sorted(ea)} != {sorted(la)}"
        )
        for k in ea:
            assert ea[k] == la[k], (
                f"[{label}] seed {s} field {k!r} VALUE/ORDER mismatch: "
                f"{ea[k]!r} != {la[k]!r}"
            )
        assert ec == lc, (
            f"[{label}] seed {s} confidence mismatch: {ec!r} != {lc!r}"
        )


# --------------------------------------------------------------------------- #
# Synthetic-data helpers
# --------------------------------------------------------------------------- #
RAW_FIELDS = [
    "pname", "kegg_ko", "kegg_ec", "kegg_pathway", "kegg_module",
    "kegg_reaction", "kegg_rclass", "kegg_brite", "kegg_tc", "kegg_cazy",
    "bigg_reaction", "pfam",
]

# A small GO vocabulary with a fixed namespace map.
GO_MF = ["GO:0003674", "GO:0003700", "GO:0016787"]
GO_BP = ["GO:0008150", "GO:0006355", "GO:0055085"]
GO_CC = ["GO:0005575", "GO:0005634", "GO:0016020"]
GO_NS_MAP = {g: "molecular_function" for g in GO_MF}
GO_NS_MAP.update({g: "biological_process" for g in GO_BP})
GO_NS_MAP.update({g: "cellular_component" for g in GO_CC})
ALL_GO = GO_MF + GO_BP + GO_CC


def meta_entry(mtype: str, depth: int, in_lineage: bool) -> dict:
    return {
        "event_id": -1,
        "ev_lca": "",
        "og_lca": "",
        "type": mtype,
        "type_tier": AnnotationEngine.TYPE_TIERS[mtype],
        "depth": depth,
        "in_seed_lineage": in_lineage,
    }


def _ortholog_fetches(db: CountingDB, seeds: List[int]) -> int:
    """Distinct *ortholog* full-fetches (seeds are always fetched — exclude)."""
    return len(db.distinct_fetched - set(seeds))


def check(
    label: str,
    prots: Dict[int, dict],
    seeds: List[int],
    meta_map: Dict[int, Dict[int, dict]],
    target: str,
    go_ns_map: Optional[Dict[str, str]],
) -> Tuple[int, int, int]:
    """Run eager, staged-lazy and masked paths; assert all three identical.

    Returns ``(eager_ortholog_fetches, staged_fetches, masked_fetches)`` — the
    number of distinct *ortholog* rows each path pulled (seeds excluded, since
    every path always fetches the seed rows).
    """
    db_e = CountingDB(prots)
    db_l = CountingDB(prots)
    db_m = CountingDB(prots)
    eng_e = make_engine(go_ns_map, db_e)
    eng_l = make_engine(go_ns_map, db_l)
    eng_m = make_engine(go_ns_map, db_m)
    eng_m.field_presence = build_masks(eng_m, prots)

    eager = run_eager(eng_e, db_e, seeds, meta_map, target)
    lazy = run_lazy(eng_l, db_l, seeds, meta_map, target)
    masked = run_masked(eng_m, db_m, seeds, meta_map, target)

    assert_identical(eager, lazy, f"{label}/staged")
    assert_identical(eager, masked, f"{label}/masked")

    e_n = _ortholog_fetches(db_e, seeds)
    l_n = _ortholog_fetches(db_l, seeds)
    m_n = _ortholog_fetches(db_m, seeds)
    # Neither lazy path may fetch MORE orthologs than eager.
    assert l_n <= e_n, f"[{label}] staged fetched more ({l_n}) than eager ({e_n})"
    assert m_n <= e_n, f"[{label}] masked fetched more ({m_n}) than eager ({e_n})"
    return e_n, l_n, m_n


# --------------------------------------------------------------------------- #
# Hand-built edge-case scenarios
# --------------------------------------------------------------------------- #
def test_field_winners_at_different_tiers() -> None:
    # seed self-donor (tier0) has kegg_ko; a tier1 ortholog has pfam; a
    # tier2 ortholog has kegg_ec. Each field wins at a different bucket.
    seed = 100
    prots = {
        seed: {"id": seed, "pname": "SEEDGENE", "kegg_ko": "K00001"},
        201: {"id": 201, "pfam": "PF00001"},
        202: {"id": 202, "kegg_ec": "1.1.1.1"},
    }
    meta = {
        seed: meta_entry("self", depth=50, in_lineage=True),
        201: meta_entry("one2many", depth=30, in_lineage=True),
        202: meta_entry("many2many", depth=10, in_lineage=True),
    }
    check("winners-different-tiers", prots, [seed], {seed: meta}, "all", None)


def test_only_bucket_is_last() -> None:
    # A field (pfam) whose only non-empty donor is the deepest/worst bucket.
    seed = 1
    prots = {
        seed: {"id": seed, "kegg_ko": "K1"},
        2: {"id": 2},  # empty ortholog row
        3: {"id": 3, "pfam": "PF9"},
    }
    meta = {
        seed: meta_entry("self", depth=50, in_lineage=True),
        2: meta_entry("one2one", depth=40, in_lineage=True),
        3: meta_entry("many2many", depth=5, in_lineage=False),
    }
    check("only-bucket-last", prots, [seed], {seed: meta}, "all", None)


def test_empty_fields_and_single_bucket() -> None:
    seed = 7
    prots = {seed: {"id": seed, "pname": "G", "kegg_ko": "K"}}
    meta = {seed: meta_entry("self", depth=12, in_lineage=True)}
    check("single-bucket-seed-only", prots, [seed], {seed: meta}, "all", None)


def test_seed_self_donor_only() -> None:
    # Seed present with all fields; no orthologs at all.
    seed = 9
    prots = {seed: {"id": seed, "pname": "X", "gos": ",".join(ALL_GO), "pfam": "PF1"}}
    meta = {seed: meta_entry("self", depth=20, in_lineage=True)}
    check("self-donor-only-nsmap", prots, [seed], {seed: meta}, "all", GO_NS_MAP)
    check("self-donor-only-legacy", prots, [seed], {seed: meta}, "all", None)


def test_go_namespaces_win_at_different_tiers() -> None:
    # MF from tier0 seed, BP from tier1, CC from tier2 — all merge into GOs
    # with best (lowest) tier = 0.
    seed = 50
    prots = {
        seed: {"id": seed, "gos": ",".join(GO_MF)},
        61: {"id": 61, "gos": ",".join(GO_BP)},
        62: {"id": 62, "gos": ",".join(GO_CC)},
    }
    meta = {
        seed: meta_entry("self", depth=40, in_lineage=True),
        61: meta_entry("one2many", depth=30, in_lineage=True),
        62: meta_entry("many2many", depth=20, in_lineage=True),
    }
    check("go-ns-different-tiers", prots, [seed], {seed: meta}, "all", GO_NS_MAP)


def test_legacy_combined_go() -> None:
    seed = 70
    prots = {
        seed: {"id": seed, "gos": ",".join(GO_MF[:1])},
        71: {"id": 71, "gos": ",".join(GO_BP)},
    }
    meta = {
        seed: meta_entry("self", depth=40, in_lineage=True),
        71: meta_entry("one2one", depth=30, in_lineage=True),
    }
    check("legacy-go", prots, [seed], {seed: meta}, "all", None)


def test_target_orthologs_floor() -> None:
    # one2one floor excludes many2many donors: pfam only on a many2many donor
    # must NOT be emitted (bucket filtered out) — same in both paths.
    seed = 80
    prots = {
        seed: {"id": seed, "kegg_ko": "K"},
        81: {"id": 81, "pfam": "PF_EXCLUDED"},
    }
    meta = {
        seed: meta_entry("self", depth=40, in_lineage=True),
        81: meta_entry("many2many", depth=30, in_lineage=True),
    }
    for tgt in ("all", "one2one", "one2many", "many2one", "many2many"):
        check(f"floor-{tgt}", prots, [seed], {seed: meta}, tgt, None)


def test_priority_key_ties() -> None:
    # Two orthologs share the exact same priority key (same in_lineage/depth/
    # tier) and both carry pname — contributor ORDER must be preserved for the
    # Counter tie-break.
    seed = 90
    prots = {
        seed: {"id": seed},  # seed has no pname of its own
        91: {"id": 91, "pname": "ALPHA", "kegg_ko": "K1"},
        92: {"id": 92, "pname": "BETA", "kegg_ko": "K2"},
    }
    meta = {
        seed: meta_entry("self", depth=99, in_lineage=True),
        91: meta_entry("one2one", depth=30, in_lineage=True),
        92: meta_entry("one2one", depth=30, in_lineage=True),
    }
    check("priority-ties", prots, [seed], {seed: meta}, "all", None)


def test_missing_prots_row_and_ortholog_is_seed() -> None:
    # oid 202 has no prots row (excluded from buckets). oid 300 is BOTH a seed
    # and an ortholog of seed 301 — its pname must be popped in both paths.
    seeds = [300, 301]
    prots = {
        300: {"id": 300, "pname": "SHARED", "kegg_ko": "K300"},
        301: {"id": 301, "pfam": "PF301"},
        # 202 deliberately absent from prots
    }
    meta_map = {
        300: {300: meta_entry("self", depth=50, in_lineage=True)},
        301: {
            301: meta_entry("self", depth=50, in_lineage=True),
            300: meta_entry("one2one", depth=30, in_lineage=True),
            202: meta_entry("many2many", depth=10, in_lineage=True),
        },
    }
    check("missing-row-and-shared-seed", prots, seeds, meta_map, "all", None)


def test_mask_derivation_is_exact() -> None:
    # The presence mask bit for a field MUST equal "field in the parsed row",
    # for both the namespace-split and legacy-combined GO configs. This is the
    # invariant the whole gate relies on.
    rng = random.Random(11)
    for go_ns_map in (GO_NS_MAP, None):
        eng = make_engine(go_ns_map, CountingDB({}))
        for _ in range(500):
            oid = rng.randint(1, 999)
            row = _rand_row(oid, rng)
            parsed = eng._pre_parse_batch({oid: row})[oid]
            mask = _presence_mask_from_parsed(parsed)
            for field, bit in ((f, 1 << i) for i, f in enumerate(
                list(AnnotationEngine.ANNOTATION_FIELDS)
                + [AnnotationEngine.LEGACY_GO_FIELD]
            )):
                assert bool(mask & bit) == (field in parsed), (
                    f"mask/parsed disagree on {field!r} for row {row!r}"
                )


def test_sparse_fields_cost_zero_fetches() -> None:
    # A seed with only pname+GOs; several deep orthologs that carry ONLY sparse
    # fields the seed also lacks would, in a naive descent, be fetched to prove
    # emptiness. Here every present field of every ortholog is ALSO absent (the
    # orthologs are empty rows), so nothing is emitted beyond the seed and the
    # masked path must fetch ZERO orthologs.
    seed = 2000
    prots = {seed: {"id": seed, "pname": "SEEDONLY", "gos": ",".join(GO_MF)}}
    meta = {seed: meta_entry("self", depth=99, in_lineage=True)}
    for oid in range(2100, 2130):
        prots[oid] = {"id": oid}  # present but no annotation fields
        meta[oid] = meta_entry("many2many", depth=oid - 2100 + 1, in_lineage=False)
    e_n, l_n, m_n = check("sparse-zero", prots, [seed], {seed: meta}, "all", GO_NS_MAP)
    assert m_n == 0, f"sparse-zero: masked should fetch 0 orthologs, got {m_n}"


def test_sparse_field_only_in_deepest_bucket() -> None:
    # CAZy present ONLY on the deepest donor; every other field either on the
    # seed or absent. Masked path fetches exactly that one deep donor (plus any
    # other winners) — proving it can reach a deep winner without fetching the
    # intervening empties, AND doesn't fetch everything.
    seed = 3000
    prots = {
        seed: {"id": seed, "gos": ",".join(GO_MF), "kegg_ko": "K1"},
        3100: {"id": 3100},  # empty intermediate
        3101: {"id": 3101},  # empty intermediate
        3200: {"id": 3200, "kegg_cazy": "GT2"},  # deepest, sole CAZy donor
    }
    meta = {
        seed: meta_entry("self", depth=99, in_lineage=True),
        3100: meta_entry("one2one", depth=60, in_lineage=True),
        3101: meta_entry("one2many", depth=40, in_lineage=True),
        3200: meta_entry("many2many", depth=5, in_lineage=False),
    }
    e_n, l_n, m_n = check(
        "sparse-deep", prots, [seed], {seed: meta}, "all", GO_NS_MAP
    )
    # Masked fetches only the CAZy donor (3200): 1 ortholog, not 3.
    assert m_n == 1, f"sparse-deep: masked should fetch 1 ortholog, got {m_n}"
    assert e_n == 3, f"sparse-deep: eager should fetch 3 orthologs, got {e_n}"


def test_pname_wins_but_cazy_absent() -> None:
    # pname wins from a tier-0 ortholog; CAZy/BiGG/etc. absent everywhere. The
    # absent sparse fields must not drag in the deep orthologs.
    seed = 4000
    # Seed (self-donor, bucket 0) already carries GOs + kegg_ko + pfam; only
    # pname is missing (and it is always popped from the seed anyway).
    prots = {
        seed: {"id": seed, "gos": ",".join(ALL_GO), "kegg_ko": "KSEED",
               "pfam": "PFSEED"},
        4100: {"id": 4100, "pname": "REALNAME"},  # sole pname donor, bucket 1
    }
    meta = {seed: meta_entry("self", depth=99, in_lineage=True),
            4100: meta_entry("one2one", depth=70, in_lineage=True)}
    # Deep orthologs carry only fields the seed already won at bucket 0, plus
    # CAZy/BiGG are absent everywhere. None of them should ever be fetched.
    for oid in range(4200, 4220):
        prots[oid] = {"id": oid, "kegg_ko": "KDEEP"}
        meta[oid] = meta_entry("many2many", depth=3, in_lineage=False)
    e_n, l_n, m_n = check(
        "pname-cazy-absent", prots, [seed], {seed: meta}, "all", GO_NS_MAP
    )
    # Masked fetches only the pname donor; the 20 deep orthologs are skipped
    # (every field they could offer is already won at bucket 0 or absent).
    assert m_n == 1, f"pname-cazy-absent: masked fetched {m_n}, expected 1"
    assert m_n < e_n


def test_early_winner_fetches_fewer() -> Tuple[int, int, int]:
    # Bucket 0 (self-donor, tier0) carries every field EXCEPT pname (the seed's
    # own pname is always popped); bucket 1 (a close one2one ortholog) supplies
    # pname. So every field decides by bucket 1, and the many deep orthologs in
    # bucket 2 are never fetched by the lazy path.
    seed = 500
    full_row = {"id": seed, "gos": ",".join(ALL_GO)}
    for f in RAW_FIELDS:
        if f != "pname":
            full_row[f] = f.upper().replace("_", "") + "1"
    prots = {seed: full_row}
    meta = {seed: meta_entry("self", depth=99, in_lineage=True)}
    # Bucket 1: closest ortholog supplies the only pname.
    prots[550] = {"id": 550, "pname": "PNAME"}
    meta[550] = meta_entry("one2one", depth=60, in_lineage=True)
    # Bucket 2 (deepest): many orthologs the cascade never needs.
    for oid in range(600, 640):
        prots[oid] = {"id": oid, "pfam": "PFX", "kegg_ko": "KX"}
        meta[oid] = meta_entry("many2many", depth=5, in_lineage=False)
    e_n, l_n, m_n = check(
        "early-winner", prots, [seed], {seed: meta}, "all", GO_NS_MAP
    )
    # The deep bucket (index 2) is only skippable by the *staged* path when the
    # round cap allows reaching bucket 1 lazily (needs >= 2 staged rounds).
    if _LAZY_MAX_ROUNDS >= 2:
        assert l_n < e_n, f"early-winner staged: eager={e_n} staged={l_n}"
    # The masked path always fetches only the one winning donor (pname), never
    # the 40 deep orthologs — independent of the round cap.
    assert m_n < e_n, f"early-winner masked: eager={e_n} masked={m_n}"
    return e_n, l_n, m_n


# --------------------------------------------------------------------------- #
# Randomized fuzz
# --------------------------------------------------------------------------- #
def _rand_value(field: str, rng: random.Random) -> str:
    return f"{field[:3].upper()}{rng.randint(0, 4)}"


def _rand_row(oid: int, rng: random.Random) -> dict:
    row: dict = {"id": oid}
    for f in RAW_FIELDS:
        if rng.random() < 0.35:
            if f == "pname":
                row[f] = f"NAME{rng.randint(0, 3)}"
            else:
                n = rng.randint(1, 3)
                row[f] = ",".join(_rand_value(f, rng) for _ in range(n))
    if rng.random() < 0.6:
        gos = rng.sample(ALL_GO, rng.randint(1, 5))
        row["gos"] = ",".join(gos)
    return row


def test_fuzz(n_seeds: int = 400, seed_rng: int = 20260731) -> Dict[str, int]:
    rng = random.Random(seed_rng)
    tot = {"eager": 0, "staged": 0, "masked": 0, "staged_fewer": 0, "masked_fewer": 0}
    for i in range(n_seeds):
        go_ns_map = GO_NS_MAP if rng.random() < 0.5 else None
        target = rng.choice(
            ["all", "one2one", "one2many", "many2one", "many2many"]
        )
        seed = 1_000_000 + i * 100
        prots: Dict[int, dict] = {}
        meta: Dict[int, dict] = {}
        # seed self-donor row (may or may not exist in prots)
        if rng.random() < 0.9:
            prots[seed] = _rand_row(seed, rng)
        meta[seed] = meta_entry("self", depth=90 + rng.randint(0, 5), in_lineage=True)
        n_orth = rng.randint(0, 8)
        types = ["one2one", "one2many", "many2one", "many2many"]
        for j in range(n_orth):
            oid = seed + 1 + j
            # ~15% of orthologs have no prots row at all.
            if rng.random() < 0.85:
                prots[oid] = _rand_row(oid, rng)
            mtype = rng.choice(types)
            depth = rng.randint(1, 80)
            in_lin = rng.random() < 0.6
            meta[oid] = meta_entry(mtype, depth, in_lin)
        e, l, m = check(f"fuzz-{i}", prots, [seed], {seed: meta}, target, go_ns_map)
        tot["eager"] += e
        tot["staged"] += l
        tot["masked"] += m
        tot["staged_fewer"] += int(l < e)
        tot["masked_fewer"] += int(m < e)
    return tot


def test_fuzz_multiseed_batch(n_batches: int = 60, seed_rng: int = 424242) -> None:
    """Fuzz whole sub-batches (cross-seed shared orthologs, staged rounds)."""
    rng = random.Random(seed_rng)
    for b in range(n_batches):
        go_ns_map = GO_NS_MAP if rng.random() < 0.5 else None
        target = rng.choice(["all", "one2one", "many2many"])
        n_seeds = rng.randint(2, 6)
        seeds: List[int] = []
        prots: Dict[int, dict] = {}
        meta_map: Dict[int, Dict[int, dict]] = {}
        # A shared ortholog pool so seeds overlap (cross-seed cache reuse).
        pool = list(range(9000, 9020))
        for oid in pool:
            if rng.random() < 0.8:
                prots[oid] = _rand_row(oid, rng)
        for k in range(n_seeds):
            seed = 8000 + b * 10 + k
            seeds.append(seed)
            if rng.random() < 0.9:
                prots[seed] = _rand_row(seed, rng)
            meta = {seed: meta_entry("self", depth=95, in_lineage=True)}
            for oid in rng.sample(pool, rng.randint(0, 8)):
                if oid == seed:
                    continue
                mtype = rng.choice(["one2one", "one2many", "many2one", "many2many"])
                meta[oid] = meta_entry(mtype, rng.randint(1, 80), rng.random() < 0.6)
            meta_map[seed] = meta
        check(f"fuzz-batch-{b}", prots, seeds, meta_map, target, go_ns_map)


# --------------------------------------------------------------------------- #
# Standalone runner
# --------------------------------------------------------------------------- #
def _main() -> None:
    print("Running lazy-cascade byte-identity verification "
          "(eager vs staged-lazy vs mask-gated)...\n")
    edge_tests = [
        test_mask_derivation_is_exact,
        test_field_winners_at_different_tiers,
        test_only_bucket_is_last,
        test_empty_fields_and_single_bucket,
        test_seed_self_donor_only,
        test_go_namespaces_win_at_different_tiers,
        test_legacy_combined_go,
        test_target_orthologs_floor,
        test_priority_key_ties,
        test_missing_prots_row_and_ortholog_is_seed,
        test_sparse_fields_cost_zero_fetches,
        test_sparse_field_only_in_deepest_bucket,
        test_pname_wins_but_cazy_absent,
    ]
    for t in edge_tests:
        t()
        print(f"  PASS  {t.__name__}")

    e_n, l_n, m_n = test_early_winner_fetches_fewer()
    print(f"  PASS  test_early_winner_fetches_fewer "
          f"(orthologs fetched: eager={e_n}, staged={l_n}, masked={m_n})")

    tot = test_fuzz()
    print(
        f"  PASS  test_fuzz (400 random seeds) identical every time.\n"
        f"          ortholog full-fetches: eager={tot['eager']} "
        f"staged={tot['staged']} masked={tot['masked']}\n"
        f"          seeds with fewer fetches: "
        f"staged={tot['staged_fewer']}/400 masked={tot['masked_fewer']}/400"
    )

    test_fuzz_multiseed_batch()
    print("  PASS  test_fuzz_multiseed_batch (60 random sub-batches)")

    print("\nALL CHECKS PASSED - staged and mask-gated cascades are "
          "byte-identical to eager.")


if __name__ == "__main__":
    _main()
