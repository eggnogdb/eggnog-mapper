"""Phase 3B cascade behavior tests.

Verifies the per-source closest-ev_lca + type-priority donor cascade
implemented in `AnnotationEngine._summarize_annotations`. Each test
constructs a minimal synthetic ortholog_meta + annot_data and asserts
the cascade output and confidence labels.

Cascade sort key per ortholog (smaller = preferred):

    (in_seed_lineage_int, -ev_lca_depth, type_tier)

with ``type_tier`` collapsing 1:many and many:1 into the medium tier:

    one2one  → 0 (high)
    one2many → 1 (medium)
    many2one → 1 (medium)
    many2many → 2 (low)
"""

from __future__ import annotations

import pytest

from eggnog_annotator.e7 import AnnotationEngine


# ---------------------------------------------------------------------------
# Lightweight engine fixture (no DB connection — _summarize_annotations only
# reads its arguments).
# ---------------------------------------------------------------------------

class _DummyDB:
    taxid_array = None


@pytest.fixture
def engine():
    return AnnotationEngine(db=_DummyDB(), lineage_filter=None)


def _meta(type_, depth=0, in_lineage=False):
    """Build one ortholog_meta entry."""
    tier = AnnotationEngine.TYPE_TIERS[type_]
    return {
        "event_id": 0,
        "ev_lca": "0",
        "type": type_,
        "type_tier": tier,
        "depth": depth,
        "in_seed_lineage": in_lineage,
    }


# ---------------------------------------------------------------------------
# Tier ordering: 1:1 wins over 1:many at the same level
# ---------------------------------------------------------------------------

def test_one2one_wins_over_one2many_same_level(engine):
    annot_data = {
        100: {"kegg_ko": "K11111"},  # 1:1 donor
        200: {"kegg_ko": "K22222"},  # 1:many donor (same depth & in_lineage)
    }
    meta = {
        100: _meta("one2one",  depth=5, in_lineage=True),
        200: _meta("one2many", depth=5, in_lineage=True),
    }
    annotations, conf = engine._summarize_annotations(annot_data, meta)
    assert annotations.get("KEGG_ko") == ["K11111"]
    assert conf.get("KEGG_ko") == "high"


# ---------------------------------------------------------------------------
# Tier ordering: medium beats low at the same level
# ---------------------------------------------------------------------------

def test_one2many_wins_over_many2many_same_level(engine):
    annot_data = {
        200: {"kegg_ko": "K22222"},
        300: {"kegg_ko": "K33333"},
    }
    meta = {
        200: _meta("one2many",  depth=5, in_lineage=True),
        300: _meta("many2many", depth=5, in_lineage=True),
    }
    annotations, conf = engine._summarize_annotations(annot_data, meta)
    assert annotations.get("KEGG_ko") == ["K22222"]
    assert conf.get("KEGG_ko") == "medium"


# ---------------------------------------------------------------------------
# Closer ev_lca wins over farther even at lower type tier
# ---------------------------------------------------------------------------

def test_closer_ev_lca_wins_over_farther_higher_tier(engine):
    # Depth 7 (closer) at tier 2 (low) BEATS depth 3 (farther) at tier 0
    # (high), because the cascade sort key is
    # (in_lineage, -depth, type_tier), and -7 < -3.
    annot_data = {
        100: {"gos": "GO:0008150"},  # 1:1 but ancient ev_lca
        200: {"gos": "GO:0009987"},  # many:many but recent ev_lca
    }
    meta = {
        100: _meta("one2one",   depth=3, in_lineage=True),
        200: _meta("many2many", depth=7, in_lineage=True),
    }
    annotations, conf = engine._summarize_annotations(annot_data, meta)
    assert annotations.get("GOs") == ["GO:0009987"]
    assert conf.get("GOs") == "low"


# ---------------------------------------------------------------------------
# In-lineage events beat out-of-lineage even when out-of-lineage has more depth
# ---------------------------------------------------------------------------

def test_in_seed_lineage_outranks_out_of_lineage(engine):
    annot_data = {
        100: {"pname": "GENEA"},  # in seed lineage, depth=2
        200: {"pname": "GENEB"},  # out of lineage, depth=10 (irrelevant)
    }
    meta = {
        100: _meta("one2one", depth=2, in_lineage=True),
        200: _meta("one2one", depth=10, in_lineage=False),
    }
    annotations, conf = engine._summarize_annotations(annot_data, meta)
    assert annotations.get("Preferred_name") == ["GENEA"]
    assert conf.get("Preferred_name") == "high"


# ---------------------------------------------------------------------------
# Per-source independence: different sources can pick different tiers
# ---------------------------------------------------------------------------

def test_per_source_winners_are_independent(engine):
    # Donor A is a 1:1 with KEGG_ko but no GOs.
    # Donor B is a many:many with both KEGG_ko and GOs.
    # Cascade should: KEGG_ko -> A (high), GOs -> B (low).
    annot_data = {
        100: {"kegg_ko": "K00001"},
        200: {"kegg_ko": "K00002", "gos": "GO:0006950"},
    }
    meta = {
        100: _meta("one2one",   depth=5, in_lineage=True),
        200: _meta("many2many", depth=5, in_lineage=True),
    }
    annotations, conf = engine._summarize_annotations(annot_data, meta)
    assert annotations.get("KEGG_ko") == ["K00001"]
    assert conf.get("KEGG_ko") == "high"
    assert annotations.get("GOs") == ["GO:0006950"]
    assert conf.get("GOs") == "low"


# ---------------------------------------------------------------------------
# Missing source → cascade silently omits it
# ---------------------------------------------------------------------------

def test_missing_source_emits_nothing(engine):
    annot_data = {
        100: {"kegg_ko": "K00001"},  # no GOs anywhere
        200: {"kegg_ko": "K00002"},
    }
    meta = {
        100: _meta("one2one", depth=5, in_lineage=True),
        200: _meta("one2one", depth=5, in_lineage=True),
    }
    annotations, conf = engine._summarize_annotations(annot_data, meta)
    assert "GOs" not in annotations
    assert "GOs" not in conf
    assert annotations.get("KEGG_ko") == ["K00001", "K00002"] or \
           annotations.get("KEGG_ko") == ["K00002", "K00001"]


# ---------------------------------------------------------------------------
# target_orthologs floor: "one2one" excludes everything else
# ---------------------------------------------------------------------------

def test_target_orthologs_one2one_only(engine):
    annot_data = {
        100: {"kegg_ko": "K11111"},
        200: {"kegg_ko": "K22222"},
    }
    meta = {
        100: _meta("one2one",   depth=5, in_lineage=True),
        200: _meta("many2many", depth=5, in_lineage=True),
    }
    annotations, conf = engine._summarize_annotations(
        annot_data, meta, target_orthologs="one2one"
    )
    assert annotations.get("KEGG_ko") == ["K11111"]
    assert conf.get("KEGG_ko") == "high"


def test_target_orthologs_one2one_drops_floor_below(engine):
    """If only many:many donors exist and floor is one2one, no source
    is emitted. The cascade has nothing to fall through to."""
    annot_data = {200: {"kegg_ko": "K22222"}}
    meta = {200: _meta("many2many", depth=5, in_lineage=True)}
    annotations, conf = engine._summarize_annotations(
        annot_data, meta, target_orthologs="one2one"
    )
    assert annotations == {}
    assert conf == {}


def test_target_orthologs_one2many_keeps_one2one_and_one2many(engine):
    annot_data = {
        100: {"kegg_ko": "K1"},
        200: {"kegg_ko": "K2"},
        300: {"kegg_ko": "K3"},
        400: {"kegg_ko": "K4"},
    }
    meta = {
        100: _meta("one2one",   depth=5, in_lineage=True),
        200: _meta("one2many",  depth=5, in_lineage=True),
        300: _meta("many2one",  depth=5, in_lineage=True),  # excluded
        400: _meta("many2many", depth=5, in_lineage=True),  # excluded
    }
    annotations, conf = engine._summarize_annotations(
        annot_data, meta, target_orthologs="one2many"
    )
    # one2one wins (cascade halts there)
    assert annotations.get("KEGG_ko") == ["K1"]
    assert conf.get("KEGG_ko") == "high"


# ---------------------------------------------------------------------------
# Empty inputs
# ---------------------------------------------------------------------------

def test_empty_annot_data_returns_empty(engine):
    annotations, conf = engine._summarize_annotations({}, {})
    assert annotations == {}
    assert conf == {}


def test_empty_meta_falls_back_to_flat(engine):
    """No metadata → legacy flat aggregation (back-compat with single
    annotate() path and pre-cascade tests). Confidence dict is empty."""
    # Three orthologs all share KEGG_ko=K1 — flat aggregation just lists it.
    annot_data = {
        100: {"kegg_ko": "K1"},
        200: {"kegg_ko": "K1"},
        300: {"kegg_ko": "K1"},
    }
    annotations, conf = engine._summarize_annotations(annot_data, None)
    assert annotations.get("KEGG_ko") == ["K1"]
    assert conf == {}


# ---------------------------------------------------------------------------
# Multiple winning donors: consensus inside the bucket only
# ---------------------------------------------------------------------------

def test_consensus_across_bucket_only(engine):
    # Two 1:1 donors disagree on kegg_ko; the cascade emits both keys.
    # A many:many donor with a third KO is NOT included.
    annot_data = {
        100: {"kegg_ko": "K1"},
        200: {"kegg_ko": "K2"},
        300: {"kegg_ko": "K3"},
    }
    meta = {
        100: _meta("one2one",   depth=5, in_lineage=True),
        200: _meta("one2one",   depth=5, in_lineage=True),
        300: _meta("many2many", depth=5, in_lineage=True),
    }
    annotations, conf = engine._summarize_annotations(annot_data, meta)
    assert set(annotations.get("KEGG_ko", [])) == {"K1", "K2"}
    assert conf.get("KEGG_ko") == "high"


# ---------------------------------------------------------------------------
# Preferred_name: cascade picks the most-common pname in the winning bucket,
# without requiring the legacy ≥2-occurrence threshold.
# ---------------------------------------------------------------------------

def test_pname_winner_takes_most_common_no_threshold(engine):
    annot_data = {
        100: {"pname": "GENEA"},   # 1 occurrence in winning bucket — wins
        200: {"pname": "GENEB"},   # in lower bucket
    }
    meta = {
        100: _meta("one2one",   depth=5, in_lineage=True),
        200: _meta("many2many", depth=5, in_lineage=True),
    }
    annotations, conf = engine._summarize_annotations(annot_data, meta)
    assert annotations.get("Preferred_name") == ["GENEA"]
    assert conf.get("Preferred_name") == "high"
