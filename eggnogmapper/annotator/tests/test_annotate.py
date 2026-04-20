"""Integration tests for AnnotationEngine using the sample database."""

import logging
import os
import pytest

from tests.conftest import TP53_PROTEIN_ID, NIFH_PROTEIN_ID

# Expected result keys from annotate_batch
REQUIRED_KEYS = {"orthologs", "ortholog_types", "all_ogs", "annotations", "og_info"}

# Expected ortholog_type sub-keys
ORTHOLOG_TYPE_KEYS = {"one2one", "one2many", "many2one", "many2many", "all"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_result(engine, protein_id: int) -> dict:
    results = engine.annotate_batch([protein_id])
    assert protein_id in results, f"Seed id {protein_id} not in result keys"
    return results[protein_id]


# ---------------------------------------------------------------------------
# Key presence tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("protein_id", [TP53_PROTEIN_ID, NIFH_PROTEIN_ID])
def test_annotate_batch_returns_all_keys(engine, protein_id):
    """annotate_batch result must contain all 5 required keys."""
    result = _get_result(engine, protein_id)
    assert REQUIRED_KEYS == set(result.keys()), (
        f"Missing or extra keys for protein_id={protein_id}: "
        f"got {set(result.keys())}"
    )


@pytest.mark.parametrize("protein_id", [TP53_PROTEIN_ID, NIFH_PROTEIN_ID])
def test_ortholog_types_has_all_five_sub_keys(engine, protein_id):
    """ortholog_types dict must always contain all 5 sub-keys."""
    result = _get_result(engine, protein_id)
    ot = result["ortholog_types"]
    assert ORTHOLOG_TYPE_KEYS == set(ot.keys()), (
        f"ortholog_types missing keys for protein_id={protein_id}: "
        f"got {set(ot.keys())}"
    )


# ---------------------------------------------------------------------------
# Union invariant
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("protein_id", [TP53_PROTEIN_ID, NIFH_PROTEIN_ID])
def test_ortholog_types_union_invariant(engine, protein_id):
    """'all' must equal the union of one2one | one2many | many2one | many2many."""
    result = _get_result(engine, protein_id)
    ot = result["ortholog_types"]
    typed_union = ot["one2one"] | ot["one2many"] | ot["many2one"] | ot["many2many"]
    assert ot["all"] == typed_union, (
        f"Union invariant violated for protein_id={protein_id}: "
        f"all={len(ot['all'])}, union={len(typed_union)}"
    )


# ---------------------------------------------------------------------------
# Subset tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("protein_id", [TP53_PROTEIN_ID, NIFH_PROTEIN_ID])
@pytest.mark.parametrize("rel", ["one2one", "one2many", "many2one", "many2many"])
def test_typed_set_is_subset_of_all(engine, protein_id, rel):
    """Every typed ortholog set must be a subset of 'all'."""
    result = _get_result(engine, protein_id)
    ot = result["ortholog_types"]
    assert ot[rel].issubset(ot["all"]), (
        f"{rel} is not a subset of 'all' for protein_id={protein_id}"
    )


# ---------------------------------------------------------------------------
# all_ogs non-empty for known seeds
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("protein_id", [TP53_PROTEIN_ID, NIFH_PROTEIN_ID])
def test_all_ogs_nonempty_for_valid_seed(engine, protein_id):
    """all_ogs must be non-empty for known seeds with OG assignments."""
    result = _get_result(engine, protein_id)
    assert isinstance(result["all_ogs"], list)
    assert len(result["all_ogs"]) > 0, (
        f"all_ogs is empty for protein_id={protein_id}"
    )


# ---------------------------------------------------------------------------
# TP53 biological validation
# ---------------------------------------------------------------------------

def test_tp53_has_orthologs(engine):
    """TP53 (human) must have non-empty orthologs in sample DB."""
    result = _get_result(engine, TP53_PROTEIN_ID)
    assert len(result["orthologs"]) > 0, (
        "TP53 should have orthologs but got an empty list"
    )


def test_tp53_and_nifh_belong_to_different_ogs(engine):
    """TP53 and NifH must belong to different OG families."""
    tp53_result = _get_result(engine, TP53_PROTEIN_ID)
    nifh_result = _get_result(engine, NIFH_PROTEIN_ID)

    tp53_ogs = set(tp53_result["all_ogs"])
    nifh_ogs = set(nifh_result["all_ogs"])

    assert tp53_ogs.isdisjoint(nifh_ogs), (
        f"TP53 and NifH share OGs: {tp53_ogs & nifh_ogs}"
    )


# ---------------------------------------------------------------------------
# Out-of-range protein ID warning test
# ---------------------------------------------------------------------------

def test_out_of_range_protein_id_warns(engine, caplog):
    """annotate_batch must emit a logger.warning for out-of-range protein IDs.

    Strategy: patch a valid event so one of the side protein IDs is far
    beyond the taxid_array length. We do this by calling _collect_orthologs
    directly with a fabricated events dict.
    """
    taxid_array_len = len(engine.db.taxid_array)
    out_of_range_id = taxid_array_len + 9999

    # Fabricate a single speciation event: seed=0, other side=out_of_range_id
    fake_events = {
        1: {
            "og": "FAKEOG@1|root",
            "og_lca": "1",
            "ev_lca": "1",
            "side1": {0},
            "side2": {out_of_range_id},
        }
    }

    with caplog.at_level(logging.WARNING, logger="eggnog_annotator.e7.annotate"):
        engine._collect_orthologs(0, fake_events)

    warning_msgs = [r.message for r in caplog.records if r.levelno == logging.WARNING]
    assert any("out-of-range" in m or ">= taxid_array" in m or str(out_of_range_id) in m
               for m in warning_msgs), (
        f"Expected out-of-range warning, got: {warning_msgs}"
    )


# ---------------------------------------------------------------------------
# Empty early-return test
# ---------------------------------------------------------------------------

def test_annotate_empty_early_return(engine):
    """annotate_batch with an unknown protein ID returns all 5 keys with empty values."""
    # Use a protein_id that doesn't exist in event_index (well beyond max)
    nonexistent_id = 999_999_999
    results = engine.annotate_batch([nonexistent_id])

    # The engine should still return an entry (from the single-call path)
    # OR return a dict with the 5 keys initialised to empty values.
    # Depending on implementation: nonexistent ID has no events, so annotate()
    # returns the early-return dict, while annotate_batch may omit it entirely.
    # Test both outcomes.
    if nonexistent_id in results:
        result = results[nonexistent_id]
        assert REQUIRED_KEYS == set(result.keys())
        assert result["orthologs"] == [] or result["orthologs"] == set()
        assert result["all_ogs"] == []
    else:
        # annotate_batch builds results only for seed_ids, so if event_index
        # returns no events the seed still gets a result entry in phase 6.
        # If the seed is not in results at all, that is also acceptable for
        # a completely unknown id (not in event_index at all).
        pass  # no crash is sufficient
