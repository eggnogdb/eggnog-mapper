"""Unit tests for eggnogmapper.annotation.batch_annotate.

Tests focus on the annotate_batch() adapter function: hit parsing,
target_orthologs filtering, warning emission, and output tuple structure.
All annotation logic is mocked — only the adapter layer is exercised here.
"""

import logging
import pytest
from unittest.mock import MagicMock, patch


# ---------------------------------------------------------------------------
# Helpers to build fake hits and mock AnnotationEngine results
# ---------------------------------------------------------------------------

def _make_hit(seed_id_str, evalue=1e-10, score=200.0):
    """Build a minimal hit tuple: (query_name, best_hit_name, evalue, score)."""
    return ("query_protein", seed_id_str, str(evalue), str(score))


def _make_engine_result(
    seed_id: int,
    one2one=None,
    one2many=None,
    many2one=None,
    many2many=None,
):
    """Build a synthetic annotate_batch() result dict for a single seed."""
    one2one = one2one or {seed_id + 100, seed_id + 101}
    one2many = one2many or {seed_id + 200}
    many2one = many2one or {seed_id + 300}
    many2many = many2many or {seed_id + 400, seed_id + 401}
    all_ids = one2one | one2many | many2one | many2many

    return {
        "orthologs": sorted(all_ids),
        "ortholog_types": {
            "one2one": one2one,
            "one2many": one2many,
            "many2one": many2one,
            "many2many": many2many,
            "all": all_ids,
        },
        "all_ogs": ["CLU_TEST@131567|root", "CLU_TEST@9606|HX-1"],
        "annotations": {"GOs": ["GO:0005515"], "KEGG_ko": ["K00001"]},
        "og_info": {
            "name": "CLU_TEST@131567|root",
            "cog_cat": "O",
            "description": "Test OG description",
            "level": "131567",
        },
    }


def _make_batch_args(
    hit,
    target_orthologs="all",
    target_taxa=None,
    excluded_taxa=None,
    tax_scope_mode=None,
    tax_scope_ids=None,
    seed_evalue=1e-3,
    seed_score=60.0,
):
    """Build the args list that annotate_batch receives."""
    return [
        (
            hit,
            None,  # annotation (not pre-annotated)
            seed_score,
            seed_evalue,
            tax_scope_mode,
            tax_scope_ids,
            target_taxa,
            target_orthologs,
            excluded_taxa,
            None,  # go_evidence
            None,  # go_excluded
            None,  # data_dir
            None,  # annotation (legacy field, unused in batch path)
        )
    ]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestTargetOrthologs:
    """annotate_batch must respect target_orthologs filtering."""

    def test_one2one_applied(self):
        """With target_orthologs='one2one', annot_orthologs must only contain
        proteins from the one2one bucket, not from many2many."""
        seed_id = 1000
        one2one_ids = {seed_id + 100, seed_id + 101}
        many2many_ids = {seed_id + 400, seed_id + 401}

        engine_result = _make_engine_result(
            seed_id,
            one2one=one2one_ids,
            many2many=many2many_ids,
        )

        mock_engine = MagicMock()
        mock_engine.annotate_batch.return_value = {seed_id: engine_result}

        hit = _make_hit(str(seed_id))

        from eggnogmapper.annotation.batch_annotate import annotate_batch
        import eggnogmapper.annotation.batch_annotate as bat_mod

        mock_db = MagicMock()
        mock_db.conn = object()
        mock_db._taxids = []

        with patch.object(bat_mod, "_get_engine", return_value=mock_engine):
            results = list(
                annotate_batch(
                    batch=_make_batch_args(hit, target_orthologs="one2one"),
                    eggnog_db=mock_db,
                    annot=True,
                    target_orthologs="one2one",
                    target_taxa=None,
                    excluded_taxa=None,
                    tax_scope_mode=None,
                    tax_scope_ids=None,
                    go_evidence=None,
                    go_excluded=None,
                    seed_ortholog_score=60.0,
                    seed_ortholog_evalue=1e-3,
                )
            )

        assert len(results) == 1
        (returned_hit, annotation), already_present = results[0]
        assert annotation is not None, "Expected a non-None annotation"

        # Element [9] = annot_orthologs
        annot_orthologs = annotation[9]
        annot_set = set(annot_orthologs)

        assert annot_set <= one2one_ids, (
            f"annot_orthologs should be a subset of one2one_ids={one2one_ids}, "
            f"but got {annot_set}"
        )
        for mid in many2many_ids:
            assert mid not in annot_set, (
                f"many2many ID {mid} should not appear in one2one result"
            )

    def test_all_orthologs_included_when_target_is_all(self):
        """With target_orthologs='all', annot_orthologs contains all typed ids."""
        seed_id = 2000
        engine_result = _make_engine_result(seed_id)
        all_ids = engine_result["ortholog_types"]["all"]

        mock_engine = MagicMock()
        mock_engine.annotate_batch.return_value = {seed_id: engine_result}

        hit = _make_hit(str(seed_id))

        from eggnogmapper.annotation.batch_annotate import annotate_batch
        import eggnogmapper.annotation.batch_annotate as bat_mod

        mock_db = MagicMock()
        mock_db.conn = object()
        mock_db._taxids = []

        with patch.object(bat_mod, "_get_engine", return_value=mock_engine):
            results = list(
                annotate_batch(
                    batch=_make_batch_args(hit, target_orthologs="all"),
                    eggnog_db=mock_db,
                    annot=True,
                    target_orthologs="all",
                    target_taxa=None,
                    excluded_taxa=None,
                    tax_scope_mode=None,
                    tax_scope_ids=None,
                    go_evidence=None,
                    go_excluded=None,
                    seed_ortholog_score=60.0,
                    seed_ortholog_evalue=1e-3,
                )
            )

        assert len(results) == 1
        (_, annotation), _ = results[0]
        assert annotation is not None
        annot_set = set(annotation[9])
        assert annot_set == all_ids


class TestInvalidSeedId:
    """Non-integer seed IDs must emit a warning and yield (hit, None)."""

    def test_non_integer_hit_name_emits_warning(self, caplog):
        """Provide a non-integer best_hit_name; check logger.warning is called."""
        hit = _make_hit("NOT_AN_INT")

        from eggnogmapper.annotation.batch_annotate import annotate_batch
        import eggnogmapper.annotation.batch_annotate as bat_mod

        mock_db = MagicMock()
        mock_db.conn = object()
        mock_db._taxids = []

        with caplog.at_level(logging.WARNING, logger="eggnogmapper.annotation.batch_annotate"):
            results = list(
                annotate_batch(
                    batch=_make_batch_args(hit),
                    eggnog_db=mock_db,
                    annot=True,
                    target_orthologs="all",
                    target_taxa=None,
                    excluded_taxa=None,
                    tax_scope_mode=None,
                    tax_scope_ids=None,
                    go_evidence=None,
                    go_excluded=None,
                    seed_ortholog_score=60.0,
                    seed_ortholog_evalue=1e-3,
                )
            )

        warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
        assert any("NOT_AN_INT" in w or "integer" in w.lower() or "seed" in w.lower()
                   for w in warnings), (
            f"Expected warning about non-integer seed ID, got: {warnings}"
        )

        # The result should be (hit, None) — annotation skipped
        assert len(results) == 1
        (returned_hit, annotation), _ = results[0]
        assert annotation is None


class TestAnnotationTupleStructure:
    """The yielded annotation tuple must conform to the 10-element contract."""

    def _get_annotation(self, target_orthologs="all"):
        seed_id = 3000
        engine_result = _make_engine_result(seed_id)

        mock_engine = MagicMock()
        mock_engine.annotate_batch.return_value = {seed_id: engine_result}

        hit = _make_hit(str(seed_id))

        from eggnogmapper.annotation.batch_annotate import annotate_batch
        import eggnogmapper.annotation.batch_annotate as bat_mod

        mock_db = MagicMock()
        mock_db.conn = object()
        mock_db._taxids = []

        with patch.object(bat_mod, "_get_engine", return_value=mock_engine):
            results = list(
                annotate_batch(
                    batch=_make_batch_args(hit, target_orthologs=target_orthologs),
                    eggnog_db=mock_db,
                    annot=True,
                    target_orthologs=target_orthologs,
                    target_taxa=None,
                    excluded_taxa=None,
                    tax_scope_mode=None,
                    tax_scope_ids=None,
                    go_evidence=None,
                    go_excluded=None,
                    seed_ortholog_score=60.0,
                    seed_ortholog_evalue=1e-3,
                )
            )

        assert len(results) == 1
        (_, annotation), _ = results[0]
        return annotation

    def test_annotation_tuple_has_10_elements(self):
        """The annotation tuple must have exactly 10 elements."""
        annotation = self._get_annotation()
        assert annotation is not None, "annotation should not be None for a valid hit"
        assert len(annotation) == 10, (
            f"Expected 10-element annotation tuple, got {len(annotation)}: {annotation}"
        )

    def test_annotation_tuple_element_5_is_3tuple(self):
        """element[5] (match_nogs_descriptions) must be a 3-tuple: (og_name, cat, desc)."""
        annotation = self._get_annotation()
        assert annotation is not None
        nogs_desc = annotation[5]
        assert isinstance(nogs_desc, tuple), (
            f"annotation[5] must be a tuple, got {type(nogs_desc).__name__}: {nogs_desc!r}"
        )
        assert len(nogs_desc) == 3, (
            f"annotation[5] must have 3 elements (og_name, cat, desc), got {len(nogs_desc)}: {nogs_desc}"
        )

    def test_annotation_element_0_is_query_name(self):
        """element[0] is the query_name string."""
        annotation = self._get_annotation()
        assert annotation is not None
        assert annotation[0] == "query_protein"

    def test_annotation_element_8_has_all_ortholog_type_keys(self):
        """element[8] (all_orthologies) must contain all 5 ortholog type keys."""
        annotation = self._get_annotation()
        assert annotation is not None
        all_orth = annotation[8]
        required = {"one2one", "one2many", "many2one", "many2many", "all"}
        assert required == set(all_orth.keys()), (
            f"annotation[8] missing keys: got {set(all_orth.keys())}"
        )

    def test_annotation_element_9_is_list(self):
        """element[9] (annot_orthologs) must be a list."""
        annotation = self._get_annotation()
        assert annotation is not None
        assert isinstance(annotation[9], list), (
            f"annotation[9] must be a list, got {type(annotation[9]).__name__}"
        )
