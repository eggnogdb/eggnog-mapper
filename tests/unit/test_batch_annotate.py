"""Unit tests for eggnogmapper.annotation.batch_annotate.

Tests focus on the annotate_batch() adapter function: hit parsing,
target_orthologs filtering, warning emission, and output tuple structure.
All annotation logic is mocked — only the adapter layer is exercised here.

Annotation tuple layout (14 elements, v3 refactor):
    [0]  query_name
    [1]  seed_id (str)
    [2]  evalue
    [3]  score
    [4]  annotations (dict)
    [5]  (og_name, og_cat, og_desc) — 3-tuple
    [6]  max_annot_lvl
    [7]  match_nog_names (list)
    [8]  all_orthologies (dict with 5 type keys)
    [9]  annot_orthologs (list filtered by target_orthologs)
    [10] annotations_confidence (dict)
    [11] tax_ceiling (str)
    [12] farthest_donor_taxid (str)
    [13] farthest_donor_lineage (str)
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
        "annotations_confidence": {"GOs": "high", "KEGG_ko": "high"},
        "tax_ceiling": "Metazoa",
        "farthest_donor_taxid": "9606",
        "farthest_donor_lineage": "Eukaryota;Metazoa;Chordata;Homo",
        "og_info": {
            "name": "CLU_TEST@131567|root",
            "cog_cat": "O",
            "description": "Test OG description",
            "level": "131567",
        },
    }


def _make_ceiling_resolver(mode: str = "auto"):
    """Build a mock TaxScopeCeilingResolver with the required interface."""
    mock_resolver = MagicMock()
    mock_resolver.mode = mode
    mock_resolver.lineage_cache = MagicMock()
    return mock_resolver


def _make_batch_args(hit, target_orthologs="all"):
    """Build the minimal args list that annotate_batch's batch param expects."""
    return [(hit,)]


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

        ceiling_resolver = _make_ceiling_resolver()

        with patch.object(bat_mod, "_get_engine", return_value=mock_engine):
            results = list(
                annotate_batch(
                    batch=_make_batch_args(hit, target_orthologs="one2one"),
                    eggnog_db=mock_db,
                    annot=True,
                    target_orthologs="one2one",
                    target_taxa=None,
                    excluded_taxa=None,
                    seed_ortholog_score=60.0,
                    seed_ortholog_evalue=1e-3,
                    ceiling_resolver=ceiling_resolver,
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

        ceiling_resolver = _make_ceiling_resolver()

        with patch.object(bat_mod, "_get_engine", return_value=mock_engine):
            results = list(
                annotate_batch(
                    batch=_make_batch_args(hit, target_orthologs="all"),
                    eggnog_db=mock_db,
                    annot=True,
                    target_orthologs="all",
                    target_taxa=None,
                    excluded_taxa=None,
                    seed_ortholog_score=60.0,
                    seed_ortholog_evalue=1e-3,
                    ceiling_resolver=ceiling_resolver,
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

        ceiling_resolver = _make_ceiling_resolver()

        with caplog.at_level(logging.WARNING, logger="eggnogmapper.annotation.batch_annotate"):
            results = list(
                annotate_batch(
                    batch=_make_batch_args(hit),
                    eggnog_db=mock_db,
                    annot=True,
                    target_orthologs="all",
                    target_taxa=None,
                    excluded_taxa=None,
                    seed_ortholog_score=60.0,
                    seed_ortholog_evalue=1e-3,
                    ceiling_resolver=ceiling_resolver,
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
    """The yielded annotation tuple must conform to the 14-element contract.

    v3 refactor: tax_scope_used replaced by tax_ceiling (element [11]),
    farthest_donor_taxid (element [12]), farthest_donor_lineage (element [13]).
    """

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

        ceiling_resolver = _make_ceiling_resolver()

        with patch.object(bat_mod, "_get_engine", return_value=mock_engine):
            results = list(
                annotate_batch(
                    batch=_make_batch_args(hit, target_orthologs=target_orthologs),
                    eggnog_db=mock_db,
                    annot=True,
                    target_orthologs=target_orthologs,
                    target_taxa=None,
                    excluded_taxa=None,
                    seed_ortholog_score=60.0,
                    seed_ortholog_evalue=1e-3,
                    ceiling_resolver=ceiling_resolver,
                )
            )

        assert len(results) == 1
        (_, annotation), _ = results[0]
        return annotation

    def test_annotation_tuple_has_14_elements(self):
        """The annotation tuple must have exactly 14 elements.

        Elements 11-13 are the three new ceiling/donor columns added in
        the v3 tax_scope + cascade refactor.
        """
        annotation = self._get_annotation()
        assert annotation is not None, "annotation should not be None for a valid hit"
        assert len(annotation) == 14, (
            f"Expected 14-element annotation tuple, got {len(annotation)}: {annotation}"
        )

    def test_annotation_element_10_is_confidence_dict(self):
        """element[10] (annotations_confidence) must be a dict mapping
        output field names to one of {high, medium, low}."""
        annotation = self._get_annotation()
        assert annotation is not None
        confidence = annotation[10]
        assert isinstance(confidence, dict), (
            f"annotation[10] must be a dict, got {type(confidence).__name__}"
        )
        for field, tier in confidence.items():
            assert tier in {"high", "medium", "low"}, (
                f"annotations_confidence[{field!r}] must be high/medium/low, "
                f"got {tier!r}"
            )

    def test_annotation_element_11_is_tax_ceiling(self):
        """element[11] (tax_ceiling) must be a non-empty string.

        Replaces the old tax_scope_used field. Values are human-readable
        clade names such as 'Metazoa', 'Fungi', 'Prokaryota', or '-'.
        """
        annotation = self._get_annotation()
        assert annotation is not None
        tax_ceiling = annotation[11]
        assert isinstance(tax_ceiling, str), (
            f"annotation[11] must be a str, got {type(tax_ceiling).__name__}"
        )
        assert tax_ceiling, "annotation[11] (tax_ceiling) must not be empty"

    def test_annotation_element_12_is_farthest_donor_taxid(self):
        """element[12] (farthest_donor_taxid) must be a non-empty string."""
        annotation = self._get_annotation()
        assert annotation is not None
        farthest_taxid = annotation[12]
        assert isinstance(farthest_taxid, str), (
            f"annotation[12] must be a str, got {type(farthest_taxid).__name__}"
        )
        assert farthest_taxid, "annotation[12] (farthest_donor_taxid) must not be empty"

    def test_annotation_element_13_is_farthest_donor_lineage(self):
        """element[13] (farthest_donor_lineage) must be a non-empty string."""
        annotation = self._get_annotation()
        assert annotation is not None
        farthest_lineage = annotation[13]
        assert isinstance(farthest_lineage, str), (
            f"annotation[13] must be a str, got {type(farthest_lineage).__name__}"
        )
        assert farthest_lineage, "annotation[13] (farthest_donor_lineage) must not be empty"

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

    def test_no_tax_scope_used_in_tuple(self):
        """The old 'tax_scope_used' field must not appear as a string in
        any tuple element — it has been replaced by tax_ceiling.

        This test catches the biologically wrong case where the old field
        name is still emitted even though its content is correct format.
        """
        annotation = self._get_annotation()
        assert annotation is not None
        # None of the string elements should be the old "tax_scope_used" header value
        # Specifically element[11] should be a clade name, not "tax_scope_used"
        assert annotation[11] != "tax_scope_used", (
            "annotation[11] is the literal string 'tax_scope_used' — "
            "old field name leaked into the data"
        )
