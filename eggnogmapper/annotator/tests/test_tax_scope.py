"""Unit tests for tax scope filtering (LineageFilter.get_valid_species_ids)."""

import pytest

from tests.conftest import TAXA_DB_PATH


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def lineage_cache(taxa_db_path):
    """LineageCache loaded from the sample taxa.db."""
    from eggnog_annotator.lineage import LineageCache
    return LineageCache(taxa_db_path=taxa_db_path)


@pytest.fixture(scope="module")
def lineage_filter(lineage_cache):
    """LineageFilter wrapping the sample lineage cache."""
    from eggnog_annotator.e7.tax_scope import LineageFilter
    return LineageFilter(lineage_cache)


# ---------------------------------------------------------------------------
# Return type tests
# ---------------------------------------------------------------------------

def test_get_valid_species_ids_returns_frozenset(lineage_filter):
    """get_valid_species_ids must return a frozenset, not a plain set."""
    # Metazoa taxid scope
    scope = {"33208"}
    result = lineage_filter.get_valid_species_ids(scope)
    assert result is not None, "Expected non-None result for Metazoa scope"
    assert isinstance(result, frozenset), (
        f"Expected frozenset, got {type(result).__name__}"
    )


def test_frozenset_immutable(lineage_filter):
    """The frozenset returned by get_valid_species_ids cannot be mutated."""
    scope = {"33208"}
    result = lineage_filter.get_valid_species_ids(scope)
    assert result is not None
    with pytest.raises(AttributeError):
        result.add("999999")  # frozenset has no .add()


def test_empty_scope_returns_none(lineage_filter):
    """Passing an empty scope should return None (no filter)."""
    result = lineage_filter.get_valid_species_ids(set())
    assert result is None


# ---------------------------------------------------------------------------
# Biological correctness: Metazoa scope
# ---------------------------------------------------------------------------

def test_metazoa_scope_contains_human(lineage_filter):
    """The Metazoa (33208) scope must include human (9606)."""
    scope = {"33208"}
    result = lineage_filter.get_valid_species_ids(scope)
    assert result is not None
    assert "9606" in result, "Human (9606) should be in Metazoa scope"


def test_metazoa_excludes_bacteria(lineage_filter):
    """The Metazoa scope must NOT include E. coli (511145), a bacterium."""
    scope = {"33208"}
    result = lineage_filter.get_valid_species_ids(scope)
    assert result is not None
    assert "511145" not in result, (
        "E. coli K-12 (511145) should NOT be in Metazoa (33208) scope"
    )


def test_metazoa_excludes_archaea(lineage_filter):
    """The Metazoa scope must NOT include Methanomicrobium (426368), an archaeon."""
    scope = {"33208"}
    result = lineage_filter.get_valid_species_ids(scope)
    assert result is not None
    assert "426368" not in result, (
        "Methanomicrobium mobile (426368) should NOT be in Metazoa scope"
    )


def test_prokaryote_scope_excludes_human(lineage_filter):
    """A Bacteria-only (2) scope must NOT include human (9606)."""
    scope = {"2"}
    result = lineage_filter.get_valid_species_ids(scope)
    if result is None:
        pytest.skip("No bacteria in sample taxa.db")
    assert "9606" not in result, "Human should NOT appear in a Bacteria-only scope"


# ---------------------------------------------------------------------------
# Caching: same call returns same object
# ---------------------------------------------------------------------------

def test_result_is_cached(lineage_filter):
    """Repeated calls with the same scope return the identical frozenset object."""
    scope = {"33208"}
    r1 = lineage_filter.get_valid_species_ids(scope)
    r2 = lineage_filter.get_valid_species_ids(scope)
    assert r1 is r2, "Expected the same frozenset object (memoized)"
