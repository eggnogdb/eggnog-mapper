"""Shared pytest fixtures for eggnog-annotator tests."""

import os
import pytest

# ---------------------------------------------------------------------------
# Sample data paths
# ---------------------------------------------------------------------------

_SAMPLE_DIR = os.environ.get("EGGNOG_SAMPLE_DATA_DIR") or os.path.normpath(
    os.path.join(
        os.path.dirname(__file__),
        "..", "..", "..", "..", "data", "e7", "sample", "final", "mapper",
    )
)

SAMPLE_DB_PATH = os.path.join(_SAMPLE_DIR, "eggnog.db")
TAXA_DB_PATH = os.path.join(_SAMPLE_DIR, "eggnog.taxa.db")
TRAVERSE_PKL_PATH = os.path.join(_SAMPLE_DIR, "eggnog.taxa.db.traverse.pkl")

# Known seed protein IDs in sample DB (from manifest.yaml)
# "9606.ENSP00000269305" -> TP53 (human) -> id=1381
# "426368.MmarC7_1147"  -> NifH (archaea) -> id=169
TP53_PROTEIN_ID = 1381
NIFH_PROTEIN_ID = 169


@pytest.fixture(scope="session")
def sample_db_path():
    """Absolute path to sample eggnog.db."""
    if not os.path.exists(SAMPLE_DB_PATH):
        pytest.skip(f"Sample DB not found: {SAMPLE_DB_PATH}")
    return SAMPLE_DB_PATH


@pytest.fixture(scope="session")
def taxa_db_path():
    """Absolute path to sample eggnog.taxa.db."""
    if not os.path.exists(TAXA_DB_PATH):
        pytest.skip(f"Taxa DB not found: {TAXA_DB_PATH}")
    return TAXA_DB_PATH


@pytest.fixture(scope="session")
def engine(sample_db_path):
    """Shared AnnotationEngine instance on the sample DB (session-scoped)."""
    from eggnogmapper.annotator.e7 import AnnotationEngine, EggnogDB

    db = EggnogDB(sample_db_path)
    eng = AnnotationEngine(db)
    yield eng
    db.close()
