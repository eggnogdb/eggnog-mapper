"""eggNOG v7 annotation module.

This module provides integer-based batch annotation logic for v7+ databases:
- Protein IDs are integers with name mapping in protein_names table
- Event side1/side2 use delta-varint encoded integer BLOBs
- Batch queries with pre-fetched taxid array for high throughput

Use eggnog_annotator.e5 for legacy v5-style databases.
"""

from .tax_scope import (
    LineageFilter,
    parse_tax_scope,
    BACTERIA,
    ARCHAEA,
    EUKARYOTA,
    METAZOA,
    FUNGI,
    VIRIDIPLANTAE,
    CELLULAR_ORGANISMS,
)
from .db import EggnogDB, EggnogDBError
from .annotate import AnnotationEngine, annotate_protein
from .orthologs import OrthologFinder, DonorTracker

__all__ = [
    # Tax scope
    "LineageFilter",
    "parse_tax_scope",
    "BACTERIA",
    "ARCHAEA",
    "EUKARYOTA",
    "METAZOA",
    "FUNGI",
    "VIRIDIPLANTAE",
    "CELLULAR_ORGANISMS",
    # Database
    "EggnogDB",
    "EggnogDBError",
    # Annotation
    "AnnotationEngine",
    "annotate_protein",
    # Orthologs
    "OrthologFinder",
    "DonorTracker",
]
