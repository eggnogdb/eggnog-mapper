"""eggNOG v7 annotation module.

Integer-based batch annotation for v7+ databases:
- Protein IDs are integers with name mapping in protein_names table
- Event side1/side2 use delta-varint encoded integer BLOBs
- Batch queries with pre-fetched taxid array for high throughput

v5-style string-encoded databases are not supported in eggnog-mapper v3.
Use eggnog-mapper v2.x for v5 databases.
"""

from .tax_scope import (
    LineageFilter,
    resolve_name_to_taxid,
    BACTERIA,
    ARCHAEA,
    EUKARYOTA,
    METAZOA,
    FUNGI,
    VIRIDIPLANTAE,
    CELLULAR_ORGANISMS,
)
from .ceiling import (
    TaxScopeCeilingResolver,
    CEILING_NAMES,
    AUTO_NARROW_PRIORITY,
    AUTO_BROAD_PRIORITY,
    OPISTHOKONTA,
    MICROSPORIDIA,
    PROKARYOTA_SYNTHETIC,
)
from .db import EggnogDB, EggnogDBError
from .annotate import AnnotationEngine, annotate_protein
from .orthologs import OrthologFinder, DonorTracker

__all__ = [
    # Tax scope
    "LineageFilter",
    "resolve_name_to_taxid",
    "BACTERIA",
    "ARCHAEA",
    "EUKARYOTA",
    "METAZOA",
    "FUNGI",
    "VIRIDIPLANTAE",
    "CELLULAR_ORGANISMS",
    # Ceiling resolver
    "TaxScopeCeilingResolver",
    "CEILING_NAMES",
    "AUTO_NARROW_PRIORITY",
    "AUTO_BROAD_PRIORITY",
    "OPISTHOKONTA",
    "MICROSPORIDIA",
    "PROKARYOTA_SYNTHETIC",
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
