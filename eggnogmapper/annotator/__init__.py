"""Annotation engine for eggnog-mapper (v3+).

Originally distributed as a separate `eggnog-annotator` package; folded
into eggnog-mapper at v3.0.0. Provides the cascade annotation engine
for v7+ integer-encoded eggnog.db files.

Subpackages:
- ``e7``: integer-based batch annotation for v7+ databases — the
  cascade engine, EggnogDB wrapper, LineageFilter and tax-scope helpers
- ``codec`` / ``_codec``: delta-varint encode/decode for integer-list
  BLOBs (Cython with pure-Python fallback)
- ``lineage``: LineageCache for taxonomy lineage lookups
- ``schema``: unified DB schema definitions, used by mapper at runtime
  and by eggnog-builder during DB construction

v5 string-based databases are not supported in v3 (the e5 stub
subpackage was removed in v3.0.0). Use eggnog-mapper v2.x for v5
databases.
"""

__version__ = "0.1.0"

from .codec import encode_intlist, decode_intlist

__all__ = ["encode_intlist", "decode_intlist", "__version__"]
