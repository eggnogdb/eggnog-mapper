"""eggnog-annotator: Unified annotation module for eggNOG ecosystem.

This package provides shared annotation logic for both eggnog-mapper
and eggnog-website, with version-specific implementations:

- eggnogmapper.annotator.e5: v5-style string-based annotation
- eggnogmapper.annotator.e7: v7-style integer-based batch annotation

Shared modules:
- codec: Delta-varint encode/decode for integer lists
- lineage: LineageCache for taxonomy lineage lookups
- schema: Unified database schema definitions
"""

__version__ = "0.1.0"

from .codec import encode_intlist, decode_intlist

__all__ = ["encode_intlist", "decode_intlist", "__version__"]
