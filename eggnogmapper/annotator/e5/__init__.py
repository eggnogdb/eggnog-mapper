"""eggNOG v5 annotation module.

This module provides legacy string-based annotation logic for v5-style databases:
- Protein IDs are stored as "TAXID.PROTEINNAME" strings
- Event side1/side2 are comma-separated string lists
- Per-hit queries (no batch optimization)

Use eggnog_annotator.e7 for v7+ integer-encoded databases.
"""
