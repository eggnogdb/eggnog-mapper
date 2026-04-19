"""Delta-varint codec for sorted integer lists.

Imported from eggnog-annotator, which is a required dependency.
"""

from eggnog_annotator import encode_intlist, decode_intlist

__all__ = ["encode_intlist", "decode_intlist"]
