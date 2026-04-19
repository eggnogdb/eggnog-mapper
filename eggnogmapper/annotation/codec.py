"""Delta-varint codec for sorted integer lists.

Imports from eggnog-annotator if available, otherwise uses local implementation
for standalone mapper installs.
"""

try:
    from eggnog_annotator import encode_intlist, decode_intlist
except ImportError:
    # Fallback for standalone mapper installs without eggnog-annotator
    def encode_intlist(ids):
        """Encode a list of non-negative integers as a sorted delta-varint BLOB."""
        if not ids:
            return b""
        ids = sorted(ids)
        out = bytearray()
        prev = 0
        for v in ids:
            delta = v - prev
            prev = v
            while delta >= 0x80:
                out.append((delta & 0x7F) | 0x80)
                delta >>= 7
            out.append(delta)
        return bytes(out)

    def decode_intlist(data):
        """Decode a delta-varint BLOB back to a sorted list of integers."""
        if not data:
            return []
        ids = []
        prev = 0
        i = 0
        n = len(data)
        while i < n:
            delta = 0
            shift = 0
            while True:
                b = data[i]
                i += 1
                delta |= (b & 0x7F) << shift
                if b < 0x80:
                    break
                shift += 7
            prev += delta
            ids.append(prev)
        return ids

__all__ = ["encode_intlist", "decode_intlist"]
