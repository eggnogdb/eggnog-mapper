"""Delta-varint codec for sorted integer lists.

Encodes a list of non-negative integers as a compact binary BLOB:
  1. Sort the list
  2. Delta-encode (store differences between consecutive values)
  3. Varint-encode each delta (1-5 bytes per value)

Typical compression: 4-5x smaller than comma-separated TEXT.
Example: "63913227,63913270,63913315" (32 bytes) -> b'...' (7 bytes)

Used for sp_events.side1, sp_events.side2 (protein ID lists) and
event_index.events (event ID lists) in eggnog.db.
"""

# encode_intlist is used by eggnog-builder to write BLOBs; not called from this package.
# decode_intlist is the hot path: called for every event in annotate_batch().


def encode_intlist(ids):
    """Encode a list of non-negative integers as a sorted delta-varint BLOB.

    Returns bytes. Empty list -> empty bytes b''.
    """
    if not ids:
        return b""
    ids = sorted(ids)
    out = bytearray()
    prev = 0
    for v in ids:
        delta = v - prev
        prev = v
        # Varint: 7 data bits per byte, high bit = continuation
        while delta >= 0x80:
            out.append((delta & 0x7F) | 0x80)
            delta >>= 7
        out.append(delta)
    return bytes(out)


def decode_intlist(data):
    """Decode a delta-varint BLOB back to a sorted list of integers.

    Accepts bytes or memoryview. Returns list[int]. Empty/None -> [].
    """
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
