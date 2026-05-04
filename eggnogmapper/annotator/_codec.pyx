# cython: boundscheck=False, wraparound=False, language_level=3
# cython: cdivision=True
"""Delta-varint codec for sorted integer lists — Cython implementation.

Same byte format and semantics as `eggnogmapper.annotator.codec_py`. Compiles to
a small `.so` that replaces the pure-Python loops with tight C; on full e7
the encoding step inside `eggnog-builder build-web --step eggnogdb` Phase A
drops from ~25 minutes to ~30 seconds.

If this extension fails to build (no `Python.h` on the host), the package
init falls back transparently to the pure-Python version. Round-trip
property tests in `eggnog-builder/tests/test_sp_events_audit.py` and
`eggnog-annotator/tests/test_codec.py` cover both paths.
"""


def encode_intlist(ids):
    """Encode a list of non-negative integers as a sorted delta-varint BLOB.

    Returns bytes. Empty list -> empty bytes b''.
    """
    if not ids:
        return b""
    cdef list sorted_ids = sorted(ids)
    cdef bytearray out = bytearray()
    cdef long long prev = 0
    cdef long long v
    cdef long long delta
    for v in sorted_ids:
        delta = v - prev
        prev = v
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
    cdef bytes buf = bytes(data)  # copy to a known type for fast indexing
    cdef list ids = []
    cdef long long prev = 0
    cdef long long delta
    cdef int shift
    cdef Py_ssize_t i = 0
    cdef Py_ssize_t n = len(buf)
    cdef unsigned char b
    while i < n:
        delta = 0
        shift = 0
        while True:
            b = buf[i]
            i += 1
            delta |= (<long long>(b & 0x7F)) << shift
            if b < 0x80:
                break
            shift += 7
        prev += delta
        ids.append(prev)
    return ids
