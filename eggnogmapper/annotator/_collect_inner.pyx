# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
# distutils: language = c++
"""Cython implementation of the `_collect_orthologs` hot inner loop.

The pure-Python `_collect_orthologs` in `annotate.py` runs a tight
inner loop over every (seed × event-with-seed × ortholog) triple, building
a `dict[oid → (sort_key_tuple, payload_dict)]` to track the highest-
priority event each ortholog participated in. On a 1000-seed plant batch
that loop fires ~6.6 M times and accounts for ~57 % of the engine's
total wall time — pure CPython interpreter overhead on small dict/tuple
ops.

This module replaces the inner loop with a C++ `unordered_map` keyed by
the int64 protein id. The map value is a 32-bit struct holding:

  * a 32-bit packed cascade priority key, and
  * a 32-bit index into a Python list of per-event payload dicts.

The cascade priority key is itself packed from `(in_lineage_int << 24) |
((DEPTH_OFFSET - depth) << 8) | type_tier` so that comparing two keys
with `<` matches the original Python tuple `(in_lineage_int, -depth,
type_tier)` ordering — smaller = better donor. 24 bits gives room for
~16 M depth values (real taxonomic depths are < 50) and the type tier
fits in 8 bits trivially.

Payload identity is preserved via the Python list — the cascade walk
in `_summarize_annotations` reads `meta["type"|"depth"|"type_tier"|
"in_seed_lineage"]` and gets back the exact same dict objects it would
under the pure-Python path.
"""

from libc.stdint cimport int64_t, int32_t, uint32_t, uint64_t
from libcpp.unordered_map cimport unordered_map
from libcpp.utility cimport pair
from cpython.unicode cimport PyUnicode_AsUTF8AndSize, PyUnicode_FromStringAndSize

# Cascade priority key packing layout (32 bits total — fits in int32):
#   bit 31           in_lineage_inv  (0 if in seed lineage, else 1)
#   bits 8..30       depth_inv       (DEPTH_OFFSET - depth)
#   bits 0..7        type_tier
cdef uint32_t _DEPTH_OFFSET = 1 << 23   # ~16 M, far above any real taxonomic depth
cdef uint32_t _PACK_DEPTH_SHIFT = 8
cdef uint32_t _PACK_LINEAGE_SHIFT = 31


cpdef uint32_t pack_sort_key(bint in_lineage, int depth, int type_tier):
    """Pack the cascade priority key tuple into a single uint32.

    Compares the same as the original `(in_lineage_int, -depth, type_tier)`
    Python tuple under `<`: smaller packed value = better donor. Uses
    `uint32_t` so the in-lineage bit at position 31 doesn't trigger
    signed-overflow weirdness — the comparison is done as unsigned.
    """
    cdef uint32_t lineage_part = 0 if in_lineage else 1
    cdef uint32_t depth_part = (_DEPTH_OFFSET - <uint32_t>depth)
    cdef uint32_t tier_part = <uint32_t>type_tier
    return (lineage_part << _PACK_LINEAGE_SHIFT) | (depth_part << _PACK_DEPTH_SHIFT) | tier_part


# C++ struct value type for the unordered_map entry. `sort_key` is the
# 32-bit packed cascade priority (unsigned, so the comparison is unsigned);
# `event_idx` is the index into the Python `_payloads` list. Total = 8
# bytes per entry.
cdef extern from *:
    """
    struct OrthEntry {
        unsigned int sort_key;
        unsigned int event_idx;
    };
    """
    cppclass OrthEntry:
        uint32_t sort_key
        uint32_t event_idx


cdef class OrthologCollector:
    """Per-seed accumulator for the inner ortholog-meta loop.

    Replaces the pure-Python `ortholog_meta_raw` dict with a C++
    unordered_map. The map value is a single int64 packing
    `(packed_sort_key << 32) | event_idx`, where event_idx points into
    `_payloads`, a Python list of the per-event payload dicts. Payload
    identity is preserved across the dict/list round-trip so downstream
    consumers see the exact same dict objects they would have under
    the pure-Python path.

    Lifecycle:
      - `__cinit__(seed_id)`: zero state.
      - `add_event(other_side, packed_sort_key, payload)`: called once
        per event the seed appears in. Updates the map for every oid in
        other_side except `seed_id` itself, keeping the smallest
        packed_sort_key seen for each oid.
      - `candidates_set()`: snapshot of all collected oids (for the
        species/taxid filter pass).
      - `export_meta(filtered_set)`: return `{oid: payload_dict}` for
        oids whose key survived filtering.
    """

    cdef unordered_map[int64_t, OrthEntry] _meta
    cdef list _payloads
    cdef list _insertion_order  # oids in first-seen order, for deterministic export
    cdef int64_t _seed_id

    def __cinit__(self, int64_t seed_id):
        self._seed_id = seed_id
        self._payloads = []
        self._insertion_order = []

    cpdef add_event(self, other_side, uint32_t packed_sort_key, payload):
        """Register one event's `other_side` proteins with the given key.

        `other_side` may be a set, list, or any iterable of int. Items
        equal to `seed_id` are skipped (the seed never enters its own
        ortholog list). For every other oid: insert if absent, otherwise
        keep whichever packed_sort_key is smaller (better donor).
        """
        cdef:
            int64_t oid
            uint32_t event_idx = <uint32_t>len(self._payloads)
            OrthEntry entry
            unordered_map[int64_t, OrthEntry].iterator it
            unordered_map[int64_t, OrthEntry].iterator end_it

        entry.sort_key = packed_sort_key
        entry.event_idx = event_idx

        self._payloads.append(payload)
        end_it = self._meta.end()
        for oid_obj in other_side:
            oid = <int64_t>oid_obj
            if oid == self._seed_id:
                continue
            it = self._meta.find(oid)
            if it == end_it:
                self._meta[oid] = entry
                # Record the first-seen position so `export_meta` can
                # return entries in insertion order. The pure-Python
                # path's downstream tie-breakers (e.g. `Counter.most_common`
                # in `_aggregate_field` for `pname`) depend on this order,
                # so preserving it keeps output byte-identical.
                self._insertion_order.append(<long>oid)
            else:
                # Strict `<`: matches the pure-Python `if sort_key < prev[0]`
                # — first event wins on ties.
                if packed_sort_key < self._meta[oid].sort_key:
                    self._meta[oid] = entry

    def candidates_set(self):
        """Return a Python set of all oids collected so far."""
        cdef set out = set()
        cdef pair[int64_t, OrthEntry] item
        for item in self._meta:
            out.add(<long>item.first)
        return out

    def export_meta(self, filtered):
        """Return `{oid: payload_dict}` for oids in `filtered`, in
        first-seen insertion order.

        Iteration order matches the pure-Python path: each oid is
        emitted in the order it was first inserted into the map, with
        oids NOT in `filtered` skipped. Downstream tie-breakers in
        `_aggregate_field` depend on this order, so preserving it
        keeps output byte-identical to the pre-Cython baseline.

        Cost is O(|all_candidates|) rather than O(|filtered|), but the
        constant factor is tiny: a Python list iteration plus a set
        membership test per item.
        """
        cdef:
            dict out = {}
            int64_t oid
            uint32_t event_idx

        for oid_obj in self._insertion_order:
            if oid_obj not in filtered:
                continue
            oid = <int64_t>oid_obj
            event_idx = self._meta[oid].event_idx
            out[<long>oid] = self._payloads[<Py_ssize_t>event_idx]
        return out

    def size(self):
        """Number of unique oids tracked. For tests / introspection."""
        return self._meta.size()


# ---------------------------------------------------------------------------
# Comma-separated string parsers
# ---------------------------------------------------------------------------
#
# The pure-Python `_pre_parse_batch` and `_parse_gos` in `annotate.py` walk
# every annotation field of every fetched ortholog (~600 k orthologs ×
# 13 fields on a 1000-seed plant batch), splitting comma-strings and
# stripping per-token whitespace. line-profiler shows this is dominated
# by `str.split(",")` loop control + per-element `str.strip()` — pure
# Python interpreter overhead on data that is overwhelmingly ASCII.
#
# These two helpers replace those inner loops with a single C-level char
# scan over the raw UTF-8 bytes of the input. They allocate exactly one
# Python `str` per emitted token via `PyUnicode_FromStringAndSize`. No
# intermediate `list` of substrings, no per-token method dispatch.
#
# Whitespace semantics match Python's `str.strip()` for the ASCII subset
# `[ \t\n\r\v\f]`. Real annotation strings are 100 % ASCII at the eggNOG
# scale, so byte-positions correspond to character-positions and slicing
# from the UTF-8 buffer is safe.


cdef inline bint _is_ascii_space(char c) noexcept nogil:
    """Return 1 iff `c` is one of the ASCII whitespace bytes that
    `str.strip()` removes (`[ \\t\\n\\r\\v\\f]`)."""
    return (c == b' ' or c == b'\t' or c == b'\n' or c == b'\r'
            or c == b'\v' or c == b'\f')


def parse_field_fast(s):
    """Comma-split, strip, and dedup-preserve-order.

    Behaviour matches the inner block in `_pre_parse_batch` for
    non-`pname`/non-`gos` fields:

        seen = []; seen_set = set()
        for v in str(raw).split(","):
            v = v.strip()
            if v and v not in seen_set:
                seen.append(v); seen_set.add(v)
        return tuple(seen)

    Returns an empty tuple for an empty / None input.
    """
    if s is None:
        return ()
    if not isinstance(s, str):
        s = str(s)

    cdef:
        const char* buf
        Py_ssize_t n
        Py_ssize_t i = 0
        Py_ssize_t start, end

    buf = PyUnicode_AsUTF8AndSize(s, &n)
    if buf == NULL or n == 0:
        return ()

    out = []
    seen = set()
    while i < n:
        # Skip leading whitespace
        while i < n and _is_ascii_space(buf[i]):
            i += 1
        start = i
        # Find next comma or end
        while i < n and buf[i] != b',':
            i += 1
        end = i
        # Trim trailing whitespace
        while end > start and _is_ascii_space(buf[end - 1]):
            end -= 1
        if end > start:
            term = PyUnicode_FromStringAndSize(buf + start, end - start)
            if term not in seen:
                out.append(term)
                seen.add(term)
        # Skip the delimiter (if any)
        if i < n:
            i += 1
    return tuple(out)


def parse_gos_fast(s):
    """Comma-split, filter to entries starting with ``GO:``, strip the
    optional ``|<evidence>`` suffix.

    Behaviour matches `_parse_gos` in `annotate.py`:

        gos = []
        for term in go_string.split(","):
            term = term.strip()
            if term.startswith("GO:"):
                gos.append(term.split("|")[0])
        return gos

    Returns an empty tuple for an empty / None input. Note: this version
    does NOT dedup — neither does the original — because the cascade
    consumes the result through `Counter(...)` which dedups implicitly.
    """
    if s is None:
        return ()
    if not isinstance(s, str):
        s = str(s)

    cdef:
        const char* buf
        Py_ssize_t n
        Py_ssize_t i = 0
        Py_ssize_t start, end, j, pipe_pos

    buf = PyUnicode_AsUTF8AndSize(s, &n)
    if buf == NULL or n == 0:
        return ()

    out = []
    while i < n:
        # Skip leading whitespace
        while i < n and _is_ascii_space(buf[i]):
            i += 1
        start = i
        # Find next comma or end
        while i < n and buf[i] != b',':
            i += 1
        end = i
        # Trim trailing whitespace
        while end > start and _is_ascii_space(buf[end - 1]):
            end -= 1
        # Need at least 3 chars for "GO:" prefix check
        if (end - start >= 3
                and buf[start] == b'G'
                and buf[start + 1] == b'O'
                and buf[start + 2] == b':'):
            # Find optional pipe to truncate the evidence code
            pipe_pos = end
            j = start + 3
            while j < end:
                if buf[j] == b'|':
                    pipe_pos = j
                    break
                j += 1
            term = PyUnicode_FromStringAndSize(buf + start, pipe_pos - start)
            out.append(term)
        # Skip the delimiter
        if i < n:
            i += 1
    return tuple(out)
