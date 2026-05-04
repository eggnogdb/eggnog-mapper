"""Unit tests for the delta-varint codec (encode_intlist / decode_intlist)."""

import pytest

from eggnog_annotator import encode_intlist, decode_intlist


# ---------------------------------------------------------------------------
# Roundtrip helpers
# ---------------------------------------------------------------------------

def _roundtrip(ids):
    """Encode then decode and return the result."""
    return decode_intlist(encode_intlist(ids))


# ---------------------------------------------------------------------------
# Parametrized roundtrip tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("ids, label", [
    ([], "empty list"),
    ([0], "single zero"),
    ([1], "single one"),
    ([0, 1_000_000], "large delta"),
    ([42], "arbitrary small"),
    ([1, 2, 3, 4, 5], "consecutive small ints"),
    ([100, 200, 300, 400, 500], "sorted multiples of 100"),
    ([63913227, 63913270, 63913315, 63913400, 63913500], "realistic ortholog IDs"),
    ([0, 1, 127, 128, 255, 256, 16383, 16384], "varint boundary values"),
])
def test_roundtrip(ids, label):
    """decode(encode(x)) == sorted(x) for all cases."""
    result = _roundtrip(ids)
    assert result == sorted(ids), f"Roundtrip failed for '{label}': {ids} -> {result}"


# ---------------------------------------------------------------------------
# Specific edge-case assertions
# ---------------------------------------------------------------------------

def test_empty_encodes_to_empty_bytes():
    """Empty list encodes to b''."""
    assert encode_intlist([]) == b""


def test_empty_bytes_decodes_to_empty_list():
    """Empty bytes decodes to []."""
    assert decode_intlist(b"") == []


def test_none_decodes_to_empty_list():
    """None decodes to [] (defensive handling)."""
    assert decode_intlist(None) == []


def test_encode_sorts_input():
    """Encode sorts the input before encoding."""
    unsorted = [5, 3, 1, 4, 2]
    result = decode_intlist(encode_intlist(unsorted))
    assert result == [1, 2, 3, 4, 5]


def test_large_delta_roundtrip():
    """A large delta (e.g. 1 million) requires multi-byte varint encoding."""
    ids = [0, 1_000_000]
    assert _roundtrip(ids) == ids


def test_typical_ortholog_list_roundtrip():
    """A realistic 5+ item ortholog list roundtrips correctly."""
    ids = [100, 1050, 2300, 45678, 63913227, 63913270, 63913315]
    assert _roundtrip(ids) == sorted(ids)


def test_decode_accepts_memoryview():
    """decode_intlist accepts memoryview in addition to bytes."""
    ids = [1, 2, 3]
    encoded = encode_intlist(ids)
    result = decode_intlist(memoryview(encoded))
    assert result == ids


def test_encode_output_is_bytes():
    """encode_intlist always returns bytes."""
    assert isinstance(encode_intlist([1, 2, 3]), bytes)
    assert isinstance(encode_intlist([]), bytes)


def test_decode_output_is_list_of_ints():
    """decode_intlist always returns list[int]."""
    ids = [10, 20, 30]
    result = decode_intlist(encode_intlist(ids))
    assert isinstance(result, list)
    assert all(isinstance(x, int) for x in result)
