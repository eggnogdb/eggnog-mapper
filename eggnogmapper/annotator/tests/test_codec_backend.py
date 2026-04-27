"""Phase 8 / v3.2 — backend selection tests for the delta-varint codec.

The public `eggnog_annotator.codec` is a shim that prefers the Cython-
compiled `_codec` extension and falls back to pure-Python `codec_py`
when the extension wasn't built (no Python.h in the install env). These
tests assert that:

1. Whichever backend is active, basic + edge-case round-trip works.
2. The pure-Python fallback is *always* available (so the shim's import
   chain is correct).
3. When both backends are present (the typical install on a real host),
   they produce byte-identical output for the same inputs.
"""

from __future__ import annotations

import importlib

import pytest

from eggnog_annotator import codec as _codec_shim
from eggnog_annotator import codec_py


def _backends():
    """Return a list of (name, module) pairs for whichever backends load."""
    backends = [("python", codec_py)]
    try:
        compiled = importlib.import_module("eggnog_annotator._codec")
    except ImportError:
        compiled = None
    if compiled is not None:
        backends.append(("cython", compiled))
    return backends


# ---------------------------------------------------------------------------
# Shim resolves to a working backend
# ---------------------------------------------------------------------------

def test_shim_exposes_a_backend():
    assert _codec_shim._BACKEND in {"python", "cython"}


def test_shim_round_trip_basic():
    blob = _codec_shim.encode_intlist([100, 50, 75, 25])
    assert _codec_shim.decode_intlist(blob) == [25, 50, 75, 100]


def test_pure_python_fallback_always_importable():
    """`codec_py` is the safety net — must always import and work."""
    assert codec_py.decode_intlist(codec_py.encode_intlist([1, 2, 3])) == [1, 2, 3]


# ---------------------------------------------------------------------------
# Per-backend round-trip
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name,mod", _backends())
@pytest.mark.parametrize("ids", [
    [],
    [0],
    [42],
    [1, 2, 3, 4, 5],
    [5, 1, 4, 2, 3],            # unsorted input → encoder sorts
    [0, 1, 0x7F, 0x80, 0x100, 0xFFFFFF],  # varint boundary stress
    [63913227, 63913270, 63913315],
    list(range(0, 1000, 7)),
])
def test_backend_round_trip(name, mod, ids):
    decoded = mod.decode_intlist(mod.encode_intlist(ids))
    assert decoded == sorted(ids), f"{name} backend failed on {ids[:5]}..."


@pytest.mark.parametrize("name,mod", _backends())
def test_backend_handles_empty_and_none(name, mod):
    assert mod.encode_intlist([]) == b""
    assert mod.decode_intlist(b"") == []
    assert mod.decode_intlist(None) == []


# ---------------------------------------------------------------------------
# Cross-backend byte-identity (only runs when both are available)
# ---------------------------------------------------------------------------

def test_backends_produce_identical_blobs():
    """If the Cython backend is built, it MUST produce the same bytes as
    the pure-Python codec for any input. A mismatch here is a corruption
    risk — the BLOBs in eggnog.db would be silently incompatible across
    install environments."""
    try:
        compiled = importlib.import_module("eggnog_annotator._codec")
    except ImportError:
        pytest.skip("Cython backend not built in this install — skipping")

    inputs = [
        [],
        [0],
        [42],
        [1, 2, 3, 4, 5],
        [0x7F, 0x80, 0x100, 0xFFFFFF, 0x7FFFFFFF],
        [63913227, 63913270, 63913315],
        list(range(0, 10_000, 13)),
    ]
    for ids in inputs:
        py_blob = codec_py.encode_intlist(ids)
        cy_blob = compiled.encode_intlist(ids)
        assert py_blob == cy_blob, (
            f"backends disagree on input {ids[:5]}...: "
            f"py={py_blob.hex()[:40]}, cy={cy_blob.hex()[:40]}"
        )
        # And both must decode back to the same sorted list
        assert codec_py.decode_intlist(cy_blob) == sorted(ids)
        assert compiled.decode_intlist(py_blob) == sorted(ids)
