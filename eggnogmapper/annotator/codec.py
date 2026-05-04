"""Delta-varint codec public entrypoint.

This module is a thin shim: it imports `encode_intlist` / `decode_intlist`
from the Cython-compiled `_codec` extension when available, and falls
back to the pure-Python implementation in `codec_py` otherwise. Callers
should always `from eggnogmapper.annotator.codec import …`; the choice of
backend is transparent.

The Cython backend (Phase 8 / v3.2) is ~20–50× faster than the pure-
Python loop on integer-heavy workloads. The fallback is byte-identical
and is used when the package is installed in an environment without
`Python.h` (e.g. minimal sandboxes).

The byte format is documented in `codec_py.py`.
"""

try:
    from ._codec import encode_intlist, decode_intlist  # type: ignore  # noqa: F401
    _BACKEND = "cython"
except ImportError:  # pragma: no cover — exercised only on no-compile installs
    from .codec_py import encode_intlist, decode_intlist  # noqa: F401
    _BACKEND = "python"


__all__ = ["encode_intlist", "decode_intlist", "_BACKEND"]
