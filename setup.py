#!/usr/bin/env python3
"""Build hook for eggnog-mapper.

The package metadata lives in setup.cfg (PEP 621-style). This file exists
to wire in the Cython extensions that ship inside eggnogmapper.annotator
(the codec and the ortholog-collector inner loop). If Cython is missing
or the C compile fails, the install will skip the extensions and the
runtime shims fall back to byte-identical pure-Python implementations
(eggnogmapper/annotator/codec.py picks the backend at import time).
"""
from setuptools import setup

try:
    from Cython.Build import cythonize
    _ext_modules = cythonize(
        [
            "eggnogmapper/annotator/_codec.pyx",
            "eggnogmapper/annotator/_collect_inner.pyx",
        ],
        language_level=3,
    )
except ImportError:  # pragma: no cover — Cython unavailable at build time
    _ext_modules = []

if __name__ == "__main__":
    setup(ext_modules=_ext_modules)
