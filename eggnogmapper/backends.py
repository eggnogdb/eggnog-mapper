##
## Database backend registry
##
## Maps backend names to their data directory paths.
## Built-in backends are always available. Additional backends are
## auto-discovered from EGGNOG_DATA_ROOT if set.

import os
from os.path import join as pjoin

BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Built-in backends (always available)
_BUILTIN_BACKENDS = {
    "eggnog5":    pjoin(BASE_PATH, "data"),
    "e5-test":    pjoin(BASE_PATH, "tests", "fixtures"),
}

DEFAULT_BACKEND = "eggnog5"

def _discover_backends():
    """Build backend registry: built-ins + auto-discovered from EGGNOG_DATA_ROOT.

    If EGGNOG_DATA_ROOT is set and points to a data/ directory with the
    ecosystem layout (data/<tag>/<dataset>/final/mapper/eggnog.db), those
    are registered as "<tag>-<dataset>" backends.
    """
    backends = dict(_BUILTIN_BACKENDS)
    data_root = os.environ.get("EGGNOG_DATA_ROOT")
    if not data_root or not os.path.isdir(data_root):
        return backends
    try:
        for tag in os.listdir(data_root):
            tag_dir = pjoin(data_root, tag)
            if not os.path.isdir(tag_dir):
                continue
            for dataset in os.listdir(tag_dir):
                mapper_dir = pjoin(tag_dir, dataset, "final", "mapper")
                if os.path.isfile(pjoin(mapper_dir, "eggnog.db")):
                    backends[f"{tag}-{dataset}"] = mapper_dir
    except OSError:
        pass
    return backends

# Lazy singleton
_backends = None

def _get_backends():
    global _backends
    if _backends is None:
        _backends = _discover_backends()
    return _backends

def get_backend_names():
    return list(_get_backends().keys())

def resolve_backend(name):
    backends = _get_backends()
    if name not in backends:
        raise ValueError(f"Unknown database backend '{name}'. Available: {', '.join(backends.keys())}")
    return backends[name]

## END
