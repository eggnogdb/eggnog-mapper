##
## Database backend registry
##
## Maps backend names to their data directory paths.
## To add a new backend, insert an entry into BACKENDS.

import os
from os.path import join as pjoin

BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

BACKENDS = {
    "eggnog5":    pjoin(BASE_PATH, "data"),
    "e5-test":    pjoin(BASE_PATH, "tests", "fixtures"),
}

DEFAULT_BACKEND = "eggnog5"

def get_backend_names():
    return list(BACKENDS.keys())

def resolve_backend(name):
    if name not in BACKENDS:
        raise ValueError(f"Unknown database backend '{name}'. Available: {', '.join(BACKENDS.keys())}")
    return BACKENDS[name]

## END
