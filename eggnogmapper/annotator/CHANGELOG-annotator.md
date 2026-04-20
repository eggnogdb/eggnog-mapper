## [Unreleased]

### Fixed

- **C7** (annotate.py): Out-of-range protein IDs now emit `logger.warning` instead of being silently discarded, making data integrity issues visible in logs.
- **C8** (tax_scope.py): `get_valid_species_ids()` now returns a `frozenset`, preventing accidental cache corruption by callers that mutate the result.
- **H4** (db.py): `EggnogDB.__init__` now raises `EggnogDBError` with the database path in the message on open failure, replacing a bare `sqlite3.OperationalError` with no context.
- **H5** (tax_scope.py): Bare `except Exception` replaced with `except sqlite3.OperationalError` combined with `logger.error()`, so unexpected exceptions are no longer swallowed silently.
