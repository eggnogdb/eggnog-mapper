"""Unified database schema definitions for eggNOG databases.

This module defines the SQLite table schemas shared between eggnog-mapper
and eggnog-website. The unified schema allows:

- eggnog-builder to create a single database with all tables
- Mapper release builds to drop website-only tables
- Both mapper and website to share annotation logic

Tables are organized into:
- CORE_TABLES: Required by both mapper and website
- WEBSITE_TABLES: Only needed by website (dropped in mapper release)
"""

# =============================================================================
# Core Tables (mapper + website)
# =============================================================================

PROTEIN_NAMES_SQL = """
CREATE TABLE IF NOT EXISTS protein_names (
    id INTEGER PRIMARY KEY,        -- integer protein ID
    name TEXT NOT NULL,            -- "TAXID.PROTEINNAME"
    taxid INTEGER NOT NULL         -- extracted taxid
);
"""

PROTEIN_NAMES_INDEX_SQL = """
CREATE UNIQUE INDEX IF NOT EXISTS idx_pn_name ON protein_names(name);
"""

SP_EVENTS_SQL = """
CREATE TABLE IF NOT EXISTS sp_events (
    i INTEGER PRIMARY KEY,         -- event ID
    name TEXT NOT NULL,            -- cluster name
    og TEXT,                       -- OG name
    og_lca TEXT,                   -- OG LCA taxid
    ev_lca TEXT,                   -- event LCA taxid
    sp_overlap REAL,               -- species Jaccard overlap between the two sides
    side1 BLOB NOT NULL,           -- delta-varint protein IDs
    side2 BLOB NOT NULL            -- delta-varint protein IDs
);
"""

EVENT_INDEX_SQL = """
CREATE TABLE IF NOT EXISTS event_index (
    protein_id INTEGER PRIMARY KEY,
    events BLOB NOT NULL           -- delta-varint event IDs
);
"""

PROTS_SQL = """
CREATE TABLE IF NOT EXISTS prots (
    id INTEGER PRIMARY KEY,
    pname VARCHAR(32),
    gos TEXT,
    kegg_ko VARCHAR(256),
    kegg_ec VARCHAR(256),
    kegg_pathway VARCHAR(256),
    kegg_module VARCHAR(256),
    kegg_reaction VARCHAR(256),
    kegg_rclass VARCHAR(256),
    kegg_brite VARCHAR(256),
    kegg_tc VARCHAR(256),
    kegg_cazy VARCHAR(256),
    kegg_cog VARCHAR(256),
    kegg_disease VARCHAR(256),
    kegg_go VARCHAR(256),
    kegg_drug VARCHAR(256),
    kegg_pubmed TEXT,
    kegg_network TEXT,
    bigg_reaction VARCHAR(32),
    pfam TEXT,
    ogs VARCHAR(256)               -- OG membership (CSV), used by mapper
);
"""

OGS_SQL = """
CREATE TABLE IF NOT EXISTS ogs (
    name TEXT PRIMARY KEY,         -- full OG name "cluster@taxid|clade"
    clust_name TEXT NOT NULL,      -- cluster name
    level TEXT NOT NULL,           -- taxonomic level (taxid)
    nm INTEGER,                    -- member count
    ns INTEGER,                    -- species count
    description TEXT,              -- aggregated from seq2bestname (mapper needs)
    cog_cat TEXT,                  -- COG categories (mapper needs)
    seqs BLOB,                     -- delta-varint protein IDs (website only)
    fprof_sum TEXT                 -- JSON functional profile (website only)
) WITHOUT ROWID;
"""

OGS_CLUST_INDEX_SQL = """
CREATE INDEX IF NOT EXISTS idx_ogs_clust ON ogs(clust_name);
"""

OGS_LEVEL_INDEX_SQL = """
CREATE INDEX IF NOT EXISTS idx_ogs_level ON ogs(level);
"""

REF_TERMS_SQL = """
CREATE TABLE IF NOT EXISTS ref_terms (
    type TEXT NOT NULL,            -- 'kegg_ko', 'go', etc.
    term TEXT NOT NULL,
    symbol TEXT,
    descr TEXT,
    PRIMARY KEY (type, term)
) WITHOUT ROWID;
"""

VERSION_SQL = """
CREATE TABLE IF NOT EXISTS version (
    version VARCHAR(16) PRIMARY KEY
) WITHOUT ROWID;
"""


# =============================================================================
# Website-only Tables
# =============================================================================

SEQINFO_SQL = """
CREATE TABLE IF NOT EXISTS seqinfo (
    id INTEGER PRIMARY KEY,
    bname TEXT,                    -- best name
    gn TEXT,                       -- gene name
    cog TEXT,
    kegg TEXT,
    pfam TEXT,
    ukb TEXT,                      -- UniProtKB ID
    ugo TEXT,                      -- UniProt GO
    cazy TEXT,
    ogs TEXT                       -- OG list as TEXT (reverse index)
);
"""

SEQINFO_UKB_INDEX_SQL = """
CREATE INDEX IF NOT EXISTS idx_si_ukb ON seqinfo(ukb);
"""

SEQUENCES_SQL = """
CREATE TABLE IF NOT EXISTS sequences (
    id INTEGER PRIMARY KEY,
    aa TEXT NOT NULL
);
"""

CLU_INFO_SQL = """
CREATE TABLE IF NOT EXISTS clu_info (
    name TEXT PRIMARY KEY,
    nm INTEGER,                    -- member count
    ns INTEGER,                    -- species count
    no INTEGER,                    -- OG count
    prots BLOB,                    -- delta-varint protein IDs
    sp TEXT                        -- species list (CSV taxids)
) WITHOUT ROWID;
"""

OG_FUNC_SQL = """
CREATE TABLE IF NOT EXISTS og_func (
    func_type TEXT NOT NULL,       -- 'kegg_ko', 'go', 'pfam_dom', etc.
    term TEXT NOT NULL,
    og_name TEXT NOT NULL,
    freq REAL NOT NULL,            -- term frequency in OG
    PRIMARY KEY (func_type, term, og_name)
) WITHOUT ROWID;
"""


# =============================================================================
# Table Groups
# =============================================================================

CORE_TABLES = [
    ("protein_names", PROTEIN_NAMES_SQL, [PROTEIN_NAMES_INDEX_SQL]),
    ("sp_events", SP_EVENTS_SQL, []),
    ("event_index", EVENT_INDEX_SQL, []),
    ("prots", PROTS_SQL, []),
    ("ogs", OGS_SQL, [OGS_CLUST_INDEX_SQL, OGS_LEVEL_INDEX_SQL]),
    ("ref_terms", REF_TERMS_SQL, []),
    ("version", VERSION_SQL, []),
]

WEBSITE_TABLES = [
    ("seqinfo", SEQINFO_SQL, [SEQINFO_UKB_INDEX_SQL]),
    ("sequences", SEQUENCES_SQL, []),
    ("clu_info", CLU_INFO_SQL, []),
    ("og_func", OG_FUNC_SQL, []),
]

ALL_TABLES = CORE_TABLES + WEBSITE_TABLES


def create_core_tables(conn):
    """Create core tables required by both mapper and website.

    Args:
        conn: sqlite3.Connection
    """
    cursor = conn.cursor()
    for name, create_sql, index_sqls in CORE_TABLES:
        cursor.execute(create_sql)
        for idx_sql in index_sqls:
            cursor.execute(idx_sql)
    conn.commit()


def create_all_tables(conn):
    """Create all tables (core + website).

    Args:
        conn: sqlite3.Connection
    """
    cursor = conn.cursor()
    for name, create_sql, index_sqls in ALL_TABLES:
        cursor.execute(create_sql)
        for idx_sql in index_sqls:
            cursor.execute(idx_sql)
    conn.commit()


def drop_website_tables(conn):
    """Drop website-only tables for mapper release build.

    Args:
        conn: sqlite3.Connection
    """
    cursor = conn.cursor()
    for name, _, _ in WEBSITE_TABLES:
        cursor.execute(f"DROP TABLE IF EXISTS {name}")
    conn.commit()


def get_table_names(conn):
    """Get list of table names in database.

    Args:
        conn: sqlite3.Connection

    Returns:
        list of table names (str)
    """
    cursor = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
    )
    return [row[0] for row in cursor.fetchall()]


def is_v7_database(db_path_or_conn):
    """Check if database uses v7 integer-encoded format.

    v7 databases have an event_index table for fast protein -> events lookup.

    Args:
        db_path_or_conn: Path to database (str) or sqlite3.Connection

    Returns:
        bool
    """
    import sqlite3
    if isinstance(db_path_or_conn, str):
        conn = sqlite3.connect(db_path_or_conn)
        tables = get_table_names(conn)
        conn.close()
    else:
        tables = get_table_names(db_path_or_conn)
    return "event_index" in tables


def has_website_tables(conn):
    """Check if database has website-only tables.

    Args:
        conn: sqlite3.Connection

    Returns:
        bool
    """
    tables = set(get_table_names(conn))
    website_table_names = {name for name, _, _ in WEBSITE_TABLES}
    return bool(tables & website_table_names)
