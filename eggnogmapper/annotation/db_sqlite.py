##

import sqlite3

from ..common import get_eggnogdb_file, existing_file
from ..utils import colorify

db = None

def get_eggnog_db(usemem = False):
    global db
    if db is None or db.conn is None:
        db = AnnotDB(usemem)
    return db

def get_fresh_eggnog_db(usemem = False):
    return AnnotDB(usemem)

class AnnotDB(object):
    conn = None
    _int_mode = False  # True for v7+ integer-encoded databases

    def __init__(self, usemem = False):
        eggnogdb_file = get_eggnogdb_file()
        existing_file(eggnogdb_file)

        if usemem == True:
            print("Loading source DB...")
            print(colorify("Warning: this can take a few minutes and load up to 45GB to RAM. "
                           "Using --dbmem is recommended to annotate a large number of sequences.", "red"))
            source = sqlite3.connect(eggnogdb_file)
            self.conn = sqlite3.connect(':memory:')
            source.backup(self.conn)
            source.close()
        else:
            existing_file(eggnogdb_file)
            self.conn = sqlite3.connect(eggnogdb_file)

        curs = self.conn.cursor()
        curs.execute("PRAGMA synchronous=OFF;")
        curs.execute("PRAGMA journal_mode=OFF;")
        curs.execute("PRAGMA cache_size=2000;")

        # Detect integer-encoded databases (v7+: has event_index table)
        curs.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='event_index'")
        self._int_mode = curs.fetchone() is not None

        # Detect unified schema (uses ogs/sp_events) vs legacy schema (uses og/event)
        curs.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='ogs'")
        self._unified_schema = curs.fetchone() is not None

        if self._int_mode:
            from .codec import decode_intlist
            import array
            self._decode_intlist = decode_intlist
            # Load compact taxid array for fast species grouping (~237 MB for 59M proteins)
            curs.execute("SELECT MAX(id) FROM protein_names")
            max_id = curs.fetchone()[0] or 0
            self._taxids = array.array('i', [0]) * (max_id + 1)
            curs.execute("SELECT id, taxid FROM protein_names")
            for pid, taxid in curs:
                self._taxids[pid] = taxid
            print(f"  {max_id + 1} protein taxids cached "
                  f"({len(self._taxids) * 4 // 1024 // 1024} MB)", flush=True)
        else:
            self._taxids = None

        curs.close()

        return

    def close(self):
        if self.conn is not None:
            self.conn.close()
            self.conn = None
        return


    def get_db_version(self):
        cmd = 'select * from version LIMIT 1'
        curs = self.conn.cursor()
        curs.execute(cmd)
        return curs.fetchone()[0]


    def get_member_ogs(self, name):
        curs = self.conn.cursor()
        if self._int_mode:
            curs.execute('SELECT ogs FROM prots WHERE id = ?;', (int(name),))
        else:
            curs.execute('SELECT ogs FROM prots WHERE name = ?;', (name,))
        return curs.fetchone()


    def get_ogs_description(self, og, level):
        curs = self.conn.cursor()
        if self._unified_schema:
            # New unified schema: `ogs` table with name, level, nm, description, cog_cat
            curs.execute('SELECT name, nm, description, cog_cat FROM ogs WHERE name = ? AND level = ?', (og, level))
        else:
            # Legacy schema: `og` table with og, level, nm, description, COG_categories
            curs.execute('SELECT og, nm, description, COG_categories FROM og WHERE og = ? AND level = ?', (og, level))
        return curs.fetchall()


    def get_annotations(self, seq_names):
        """Get annotations for proteins. seq_names is comma-separated.
        In int_mode, values are integer protein IDs — batched into one IN query.
        In legacy mode, values are protein name strings — queried one by one.
        """
        curs = self.conn.cursor()
        _cols = ("pname, gos, kegg_ec, kegg_ko, kegg_pathway, kegg_module, kegg_reaction, "
                 "kegg_rclass, kegg_brite, kegg_tc, kegg_cazy, bigg_reaction, pfam")

        if self._int_mode:
            ids = [int(s.replace('"', '').strip()) for s in seq_names.split(",")
                   if s.replace('"', '').strip()]
            if not ids:
                return
            # Batch query: one round-trip for all orthologs
            BATCH = 900
            for i in range(0, len(ids), BATCH):
                chunk = ids[i:i + BATCH]
                ph = ",".join(["?"] * len(chunk))
                curs.execute(f"SELECT {_cols} FROM prots WHERE id IN ({ph})", chunk)
                yield from curs.fetchall()
        else:
            for seq in seq_names.split(","):
                seq = seq.replace('"', '').strip()
                if not seq:
                    continue
                curs.execute(f"SELECT {_cols} FROM prots WHERE name = ?", (seq,))
                prot_data = curs.fetchone()
                if prot_data is not None:
                    yield prot_data

        return

    def get_pfam_annotations(self, seq_names):
        """Get Pfam annotations. seq_names is comma-separated (int IDs or names)."""
        names = [s.replace('"', '').strip() for s in seq_names.split(",")]
        placeholders = ','.join(['?' for _ in names])
        if self._int_mode:
            params = [int(n) for n in names if n]
            placeholders = ','.join(['?' for _ in params])
            cmd = f"SELECT pfam FROM prots WHERE id IN ({placeholders})"
        else:
            params = names
            cmd = f"SELECT pfam FROM prots WHERE name IN ({placeholders})"

        curs = self.conn.cursor()
        curs.execute(cmd, params)
        return curs.fetchall()

    def get_taxid(self, protein_id):
        """Get taxid for a protein ID. Uses in-memory array (O(1))."""
        if self._taxids is not None and protein_id < len(self._taxids):
            return self._taxids[protein_id]
        return 0

    def decode_protein_ids(self, int_ids):
        """Batch decode integer IDs to protein names. Returns {id: name}.
        Only needed for output display — NOT for ortholog resolution.
        """
        if not int_ids:
            return {}
        id_to_name = {}
        curs = self.conn.cursor()
        BATCH = 900
        ids = list(int_ids)
        for i in range(0, len(ids), BATCH):
            chunk = ids[i:i + BATCH]
            ph = ",".join(["?"] * len(chunk))
            curs.execute(f"SELECT id, name FROM protein_names WHERE id IN ({ph})", chunk)
            id_to_name.update(curs.fetchall())
        return id_to_name

    def get_protein_name(self, protein_id):
        """Decode a single integer protein ID to its name."""
        curs = self.conn.cursor()
        curs.execute("SELECT name FROM protein_names WHERE id = ?", (int(protein_id),))
        row = curs.fetchone()
        return row[0] if row else str(protein_id)

    def get_protein_id(self, protein_name):
        """Look up integer protein ID from name. Returns None if not found."""
        curs = self.conn.cursor()
        curs.execute("SELECT id FROM protein_names WHERE name = ?", (protein_name,))
        row = curs.fetchone()
        return row[0] if row else None

    def bulk_get_protein_ids(self, names):
        """Look up integer protein IDs for a batch of protein names.

        Args:
            names: iterable of protein name strings

        Returns:
            dict mapping name -> integer ID (missing names not included)
        """
        if not names:
            return {}
        name_to_id = {}
        curs = self.conn.cursor()
        BATCH = 900
        name_list = list(names)
        for i in range(0, len(name_list), BATCH):
            chunk = name_list[i:i + BATCH]
            ph = ",".join(["?"] * len(chunk))
            curs.execute(f"SELECT name, id FROM protein_names WHERE name IN ({ph})", chunk)
            name_to_id.update(curs.fetchall())
        return name_to_id

    def get_member_events(self, member, target_levels):
        """Get orthology events for a protein.

        In int_mode: member is an integer protein ID, queries event_index table.
        side1/side2 are delta-varint BLOBs, decoded to Python int lists.
        Yields (level, side1_ids, side2_ids) where sides are list[int].

        In legacy mode: member is a protein name string, queries prots.orthoindex.
        Yields (level, side1_csv, side2_csv) where sides are "TAXID.NAME" CSV strings.
        """
        curs = self.conn.cursor()

        if self._int_mode:
            curs.execute('SELECT events FROM event_index WHERE protein_id = ?',
                         (int(member),))
            row = curs.fetchone()
            if row is None or row[0] is None:
                return
            idx_list = self._decode_intlist(row[0])
        else:
            curs.execute('SELECT orthoindex FROM prots WHERE name = ?',
                         (member.strip(),))
            event_indexes = curs.fetchone()
            if event_indexes is None or len(event_indexes) == 0:
                return
            if event_indexes[0] is None:
                return
            idx_list = [int(x) for x in str(event_indexes[0]).split(",")]

        level_list = list(target_levels)
        if idx_list and level_list:
            idx_ph = ','.join(['?' for _ in idx_list])
            lvl_ph = ','.join(['?' for _ in level_list])
            if self._unified_schema:
                # New unified schema: `sp_events` table with og_lca column
                query = f'SELECT og_lca, side1, side2 FROM sp_events WHERE i IN ({idx_ph}) AND og_lca IN ({lvl_ph})'
            else:
                # Legacy schema: `event` table with level column
                query = f'SELECT level, side1, side2 FROM event WHERE i IN ({idx_ph}) AND level IN ({lvl_ph})'
            params = idx_list + level_list
            curs.execute(query, params)
            if self._int_mode:
                for level, _side1, _side2 in curs.fetchall():
                    yield level, self._decode_intlist(_side1), self._decode_intlist(_side2)
            else:
                for level, _side1, _side2 in curs.fetchall():
                    yield level, _side1, _side2

        return



    # --- Bulk query methods for batch pre-fetch ---

    BATCH = 900  # SQLite variable limit safe threshold

    def bulk_get_ogs(self, seed_ids):
        """Fetch OGs for multiple seeds. Returns {id: ogs_string}."""
        result = {}
        curs = self.conn.cursor()
        for i in range(0, len(seed_ids), self.BATCH):
            chunk = seed_ids[i:i + self.BATCH]
            ph = ",".join(["?"] * len(chunk))
            curs.execute(f"SELECT id, ogs FROM prots WHERE id IN ({ph})", chunk)
            result.update(curs.fetchall())
        return result

    def bulk_get_event_indices(self, seed_ids):
        """Fetch event indices for multiple seeds. Returns {protein_id: [event_ids]}."""
        result = {}
        curs = self.conn.cursor()
        for i in range(0, len(seed_ids), self.BATCH):
            chunk = seed_ids[i:i + self.BATCH]
            ph = ",".join(["?"] * len(chunk))
            curs.execute(f"SELECT protein_id, events FROM event_index WHERE protein_id IN ({ph})", chunk)
            for pid, blob in curs.fetchall():
                result[pid] = self._decode_intlist(blob)
        return result

    def bulk_get_events(self, event_ids):
        """Fetch events by ID.

        Returns {i: (level, side1_ids, side2_ids)} where side1/side2 are
        sets of protein integer IDs (converted from delta-varint BLOBs).

        Uses sets for fast membership testing in ortholog collection.
        """
        result = {}
        curs = self.conn.cursor()
        ids = list(event_ids)

        for i in range(0, len(ids), self.BATCH):
            chunk = ids[i:i + self.BATCH]
            ph = ",".join(["?"] * len(chunk))
            if self._unified_schema:
                curs.execute(f"SELECT i, og_lca, side1, side2 FROM sp_events WHERE i IN ({ph})", chunk)
            else:
                curs.execute(f"SELECT i, level, side1, side2 FROM event WHERE i IN ({ph})", chunk)

            for eid, level, s1, s2 in curs.fetchall():
                side1_ids = set(self._decode_intlist(s1))
                side2_ids = set(self._decode_intlist(s2))
                result[eid] = (level, side1_ids, side2_ids)
        return result

    def bulk_get_annotations(self, protein_ids):
        """Fetch annotations for multiple proteins. Returns {id: (pname, gos, ...)}."""
        _cols = ("pname, gos, kegg_ec, kegg_ko, kegg_pathway, kegg_module, kegg_reaction, "
                 "kegg_rclass, kegg_brite, kegg_tc, kegg_cazy, bigg_reaction, pfam")
        result = {}
        curs = self.conn.cursor()
        ids = list(protein_ids)
        for i in range(0, len(ids), self.BATCH):
            chunk = ids[i:i + self.BATCH]
            ph = ",".join(["?"] * len(chunk))
            curs.execute(f"SELECT id, {_cols} FROM prots WHERE id IN ({ph})", chunk)
            for row in curs.fetchall():
                result[row[0]] = row[1:]  # id → (pname, gos, ...)
        return result

    def bulk_get_og_descriptions(self, og_level_pairs):
        """Fetch OG descriptions. Returns {(og, level): (nm, cat, desc)}."""
        result = {}
        curs = self.conn.cursor()
        pairs = list(og_level_pairs)
        for i in range(0, len(pairs), self.BATCH // 2):
            chunk = pairs[i:i + self.BATCH // 2]
            if self._unified_schema:
                # New unified schema: `ogs` table with name, level, nm, description, cog_cat
                conditions = " OR ".join(["(name = ? AND level = ?)"] * len(chunk))
                params = [v for pair in chunk for v in pair]
                curs.execute(f"SELECT name, level, nm, description, cog_cat FROM ogs WHERE {conditions}", params)
            else:
                # Legacy schema: `og` table with og, level, nm, description, COG_categories
                conditions = " OR ".join(["(og = ? AND level = ?)"] * len(chunk))
                params = [v for pair in chunk for v in pair]
                curs.execute(f"SELECT og, level, nm, description, COG_categories FROM og WHERE {conditions}", params)
            for og, level, nm, desc, cat in curs.fetchall():
                result[(og, level)] = (nm, cat, desc)
        return result


## END
