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
        # curs.execute("PRAGMA temp_store=MEMORY;")
        # curs.execute("PRAGMA locking_mode=EXCLUSIVE;")
        # curs.execute("PRAGMA mmap_size=0;")
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
        curs.execute('SELECT ogs FROM prots WHERE name = ?;', (name,))
        return curs.fetchone()


    def get_ogs_description(self, og, level):
        curs = self.conn.cursor()
        curs.execute('SELECT og.og, nm, description, COG_categories FROM og WHERE og.og = ? AND og.level = ?', (og, level))
        return curs.fetchall()


    def get_annotations(self, seq_names):
        curs = self.conn.cursor()

        for seq in seq_names.split(","):
            seq = seq.replace('"', '')
            curs.execute(("SELECT pname, gos, kegg_ec, kegg_ko, kegg_pathway, kegg_module, kegg_reaction, "
                          "kegg_rclass, kegg_brite, kegg_tc, kegg_cazy, bigg_reaction, pfam FROM prots WHERE name = ?")
                         , (seq,))

            prot_data = curs.fetchone()
            if prot_data is not None:
                yield prot_data

        return

    def get_pfam_annotations(self, seq_names):
        names = [s.replace('"', '').strip() for s in seq_names.split(",")]
        placeholders = ','.join(['?' for _ in names])
        cmd = f"SELECT pfam FROM prots WHERE name IN ({placeholders})"

        curs = self.conn.cursor()
        curs.execute(cmd, names)
        return curs.fetchall()

    def get_member_events(self, member, target_levels):
        curs = self.conn.cursor()
        curs.execute('SELECT orthoindex FROM prots WHERE name = ?', (member.strip(),))
        event_indexes = curs.fetchone()
        if event_indexes is not None and len(event_indexes) > 0:
            if event_indexes[0] is not None:
                idx_list = str(event_indexes[0]).split(",")
                level_list = list(target_levels)
                if idx_list and level_list:
                    idx_ph = ','.join(['?' for _ in idx_list])
                    lvl_ph = ','.join(['?' for _ in level_list])
                    query = f'SELECT level, side1, side2 FROM event WHERE i IN ({idx_ph}) AND level IN ({lvl_ph})'
                    params = [int(x) for x in idx_list] + level_list
                    curs.execute(query, params)
                    for level, _side1, _side2 in curs.fetchall():
                        yield level, _side1, _side2

        return



## END
