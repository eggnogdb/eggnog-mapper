##

import sqlite3

from ..common import get_eggnogdb_file, existing_file
from ..utils import colorify

db = None

def get_eggnog_db(usemem = False):
    global db
    if db is None or db.conn is None:
        print(f"Creating new AnnotDB with usemem {usemem}")
        db = AnnotDB(usemem)
    else:
        print(f"AnnotDB already exists")        
    return db

def get_fresh_eggnog_db(usemem = False):
    print(f"Creating fresh AnnotDB with usemem {usemem}")
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
        curs.execute('SELECT groups FROM eggnog WHERE name == ?;', (name,))
        return curs.fetchone()


    def get_ogs_description(self, og, level):
        curs = self.conn.cursor()
        curs.execute('SELECT og.og, nm, description, COG_categories FROM og WHERE og.og == ? AND og.level == ?', (og, level))
        return curs.fetchall()


    def get_annotations(self, seq_names):
        curs = self.conn.cursor()

        for seq in seq_names.split(","):
            seq = seq.replace('"', '')
            curs.execute("SELECT name FROM eggnog WHERE name == ?", (seq,))

            eggnog = curs.fetchone()
            if eggnog is not None:
                curs.execute("SELECT pname FROM seq WHERE name == ?", (seq,))
                pname = curs.fetchone()
                pname = pname[0] if pname is not None else None

                curs.execute("SELECT gos FROM gene_ontology WHERE name == ?", (seq,))
                gos = curs.fetchone()
                gos = gos[0] if gos is not None else None

                curs.execute("SELECT ec, ko, pathway, module, reaction, rclass, brite, tc, cazy FROM kegg WHERE name == ?", (seq,))
                kegg_fields = curs.fetchone()
                kegg_fields = list(kegg_fields) if kegg_fields is not None else [None] * 9

                curs.execute("SELECT reaction FROM bigg WHERE name == ?", (seq,))
                bigg = curs.fetchone()
                bigg = bigg[0] if bigg is not None else None

                curs.execute("SELECT pfam FROM pfam WHERE name == ?", (seq,))
                pfam = curs.fetchone()
                pfam = pfam[0] if pfam is not None else None

                yield [pname, gos] + kegg_fields + [bigg, pfam]

        return

    def get_pfam_annotations(self, seq_names):
        cmd = """SELECT pfam.pfam
            FROM pfam
            WHERE pfam.name in (%s)
            """ % seq_names

        curs = self.conn.cursor()
        s = curs.execute(cmd)
        return curs.fetchall()

    def get_member_events(self, member, target_levels):
        curs = self.conn.cursor()
        curs.execute('SELECT orthoindex FROM orthologs WHERE name == ?', (member.strip(),))
        event_indexes = curs.fetchone()
        if event_indexes is not None and len(event_indexes) > 0:
            event_indexes = str(event_indexes[0])

            # if target_levels is not None and len(target_levels) > 0:
            #     levels = ",".join([f'"{x}"' for x in target_levels])
            #     db.execute(f'SELECT level, side1, side2 FROM event WHERE i IN ({event_indexes}) AND level IN ({levels})')
            # else:
            #     db.execute(f'SELECT level, side1, side2 FROM event WHERE i IN ({event_indexes})')

            # return db.fetchall()

            # this one looks like faster than the code above
            all_queries = [(int(x),y) for x in event_indexes.split(",") for y in target_levels]

            for query in all_queries:
                curs.execute('SELECT level, side1, side2 FROM event WHERE i == ? AND level == ?', query)
                all_levels = curs.fetchall()
                for level, _side1, _side2 in all_levels:
                    yield level, _side1, _side2

        return



## END
