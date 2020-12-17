##

import sqlite3

from ..common import get_eggnogdb_file, existing_file
from ..utils import colorify

conn = None
db = None

def connect(usemem = False):
    global conn, db
    if not conn:
        eggnogdb_file = get_eggnogdb_file()
        existing_file(eggnogdb_file)
        
        if usemem == True:
            print("Loading source DB...")
            print(colorify("Warning: this can take a few minutes and load up to 45GB to RAM. "
                           "Using --dbmem is recommended to annotate a large number of sequences.", "red"))
            source = sqlite3.connect(eggnogdb_file)
            conn = sqlite3.connect(':memory:')
            source.backup(conn)
            source.close()
        else:
            existing_file(eggnogdb_file)
            conn = sqlite3.connect(eggnogdb_file)

        db = conn.cursor()
        db.execute("PRAGMA synchronous=OFF;")
        db.execute("PRAGMA journal_mode=OFF;")
        db.execute("PRAGMA cache_size=2000;")
        # db.execute("PRAGMA temp_store=MEMORY;")
        # db.execute("PRAGMA locking_mode=EXCLUSIVE;")
        # db.execute("PRAGMA mmap_size=0;")

        
def close():
    global conn
    if conn:
        conn.close()
        conn = None
    return


def get_db_version():
    cmd = 'select * from version LIMIT 1'
    db.execute(cmd)
    return db.fetchone()[0]


def get_member_ogs(name):
    db.execute('SELECT groups FROM eggnog WHERE name == ?;', (name,))
    return db.fetchone()


def get_ogs_description(og, level):
    db.execute('SELECT og.og, nm, description, COG_categories FROM og WHERE og.og == ? AND og.level == ?', (og, level))
    return db.fetchall()


def get_annotations(seq_names):
    
    for seq in seq_names.split(","):
        seq = seq.replace('"', '')
        db.execute("SELECT name FROM eggnog WHERE name == ?", (seq,))

        eggnog = db.fetchone()
        if eggnog is not None:
            db.execute("SELECT pname FROM seq WHERE name == ?", (seq,))
            pname = db.fetchone()
            pname = pname[0] if pname is not None else None
            
            db.execute("SELECT gos FROM gene_ontology WHERE name == ?", (seq,))
            gos = db.fetchone()
            gos = gos[0] if gos is not None else None
            
            db.execute("SELECT ec, ko, pathway, module, reaction, rclass, brite, tc, cazy FROM kegg WHERE name == ?", (seq,))
            kegg_fields = db.fetchone()
            kegg_fields = list(kegg_fields) if kegg_fields is not None else [None] * 9
            
            db.execute("SELECT reaction FROM bigg WHERE name == ?", (seq,))
            bigg = db.fetchone()
            bigg = bigg[0] if bigg is not None else None
            
            db.execute("SELECT pfam FROM pfam WHERE name == ?", (seq,))
            pfam = db.fetchone()
            pfam = pfam[0] if pfam is not None else None
            
            yield [pname, gos] + kegg_fields + [bigg, pfam]
            
    return

    # The next 2 are slower than the one above
    # for seq in seq_names.split(","):
    #     seq = seq.replace('"', '')
    #     db.execute("""SELECT seq.pname, gene_ontology.gos,
    #     kegg.ec, kegg.ko, kegg.pathway, kegg.module, kegg.reaction, kegg.rclass, kegg.brite, kegg.tc, kegg.cazy, 
    #     bigg.reaction, pfam.pfam
    #     FROM eggnog
    #     LEFT JOIN seq on seq.name = eggnog.name
    #     LEFT JOIN gene_ontology on gene_ontology.name = eggnog.name
    #     LEFT JOIN kegg on kegg.name = eggnog.name
    #     LEFT JOIN bigg on bigg.name = eggnog.name
    #     LEFT JOIN pfam on pfam.name = eggnog.name
    #     WHERE eggnog.name == ?
    #     """, (seq,))
        
    #     for field in db.fetchall():
    #         yield field
            
    # return
        
    # cmd = """SELECT seq.pname, gene_ontology.gos,
    # kegg.ec, kegg.ko, kegg.pathway, kegg.module, kegg.reaction, kegg.rclass, kegg.brite, kegg.tc, kegg.cazy, 
    # bigg.reaction, pfam.pfam
    #     FROM eggnog
    #     LEFT JOIN seq on seq.name = eggnog.name
    #     LEFT JOIN gene_ontology on gene_ontology.name = eggnog.name
    #     LEFT JOIN kegg on kegg.name = eggnog.name
    #     LEFT JOIN bigg on bigg.name = eggnog.name
    #     LEFT JOIN pfam on pfam.name = eggnog.name
    #     WHERE eggnog.name in (%s)
    #     """ % seq_names
    
    # s = db.execute(cmd)
    # return db.fetchall()

def get_pfam_annotations(seq_names):
    cmd = """SELECT pfam.pfam
        FROM pfam
        WHERE pfam.name in (%s)
        """ % seq_names
    
    s = db.execute(cmd)
    return db.fetchall()

def get_member_events(member, target_levels):
    db.execute('SELECT orthoindex FROM orthologs WHERE name == ?', (member.strip(),))
    event_indexes = db.fetchone()
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
            db.execute('SELECT level, side1, side2 FROM event WHERE i == ? AND level == ?', query)
            all_levels = db.fetchall()
            for level, _side1, _side2 in all_levels:
                yield level, _side1, _side2

    return



## END
