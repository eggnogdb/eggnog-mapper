##

import sqlite3

from ..common import get_eggnogdb_file

conn = None
db = None


def connect():
    global conn, db
    if not conn:
        conn = sqlite3.connect(get_eggnogdb_file())
        db = conn.cursor()

        
def close():
    global conn
    if conn:
        conn.close()
        conn = None
    return


def get_db_version():
    cmd = 'select * from version'
    db.execute(cmd)
    return db.fetchone()[0]


def get_member_ogs(name):
    cmd = 'SELECT groups FROM eggnog WHERE name == "%s";' % (name)
    db.execute(cmd)
    return db.fetchone()


def get_ogs_description(ogs):
    cmd = 'SELECT og.og, nm, description, COG_categories FROM og WHERE og.og IN ("%s")' % ogs
    db.execute(cmd)
    return db.fetchall()


def get_annotations(seq_names):
    cmd = """SELECT seq.pname, gene_ontology.gos,
    kegg.ec, kegg.ko, kegg.pathway, kegg.module, kegg.reaction, kegg.rclass, kegg.brite, kegg.tc, kegg.cazy, 
    bigg.reaction, pfam.pfam
        FROM eggnog
        LEFT JOIN seq on seq.name = eggnog.name
        LEFT JOIN gene_ontology on gene_ontology.name = eggnog.name
        LEFT JOIN kegg on kegg.name = eggnog.name
        LEFT JOIN bigg on bigg.name = eggnog.name
        LEFT JOIN pfam on pfam.name = eggnog.name
        WHERE eggnog.name in (%s)
        """ % seq_names
    
    s = db.execute(cmd)
    return db.fetchall()


def get_member_events(member, target_levels):
    
    cmd = 'SELECT orthoindex FROM orthologs WHERE name = "%s"' % member.strip()
    db.execute(cmd)
    event_indexes = str(db.fetchone()[0])
    
    cmd2 = 'SELECT level, side1, side2 FROM event WHERE i IN (%s)' % event_indexes
    if target_levels:
        cmd2 += " AND level IN (%s)" % (','.join(['"%s"' %x for x in target_levels]))
    db.execute(cmd2)

    return db.fetchall()


## END
