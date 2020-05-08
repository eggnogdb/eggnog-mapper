#
## JHCepas


from collections import Counter, defaultdict
import sqlite3
import re
import time
import multiprocessing

from ..common import get_eggnogdb_file, ANNOTATIONS_HEADER
from ..utils import timeit

conn = None
db = None

cog_cat_cleaner = re.compile('[\[\]u\'\"]+')


def connect():
    global conn, db
    if not conn:
        conn = sqlite3.connect(get_eggnogdb_file())
        db = conn.cursor()

        
def close():
    global conn
    if conn:
        conn.close()

        
def get_member_ogs(name):
    cmd = 'SELECT groups FROM eggnog WHERE name == "%s";' % (name)
    db.execute(cmd)
    match = db.fetchone()
    ogs = None
    if match:
        ogs = [str(x).strip() for x in match[0].split(',')]
    return ogs


def get_deepest_og_description(deeper_og):
    cmd = 'SELECT og.og, nm, description, COG_categories FROM og WHERE og.og IN ("%s")' % deeper_og
    best = [None, '', '']
    if db.execute(cmd):
        for og, nm, desc, cat in db.fetchall():
            desc = desc.strip()
            if desc and desc != 'N/A' and desc != 'NA':
                best = [nm, cat, desc]
                break
    
    return best[1], best[2]


def summarize_annotations(seq_names, target_go_ev, excluded_go_ev):
    in_clause = ','.join(['"%s"' % n for n in seq_names])
    cmd = """SELECT seq.pname, gene_ontology.gos,
    kegg.ec, kegg.ko, kegg.pathway, kegg.module, kegg.reaction, kegg.rclass, kegg.brite, kegg.tc, kegg.cazy, 
    bigg.reaction
        FROM eggnog
        LEFT JOIN seq on seq.name = eggnog.name
        LEFT JOIN gene_ontology on gene_ontology.name = eggnog.name
        LEFT JOIN kegg on kegg.name = eggnog.name
        LEFT JOIN bigg on bigg.name = eggnog.name
        WHERE eggnog.name in (%s)
        """ %in_clause

    annotations = defaultdict(Counter)
    s = db.execute(cmd)
    results = db.fetchall()

    for fields in results:
        for i, h in enumerate(ANNOTATIONS_HEADER):
            if not fields[i]:
                continue
            if h == 'GOs':
                gos = fields[i]
                annotations[h].update(parse_gos(gos, target_go_ev, excluded_go_ev))
            elif h == 'Preferred_name':
                annotations[h].update([fields[i].strip()])
            else:
                values = [str(x).strip() for x in fields[i].split(',')]
                annotations[h].update(values)

    for h in annotations:
        del annotations[h]['']

    if annotations:
        try:
            pname = annotations['Preferred_name'].most_common(1)
            if pname: 
                name_candidate, freq = annotations['Preferred_name'].most_common(1)[0]
            else:
                freq =  0
        except:
            print(annotations)
            raise 
        if freq >= 2:
            annotations['Preferred_name'] = [name_candidate]
        else:
            annotations['Preferred_name'] = ['']

    return annotations


def parse_gos(gos, target_go_ev, excluded_go_ev):
    selected_gos = set()
    for g in gos.strip().split(','):
        if not g:
            continue
        gocat, gid, gevidence = list(map(str, g.strip().split('|')))
        if not target_go_ev or gevidence in target_go_ev:
            if not excluded_go_ev or gevidence not in excluded_go_ev:
                selected_gos.add(gid)
    return selected_gos


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
