from __future__ import absolute_import
from collections import Counter
import sqlite3
from common import EGGNOGDB_FILE

conn = None
db = None

def connect():
    global conn, db
    conn = sqlite3.connect(EGGNOGDB_FILE)
    db = conn.cursor()

def get_member_annotations(names, excluded_gos):
    in_clause = ','.join(['"%s"' %n for n in names])    
    cmd = 'SELECT name, pname, go, kegg FROM member WHERE name in (%s);' %in_clause
    all_gos = Counter()
    all_kegg = Counter()
    all_pnames = Counter()
    db.execute(cmd)
    for name, pname, gos, kegg, in db.fetchall():
        all_gos.update([g.split('|')[1] for g in gos.strip().split(',') if g and g.split('|')[2] not in excluded_gos])
        all_kegg.update(map(lambda x: str(x).strip(), kegg.strip().split(',')))
        all_pnames.update([pname.strip()])
    del all_kegg['']
    del all_gos['']
    
    return all_pnames, all_gos, all_kegg

def get_member_orthologs(member, target_taxa=None):
    target_members = set([member])
    cmd = 'SELECT orthoindex FROM member WHERE name = "%s"' %member.strip()
    db.execute(cmd)
    event_indexes = str(db.fetchone()[0])
    cmd2 = 'SELECT level, side1, side2 FROM event WHERE i IN (%s)' %event_indexes
    db.execute(cmd2)
    orthology = {}
    for level, _side1, _side2 in db.fetchall():
        side1 = [m.split('.', 1) for m in _side1.split(',')]
        side2 = [m.split('.', 1) for m in _side2.split(',')]
                
        by_sp1, by_sp2 = {}, {}
        for _sp, _side in [(by_sp1, side1),
                           (by_sp2, side2)]:
            for t, s in _side:
                if not target_taxa or t in target_taxa or t in query_taxa:
                    mid = "%s.%s" %(t, s)
                    _sp.setdefault(t, set()).add(mid)

        # merge by side1 coorthologs
        targets = target_taxa or by_sp2.keys()
        for sp1, co1 in by_sp1.iteritems():
            if target_members & co1:
                key1 = (sp1, tuple(sorted((co1))))
                for sp2 in targets:
                    if sp2 not in by_sp2: continue
                    co2 = by_sp2[sp2]
                    key2 = (sp2, tuple(sorted(co2)))
                    orthology.setdefault(key1, set()).add(key2)
                    
        # merge by side2 coorthologs
        targets = target_taxa or by_sp1.keys()
        for sp1, co1 in by_sp2.iteritems():
            if target_members & co1:
                key1 = (sp1, tuple(sorted((co1))))
                for sp2 in targets:
                    if sp2 not in by_sp1: continue
                    co2 = by_sp1[sp2]
                    key2 = (sp2, tuple(sorted(co2)))
                    orthology.setdefault(key1, set()).add(key2)

    all_orthologs = {
        "one2one": set(),
        "one2many": set(),
        "many2many": set(),
        "many2one": set(),
        "all": set(),
    }

    for k, v in orthology.iteritems():
        if len(k[1]) == 1:
            otype_prefix = "one2"
        else:
            otype_prefix = "many2"            
        all_orthologs['all'].update(k[1])            
        for t2, co2 in v:
            if len(co2) == 1:
                otype = otype_prefix + "one"
            else:
                otype = otype_prefix + "many"                
            all_orthologs[otype].update(k[1])
            all_orthologs[otype].update(co2)
            all_orthologs['all'].update(co2)
    return all_orthologs
    
    

