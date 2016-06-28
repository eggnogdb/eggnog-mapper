from __future__ import absolute_import
from collections import Counter
import sqlite3

from pymongo import MongoClient
mongoCli = MongoClient()
mongoDB = mongoCli.eggnog4_1
db_speciation = mongoDB.sp_events
db_members = mongoDB.members

from common import EGGNOGDB_PATH

conn = sqlite3.connect(EGGNOGDB_PATH)
db = conn.cursor()

conn2 = sqlite3.connect('/home/huerta/eggnog-mapper-dev/test.db')
db2 = conn2.cursor()


def get_nogname(name):
    if len(name) == 5:
        return "ENOG41%s" %name
    return name

def get_og_data(ogs):
    in_clause = ','.join(['"%s"' %oname for oname in set(ogs)])    
    cmd = 'SELECT og, level, description, COG_categories, nm, GO_freq, KEGG_freq, SMART_freq from annotations WHERE og in (%s);' %in_clause
    db.execute(cmd)
    og2data = {}
    for og, level, description, cat, nm, GO_freq, KEGG_freq, SMART_freq in db.fetchall():
        cat = re.sub('[^A-Z]', '', cat)
        og2data[og] = [level, description, cat, nm, json.loads(GO_freq), json.loads(KEGG_freq), json.loads(SMART_freq)]
    return og2data
    
def get_preferred_names_dict(names):
    query = {'$or': [{"n":n.split('.', 1)[1], "t":int(n.split('.', 1)[0])} for n in names]}
    return dict([("%s.%s" %(e['t'], e['n']), e['p']) for e in db_members.find(query, {'n':1, 't':1, 'p':1})])

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

    
def refine_orthologs_by_member(target_members, target_taxa=None, target_nogs=None, target_level=None):        
    query = {}
    
    target_seqids = set()
    query_taxa = set()
    for m in target_members: 
        taxid, seqid = m.split(".", 1)
        query_taxa.add(taxid)
        target_seqids.add(seqid)

    if len(target_members) >1 :
        query = {'m': {"$in":target_members}}
    else:
        query = {'m': target_members[0]}

    if target_level:
        query['l'] = target_level
        
    if target_taxa:
        target_taxa = set(map(str, target_taxa))
    else:
        target_taxa = set()
    
    orthology = {}
    target_members = set(target_members)
    for event in db_speciation.find(query, {'z':1, 'm':1, 'n':1, 'l':1}):
        if target_nogs and event['n'] not in target_nogs:
            continue
        
        all_seqs = [m.split(".", 1) for m in event['m']]
        by_sp1, by_sp2 = {}, {}
        for _sp, _side in [(by_sp1, all_seqs[:event['z']]),
                           (by_sp2, all_seqs[event['z']:])]:
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
                    
    #as_list = [((k[0], get_sp(k[0]), k[1]), [(t2, get_sp(t2), co2) for t2,co2 in v])
    #            for k,v in orthology.iteritems()]

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
    
def get_gos(seqnames, ignore_type=None):
    query = ','.join(['"%s"'%n for n in set(seqnames)])    
    cmd = 'SELECT seqname, terms from seq2go WHERE seqname IN (%s);' %query
    db2.execute(cmd)
    by_seq = {}
    all_gos = set()
    for reg in db2.fetchall():
        gos = set()
        for go in reg[1].split(','):
            cat, term, atype, flag = go.split('|')
            if (not ignore_type or atype not in ignore_type):
                gos.add(term)
        
        #gos = [term for go in reg[1].split(',') for cat,term,atype,_ in go.split('|') if (not restrict_type or atype not in ignore_type)]
        by_seq[reg[0]] = set(gos)
        all_gos.update(gos)
       
    return by_seq, all_gos

