
from collections import Counter
import sqlite3
import os

from pymongo import MongoClient
mongoCli = MongoClient()
mongoDB = mongoCli.eggnog4_1
db_speciation = mongoDB.sp_events
db_members = mongoDB.members

from .common import BASE_PATH

conn2 = sqlite3.connect(os.path.join(BASE_PATH, 'db', 'test.db'))
db2 = conn2.cursor()

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
       
        by_seq[reg[0]] = set(gos)
        all_gos.update(gos)
      
    return by_seq, all_gos
def get_nogname(name):
    if len(name) == 5:
        return "ENOG41%s" %name
    return name

def get_preferred_names_dict(names):
    query = {'$or': [{"n":n.split('.', 1)[1], "t":int(n.split('.', 1)[0])} for n in names]}
    return dict([("%s.%s" %(e['t'], e['n']), e['p']) for e in db_members.find(query, {'n':1, 't':1, 'p':1})])
    
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
        targets = target_taxa or list(by_sp2.keys())
        for sp1, co1 in by_sp1.items():
            if target_members & co1:
                key1 = (sp1, tuple(sorted((co1))))
                for sp2 in targets:
                    if sp2 not in by_sp2: continue
                    co2 = by_sp2[sp2]
                    key2 = (sp2, tuple(sorted(co2)))
                    orthology.setdefault(key1, set()).add(key2)
                    
        # merge by side2 coorthologs
        targets = target_taxa or list(by_sp1.keys())
        for sp1, co1 in by_sp2.items():
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

    for k, v in orthology.items():
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



    
