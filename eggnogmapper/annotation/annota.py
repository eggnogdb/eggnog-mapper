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

def get_ogs_annotations(ognames):
    # og VARCHAR(16) PRIMARY KEY,
    # level VARCHAR(16),
    # nm INTEGER,
    # description TEXT,
    # COG_categories VARCHAR(8),
    # GO_freq TEXT,
    # KEGG_freq TEXT,
    # SMART_freq TEXT,
    # proteins TEXT);
    query = ','.join(['"%s"'%x for x in ognames])
    cmd = 'SELECT og.og, description, COG_categories FROM og WHERE og.og IN (%s)' % query
    og2desc = {}
    if db.execute(cmd):
        for og, desc, cat in db.fetchall():
            cat = re.sub(cog_cat_cleaner, '', cat)
            og2desc[og] = [cat, desc]
    
    return og2desc

def get_best_og_description(ognames):
    query = ','.join(['"%s"'%x.split('@')[0] for x in ognames])
    cmd = 'SELECT og.og, nm, description, COG_categories FROM og WHERE og.og IN (%s)' % query
    best = [None, '', '']
    if db.execute(cmd):
        for og, nm, desc, cat in db.fetchall():
            nm = int(nm) # number of members (proteins) in this OG
            desc = desc.strip()
            if desc and desc != 'N/A' and desc != 'NA':
                if not best[0] or nm <= best[0]:
                    cat = re.sub(cog_cat_cleaner, '', cat)
                    best = [nm, cat, desc]
    
    return best[1], best[2]

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


def get_member_ogs(name):
    cmd = 'SELECT groups FROM eggnog WHERE name == "%s";' % (name)
    db.execute(cmd)
    match = db.fetchone()
    ogs = None
    if match:
        ogs = [str(x).strip() for x in match[0].split(',')]
    return ogs


# def set_coorthologs(by_sp1, by_sp2, target_members, target_taxa, orthology):
def set_coorthologs(by_sp1, by_sp2, target_members, orthology):
    # targets = target_taxa or list(by_sp2.keys())

    for sp1, co1 in by_sp1.items():
        if target_members & co1:
            # print("target member in sp1")
            key1 = (sp1, tuple(sorted((co1))))

            # for sp2 in targets:
            #     if sp2 not in by_sp2:
            #         continue
            #     co2 = by_sp2[sp2]
            #     key2 = (sp2, tuple(sorted(co2)))
            #     orthology.setdefault(key1, set()).add(key2)
            for sp2, co2 in by_sp2.items():
                # if sp2 not in by_sp2:
                #     continue
                # co2 = by_sp2[sp2]
                key2 = (sp2, tuple(sorted(co2)))
                orthology.setdefault(key1, set()).add(key2)
                

    return


def by_species(side, target_taxa, query_taxa):
    by_sp = {}
    for t, s in side:            
        if not target_taxa or t in target_taxa or t == query_taxa:
            mid = "%s.%s" % (t, s)
            by_sp.setdefault(t, set()).add(mid)
            
    return by_sp


def setup_orthology(member, target_taxa, target_levels):
    orthology = {}
    
    query_taxa = member.split('.', 1)[0]
    if target_taxa:
        target_taxa = list(map(str, target_taxa))
    member_as_set = set([member])
    
    cmd = 'SELECT orthoindex FROM orthologs WHERE name = "%s"' % member.strip()
    db.execute(cmd)
    event_indexes = str(db.fetchone()[0])
    
    cmd2 = 'SELECT level, side1, side2 FROM event WHERE i IN (%s)' % event_indexes
    if target_levels:
        cmd2 += " AND level IN (%s)" % (','.join(['"%s"' %x for x in target_levels]))
    db.execute(cmd2)

    for level, _side1, _side2 in db.fetchall():

        # print()
        # print("annota: event result")
        # print(" - ".join([str(x) for x in [level, _side1, _side2]]))
        # orthology = {} # REMOVE THIS

        side1 = [m.split('.', 1) for m in _side1.split(',')]
        side2 = [m.split('.', 1) for m in _side2.split(',')]

        # filter by taxa (by species)
        by_sp1 = by_species(side1, target_taxa, query_taxa)
        by_sp2 = by_species(side2, target_taxa, query_taxa)

        # print()
        # print(by_sp1)
        # print("-")
        # print(by_sp2)
        
        # merge by coorthologs
        # set_coorthologs(by_sp1, by_sp2, member_as_set, target_taxa, orthology)
        # set_coorthologs(by_sp2, by_sp1, member_as_set, target_taxa, orthology)
        set_coorthologs(by_sp1, by_sp2, member_as_set, orthology)
        set_coorthologs(by_sp2, by_sp1, member_as_set, orthology)        
        
        # print(orthology)
    
    return orthology

def get_member_orthologs(member, target_taxa=None, target_levels=None):

    all_orthologs = {
        "one2one": set(),
        "one2many": set(),
        "many2many": set(),
        "many2one": set(),
        "all": set(),
    }

    # print("annota: get_member_orthologs")
    # print(member)
    # print(target_taxa)
    # print(target_levels)
    
    orthology = setup_orthology(member, target_taxa, target_levels)

    for k, v in orthology.items():

        all_orthologs['all'].update(k[1])
        
        if len(k[1]) == 1:
            otype_prefix = "one2"
        else:
            otype_prefix = "many2"
        
        for t2, co2 in v:

            all_orthologs['all'].update(co2)
            
            if len(co2) == 1:
                otype = otype_prefix + "one"
            else:
                otype = otype_prefix + "many"
                
            all_orthologs[otype].update(k[1])
            all_orthologs[otype].update(co2)


    return all_orthologs

## END
