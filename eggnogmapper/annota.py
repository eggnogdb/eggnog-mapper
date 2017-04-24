from __future__ import absolute_import
from collections import Counter, defaultdict
import sqlite3
import re
import time
import multiprocessing

from .common import EGGNOGDB_FILE
from .utils import timeit

conn = None
db = None

cog_cat_cleaner = re.compile('[\[\]u\'\"]+')

def connect():
    global conn, db
    conn = sqlite3.connect(EGGNOGDB_FILE)
    db = conn.cursor()

def get_og_annotations(ogname):
    # og VARCHAR(16) PRIMARY KEY,
    # level VARCHAR(16),
    # nm INTEGER,
    # description TEXT,
    # COG_categories VARCHAR(8),
    # GO_freq TEXT,
    # KEGG_freq TEXT,
    # SMART_freq TEXT,
    # proteins TEXT);
    cmd = 'SELECT level, nm, description, COG_categories FROM og WHERE og.og = "%s"' % ogname
    level, nm, desc, cat = None, None, None, None
    if db.execute(cmd):
        try:
            level, nm, desc, cat = db.fetchone()
        except TypeError:                       # avoids crashing ig no annotation for OG viral
            desc, cat = '', ''
        else:
            cat = re.sub(cog_cat_cleaner, '', cat)
    return level, nm, desc, cat


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
    query = ','.join(map(lambda x: '"%s"'%x, ognames))
    cmd = 'SELECT og.og, description, COG_categories FROM og WHERE og.og IN (%s)' % query
    og2desc = {}
    if db.execute(cmd):
        for og, desc, cat in db.fetchall():
            cat = re.sub(cog_cat_cleaner, '', cat)
            og2desc[og] = [cat, desc]
    return og2desc

def get_annotated_orthologs(target_members, orthotype, excluded_gos, cpu):
    # get speciation events associated to target_members
    in_clause = ','.join(['"%s"' %name.strip() for name in target_members])
    cmd = 'SELECT name, orthoindex FROM member WHERE name IN (%s)' % in_clause
    db.execute(cmd)
    name2index = defaultdict(set)
    event_indexes = set()
    for name, index in db.fetchall():
        indexes = map(int, index.split(','))
        name2index[name].update(indexes)
        event_indexes.update(indexes)

    # get referenced events data
    cmd2 = 'SELECT i, level, side1, side2 FROM event WHERE i IN (%s)' % ','.join(map(str, event_indexes))
    db.execute(cmd2)
    events = {e[0]: e[1:] for e in db.fetchall()}

    # compute orthologs for each target_member
    all_orthologs = set()
    m2or = {}
    # for member in target_members:
    #     orthologs = build_orthologs(set([member]),
    #                                 [events[eid] for eid in name2index[member]])[orthotype]
    #     all_orthologs.update(orthologs)
    #     m2or[member] = orthologs
    cmds = []
    for member in target_members:
        cmds.append((member, [events[eid] for eid in name2index[member]]))

    pool = multiprocessing.Pool(cpu)
    for member, orthologs in pool.imap_unordered(_parallel_orthologs, cmds):
        all_orthologs.update(orthologs[orthotype])
        m2or[member] = orthologs[orthotype]

    # preloads functional annotations
    in_clause = ','.join(['"%s"' % n for n in all_orthologs])
    cmd3 = 'SELECT name, pname, go, kegg FROM member WHERE name in (%s);' % in_clause
    t1 = time.time()
    db.execute(cmd3)
    print time.time() - t1
    functions = {e[0]: e[1:] for e in db.fetchall()}
    print time.time() - t1

    #print len(all_orthologs), list(all_orthologs)[:10], in_clause[:10]
    #print len([b for a,b in functions.iteritems() if b[0] ])
    annotations = {}
    for m, orthologs in m2or.iteritems():
        all_gos = set()
        all_kegg = set()
        all_pnames = Counter()
        for o in orthologs:
            pname, gos, kegg = functions.get(o, ['', '', ''])
            all_gos.update([g.split('|')[1] for g in gos.strip().split(
                ',') if g and g.split('|')[2] not in excluded_gos])
            all_kegg.update(map(lambda x: str(x).strip(), kegg.strip().split(',')))
            all_pnames.update([pname.strip()])
            del all_pnames['']
            all_kegg.discard('')
            all_gos.discard('')
        annotations[m] = [orthologs, all_pnames, all_gos, all_kegg]
    return annotations

def _parallel_orthologs(values):
    _orthologs = build_orthologs(set([values[0]]), values[1])
    return (values[0], _orthologs)

def build_orthologs(target_members, events, target_taxa=None):
    orthology = {}
    for level, _side1, _side2 in events:
        side1 = [m.split('.', 1) for m in _side1.split(',')]
        side2 = [m.split('.', 1) for m in _side2.split(',')]

        by_sp1, by_sp2 = {}, {}
        for _sp, _side in [(by_sp1, side1),
                           (by_sp2, side2)]:
            for t, s in _side:
                if not target_taxa or t in target_taxa or t in query_taxa:
                    mid = "%s.%s" % (t, s)
                    _sp.setdefault(t, set()).add(mid)

        # merge by side1 coorthologs
        targets = target_taxa or by_sp2.keys()
        for sp1, co1 in by_sp1.iteritems():
            if target_members & co1:
                key1 = (sp1, tuple(sorted((co1))))
                for sp2 in targets:
                    if sp2 not in by_sp2:
                        continue
                    co2 = by_sp2[sp2]
                    key2 = (sp2, tuple(sorted(co2)))
                    orthology.setdefault(key1, set()).add(key2)

        # merge by side2 coorthologs
        targets = target_taxa or by_sp1.keys()
        for sp1, co1 in by_sp2.iteritems():
            if target_members & co1:
                key1 = (sp1, tuple(sorted((co1))))
                for sp2 in targets:
                    if sp2 not in by_sp1:
                        continue
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

def parse_gos(gos, target_go_ev, excluded_go_ev):
    selected_gos = set()
    for g in gos.strip().split(','):
        if not g:
            continue
        gocat, gid, gevidence = map(str, g.strip().split('|'))
        if not target_go_ev or gevidence in target_go_ev:
            if not excluded_go_ev or gevidence not in excluded_go_ev:
                selected_gos.add(gid)
    return selected_gos


def get_member_annotations(names, target_go_ev, excluded_go_ev):
    in_clause = ','.join(['"%s"' % n for n in names])
    cmd = 'SELECT name, pname, go, kegg FROM member WHERE name in (%s);' % in_clause
    all_gos = set()
    all_kegg = set()
    all_pnames = Counter()
    db.execute(cmd)
    for name, pname, gos, kegg, in db.fetchall():
        all_gos.update(parse_gos(gos, target_go_ev, excluded_go_ev))
        all_kegg.update(map(lambda x: str(x).strip(), kegg.strip().split(',')))
        all_pnames.update([pname.strip()])
    all_kegg.discard('')
    all_gos.discard('')
    del all_pnames['']

    return all_pnames, all_gos, all_kegg

def get_by_member_annotations(names, target_go_ev, excluded_go_ev):
    in_clause = ','.join(['"%s"' % n for n in names])
    cmd = 'SELECT name, pname, go, kegg FROM member WHERE name in (%s);' % in_clause
    db.execute(cmd)
    by_member = {n:[set(), set(), None] for n in names}
    for name, pname, gos, kegg, in db.fetchall():
        selected_gos = parse_gos(gos, target_go_ev, excluded_go_ev)
        keggs = set(map(lambda x: str(x).strip(), kegg.strip().split(',')))
        by_member[str(name)] = [selected_gos, keggs, str(pname)]
    return by_member

def get_by_member_gos(names, target_go_ev, excluded_go_ev):
    in_clause = ','.join(['"%s"' % n for n in names])
    cmd = 'SELECT name, go FROM member WHERE name in (%s);' % in_clause
    db.execute(cmd)
    by_member = {n:set() for n in names}
    for name, gos in db.fetchall():
        selected_gos = parse_gos(gos, target_go_ev, excluded_go_ev)
        by_member[str(name)] = selected_gos
    return by_member

def get_member_ogs(name):
    cmd = 'SELECT groups FROM member WHERE name == "%s";' % (name)
    db.execute(cmd)
    match = db.fetchone()
    ogs = None
    if match:
        ogs = map(lambda x: str(x).strip(), match[0].split(','))
    return ogs

def get_member_orthologs(member, target_taxa=None, target_levels=None):
    query_taxa = member.split('.', 1)[0]
    target_members = set([member])
    cmd = 'SELECT orthoindex FROM member WHERE name = "%s"' % member.strip()
    db.execute(cmd)
    event_indexes = str(db.fetchone()[0])
    cmd2 = 'SELECT level, side1, side2 FROM event WHERE i IN (%s)' % event_indexes
    if target_levels:
        cmd2 += " AND level IN (%s)" % (','.join(map(lambda x: '"%s"' %x, target_levels)))
    db.execute(cmd2)
    orthology = {}
    for level, _side1, _side2 in db.fetchall():
        side1 = [m.split('.', 1) for m in _side1.split(',')]
        side2 = [m.split('.', 1) for m in _side2.split(',')]

        by_sp1, by_sp2 = {}, {}
        for _sp, _side in [(by_sp1, side1),
                           (by_sp2, side2)]:
            for t, s in _side:
                if not target_taxa or t in target_taxa or t == query_taxa:
                    mid = "%s.%s" % (t, s)
                    _sp.setdefault(t, set()).add(mid)

        # merge by side1 coorthologs
        targets = target_taxa or by_sp2.keys()
        for sp1, co1 in by_sp1.iteritems():
            if target_members & co1:
                key1 = (sp1, tuple(sorted((co1))))
                for sp2 in targets:
                    if sp2 not in by_sp2:
                        continue
                    co2 = by_sp2[sp2]
                    key2 = (sp2, tuple(sorted(co2)))
                    orthology.setdefault(key1, set()).add(key2)

        # merge by side2 coorthologs
        targets = target_taxa or by_sp1.keys()
        for sp1, co1 in by_sp2.iteritems():
            if target_members & co1:
                key1 = (sp1, tuple(sorted((co1))))
                for sp2 in targets:
                    if sp2 not in by_sp1:
                        continue
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
