import sys
from collections import defaultdict
import cPickle
import re

from pymongo import MongoClient
import sqlite3

mongoCli = MongoClient()
mongoDB = mongoCli.eggnog4_1
db_speciation = mongoDB.sp_events
db_members = mongoDB.members

TARGETS = set(sys.argv[1:])
assert TARGETS

#og2group
#og2group = dict([tuple(map(str.strip, line.split())) for line in open('og2groups.tsv')])
#with open("og2group.pkl", "wb") as GROUPS:
#    cPickle.dump(og2group, GROUPS, 2)
og2group = cPickle.load(open("og2group.pkl"))

# Load new names
print "loading name convesion"
old2new = {}
for line in open("group_name_conversion.tsv"):
    old, new = map(str.strip, line.split('\t'))
    if old.startswith("COG") or old.startswith("KOG") or old.startswith("arCOG"):
        old2new[old] = old
    else:
        old2new[old] = new


if "events" in TARGETS:
    # load/dump events
    seq2index = defaultdict(list)
    EVENTS = open('events.tsv', 'w')
    nevents = db_speciation.count()
    for index, e in enumerate(db_speciation.find({}, {'z':1, 'm':1, 'n':1, 'l':1})):
        for member in e["m"]:
            seq2index[member].append(index)
        if index % 10000 == 0:
            print >>sys.stderr, "\r",index,nevents,
            sys.stdout.flush()
        s2 = e['m'][:e['z']]
        s1 = e['m'][e['z']:]
        print >>EVENTS, '\t'.join(map(str, [index, e["l"], old2new[e["n"]], ','.join(s1), ','.join(s2)]))
    EVENTS.close()
    with open('seq2index.pkl', 'w') as OUT:
        cPickle.dump(seq2index, OUT)


if 'members' in TARGETS:
    print 'loading seq2index'
    seq2index = cPickle.load(open('seq2index.pkl'))
    # dump member annotations
    member2nogs = {}
    print "loading member nogs"
    for line in open("member2nogs.tsv"):
        member, nogs = map(str.strip, line.split('\t'))
        level_groups = ["%s@%s" %(re.sub('^ENOG41', '', old2new[n]),
                                  og2group[re.sub('^ENOG41', '', old2new[n])]) for n in nogs.split(',')]
        level_groups.sort()
        member2nogs[member] = ','.join(level_groups)

    print 'dumping member annotations'
    ANN = open('members.tsv', 'w')
    print >>sys.stderr, "Dumping..."
    for m in db_members.find({}, {'n':1, 't':1, 'p':1, 'go':1, 'kg':1}):
        gos = ["%s|%s|%s" %(g["c"], g["n"], g["e"]) for g in m.get("go", [])]
        kegg = m.get("kg", [])
        name = "%s.%s" %(m["t"], m["n"])
        print >>ANN, '\t'.join([name, m["p"], member2nogs[name], ','.join(gos), ','.join(map(str, kegg)), ','.join(map(str, seq2index.get(name, [])))])
    ANN.close()
