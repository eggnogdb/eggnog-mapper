#!/usr/bin/env python

import argparse
import sqlite3
import gzip
import json
import re
import os
import sys
from collections import defaultdict

BASEPATH = os.path.split(os.path.abspath(__file__))[0]
DBFILE = os.path.join(BASEPATH, "annotations.db")
print DBFILE 
conn = sqlite3.connect(DBFILE)
db = conn.cursor()

def tabprint(*args):
    print '\t'.join(map(str, args))

def dbget(og):
    cmd = 'SELECT level, description, COG_categories, nm, GO_freq, KEGG_freq, SMART_freq from annotations WHERE og="%s";' %og
    db.execute(cmd)    
    level, description, cat, nm, GO_freq, KEGG_freq, SMART_freq = db.fetchone()
    cat = re.sub('[^A-Z]', '', cat)
    return [level, description, cat, nm, json.loads(GO_freq), json.loads(KEGG_freq), json.loads(SMART_freq)]
        
def main(args):
    if args.hitsfile[0].endswith('.gz'):
        INPUT = gzip.open(args.hitsfile[0], "r:gz")
    else:
        INPUT = open(args.hitsfile[0], "rU")

    if args.go:
        tabprint("# query", "OG", "level", "evalue", "score", "GO category", "GO id", "GO desc", "evindence", "nseqs", "freq in OG (%)")
    if args.kegg:
        tabprint("# query", "OG", "level", "evalue", "score", "KEGG pathway", "nseqs", "freq in OG")
    if args.smart:
        tabprint("# query", "OG", "level", "evalue", "score", "Domain source", "domain name", "nseqs", "freq in OG (%)")
        
    visited = defaultdict(int)
    with INPUT:        
        for line in INPUT:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            query, og, evalue, score = line.split('\t')
            if og != '-' and og != 'ERROR':
                try:
                    level, desc, cats, nm, gos, kegg, domain = dbget(og)
                except Exception: 
                    print >>sys.stderr, "Problem processing line:\n%s" % line
                    raise 
                if args.maxhits and visited[query] == args.maxhits:
                    continue
                
                if args.level and level not in args.level:
                    continue
                
                if args.maxhits: 
                    visited[query] += 1

                    
                if args.desc:
                    tabprint(query, og, level, evalue, score, "Description", desc)
                    tabprint(query, og, level, evalue, score, "COG Categories", cats)
                if args.go: 
                    for go_cat, terms in gos.iteritems():
                        for goid, goname, evidence, nseqs, freq, _ in terms:
                            tabprint(query, og, level, evalue, score, "GO_"+go_cat, goid, goname, evidence, nseqs, freq)
                if args.kegg: 
                    for pathway, nseqs, freq, _ in kegg:
                        tabprint(query, og, level, evalue, score, pathway, nseqs, freq)                        
                if args.smart:
                    for dom_source, terms in domain.iteritems():
                        for dom_name, nseqs, freq, _ in terms:
                            tabprint(query, og, level, evalue, score, dom_source, dom_name, nseqs, freq)
            else:
                tabprint(query, "-")
                
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('hitsfile', metavar="hitsfile", nargs=1, help='query file containng a list of hits returned by eggnog_mapper.py')
    parser.add_argument('--go', dest='go', action='store_true')
    parser.add_argument('--kegg', dest='kegg', action='store_true')
    parser.add_argument('--desc', dest='desc', action='store_true')
    parser.add_argument('--smart', dest='smart', action='store_true')
    parser.add_argument('--restrict_level', dest='level', type=str, nargs="+", help="report only hits from the provided taxonomic level")
    parser.add_argument('--maxhits', dest='maxhits', type=int, help="report only the first `maxhits` hits")

    args = parser.parse_args()
    if args.level:
        args.level = set(args.level)
        
    main(args)
