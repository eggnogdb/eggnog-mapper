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

get_setdict = lambda: defaultdict(set)
get_listdict = lambda: defaultdict(list)
    
def tabprint(*args):
    print '\t'.join(map(str, args))

def dbget(og):
    cmd = 'SELECT level, description, COG_categories, nm, GO_freq, KEGG_freq, SMART_freq from annotations WHERE og="%s";' %og
    db.execute(cmd)    
    level, description, cat, nm, GO_freq, KEGG_freq, SMART_freq = db.fetchone()
    cat = re.sub('[^A-Z]', '', cat)
    return [level, description, cat, nm, json.loads(GO_freq), json.loads(KEGG_freq), json.loads(SMART_freq)]

def get_og_data(ogs):
    in_clause = ','.join(['"%s"' %oname for oname in set(ogs)])    
    cmd = 'SELECT og, level, description, COG_categories, nm, GO_freq, KEGG_freq, SMART_freq from annotations WHERE og in (%s);' %in_clause
    db.execute(cmd)
    og2data = {}
    for og, level, description, cat, nm, GO_freq, KEGG_freq, SMART_freq in db.fetchall():
        cat = re.sub('[^A-Z]', '', cat)
        og2data[og] = [level, description, cat, nm, json.loads(GO_freq), json.loads(KEGG_freq), json.loads(SMART_freq)]
    return og2data
    
    
def main(args):
    if args.hitsfile[0].endswith('.gz'):
        INPUT = gzip.open(args.hitsfile[0], "r:gz")
    else:
        INPUT = open(args.hitsfile[0], "rU")

    # if args.go:
    #     tabprint("# query", "OG", "level", "evalue", "score", "GO category", "GO id", "GO desc", "evindence", "nseqs", "freq in OG (%)")
    # if args.kegg:
    #     tabprint("# query", "OG", "level", "evalue", "score", "KEGG pathway", "nseqs", "freq in OG")
    # if args.smart:
    #     tabprint("# query", "OG", "level", "evalue", "score", "Domain source", "domain name", "nseqs", "freq in OG (%)")

    visited = defaultdict(int)
    query_list = []
    lines = []

    ogs = set()
    for line in INPUT:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        query, og, evalue, score, hmmfrom, hmmto, qfrom, qto, devalue = line.split('\t')       
        if query not in visited: 
            query_list.append(query)
            visited[query] = 0

        if og != '-' and og != 'ERROR':
            if args.maxhits and visited[query] == args.maxhits:
                continue

            if args.level and level not in args.level:
                continue

            if args.maxhits: 
                visited[query] += 1

            ogs.add(og)            
            lines.append((query, og, evalue, score, hmmfrom, hmmto, qfrom, qto, devalue))
            
    print >>sys.stderr, "Loading OG meta data for %s OGs..." %len(ogs)
    og_data = get_og_data(ogs)
    INPUT.close()
    
    print >>sys.stderr, "Annotating %s queries..." %len(query_list)
    per_query = defaultdict(get_listdict)
    
    for query, og, evalue, score, hmmfrom, hmmto, qfrom, qto, devalue in lines:
        try:
            level, desc, cats, nm, gos, kegg, domain = og_data[og]
        except KeyError: 
            print >>sys.stderr, "target OG not found: %s" % line
            continue
        hitinfo = [og, level, nm, evalue, hmmfrom, hmmto, qfrom, qto]            
        if args.desc:
            #tabprint(query, og, level, evalue, score, "Description", desc)
            #tabprint(query, og, level, evalue, score, "COG Categories", cats)
            per_query[query][("desc", desc)].append(hitinfo + [0, 100.0])
            per_query[query][("cats", desc)].append(hitinfo + [0, 100.0])

        if args.go: 
            for go_cat, terms in gos.iteritems():
                for goid, goname, evidence, nseqs, freq, _ in terms:
                    #tabprint(query, og, level, evalue, score, "GO_"+go_cat, goid, goname, evidence, nseqs, freq)
                    per_query[query][("go", goid, goname)].append(hitinfo+[nseqs, freq, evidence])

        if args.kegg: 
            for pathway, nseqs, freq, _ in kegg:
                #tabprint(query, og, level, evalue, score, pathway, nseqs, freq)
                per_query[query][("kegg", pathway)].append(hitinfo+[nseqs, freq])

        if args.smart:
            for dom_source, terms in domain.iteritems():
                for dom_name, nseqs, freq, _ in terms:
                    #tabprint(query, og, level, evalue, score, dom_source, dom_name, nseqs, freq)
                    per_query[query][(dom_source, dom_name)].append(hitinfo+[nseqs, freq])
        #else:
        #    tabprint(query, "-")


    #print tabprint("#query", "nhits", "nOGs", "max freq")
    for qname in query_list:
        for k in sorted(per_query[qname]):
            hit_details = per_query[qname][k]
            n_ogs = len(set([v[0] for v in hit_details]))
            n_hits = len(hit_details)            
            tabprint(qname, n_ogs, n_hits, *(list(k) + hit_details[0]))
            for i in xrange(1, len(hit_details)):
                tabprint(" "*len(qname), "", "", *(list(k) + hit_details[i]))
            
                
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('hitsfile', metavar="hitsfile", nargs=1, help='query file containng a list of hits returned by eggnog_mapper.py')
    parser.add_argument('--go', dest='go', action='store_true')
    parser.add_argument('--kegg', dest='kegg', action='store_true')
    parser.add_argument('--desc', dest='desc', action='store_true')
    parser.add_argument('--smart', dest='smart', action='store_true')
    parser.add_argument('--restrict_level', dest='level', type=str, nargs="+", help="report only hits from the provided taxonomic level")
    parser.add_argument('--nr_output', dest="nr", type=str, help="non redundant outgroup")
    parser.add_argument('--maxhits', dest='maxhits', type=int, help="report only the first `maxhits` hits")

    args = parser.parse_args()
    if args.level:
        args.level = set(args.level)
        
    main(args)
