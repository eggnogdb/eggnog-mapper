#!/usr/bin/env python

import argparse
import sqlite3
import gzip
import json
import re
import os
import sys
import gzip
from collections import defaultdict

BASEPATH = os.path.split(os.path.abspath(__file__))[0]
DBFILE = os.path.join(BASEPATH, "annotations.db")

conn = sqlite3.connect(DBFILE)
db = conn.cursor()

get_setdict = lambda: defaultdict(set)
get_listdict = lambda: defaultdict(list)
    
def tabprint(buff, *args):
    print >>buff, '\t'.join(map(str, args))

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
        query, og, evalue, score, querylen, hmmfrom, hmmto, qfrom, qto  = line.split('\t')       
        if query not in visited: 
            query_list.append(query)
            visited[query] = 0

        if og != '-' and og != 'ERROR':            
            if args.maxhits and visited[query] == args.maxhits:
                continue

            if args.level and level not in args.level:
                continue

            q_cov = abs(float(qto)-float(qfrom))/float(querylen)
            if q_cov < args.min_q_coverage:
                continue
                
            if args.maxhits: 
                visited[query] += 1

            ogs.add(og)            
            lines.append((query, og, evalue, score, hmmfrom, hmmto, qfrom, qto, q_cov))
            
    print >>sys.stderr, "Loading OG meta data for %s OGs..." %len(ogs)
    og_data = get_og_data(ogs)
    INPUT.close()
    
    print >>sys.stderr, "Annotating %s queries..." %len(query_list)
    per_query = defaultdict(get_listdict)
    
    for query, og, evalue, score, hmmfrom, hmmto, qfrom, qto, q_cov in lines:
        try:
            level, desc, cats, nm, gos, kegg, domain = og_data[og]
        except KeyError: 
            print >>sys.stderr, "target OG not found: %s" % line
            continue
        hitinfo = [og, level, nm, evalue, score, hmmfrom, hmmto, qfrom, qto, q_cov]            

        #tabprint(query, og, level, evalue, score, "Description", desc)
        #tabprint(query, og, level, evalue, score, "COG Categories", cats)
        per_query[query][("desc", desc)].append(hitinfo + ['NA', 'NA'])
        per_query[query][("cats", cats)].append(hitinfo + ['NA', 'NA'])
        for go_cat, terms in gos.iteritems():
            for goid, goname, evidence, nseqs, freq, _ in terms:
                #tabprint(query, og, level, evalue, score, "GO_"+go_cat, goid, goname, evidence, nseqs, freq)
                if float(freq) >= args.min_prevalence:
                    per_query[query][("go", goid, goname)].append(hitinfo + [nseqs, freq])

        for pathway, nseqs, freq, _ in kegg:
            #tabprint(query, og, level, evalue, score, pathway, nseqs, freq)
            if float(freq) >= args.min_prevalence:
                per_query[query][("kegg", pathway)].append(hitinfo+[nseqs, freq])

        for dom_source, terms in domain.iteritems():
            for dom_name, nseqs, freq, _ in terms:
                #tabprint(query, og, level, evalue, score, dom_source, dom_name, nseqs, freq)
                if float(freq) >= args.min_prevalence:
                    per_query[query][(dom_source, dom_name)].append(hitinfo+[nseqs, freq])

    OUT = {
        'kegg': gzip.open(args.output+'.KEGG_pathways.txt.gz', "w:gz"),
        'go':   gzip.open(args.output+'.GeneOntology.txt.gz', "w:gz"),
        'desc': gzip.open(args.output+'.description.txt.gz', "w:gz"),
        'cats': gzip.open(args.output+'.COG_categories.txt.gz', "w:gz"),
        'PFAM': gzip.open(args.output+'.PFAM.txt.gz', "w:gz"),
        'SMART': gzip.open(args.output+'.SMART.txt.gz', "w:gz"),
    }
    
    header = map(str.strip, "query_name, eggNOG_OG, OG_taxonomic_level, OG_size, evalue, score, hmmfrom, hmmto, seqfrom, seqto, query_coverage, nseqs_in_OG_with_term, term_prevalence_in_OG, term_type, term_info".split(','))
    for F in OUT.values():
        print >>F, '#'+'\t'.join(header)
        
    for qname in query_list:
        for k in sorted(per_query[qname]):
            hit_details = per_query[qname][k]
            n_ogs = len(set([v[0] for v in hit_details]))
            n_hits = len(hit_details)

            tabprint(OUT[k[0]], qname, *(hit_details[0]+list(k)))
            if args.full_report:
                for i in xrange(1, len(hit_details)):
                    tabprint(OUT[k[1]], " "*len(qname), *(hit_details[i]+list(k)))

                    
    for v in OUT.values():
        v.close()        
    print >>sys.stderr, "Done!"
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('hitsfile', metavar="hitsfile", nargs=1, help='query file containng a list of hits returned by eggnog_mapper.py')
    parser.add_argument('--full_report', dest='full_report', action='store_true')
    #parser.add_argument('--go', dest='go', action='store_true')    
    #parser.add_argument('--kegg', dest='kegg', action='store_true')
    #parser.add_argument('--desc', dest='desc', action='store_true')
    #parser.add_argument('--smart', dest='smart', action='store_true')
    parser.add_argument('--restrict_level', dest='level', type=str, nargs="+", help="report only hits from the provided taxonomic level")
    parser.add_argument('--min_prevalence', dest='min_prevalence', type=float, default=50.0)
    parser.add_argument('--min_query_coverage', dest='min_q_coverage', type=float, default=0.66)
    parser.add_argument('--output', dest="output", type=str, required=True)
    parser.add_argument('--maxhits', dest='maxhits', type=int, help="report only the first `maxhits` hits")

    args = parser.parse_args()
    if args.level:
        args.level = set(args.level)
        
    main(args)
