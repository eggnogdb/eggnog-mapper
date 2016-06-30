#!/usr/bin/env python

__description__ = 'Run protein sequences searches against eggNOG precompiled hmm databases or custom hmm files, allowing to run computaions on memmory.'
__author__ = 'Jaime Huerta Cepas'
__license__ = "GPL v2"

import sys
import os
import time
import cPickle
from collections import defaultdict, Counter
import multiprocessing
import argparse
import re

SCRIPT_PATH = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.insert(0, SCRIPT_PATH)

from eggnogmapper.common import *
from eggnogmapper import search
from eggnogmapper import annota
#from eggnogmapper import annota_mongo
from eggnogmapper import seqio
from eggnogmapper.server import server_functional, load_server

def main(args):
    host = 'localhost'
    if args.db in EGGNOG_DATABASES:
        dbpath, port = get_db_info(args.db)
        if not pexists(dbpath+".h3f"):
            print dbpath
            print download_database(args.db)
            raw_input("press")
        try:
            idmap = cPickle.load(open(dbpath.replace('.hmm', '.pkl'), 'rb'))
        except IOError:
            idmap = None
            
        if port: 
            end_port = port+1
            
        if args.usemem or server_functional(host, port, args.dbtype):
            scantype = 'mem'
        else:
            port = None
            scantype = 'disk'

    elif os.path.isfile(args.db+'.h3f'):
        if args.usemem:
            scantype = 'mem'
            if not args.idmap:
                if server.generate_idmap(args.db):
                    args.idmap = args.db+".idmap"
                    print >>sys.stderr, "idmap succesfully created!"
                else:
                    print >>sys.stderr, "idmap could not be created!"
                    idmap = None
            else:
                idmap = None
            port = 53000
            end_port = 53200
        else:
            idmap = None
            port = None
            dbpath = args.db
            scantype = 'disk'
        
    elif ":" in args.db:
        host, port = args.db.split(":")
        host = host.strip()
        port = int(port)
        idmap = None
        dbpath = None
        scantype = 'mem'
    else:
        ValueError('Invalid database name/server')

    # If memory based searches requested, load server if necessary
    if args.usemem and not server_functional(host, port, args.dbtype):
        master_db, worker_db = None, None
        for try_port in range(port, end_port, 2):
            print >>sys.stderr, "Loading server at localhost, port", try_port, try_port+1
            dbpath, master_db, worker_db = load_server(dbpath, try_port, try_port+1, args.cpu)
            port = try_port

            ready = False
            while 1:
                print >>sys.stderr, "Waiting for server to become ready...", host, try_port
                time.sleep(1)
                if not master_db.is_alive() or not worker_db.is_alive():
                    master_db.terminate()
                    master_db.join()
                    worker_db.terminate()
                    worker_db.join()
                    break
                elif server_functional(host, port, args.dbtype):
                    ready = True
                    break
            if ready:
                break
    if port:
        scantype = 'mem'        
    else:
        scantype = 'disk'
        host = dbpath
    
    # Exits if trying to connect to a dead server
    if port and not server_functional(host, port, args.dbtype):
        print >>sys.stderr, "hmmpgmd Server not found at %s:%s" %(host, port)
        exit(-1)

    # Loads idmaps if file is provided 
    if args.idmap:
        print >>sys.stderr, "Reading idmap"
        idmap = {}
        for _lnum, _line in enumerate(open(args.idmap)):
            if not _line.strip():
                continue
            try:
                _seqid, _seqname = map(str, _line.strip().split(' '))
            except ValueError:
                if _lnum == 0:
                    continue # idmap generate by esl_reformat has an info line at begining
                else:
                    raise
            _seqid = int(_seqid)
            idmap[_seqid] = [_seqname]
        print >>sys.stderr, len(idmap), "names loaded"        

    hits_file = "%s.hits" %args.output
    annot_file = "%s.annot" %args.output

    if pexists(hits_file) and not args.resume and not args.override:
        print "Output files already present. User --resume or --override to continue"
        sys.exit(1)
            
    # Start scanning sequences 
    if not args.annotate_only:
        VISITED = set()
        if args.resume:
            print "Resuming previous run. Reading computed output from", hits_file
            VISITED = set([line.split('\t')[0].strip() for line in open(hits_file) if not line.startswith('#')])
            print len(VISITED), 'processed queries skipped'
            OUT = open(hits_file, 'a')
        else:
            OUT = open(hits_file, 'w')

        print >>sys.stderr, "Analysis starts now."

        HEADER = ['query', 'hit', 'e-value', 'sum_score', 'query_length', 'hmmfrom', 'hmmto', 'seqfrom', 'seqto', 'q_coverage']
        print >>OUT, '# ' + time.ctime()
        print >>OUT, '# ' + ' '.join(sys.argv)
        print >>OUT, '# ' + '\t'.join(HEADER)
        total_time = 0
        last_time = time.time()
        start_time = time.time()
        qn = 0

        for qn, (name, elapsed, hits, querylen, seq) in enumerate(search.iter_hits(args.input,
                                                                                   args.translate,
                                                                                   args.qtype,
                                                                                   args.dbtype,
                                                                                   scantype,
                                                                                   host,
                                                                                   port,
                                                                                   evalue_thr=args.evalue,
                                                                                   score_thr=args.score,                                                                               
                                                                                   qcov_thr=args.qcov,
                                                                                   fixed_Z=args.Z,
                                                                                   max_hits=args.maxhits,
                                                                                   skip=VISITED,
                                                                                   maxseqlen=args.maxseqlen,
                                                                                   cpus=args.cpu)):
            if elapsed == -1:
                # error occured
                print >>OUT, '\t'.join([name] + ['ERROR'] * len(HEADER))
            elif not hits:
                print >>OUT, '\t'.join([name] + ['-'] * len(HEADER))
            else:
                for hitindex, (hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore) in enumerate(hits):
                    hitname = hid
                    if idmap:
                        hitname = idmap[hid][0]
                                                    
                    print >>OUT, '\t'.join(map(str, [name, hitname, heval, hscore, int(querylen), int(hmmfrom), int(hmmto), int(sqfrom), int(sqto), float(sqto-sqfrom)/querylen]))
            OUT.flush()

            # monitoring
            total_time += time.time() - last_time
            last_time = time.time()
            if qn and (qn % 25 == 0):
                print >>sys.stderr, qn, total_time, "%0.2f q/s" %((float(qn)/total_time))
                sys.stderr.flush()

        # finish
        ellapsed_time = time.time()-start_time
        print >>sys.stderr, qn, total_time, "%0.2f q/s" %((float(qn)/ellapsed_time))
        sys.stderr.flush()
        print >>OUT, '# %d queries scanned' %(qn + 1)
        print >>OUT, '# Total time (seconds):', ellapsed_time
        OUT.close()

    start_time = time.time()
    if not args.hits_only and args.db in EGGNOG_DATABASES:
        OUT = open(annot_file, "w")
        print >>OUT, '\t'.join("#query_seq, best_hit_eggNOG_ortholog, best_hit_evalue, best_hit_score, predicted_name, strict_orthologs, GO, KEGG(pathway)".split(','))
        for qn, r in enumerate(process_hits_file(hits_file, args.input, args.db, translate=args.translate, cpu=args.cpu)):
            if qn and (qn % 25 == 0):
                total_time = time.time() - start_time
                print >>sys.stderr, qn, total_time, "%0.2f q/s" %((float(qn)/total_time))
                sys.stderr.flush()
            
            query_name = r[0]
            best_hit_name = r[1]
            best_hit_evalue = r[2]
            best_hit_score = r[3]
            if best_hit_name != '-' and float(best_hit_score) >= 20: 
                #_orthologs = sorted(annota_mongo.refine_orthologs_by_member([best_hit_name])['one2one'])
                orthologs = sorted(annota.get_member_orthologs(best_hit_name)[args.orthotype])
                
                if orthologs:
                    pname, gos, keggs = annota.get_member_annotations(orthologs, excluded_gos=set(["IEA", "ND"]))
                    #_pname = Counter(annota_mongo.get_preferred_names_dict(orthologs).values())
                    name_ranking = sorted(pname.items(), key=lambda x:x[1], reverse=True)
                else:
                    name_ranking = [[0, 0, 0]]
            
                if name_ranking[0][1] > 2:
                    best_name = name_ranking[0][0]
                else:
                    best_name = '-'
                    
                # TEST
                #by_seq, _gos = annota_mongo.get_gos(orthologs, set(["IEA", "ND"]))
                # assert sorted(orthologs) == sorted(_orthologs)
                # assert sorted(pname.items()) == sorted(_pname.items())
                # print sorted(orthologs) == sorted(_orthologs)
                # print sorted(pname.items()) == sorted(_pname.items())
                
                print >>OUT, '\t'.join(map(str, (query_name, best_hit_name, best_hit_evalue, best_hit_score, best_name,
                                                 ','.join(orthologs),
                                                 ','.join(sorted(gos)),
                                                 ','.join(sorted(keggs))
                                             )))

        print >>OUT, '# Total time (seconds):', time.time()-start_time
        OUT.close()


    for p in multiprocessing.active_children():
        p.terminate()
    print 'finish'

def process_hits_file(hits_file, query_fasta, eggnogdb, skip_queries=None, translate=False, cpu=1):        
    if eggnogdb in set(['euk', 'bact', 'arch']):
        print 'loading OG levels'
        og2level = cPickle.load(open(pjoin(DATA_PATH, 'og2level.pkl')))
        level = None
    else:
        og2level = None
        level = eggnogdb
        
    sequences = {name:seq for name, seq in seqio.iter_fasta_seqs(query_fasta, translate=translate)}    
    cmds = []
    visited_queries = set()

    if skip_queries:
        visited_queries.update(skip_queries)
    
    for line in gopen(hits_file):
        if line.startswith('#'):
            continue        

        fields = map(str.strip, line.split('\t'))
        seqname = fields[0]
        hitname = fields[1]
        
        if hitname == '-' or hitname == 'ERROR':
            continue
            
        if seqname in visited_queries:
            continue

        # names in hmm database files are not clean eggnog OG names
        m = re.search('\w+\.((ENOG41|COG|KOG|arCOG)\w+)\.', hitname)
        if m: 
            hitname = re.sub("^ENOG41", "", m.groups()[0])

        if eggnogdb != 'viruses':
            seq = sequences[seqname] 
            visited_queries.add(seqname)
            if og2level:
                level = og2level[hitname]            
            target_fasta = os.path.join(FASTA_PATH, level, "%s.fa" %hitname)
            cmds.append([seqname, seq, target_fasta])

    if cmds:
        pool = multiprocessing.Pool(cpu)
        print 'Predicting orthologs...'
        for r in pool.imap(search.refine_hit, cmds):
            yield r
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # server
    g1 = parser.add_argument_group('Target database options')
    g1.add_argument('--db',  dest='db', help='specify the target database for sequence searches. Choose among: euk,bact,arch, host:port, or local hmmpressed database')    
    g1.add_argument('--dbtype', dest="dbtype", choices=["hmmdb", "seqdb"], default="hmmdb")
    g1.add_argument('--qtype',  choices=["hmm", "seq"], default="seq")
    g1.add_argument('--idmap', dest='idmap', type=str)
    
    g2 = parser.add_argument_group('Sequence search options')
    g2.add_argument('--maxhits', dest='maxhits', type=int, default=25,
                    help="Max number of hits to report. Default=25")
    g2.add_argument('--evalue', dest='evalue', default=0.001, type=float,
                    help="E-value threshold. Default=0.001")
    g2.add_argument('--score', dest='score', default=20, type=float,
                    help="Bit score threshold. Default=20")
    g2.add_argument('--maxseqlen', type=int, default=5000,
                    help="Ignore query sequences larger than `maxseqlen`. Default=5000")
    g2.add_argument('--qcov', type=float,
                    help="min query coverage (from 0 to 1). Default=(disabled)")
    g2.add_argument('--Z', dest='Z', type=float, default=40000000,
                    help='fix database size (allows comparing e-values among databases). Default=40,000,000')

    g2.add_argument('--othologs_evalue', dest='ortho_evalue', default=0.001, type=float,
                    help="E-value threshold for ortholog dectection. Default=0.001")
    g2.add_argument('--orthologs_score', dest='ortho_score', default=60, type=float,
                    help="Bit score threshold for ortholog detection. Default=60")
    
    
    g3 = parser.add_argument_group('Output options')
    g3.add_argument('--output', type=str, help="base name for output files", required=True)
    g3.add_argument('--resume', action="store_true",
                    help="Resumes a previous execution skipping reported hits in the output file.")
    g3.add_argument('--override', action="store_true",
                    help="Overwrites output files if they exist.")
    g3.add_argument("--hits_only", action="store_true",
                    help="Skip fine-grained orthology basedannotation, reporting only HMM hits.")    
    g3.add_argument("--annotate_only", action="store_true",
                    help="Skip mapping. Use existing hits file")
    
    # exec mode 
    g4 = parser.add_argument_group('Exection options') 
    g4.add_argument('-i', dest="input",
                    help='Computes annotations for the provided FASTA file')

    g4.add_argument('--translate', action="store_true",
                    help='Assume sequences are genes instead of proteins')

    g4.add_argument("--servermode", action="store_true",
                    help='Loads target database in memory and keeps running in server mode,'
                    ' so another instance of eggnog-mapper can connect to this sever.'
                    ' Auto turns on the --usemem flag')
    
    g4.add_argument('--usemem', action="store_true",
                    help="""If a local hmmpressed database is provided as target using --db,
                    this flag will allocate the whole database in memory using hmmpgmd.
                    Database will be unloaded after execution.""")
    
    g4.add_argument('--cpu', type=int, default=2)

    parser.add_argument('--orthotype', choices=["one2one", "many2one", "one2many", "many2many", "all"],
                      default="one2one")
    
    args = parser.parse_args()

    if args.servermode and args.input:
        parser.error('Incompatible execution modes. Choose between -i or --servermode')
    if not (args.servermode or args.input):
        parser.error('Execution must be specified. Choose between: -i [fastafile]  or --servermode')
    if not args.db:
        parser.error('A target databse must be specified with --db')

    if args.db in EGGNOG_DATABASES and not args.hits_only:
        annota.connect()
        
    main(args)








    
