#!/usr/bin/env python
import sys
import os
import time
import cPickle
from collections import defaultdict, Counter
import multiprocessing
import argparse
import re
import shutil

SCRIPT_PATH = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.insert(0, SCRIPT_PATH)

from eggnogmapper.common import *
from eggnogmapper import search
from eggnogmapper import annota
#from eggnogmapper import annota_mongo
from eggnogmapper import seqio
from eggnogmapper.server import (
    server_functional, load_server, generate_idmap, shutdown_server)
from eggnogmapper.utils import colorify

__description__ = ('A program for bulk functional annotation of novel '
                    'sequences using the EggNOG orthology database')
__author__ = 'Jaime Huerta Cepas'
__license__ = "GPL v2"

def cleanup_og_name(name):
    # names in hmm database files are not clean eggnog OG names
    m = re.search('\w+\.((ENOG41|COG|KOG|arCOG)\w+)\.', name)
    if m:
        name = m.groups()[0]
    name = re.sub("^ENOG41", "", name)
    return name


def main(args):
    host = 'localhost'
    idmap = None
    if args.usemem:
        scantype = 'mem'
    else:
        scantype = 'disk'
    connecting_to_server = False
    if args.db in EGGNOG_DATABASES:
        dbpath, port = get_db_info(args.db)
        db_present = [pexists(dbpath + "." + ext)
                      for ext in 'h3f h3i h3m h3p idmap'.split()]
        if False in db_present:
            print db_present
            print colorify('Database %s not present. Use download_eggnog_database.py to fetch it' % (args.db), 'red')
            raise ValueError('Database not found')

        if not args.no_annot:
            if not pexists(pjoin(DATA_PATH, 'eggnog.db')):
                print colorify('Database eggnog.db not present. Use download_eggnog_database.py to fetch it', 'red')
                raise ValueError('Database not found')
        if not args.no_refine:
            if not pexists(pjoin(DATA_PATH, 'OG_fasta')):
                print colorify('Database OG_fasta not present. Use download_eggnog_database.py to fetch it', 'red')
                raise ValueError('Database not found')

        if scantype == 'mem':
            idmap_file = dbpath + '.idmap'
            end_port = port + 200

    elif os.path.isfile(args.db + '.h3f'):
        dbpath = args.db
        if scantype == 'mem':
            idmap_file = args.db + ".idmap"
            if not pexists(idmap_file):
                if generate_idmap(args.db):
                    idmap_file = args.db + ".idmap"
                    print >>sys.stderr, "idmap succesfully created!"
                else:
                    raise ValueError("idmap could not be created!")
            port = 53000
            end_port = 53200
        else:
            idmap_file = None
            port = None

    elif ":" in args.db:
        dbname, host, port = map(str.strip, args.db.split(":"))
        scantype = 'mem'
        port = int(port)
        if dbname in EGGNOG_DATABASES:
            dbfile, port = get_db_info(dbname)
            args.db = dbname
        else:
            dbfile = dbname

        idmap_file = dbfile + '.idmap'
        if not pexists(idmap_file):
            raise ValueError("idmap file not found: %s" % idmap_file)
        dbpath = host

        if not server_functional(host, port, args.dbtype):
            print colorify("eggnog-mapper server not found at %s:%s" %
                           (host, port), 'red')
            exit(-1)
        connecting_to_server = True
    else:
        ValueError('Invalid database name/server')

    # If memory based searches requested, load server if necessary
    if scantype == "mem" and not connecting_to_server and not args.annotate_only:
        master_db, worker_db = None, None
        for try_port in range(port, end_port, 2):
            print colorify("Loading server at localhost, port %s-%s" %
                           (try_port, try_port + 1), 'lblue')
            dbpath, master_db, worker_db = load_server(
                dbpath, try_port, try_port + 1, args.cpu)
            port = try_port
            ready = False
            for _ in xrange(TIMEOUT_LOAD_SERVER):
                print "Waiting for server to become ready...", host, try_port
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
                dbpath = host
                break
    elif scantype == "mem":
        print colorify("DB Server already running or not needed!", 'yellow')
        dbpath = host

    if scantype == "mem" and not args.annotate_only:
        print colorify("Reading idmap %s" % idmap_file, color='lblue')
        idmap = {}
        for _lnum, _line in enumerate(open(idmap_file)):
            if not _line.strip():
                continue
            try:
                _seqid, _seqname = map(str, _line.strip().split(' '))
            except ValueError:
                if _lnum == 0:
                    continue  # idmap generated by esl_reformat has an info line at beginning
                else:
                    raise
            _seqid = int(_seqid)
            idmap[_seqid] = [_seqname]
        print >>sys.stderr, len(idmap), "names loaded"

    # If server mode, just listen for connections
    if args.servermode:
        while True:
            print colorify("Server ready listening at %s:%s and using %d CPU cores" % (host, port, args.cpu), 'green')
            print colorify("Use `emapper.py -d %s:%s:%s (...)` to search against this server" % (args.db, host, port), 'lblue')
            time.sleep(10)
        sys.exit(0)

    hits_file = "%s.hits" % args.output
    refine_file = "%s.refined_hits" % args.output
    hits_annot_file = "%s.hits.annot" % args.output
    annot_file = "%s.refined_hits.annot" % args.output
    os.chdir(args.output_dir)

    if pexists(hits_file) and not args.resume and not args.override:
        print "Output files already present. User --resume or --override to continue"
        sys.exit(1)

    if args.resume and args.scratch_dir:
        for f in [hits_file, refine_file, hits_annot_file, annot_file]:
            if pexists(f):
                print "   Copying input file %s to scratch dir %s" % (f, args.scratch_dir)
                shutil.copy(f, args.scratch_dir)

    if args.scratch_dir:
        os.chdir(args.scratch_dir)

    hits_header = map(
        str.strip, "#query_name, hit, evalue, sum_score, query_length, hmmfrom, hmmto, seqfrom, seqto, query_coverage".split(','))
    refine_header = map(
        str.strip, "#query_name, best_hit_eggNOG_ortholog, best_hit_evalue, best_hit_score".split(','))
    hits_annot_header = map(
        str.strip, "#query_name, hit, level, evalue, sum_score, query_length, hmmfrom, hmmto, seqfrom, seqto, query_coverage, members_in_og, og_description, og_COG_categories".split(','))
    annot_header = ("query_name", "best_hit_eggNOG_ortholog", "best_hit_evalue",
                    "best_hit_score", "predicted_name", "strict_orthologs", "GO",
                    "KEGG(pathway)")

    # Start scanning sequences
    if not args.annotate_only:
        VISITED = set()
        if args.resume and pexists(hits_file):
            print colorify("Resuming previous run. Reading computed output from %s" % hits_file, 'yellow')
            VISITED = set([line.split('\t')[0].strip()
                           for line in open(hits_file) if not line.startswith('#')])
            print len(VISITED), 'queries skipped'
            OUT = open(hits_file, 'a')
        else:
            OUT = open(hits_file, 'w')

        print colorify("Sequence mapping starts now!", 'green')

        print >>OUT, '# ' + time.ctime()
        print >>OUT, '# ' + ' '.join(sys.argv)
        print >>OUT, '# ' + '\t'.join(hits_header)
        total_time = 0
        last_time = time.time()
        start_time = time.time()
        qn = 0

        for qn, (name, elapsed, hits, querylen, seq) in enumerate(
            search.iter_hits(
                args.input,
                args.translate,
                args.qtype,
                args.dbtype,
                scantype,
                dbpath,
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
                print >>OUT, '\t'.join(
                    [name] + ['ERROR'] * (len(hits_header) - 1))
            elif not hits:
                print >>OUT, '\t'.join([name] + ['-'] * (len(hits_header) - 1))
            else:
                for hitindex, (hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore) in enumerate(hits):
                    hitname = hid
                    if idmap:
                        hitname = idmap[hid][0]

                    print >>OUT, '\t'.join(map(str, [name, hitname, heval, hscore, int(querylen),
                                                     int(hmmfrom), int(hmmto), int(sqfrom), int(sqto), float(sqto - sqfrom) / querylen]))
            OUT.flush()

            # monitoring
            total_time += time.time() - last_time
            last_time = time.time()
            if qn and (qn % 25 == 0):
                print >>sys.stderr, qn + \
                    1, total_time, "%0.2f q/s" % ((float(qn + 1) / total_time))
                sys.stderr.flush()

        # finish
        ellapsed_time = time.time() - start_time
        print colorify("processed queries:%s total_time:%s rate:%s" %
                       (qn + 1, total_time,
                        "%0.2f q/s" % ((float(qn + 1) / ellapsed_time))), 'lblue')

        print >>OUT, '# %d queries scanned' % (qn + 1)
        print >>OUT, '# Total time (seconds):', ellapsed_time
        print >>OUT, '# Rate:', "%0.2f q/s" % ((float(qn + 1) / ellapsed_time))
        OUT.close()
        if args.scratch_dir:
            print "   Copying result file %s from scratch to %s" % (hits_file, args.output_dir)
            shutil.copy(hits_file, args.output_dir)

    if args.db in EGGNOG_DATABASES:
        if not args.no_refine:
            print colorify("Hit refinement starts now", 'green')
            start_time = time.time()
            og2level = dict([tuple(map(str.strip, line.split('\t')))
                             for line in gopen(OGLEVELS_FILE)])
            OUT = open(refine_file, "w")
            print >>OUT, '\t'.join(refine_header)
            qn = 0
            for qn, r in enumerate(process_hits_file(hits_file, args.input, og2level, translate=args.translate, cpu=args.cpu)):
                if qn and (qn % 25 == 0):
                    total_time = time.time() - start_time
                    print >>sys.stderr, qn + \
                        1, total_time, "%0.2f q/s (refinement)" % (
                            (float(qn + 1) / total_time))
                    sys.stderr.flush()
                query_name = r[0]
                best_hit_name = r[1]
                best_hit_evalue = r[2]
                best_hit_score = r[3]
                if best_hit_name != '-' and float(best_hit_score) >= 20:
                    print >>OUT, '\t'.join(
                        map(str, (query_name, best_hit_name, best_hit_evalue, best_hit_score)))
                    OUT.flush()
            ellapsed_time = time.time() - start_time
            print >>OUT, '# %d queries scanned' % (qn + 1)
            print >>OUT, '# Total time (seconds):', ellapsed_time
            print >>OUT, '# Rate:', "%0.2f q/s" % (
                (float(qn + 1) / ellapsed_time))
            OUT.close()
            print colorify("processed queries:%s total_time:%s rate:%s" % (qn, total_time, "%0.2f q/s" % ((float(qn) / ellapsed_time))), 'lblue')
            if args.scratch_dir:
                print "   Copying result file %s from scratch to %s" % (refine_file, args.output_dir)
                shutil.copy(refine_file, args.output_dir)

        if not args.no_annot:
            annota.connect()
            print colorify("Functional annotation of hits starts now", 'green')
            start_time = time.time()
            if pexists(hits_file):
                OUT = open(hits_annot_file, "w")
                print >>OUT, '\t'.join(hits_annot_header)
                qn = 0
                t1 = time.time()
                for line in open(hits_file):
                    if not line.strip() or line.startswith('#'):
                        continue
                    if qn and (qn % 10000 == 0):
                        total_time = time.time() - start_time
                        print >>sys.stderr, qn+1, total_time, "%0.2f q/s (refinement)" % (
                            (float(qn + 1) / total_time))
                        sys.stderr.flush()
                    qn += 1
                    (query, hit, evalue, sum_score, query_length, hmmfrom, hmmto,
                     seqfrom, seqto, q_coverage) = map(str.strip, line.split('\t'))
                    if hit not in ['ERROR', '-']:
                        hitname = cleanup_og_name(hit)
                        level, nm, desc, cats = annota.get_og_annotations(hitname)
                        print >>OUT, '\t'.join(map(
                            str, [query, hitname, level, evalue, sum_score, query_length, hmmfrom, hmmto, seqfrom, seqto, q_coverage, nm, desc, cats]))
                    else:
                        print >>OUT, '\t'.join(
                            [query] + [hit] * (len(hits_annot_header) - 1))
                elapsed_time = time.time() - t1
                print >>OUT, '# %d queries scanned' % (qn + 1)
                print >>OUT, '# Total time (seconds):', elapsed_time
                print >>OUT, '# Rate:', "%0.2f q/s" % (
                    (float(qn + 1) / elapsed_time))
                OUT.close()
                if args.scratch_dir:
                    print "   Copying result file %s from scratch to %s" % (hits_annot_file, args.output_dir)
                    shutil.copy(hits_annot_file, args.output_dir)

            if args.db != 'viruses' and pexists(refine_file):
                annotate_refined_hits_sequential(refine_file, annot_file, annot_header)
                #annotate_refined_hits(refine_file, annot_file+'_new', args.orthotype)

                if args.scratch_dir:
                    print "   Copying result file %s from scratch to %s" % (annot_file, args.output_dir)
                    shutil.copy(annot_file, args.output_dir)

    if args.scratch_dir:
        for f in [hits_file, refine_file, hits_annot_file, annot_file]:
            print " Cleaning", f
            silent_rm(f)

    # for p in multiprocessing.active_children():
    #    p.terminate()
    #    p.join()
    print colorify('Done', 'green')
    for f in [hits_file, refine_file, hits_annot_file, annot_file]:
        colorify('Result files:', 'yellow')
        if pexists(f):
            print "   %s" % (f)

    shutdown_server()


def annotate_refined_hits_sequential(refine_file, annot_file, annot_header):
    start_time = time.time()
    print colorify("Functional annotation of refined hits starts now", 'green')
    OUT = open(annot_file, "w")
    print >>OUT, '\t'.join(annot_header)
    qn = 0
    for line in open(refine_file):
        if not line.strip() or line.startswith('#'):
            continue
        if qn and (qn % 25 == 0):
            total_time = time.time() - start_time
            print >>sys.stderr, qn+1, total_time, "%0.2f q/s (refinement)" % (
                (float(qn + 1) / total_time))
            sys.stderr.flush()
        qn += 1
        r = map(str.strip, line.split('\t'))
        query_name = r[0]
        best_hit_name = r[1]
        best_hit_evalue = r[2]
        best_hit_score = r[3]
        if best_hit_name != '-' and float(best_hit_score) >= 20:
            #_orthologs = sorted(annota_mongo.refine_orthologs_by_member([best_hit_name])['one2one'])
            orthologs = sorted(annota.get_member_orthologs(
                best_hit_name)[args.orthotype])
            if orthologs:
                pname, gos, keggs = annota.get_member_annotations(
                    orthologs, excluded_gos=set(["IEA", "ND"]))
                #_pname = Counter(annota_mongo.get_preferred_names_dict(orthologs).values())
                name_ranking = sorted(
                    pname.items(), key=lambda x: x[1], reverse=True)
            else:
                name_ranking = [[0, 0, 0]]
                gos, keggs = '', ''

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
            OUT.flush()
    print >>OUT, '# Total time (seconds):', time.time() - start_time
    OUT.close()


def annotate_refined_hits(refine_file, annot_file, orthotype):
    start_time = time.time()
    annot_header = ("query_name", "best_hit_eggNOG_ortholog", "best_hit_evalue",
                    "best_hit_score", "predicted_name", "strict_orthologs", "GO",
                    "KEGG(pathway)2")
    print colorify("Functional annotation of refined hits starts now", 'green')
    # Preloads all raw best entries
    entries = []
    for line in open(refine_file):
        if not line.strip() or line.startswith('#'):
            continue
        fields = map(str.strip, line.split('\t'))
        query = fields[0]
        best_hit_name = fields[1]
        best_hit_evalue = fields[2]
        best_hit_score = fields[3]
        if best_hit_name != '-' and float(best_hit_score) >= 20:
            entries.append(fields)

    # Get annotations for all unique best hit names
    all_best_hits = set([e[1] for e in entries])
    annotations = annota.get_annotated_orthologs(all_best_hits,
                                                 orthotype,
                                                 set(["IEA", "ND"]),
                                                 args.cpu)

    # augment each entry with the annotations available
    OUT = open(annot_file, "w")
    print >>OUT, '\t'.join(annot_header)
    for entry in entries:
        (query_name, best_hit_name, best_hit_evalue, best_hit_score) = entry
        orthologs, pname, gos, keggs = annotations[best_hit_name]
        best_name = '-'
        if pname:
            name_ranking = sorted(pname.items(), key=lambda x: x[1], reverse=True)
            if name_ranking[0][1] > 2:
                best_name = name_ranking[0][0]

        print >>OUT, '\t'.join(map(str, (query_name, best_hit_name,
                                         best_hit_evalue, best_hit_score,
                                         best_name,
                                         ','.join(sorted(orthologs)),
                                         ','.join(sorted(gos)),
                                         ','.join(sorted(keggs))
                                         )))
        # print '\t'.join(map(str, (query_name, best_hit_name,
        #                          best_hit_evalue, best_hit_score,
        #                          best_name,
        #                          ','.join(sorted(orthologs)),
        #                          ','.join(sorted(gos)),
        #                          ','.join(sorted(keggs))
        #                          )))

        #OUT.flush()

    print >>OUT, '# Total time (seconds):', time.time() - start_time
    OUT.close()



def process_hits_file(hits_file, query_fasta, og2level, skip_queries=None, translate=False, cpu=1):
    sequences = {name: seq for name, seq in seqio.iter_fasta_seqs(
        query_fasta, translate=translate)}
    cmds = []
    visited_queries = set()

    if skip_queries:
        visited_queries.update(skip_queries)

    for line in gopen(hits_file):
        if line.startswith('#'):
            continue

        fields = map(str.strip, line.split('\t'))
        seqname = fields[0]

        if fields[1] == '-' or fields[1] == 'ERROR':
            continue

        if seqname in visited_queries:
            continue

        hitname = cleanup_og_name(fields[1])
        level = og2level[hitname]

        seq = sequences[seqname]
        visited_queries.add(seqname)
        target_fasta = os.path.join(FASTA_PATH, level, "%s.fa" % hitname)
        cmds.append([seqname, seq, target_fasta])

    if cmds:
        pool = multiprocessing.Pool(cpu)
        for r in pool.imap(search.refine_hit, cmds):
            yield r
        pool.terminate()
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # server
    g1 = parser.add_argument_group('Target database options')
    g1.add_argument(
        '--guessdb', help='guess eggnog db based on the provided taxid', type=int)

    g1.add_argument('-d', dest='db', help=('specify the target database for sequence searches'
                                           '. Choose among: euk,bact,arch, host:port, or a local hmmpressed database'))
    g1.add_argument('--dbtype', dest="dbtype",
                    choices=["hmmdb", "seqdb"], default="hmmdb")
    g1.add_argument('--qtype',  choices=["hmm", "seq"], default="seq")

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

    # g2.add_argument('--othologs_evalue', dest='ortho_evalue', default=0.001, type=float,
    #                help="E-value threshold for ortholog dectection. Default=0.001")
    # g2.add_argument('--orthologs_score', dest='ortho_score', default=60, type=float,
    # help="Bit score threshold for ortholog detection. Default=60")

    g3 = parser.add_argument_group('Output options')
    g3.add_argument('--output', '-o', type=str,
                    help="base name for output files")
    g3.add_argument('--resume', action="store_true",
                    help="Resumes a previous execution skipping reported hits in the output file.")
    g3.add_argument('--override', action="store_true",
                    help="Overwrites output files if they exist.")
    g3.add_argument("--no_refine", action="store_true",
                    help="Skip hit refinement, reporting only HMM hits.")
    g3.add_argument("--no_annot", action="store_true",
                    help="Skip functional annotation, reporting only hits")
    g3.add_argument("--annotate_only", action="store_true",
                    help="Skip mapping. Use existing hits file")

    g3.add_argument("--scratch_dir", type=str,
                    help="Write output files in a temporary scratch dir, move them to final the final output dir when finished. Speed up large computations using network file systems.")

    g3.add_argument("--output_dir", default=os.getcwd(), type=str,
                    help="Where output files should be written")

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
    if args.input:
        args.input = os.path.realpath(args.input)

    if args.servermode and args.input:
        parser.error(
            'Incompatible execution modes. Choose between -i or --servermode')
    if not (args.servermode or args.input):
        parser.error(
            'Execution must be specified. Choose between: -i [fastafile]  or --servermode')
    if not args.db:
        if args.guessdb:
            from ete3 import NCBITaxa
            ncbi = NCBITaxa()
            lineage = ncbi.get_lineage(args.guessdb)
            for tid in reversed(lineage):
                if tid in TAXID2LEVEL:
                    print tid, TAXID2LEVEL[tid]
                    args.db = TAXID2LEVEL[tid]
                    break
        if not args.db:
            parser.error('A target databse must be specified with --db')

    if args.servermode:
        args.usemem = True

    if not args.output and not args.servermode:
        parser.error(
            "a base name for output files has to be provided with --output")

    if args.scratch_dir and not os.path.isdir(args.scratch_dir):
        parser.error("--scratch_dir should point to an existing directory")

    if args.output_dir and not os.path.isdir(args.output_dir):
        print args.output_dir
        parser.error("--output_dir should point to an existing directory")

    _total_time = time.time()
    try:
        main(args)
    except:
        #shutdown_server()
        raise
        sys.exit(1)
    print 'Total time: %g secs' % (time.time()-_total_time)
