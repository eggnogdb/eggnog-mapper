#!/usr/bin/env python2
import sys
import os
import errno
import time
import cPickle
import multiprocessing
import argparse
import re
import shutil
from collections import defaultdict, Counter
from tempfile import mkdtemp
import uuid
import shutil
import subprocess
import json

SCRIPT_PATH = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.insert(0, SCRIPT_PATH)

from eggnogmapper.common import *
from eggnogmapper.vars import LEVEL_PARENTS, LEVEL_NAMES, LEVEL_DEPTH

from eggnogmapper import search
from eggnogmapper import annota
from eggnogmapper import seqio
from eggnogmapper import orthology
from eggnogmapper.utils import colorify
from eggnogmapper.server import (server_functional, load_server,
                                 generate_idmap, shutdown_server)


__description__ = ('A program for bulk functional annotation of novel '
                    'sequences using EggNOG database orthology assignments')
__author__ = 'Jaime Huerta Cepas'
__license__ = "GPL v2"

class emapperException(Exception):
    def __init__(self, *args, **kargs):
        sys.excepthook = lambda exctype,exc,traceback: ""
        super(emapperException, self).__init__(*args, **kargs)


def cleanup_og_name(name):
    # names in the hmm databases are sometiemes not clean eggnog OG names
    m = re.search('\w+\.((ENOG41|COG|KOG|arCOG)\w+)\.', name)
    if m:
        name = m.groups()[0]
    name = re.sub("^ENOG41", "", name)
    return name

def setup_hmm_search(args):
    host = 'localhost'
    idmap = None
    if args.usemem:
        scantype = 'mem'
    else:
       scantype = 'disk'

    connecting_to_server = False
    # If searching against a predefined database name
    if args.db in EGGNOG_DATABASES:
        dbpath, port = get_db_info(args.db)
        print dbpath
        db_present = [pexists(dbpath + "." + ext)
                      for ext in 'h3f h3i h3m h3p idmap'.split()]

        if False in db_present:
            print db_present
            print colorify('Database %s not present. Use download_eggnog_database.py to fetch it' % (args.db), 'red')
            raise ValueError('Database not found')

        if not args.no_refine:
            if not pexists(pjoin(get_data_path(), 'OG_fasta')):
                print colorify('Database data/OG_fasta/ not present. Use download_eggnog_database.py to fetch it', 'red')
                raise ValueError('Database not found')

        if scantype == 'mem':
            idmap_file = dbpath + '.idmap'
            end_port = 53200

    # If searching against a custom hmm database
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

    # If searching against a emapper hmm server
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
            print colorify("eggnog-mapper server not found at %s:%s" % (host, port), 'red')
            exit(1)
        connecting_to_server = True
    else:
        raise ValueError('Invalid database name/server')


    # If memory based searches requested, start server
    if scantype == "mem" and not connecting_to_server:
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

    # Preload seqid map to translate hits from hmmpgmd
    if scantype == "mem":
        print colorify("Reading idmap %s" % idmap_file, color='lblue')
        idmap = {}
        for _lnum, _line in enumerate(open(idmap_file)):
            if not _line.strip():
                continue
            try:
                _seqid, _seqname = map(str, _line.strip().split(' '))
            except ValueError:
                if _lnum == 0:
                    # idmap generated by esl_reformat has info line at beginning
                    continue  
                else:
                    raise
            _seqid = int(_seqid)
            idmap[_seqid] = [_seqname]
        print len(idmap), "names loaded"

    # If server mode, just listen for connections and exit when interrupted
    if args.servermode:
        while True:
            print colorify("Server ready listening at %s:%s and using %d CPU cores" % (host, port, args.cpu), 'green')
            print colorify("Use `emapper.py -d %s:%s:%s (...)` to search against this server" % (args.db, host, port), 'lblue')
            time.sleep(10)
        raise emapperException()
    else:
        return host, port, dbpath, scantype, idmap

def main(args):
    # Output and intermediate files
    hmm_hits_file = "%s.emapper.hmm_hits" % args.output
    seed_orthologs_file = "%s.emapper.seed_orthologs" % args.output
    annot_file = "%s.emapper.annotations" % args.output
    orthologs_file = "%s.emapper.predict_orthologs" %args.output

    if args.no_search:
        output_files = [annot_file]
    elif args.no_annot:
        output_files = [hmm_hits_file, seed_orthologs_file]
    else:
        output_files = [hmm_hits_file, seed_orthologs_file, annot_file]

    # convert to absolute path before changing directory
    if args.annotate_hits_table:
        args.annotate_hits_table = os.path.abspath(args.annotate_hits_table)
    # force user to decide what to do with existing files
    os.chdir(args.output_dir)
    files_present = set([pexists(fname) for fname in output_files])
    if True in files_present and not args.resume and not args.override:
        print "Output files detected in disk. Use --resume or --override to continue"
        raise emapperException()

    if args.override:
        for outf in output_files:
            silent_rm(outf)

    print '# ', get_version()
    print '# ./emapper.py ', ' '.join(sys.argv[1:])

    if args.scratch_dir:
        # If resuming in and using --scratch_dir, transfer existing files.
        if args.resume and args.scratch_dir:
            for f in output_files:
                if pexists(f):
                    print "   Copying input file %s to scratch dir %s" % (f, args.scratch_dir)
                    shutil.copy(f, args.scratch_dir)

        # Change working dir
        os.chdir(args.scratch_dir)

    # Step 1. Sequence search
    if not args.no_search:
        if args.mode == 'diamond' and not args.no_search:
            dump_diamond_matches(args.input, seed_orthologs_file, args)

        elif args.mode == 'hmmer' and not args.no_search:
            host, port, dbpath, scantype, idmap = setup_hmm_search(args)
            # Start HMM SCANNING sequences (if requested)
            if not pexists(hmm_hits_file) or args.override:
                dump_hmm_matches(args.input, hmm_hits_file, dbpath, port, scantype, idmap, args)

            if not args.no_refine and (not pexists(seed_orthologs_file) or args.override):
                if args.db == 'viruses':
                    print 'Skipping seed ortholog detection in "viruses" database'
                elif args.db in EGGNOG_DATABASES:
                    refine_matches(args.input, seed_orthologs_file, hmm_hits_file, args)
                else:
                    print 'refined hits not available for custom hmm databases.'

    # Step 2. Annotation
    if not args.no_annot:
        annota.connect()
        if args.annotate_hits_table:
            if not os.path.exists(args.annotate_hits_table):
                raise IOError(errno.ENOENT,
                              os.strerror(errno.ENOENT),
                              args.annotate_hits_table)
            annotate_hits_file(args.annotate_hits_table, annot_file, hmm_hits_file, args)
        elif args.db == 'viruses':
            annotate_hmm_matches(hmm_hits_file, hmm_hits_file+'.annotations', args)
            OUT = open(annot_file, 'w')
            for line in open(hmm_hits_file+'.annotations'):
                if line.startswith('#') or not line.strip():
                    continue
                (query, hitname, level, evalue, sum_score, query_length,
                 hmmfrom, hmmto, seqfrom, seqto, q_coverage, nm, desc, cats) = line.split("\t")

                if hitname != '-' and hitname != 'ERROR':
                    print >>OUT, '\t'.join(map(str, (query,
                                                     hitname,
                                                     evalue,
                                                     sum_score,
                                                     '',
                                                     '',
                                                     '',
                                                     'viruses',
                                                     hitname+"@viruses",
                                                     "%s|%s|%s" %(hitname, evalue, sum_score),
                                                     cats.replace('\n', ''),
                                                     desc.replace('\n', ' '))))
            OUT.close()
        else:
            annotate_hits_file(seed_orthologs_file, annot_file, hmm_hits_file, args)

    if args.predict_ortho:
        orthology.connect()
        dump_orthologs(seed_orthologs_file, orthologs_file, args)
            
    # If running in scratch, move files to real output dir and clean up

    if args.scratch_dir:
        for fname in output_files:
            if pexists(fname):
                print " Copying result file %s from scratch to %s" % (fname, args.output_dir)
                shutil.copy(annot_file, args.output_dir)
                print "  Cleaning result file %s from scratch dir" %(fname)

    # Finalize and exit
    print colorify('Done', 'green')
    for f in output_files:
        colorify('Result files:', 'yellow')
        if pexists(f):
            print "   %s" % (f)

    print 'Total time: %g secs' % (time.time()-_total_time)

    if args.mode == 'hmmer':
        print get_citation(['hmmer'])
    elif args.mode == 'diamond':
        print get_citation(['diamond'])

    shutdown_server()

def dump_diamond_matches(fasta_file, seed_orthologs_file, args):
    cpu = args.cpu
    score_thr = args.seed_ortholog_score
    evalue_thr = args.seed_ortholog_evalue
    excluded_taxa = args.excluded_taxa if args.excluded_taxa else None
    
    if args.translate:
        tool = 'blastx'
    else:
        tool = 'blastp'
    
    dmnd_db = args.dmnd_db if args.dmnd_db else get_eggnog_dmnd_db()
    query_cov = args.query_cover
    subject_cov = args.subject_cover
    dmnd_opts = ''
    
    if args.matrix is not None:
        dmnd_opts += ' --matrix %s' % args.matrix
    if args.gapopen is not None:
        dmnd_opts += ' --gapopen %d' % args.gapopen
    if args.gapextend is not None:
        dmnd_opts += ' --gapextend %d' % args.gapextend

    if not DIAMOND:
        raise ValueError("diamond not found in path")

    tempdir = mkdtemp(prefix='emappertmp_dmdn_', dir=args.temp_dir)

    raw_output_file = pjoin(tempdir, uuid.uuid4().hex)
    
    if excluded_taxa:
        cmd = '%s %s -d %s -q %s --more-sensitive --threads %s -e %f -o %s --max-target-seqs 25 --query-cover %s --subject-cover %s' %\
          (DIAMOND, tool, dmnd_db, fasta_file, cpu, evalue_thr, raw_output_file, query_cov, subject_cov)
    else:
        cmd = '%s %s -d %s -q %s --more-sensitive --threads %s -e %f -o %s --top 3 --query-cover %s --subject-cover %s' %\
          (DIAMOND, tool, dmnd_db, fasta_file, cpu, evalue_thr, raw_output_file, query_cov, subject_cov)


    print colorify('  '+cmd, 'yellow')

    try:
        with open(raw_output_file+'.stdout', 'w') as STDOUT:
            subprocess.check_call(cmd, shell=True, stdout=STDOUT)
            
        OUT = open('%s' %seed_orthologs_file, 'w')

        if not args.no_file_comments:
            print >>OUT, get_call_info()
            print >>OUT, '#', cmd

        visited = set()
        for line in open(raw_output_file):
            if not line.strip() or line.startswith('#'):
                continue
            fields = map(str.strip, line.split('\t'))
            query = fields[0]
            hit = fields[1]
            evalue = float(fields[10])
            score = float(fields[11])

            if query in visited:
                continue

            if evalue > evalue_thr or score < score_thr:
                continue

            if excluded_taxa and hit.startswith("%s." % excluded_taxa):
                continue

            visited.add(query)
            print >>OUT, '\t'.join(map(str, [query, hit, evalue, score]))
        OUT.close()

    except subprocess.CalledProcessError as e:
        raise e
    finally:
        shutil.rmtree(tempdir)


def dump_hmm_matches(fasta_file, hits_file, dbpath, port, scantype, idmap, args):
    hits_header = ("#query_name", "hit", "evalue", "sum_score", "query_length",
                   "hmmfrom", "hmmto", "seqfrom", "seqto", "query_coverage")

    # Cache previous results if resuming is enabled
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
    if not args.no_file_comments:
        print >>OUT, get_call_info()
        print >>OUT, '# ' + '\t'.join(hits_header)
    total_time = 0
    last_time = time.time()
    start_time = time.time()
    qn = 0 # in case nothing to loop bellow
    for qn, (name, elapsed, hits, querylen, seq) in enumerate(search.iter_hits(
                                                        fasta_file,
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
                                                        cpus=args.cpu,
                                                        base_tempdir=args.temp_dir)):

        if elapsed == -1:
            # error occurred
            print >>OUT, '\t'.join(
                [name] + ['ERROR'] * (len(hits_header) - 1))
        elif not hits:
            print >>OUT, '\t'.join([name] + ['-'] * (len(hits_header) - 1))
        else:
            for hitindex, (hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore) in enumerate(hits):
                hitname = hid
                if idmap:
                    hitname = idmap[hid][0]

                print >>OUT, '\t'.join(map(str, [name, hitname, heval, hscore,
                                                 int(querylen), int(hmmfrom),
                                                 int(hmmto), int(sqfrom),
                                                 int(sqto),
                                                 float(sqto - sqfrom) / querylen]))
        OUT.flush()

        # monitoring
        total_time += time.time() - last_time
        last_time = time.time()
        if qn and (qn % 25 == 0):
            print >>sys.stderr, qn + \
                1, total_time, "%0.2f q/s" % ((float(qn + 1) / total_time))
            sys.stderr.flush()

    # Writes final stats
    elapsed_time = time.time() - start_time
    if not args.no_file_comments:
        print >>OUT, '# %d queries scanned' % (qn + 1)
        print >>OUT, '# Total time (seconds):', elapsed_time
        print >>OUT, '# Rate:', "%0.2f q/s" % ((float(qn + 1) / elapsed_time))
    OUT.close()
    print colorify(" Processed queries:%s total_time:%s rate:%s" %\
                   (qn+1, elapsed_time, "%0.2f q/s" % ((float(qn+1) / elapsed_time))), 'lblue')


def annotate_hmm_matches(hits_file, hits_annot_file, args):
    hits_annot_header = map(str.strip, '''#query_name, hit, level, evalue,
                         sum_score, query_length, hmmfrom, hmmto, seqfrom, seqto, query_coverage,
                         members_in_og, og_description, og_COG_categories'''.split(','))

    annota.connect()
    print colorify("Functional annotation of hits starts now", 'green')
    start_time = time.time()
    if pexists(hits_file):
        OUT = open(hits_annot_file, "w")
        if not args.no_file_comments:
            print >>OUT, get_call_info()
            print >>OUT, '\t'.join(hits_annot_header)
        qn = 0
        t1 = time.time()
        for line in open(hits_file):
            if not line.strip() or line.startswith('#'):
                continue
            qn += 1
            if qn and (qn % 10000 == 0):
                total_time = time.time() - start_time
                print >>sys.stderr, qn, total_time, "%0.2f q/s (refinement)" %\
                    ((float(qn) / total_time))
                sys.stderr.flush()

            (query, hit, evalue, sum_score, query_length, hmmfrom, hmmto,
             seqfrom, seqto, q_coverage) = map(str.strip, line.split('\t'))
            if hit not in ['ERROR', '-']:
                hitname = cleanup_og_name(hit)
                level, nm, desc, cats = annota.get_og_annotations(hitname)
                print >>OUT, '\t'.join(map( str, [query, hitname, level, evalue,
                                                  sum_score, query_length,
                                                  hmmfrom, hmmto, seqfrom,
                                                  seqto, q_coverage, nm, desc,
                                                  cats]))
            else:
                print >>OUT, '\t'.join(
                    [query] + [hit] * (len(hits_annot_header) - 1))
        elapsed_time = time.time() - t1
        if not args.no_file_comments:
            print >>OUT, '# %d queries scanned' % (qn)
            print >>OUT, '# Total time (seconds):', elapsed_time
            print >>OUT, '# Rate:', "%0.2f q/s" % ((float(qn) / elapsed_time))
        OUT.close()
        print colorify(" Processed queries:%s total_time:%s rate:%s" %\
                       (qn, elapsed_time, "%0.2f q/s" % ((float(qn) / elapsed_time))), 'lblue')


def get_seq_hmm_matches(hits_file):
    annota.connect()
    print colorify("Reading HMM matches", 'green')
    seq2oginfo = {}
    start_time = time.time()
    hitnames = set()
    if pexists(hits_file):
        for line in open(hits_file):
            if not line.strip() or line.startswith('#'):
                continue

            (query, hit, evalue, sum_score, query_length, hmmfrom, hmmto,
             seqfrom, seqto, q_coverage) = map(str.strip, line.split('\t'))

            if query not in seq2oginfo and hit not in ['ERROR', '-']:
                hitname = cleanup_og_name(hit)
                seq2oginfo[query] = [hitname, evalue, sum_score, query_length,
                                     hmmfrom, hmmto, seqfrom, seqto,
                                     q_coverage]
    return seq2oginfo

def refine_matches(fasta_file, refine_file, hits_file, args):
    refine_header = map(str.strip, '''#query_name, best_hit_eggNOG_ortholog,
                        best_hit_evalue, best_hit_score'''.split(','))

    print colorify("Hit refinement starts now", 'green')
    start_time = time.time()
    og2level = dict([tuple(map(str.strip, line.split('\t')))
                     for line in gopen(get_oglevels_file())])
    OUT = open(refine_file, "w")

    if not args.no_file_comments:
        print >>OUT, get_call_info()
        print >>OUT, '\t'.join(refine_header)

    qn = 0 # in case no hits in loop bellow
    for qn, r in enumerate(process_nog_hits_file(hits_file, fasta_file, og2level,
                                                 translate=args.translate,
                                                 cpu=args.cpu,
                                                 excluded_taxa=args.excluded_taxa,
                                                 base_tempdir=args.temp_dir)):
        if qn and (qn % 25 == 0):
            total_time = time.time() - start_time
            print >>sys.stderr, qn + 1, total_time, "%0.2f q/s (refinement)" % ((float(qn + 1) / total_time))
            sys.stderr.flush()
        query_name = r[0]
        best_hit_name = r[1]
        if best_hit_name == '-' or best_hit_name == 'ERROR':
            continue
        best_hit_evalue = float(r[2])
        best_hit_score = float(r[3])
        print >>OUT, '\t'.join(map(str, (query_name, best_hit_name,
                                         best_hit_evalue, best_hit_score)))
        #OUT.flush()

    elapsed_time = time.time() - start_time
    if not args.no_file_comments:
        print >>OUT, '# %d queries scanned' % (qn + 1)
        print >>OUT, '# Total time (seconds):', elapsed_time
        print >>OUT, '# Rate:', "%0.2f q/s" % ((float(qn + 1) / elapsed_time))
    OUT.close()
    print colorify(" Processed queries:%s total_time:%s rate:%s" %\
                   (qn+1, elapsed_time, "%0.2f q/s" % ((float(qn+1) / elapsed_time))), 'lblue')


def process_nog_hits_file(hits_file, query_fasta, og2level, skip_queries=None,
                          translate=False, cpu=1, excluded_taxa=None, base_tempdir=None):
    sequences = {name: seq for name, seq in seqio.iter_fasta_seqs(
        query_fasta, translate=translate)}
    cmds = []
    visited_queries = set()

    if skip_queries:
        visited_queries.update(skip_queries)

    tempdir = mkdtemp(prefix='emappertmp_phmmer_', dir=base_tempdir)

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
        target_fasta = os.path.join(get_fasta_path(), level, "%s.fa" % hitname)
        cmds.append([seqname, seq, target_fasta, excluded_taxa, tempdir])

    if cmds:
        pool = multiprocessing.Pool(cpu)
        for r in pool.imap(search.refine_hit, cmds):
            yield r
        pool.terminate()

    shutil.rmtree(tempdir)

def annotate_hit_line(arguments):
    try:
        return _annotate_hit_line(arguments)
    except:
        import traceback
        traceback.print_exc(file=sys.stdout)
        raise

def _annotate_hit_line(arguments):
    annota.connect()
    line, args = arguments

    if not line.strip() or line.startswith('#'):
        return None
    r = map(str.strip, line.split('\t'))

    query_name = r[0]
    best_hit_name = r[1]
    if best_hit_name == '-' or best_hit_name == 'ERROR':
        return None

    best_hit_evalue = float(r[2])
    best_hit_score = float(r[3])
    if best_hit_score < args.seed_ortholog_score or best_hit_evalue > args.seed_ortholog_evalue:
        return None

    match_nogs = annota.get_member_ogs(best_hit_name)
    if not match_nogs:
        return None

    match_levels = set()
    for nog in match_nogs:
        match_levels.update(LEVEL_PARENTS[nog.split("@")[1]])

    swallowest_level = sorted(match_levels & set(LEVEL_DEPTH.keys()),
                              key=lambda x: LEVEL_DEPTH[x], reverse=True)[0]

    annot_levels = set()
    if args.tax_scope == "auto":
        for level in TAXONOMIC_RESOLUTION:
            if level in match_levels:
                annot_levels.add(level)
                annot_level_max = LEVEL_NAMES.get(level, level)
                break
    else:
        annot_levels.add(args.tax_scope)
        annot_level_max = LEVEL_NAMES.get(args.tax_scope, args.tax_scope)

    if args.target_taxa != 'all':
        target_taxa = orthology.normalize_target_taxa(args.target_taxa)
    else:
        target_taxa = None

    try:
        all_orthologies = annota.get_member_orthologs(best_hit_name, target_taxa=target_taxa, target_levels=annot_levels)
    except Exception:
        orthologs = None
        status = 'Error'
    else:
        orthologs = sorted(all_orthologies[args.target_orthologs])
        if args.excluded_taxa:
            orthologs = [o for o in orthologs if not o.startswith("%s." %args.excluded_taxa)]
        status = 'OK'
        
    if orthologs:
        annotations = annota.summarize_annotations(orthologs,
                                                   target_go_ev=args.go_evidence,
                                                   excluded_go_ev=args.go_excluded)
    else:
        annotations = {}

    return (query_name, best_hit_name, best_hit_evalue, best_hit_score,
            annotations, annot_level_max, swallowest_level, match_nogs, orthologs)


def iter_hit_lines(filename, args):
    for line in open(filename):
        if line.startswith('#') or not line.strip():
            continue
        yield (line, args)

def annotate_hits_file(seed_orthologs_file, annot_file, hmm_hits_file, args):
    HIT_HEADER = ["#query_name",
                          "seed_eggNOG_ortholog",
                          "seed_ortholog_evalue",
                          "seed_ortholog_score",
                          "best_tax_level", ]

    HIT_OG_HEADER = ["taxonomic scope", "eggNOG OGs", "best eggNOG OG",
                     "COG Functional cat.", "eggNOG free text desc."]

    start_time = time.time()
    seq2bestOG = {}
    if pexists(hmm_hits_file):
        seq2bestOG = get_seq_hmm_matches(hmm_hits_file)

    seq2annotOG = annota.get_ogs_annotations(set([v[0] for v in seq2bestOG.itervalues()]))

    print colorify("Functional annotation of refined hits starts now", 'green')

    OUT = open(annot_file, "w")

    if args.report_orthologs:
        ORTHOLOGS = open(annot_file+".orthologs", "w")

    if not args.no_file_comments:
        print >>OUT, '# emapper version:', get_version(), 'emapper DB:', get_db_version()
        print >>OUT, '# command: ./emapper.py ', ' '.join(sys.argv[1:])
        print >>OUT, '# time: ' + time.ctime()
        print >>OUT, '\t'.join(HIT_HEADER + ANNOTATIONS_HEADER + HIT_OG_HEADER)
    qn = 0

    pool = multiprocessing.Pool(args.cpu)

    for result in pool.imap(annotate_hit_line, iter_hit_lines(seed_orthologs_file, args)):
        qn += 1
        if qn and (qn % 500 == 0):
            total_time = time.time() - start_time
            print >>sys.stderr, qn, total_time, "%0.2f q/s (func. annotation)" % (
                (float(qn) / total_time))
            sys.stderr.flush()

        if result:
            (query_name, best_hit_name, best_hit_evalue, best_hit_score,
             annotations, annot_level_max, swallowest_level, match_nogs, orthologs) = result
            if query_name in seq2bestOG:
                (hitname, evalue, score, qlength, hmmfrom, hmmto, seqfrom,
                 seqto, q_coverage) = seq2bestOG[query_name]
                bestOG = '%s|%s|%s' %(hitname, evalue, score)
                og_cat, og_desc = seq2annotOG.get(hitname, ['', ''])
            else:
                bestOG = 'NA|NA|NA'
                og_cat, og_desc = annota.get_best_og_description(match_nogs)

            if args.report_orthologs:
                print >>ORTHOLOGS, '\t'.join(map(str, (query_name, ','.join(orthologs))))

            # prepare annotations for printing
            annot_columns = [query_name,
                             best_hit_name,
                             str(best_hit_evalue),
                             str(best_hit_score),
                             LEVEL_NAMES[swallowest_level]]

            for h in ANNOTATIONS_HEADER:
                if h in annotations:
                    annot_columns.append(','.join(sorted(annotations[h])))
                else:
                    annot_columns.append('')

            annot_columns.extend([annot_level_max,
                                    ','.join(match_nogs),
                                    bestOG,
                                    og_cat.replace('\n', ''),
                                    og_desc.replace('\n', ' ')])

            print >>OUT, '\t'.join(annot_columns)

        #OUT.flush()

    pool.terminate()

    elapsed_time = time.time() - start_time
    if not args.no_file_comments:
        print >>OUT, '# %d queries scanned' % (qn)
        print >>OUT, '# Total time (seconds):', elapsed_time
        print >>OUT, '# Rate:', "%0.2f q/s" % ((float(qn) / elapsed_time))
    OUT.close()

    if args.report_orthologs:
        ORTHOLOGS.close()

    print colorify(" Processed queries:%s total_time:%s rate:%s" %\
                   (qn, elapsed_time, "%0.2f q/s" % ((float(qn) / elapsed_time))), 'lblue')


def dump_orthologs(seed_orthologs_file, orthologs_file, args):
#Copy from predict_orthologs.py
    OUT = open(orthologs_file, "w")

    if args.predict_output_format == "per_query":
        ortholog_header = ("#Query", "Orthologs")
    elif args.predict_output_format == "per_species":
        ortholog_header = ("#Query", "Species", "Orthologs")

    print >> OUT, "\t".join(ortholog_header)

    if args.target_taxa != 'all':
        args._expanded_target_taxa = orthology.normalize_target_taxa(args.target_taxa)
    else:
        # report orthologs from any species by default
        args._expanded_target_taxa = None

    pool = multiprocessing.Pool(args.cpu)
    for result in pool.imap(find_orthologs_per_hit, iter_hit_lines(seed_orthologs_file, args)):
        if result:
            write_orthologs_in_file(result, OUT, args)

    pool.terminate()


def find_orthologs_per_hit(arguments):
#Copy from predict_orthologs.py

    orthology.connect()
    line, args = arguments

    if not line.strip() or line.startswith('#'):
        return None
    r = map(str.strip, line.split('\t'))

    query_name = r[0]
    best_hit_name = r[1]
    if best_hit_name == '-' or best_hit_name == 'ERROR':
        return None

    best_hit_evalue = float(r[2])
    best_hit_score = float(r[3])

    if best_hit_score < args.seed_ortholog_score or best_hit_evalue > args.seed_ortholog_evalue:
        return None

    target_taxa = args._expanded_target_taxa
    
    orthologs_pred = orthology.predict_orthologs_by_seed(best_hit_name, target_taxa=target_taxa, target_levels = None)
    return (query_name, best_hit_name, orthologs_pred)


def write_orthologs_in_file(result_line, ORTHOLOGS, args):
#Copy from predict_orthologs.py

    """
    Writes orthologs in file for all output formats except json
    """
    query_name, best_hit_name, orthologs_pred = result_line
    
    if args._expanded_target_taxa:  
        target_taxa = list(args._expanded_target_taxa)
    else:
        target_taxa = None
    
    if args.predict_output_format == "per_query":
        orthologs = []
        for key in orthologs_pred:
            if target_taxa is not None:
                    if key in target_taxa:
                        members = (','.join(orthologs_pred[key]))
                        orthologs.append(members)
                        print >> ORTHOLOGS, '\t'.join(map(str, (query_name, ','.join(orthologs))))

            else:
                members = (','.join(orthologs_pred[key]))
                orthologs.append(members)
                print >> ORTHOLOGS, '\t'.join(map(str, (query_name, ','.join(orthologs))))

    elif args.predict_output_format == "per_species":
        for key in orthologs_pred:
            sp_taxid = int(key)
            if target_taxa is not None: 
                if sp_taxid in target_taxa:
                    print >> ORTHOLOGS, '\t'.join(map(str, (query_name, key,
                                                    ','.join(orthologs_pred[key])))) 

            else:
                print >> ORTHOLOGS, '\t'.join(map(str, (query_name, key,
                                                    ','.join(orthologs_pred[key]))))
        '''
        sorted_orthologs = orthology.sort_orthologs_by_species(predict_ortho, best_hit_name)
        for (sp, _, ortho_type), ortho_list in sorted_orthologs.items():
            if ortho_type == 'all':
                print >>ORTHOLOGS, '\t'.join(map(str, [query_name, sp, ','.join(sorted(ortho_list))]))
'''
    ORTHOLOGS.flush()
    


'''               
    elif args.output_format == "per_species_and_type" :
        sorted_orthologs = orthology.sort_orthologs_by_species(all_orthologs, best_hit_name)
        seed_ortholog_sp = best_hit_name.split(".", 1)[0]
        for (sp, inparalogs, ortho_type), ortho_list in sorted_orthologs.items():
            if ortho_type == 'all':
                continue
            if sp == seed_ortholog_sp:
                if len(inparalogs) > 1:
                    ortho_type_temp = 'one2many'
                else:
                    ortho_type_temp = 'one2one'
                print >>ORTHOLOGS, '\t'.join(map(str, (query_name, sp, ortho_type_temp,
                                                       best_hit_name,
                                                       ','.join(inparalogs))))
            else:
                inparalogs = tuple(sorted(inparalogs - set([best_hit_name])))
                print >>ORTHOLOGS, '\t'.join(map(str, (query_name, sp, ortho_type,
                                                       ",".join(inparalogs),
                                                       ','.join(ortho_list))))


def build_json_format(result_line, json_dict):
    #Copy from predict_orthologs.py

    query_name, all_orthologs, best_hit_name = result_line
    json_dict[query_name] = {}

    sorted_orthologs = orthology.sort_orthologs_by_species(all_orthologs, best_hit_name)
    seed_ortholog_sp = best_hit_name.split(".", 1)[0]
    for (sp, inparalogs, ortho_type), ortho_list in sorted_orthologs.items():
        if sp not in json_dict[query_name].keys():
            json_dict[query_name][sp]= {}
        if ortho_type == 'all':
            continue
        if sp == seed_ortholog_sp:
            if len(inparalogs) > 1:
                ortho_type_temp = 'one2many'
            else:
                ortho_type_temp = 'one2one'

            json_dict[query_name][sp][ortho_type_temp]= [best_hit_name, list(inparalogs)]

        else:
            inparalogs = tuple(sorted(inparalogs - set([best_hit_name])))
            json_dict[query_name][sp][ortho_type]= [list(inparalogs), list(ortho_list)]

    return json_dict
'''


def parse_args(parser):
    args = parser.parse_args()

    if args.version:
        print get_version()
        sys.exit(0)

    if "EGGNOG_DATA_DIR" in os.environ:
        set_data_path(os.environ["EGGNOG_DATA_DIR"])

    if args.data_dir:
        set_data_path(args.data_dir)

    if not args.no_annot and not pexists(get_eggnogdb_file()):
        print colorify('Annotation database data/eggnog.db not present. Use download_eggnog_database.py to fetch it', 'red')
        raise emapperException()

    if args.mode == 'diamond':
        dmnd_db = args.dmnd_db if args.dmnd_db else get_eggnog_dmnd_db()
        if not pexists(dmnd_db):
            print colorify('DIAMOND database %s not present. Use download_eggnog_database.py to fetch it' % dmnd_db, 'red')
            raise emapperException()

    if args.cpu == 0:
        args.cpu = multiprocessing.cpu_count()

    # No --servermode available for diamond
    if args.mode == 'diamond' and args.servermode:
        parser.error('--mode [diamond] and --servermode are mutually exclusive')

    # Output file required unless running in servermode
    if not args.servermode and not args.output:
        parser.error('An output project name is required (-o)')

    # Servermode implies using mem-based databases
    if args.servermode:
        args.usemem = True

    # Direct annotation implies no searches
    if args.annotate_hits_table:
        args.no_search = True
        args.no_annot = False


    # Sets GO evidence bases
    if args.go_evidence == 'experimental':
        args.go_evidence = set(["EXP","IDA","IPI","IMP","IGI","IEP"])
        args.go_excluded = set(["ND", "IEA"])

    elif args.go_evidence == 'non-electronic':
        args.go_evidence = None
        args.go_excluded = set(["ND", "IEA"])
    else:
        raise ValueError('Invalid --go_evidence value')

    # Check inputs for running sequence searches
    if not args.no_search and not args.servermode:
        if not args.input:
            parser.error('An input fasta file is required (-i)')

        # HMM
        if args.mode == 'hmmer':
            if not args.db and not args.guessdb:
                parser.error('HMMER mode requires specifying a target database (i.e. -d, --guessdb ))')
            if args.db and args.guessdb:
                parser.error('-d and --guessdb options are mutually exclusive')

            if args.guessdb:
                from ete3 import NCBITaxa
                ncbi = NCBITaxa()
                lineage = ncbi.get_lineage(args.guessdb)
                for tid in reversed(lineage):
                    if tid in TAXID2LEVEL:
                        print tid, TAXID2LEVEL[tid]
                        args.db = TAXID2LEVEL[tid]
                        break
        # DIAMOND
        elif args.mode == 'diamond':
            #if args.db or args.guessdb:
            #    parser.error('diamond mode does not require -d or --guessdb options')
            pass

    return args


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # server
    pg_db = parser.add_argument_group('Target HMM Database Options')

    pg_db.add_argument('--guessdb', type=int, metavar='',
                       help='guess eggnog db based on the provided taxid')

    pg_db.add_argument('--database', '-d', dest='db', metavar='',
                       help=('specify the target database for sequence searches'
                             '. Choose among: euk,bact,arch, host:port, or a local hmmpressed database'))

    pg_db.add_argument('--dbtype', dest="dbtype",
                    choices=["hmmdb", "seqdb"], default="hmmdb")

    pg_db.add_argument("--data_dir", metavar='', type=existing_dir,
                    help='Directory to use for DATA_PATH.')

    pg_db.add_argument('--qtype',  choices=["hmm", "seq"], default="seq")


    pg_annot = parser.add_argument_group('Annotation Options')

    pg_annot.add_argument("--tax_scope", type=str, choices=LEVEL_NAMES.keys()+["auto"],
                    default='auto', metavar='',
                    help=("Fix the taxonomic scope used for annotation, so only orthologs from a "
                          "particular clade are used for functional transfer. "
                          "By default, this is automatically adjusted for every query sequence."))

    pg_annot.add_argument('--target_orthologs', choices=["one2one", "many2one",
                                                         "one2many","many2many", "all"],
                          default="all",
                          help='defines what type of orthologs should be used for functional transfer')

    pg_annot.add_argument('--excluded_taxa', type=int, metavar='',
                          help='(for debugging and benchmark purposes)')

    pg_annot.add_argument('--go_evidence', type=str, choices=('experimental', 'non-electronic'),
                          default='non-electronic',
                          help='Defines what type of GO terms should be used for annotation:'
                          'experimental = Use only terms inferred from experimental evidence'
                          'non-electronic = Use only non-electronically curated terms')

    pg_hmm = parser.add_argument_group('HMM search_options')

    pg_hmm.add_argument('--hmm_maxhits', dest='maxhits', type=int, default=1, metavar='',
                    help="Max number of hits to report. Default=1")

    pg_hmm.add_argument('--hmm_evalue', dest='evalue', default=0.001, type=float, metavar='',
                    help="E-value threshold. Default=0.001")

    pg_hmm.add_argument('--hmm_score', dest='score', default=20, type=float, metavar='',
                    help="Bit score threshold. Default=20")

    pg_hmm.add_argument('--hmm_maxseqlen', dest='maxseqlen', type=int, default=5000, metavar='',
                    help="Ignore query sequences larger than `maxseqlen`. Default=5000")

    pg_hmm.add_argument('--hmm_qcov', dest='qcov', type=float, metavar='',
                    help="min query coverage (from 0 to 1). Default=(disabled)")

    pg_hmm.add_argument('--Z', dest='Z', type=float, default=40000000, metavar='',
                    help='Fixed database size used in phmmer/hmmscan'
                        ' (allows comparing e-values among databases). Default=40,000,000')

    pg_diamond = parser.add_argument_group('diamond search_options')
	
    pg_diamond.add_argument('--dmnd_db',
		    help="Path to DIAMOND-compatible database")

    pg_diamond.add_argument('--matrix', dest='matrix', 
                    choices = ['BLOSUM62', 'BLOSUM90','BLOSUM80','BLOSUM50','BLOSUM45','PAM250','PAM70','PAM30'], 
                    default=None, help='Scoring matrix')

    pg_diamond.add_argument('--gapopen', dest='gapopen', type=int, default=None, 
                    help='Gap open penalty')

    pg_diamond.add_argument('--gapextend', dest='gapextend', type=int, default=None, 
                    help='Gap extend  penalty')

    pg_diamond.add_argument('--query-cover', dest='query_cover', type=float, default=0,
                    help='Report only alignments above the given percentage of query cover. Default=0')

    pg_diamond.add_argument('--subject-cover', dest='subject_cover', type=float, default=0,
                    help='Report only alignments above the given percentage of subject cover. Default=0')

    pg_seed = parser.add_argument_group('Seed ortholog search option')

    pg_seed.add_argument('--seed_ortholog_evalue', default=0.001, type=float, metavar='',
                    help='Min E-value expected when searching for seed eggNOG ortholog.'
                         ' Applies to phmmer/diamond searches. Queries not having a significant'
                         ' seed orthologs will not be annotated. Default=0.001')

    pg_seed.add_argument('--seed_ortholog_score', default=60, type=float, metavar='',
                    help='Min bit score expected when searching for seed eggNOG ortholog.'
                         ' Applies to phmmer/diamond searches. Queries not having a significant'
                         ' seed orthologs will not be annotated. Default=60')


    pg_out = parser.add_argument_group('Output options')

    pg_out.add_argument('--output', '-o', type=str, metavar='',
                    help="base name for output files")

    pg_out.add_argument('--resume', action="store_true",
                    help="Resumes a previous execution skipping reported hits in the output file.")

    pg_out.add_argument('--override', action="store_true",
                    help="Overwrites output files if they exist.")

    pg_out.add_argument("--no_refine", action="store_true",
                    help="Skip hit refinement, reporting only HMM hits.")

    pg_out.add_argument("--no_annot", action="store_true",
                    help="Skip functional annotation, reporting only hits")

    pg_out.add_argument("--no_search", action="store_true",
                    help="Skip HMM search mapping. Use existing hits file")

    pg_out.add_argument("--predict_ortho", action="store_true", help="The list of predicted orthologs")
    
    pg_out.add_argument("--report_orthologs", action="store_true",
                    help="The list of orthologs used for functional transferred are dumped into a separate file")

    pg_out.add_argument("--scratch_dir", metavar='', type=existing_dir,
                    help='Write output files in a temporary scratch dir, move them to final the final'
                        ' output dir when finished. Speed up large computations using network file'
                        ' systems.')

    pg_out.add_argument("--output_dir", default=os.getcwd(), type=existing_dir, metavar='',
                    help="Where output files should be written")

    pg_out.add_argument("--temp_dir", default=os.getcwd(), type=existing_dir, metavar='',
                    help="Where temporary files are created. Better if this is a local disk.")

    pg_out.add_argument('--no_file_comments', action="store_true",
                        help="No header lines nor stats are included in the output files")

    pg_out.add_argument('--keep_mapping_files', action='store_true',
                        help='Do not delete temporary mapping files used for annotation (i.e. HMMER and'
                        ' DIAMOND search outputs)')

    pg_predict = parser.add_argument_group('Predict orthologs options')

    pg_predict.add_argument('--target_taxa', type=str,
                          default= "all", nargs="+",
                            help='taxa that will be searched for orthologs')

    pg_predict.add_argument('--predict_output_format', choices=["per_query", "per_species"],
                            default= "per_species", help="Choose the output format among: per_query, per_species .Default = per_species")
    
    # exec mode
    g4 = parser.add_argument_group('Execution options')
    g4.add_argument('-m', dest='mode', choices = ['hmmer', 'diamond'], default='hmmer',
                    help='Default:hmmer')


    g4.add_argument('-i', dest="input", metavar='', type=existing_file,
                    help='Input FASTA file containing query sequences')

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

    g4.add_argument('--cpu', type=int, default=2, metavar='')

    g4.add_argument('--annotate_hits_table', type=str, metavar='',
                    help='Annotatate TSV formatted table of query->hits. 4 fields required:'
                    ' query, hit, evalue, score. Implies --no_search and --no_refine.')


    parser.add_argument('--version', action='store_true')

    args = parse_args(parser)

    _total_time = time.time()
    try:
        main(args)
    except emapperException:
        sys.exit(1)
    except:
        raise
    else:
        sys.exit(0)

