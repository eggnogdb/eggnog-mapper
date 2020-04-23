##
## CPCantalapiedra 2020

import os, sys
import time
from tempfile import mkdtemp
import multiprocessing
import shutil

from os.path import exists as pexists
from os.path import join as pjoin

from ..common import EGGNOG_DATABASES, get_db_info, TIMEOUT_LOAD_SERVER, get_oglevels_file, get_fasta_path, get_call_info, cleanup_og_name, get_data_path, gopen
from ..utils import colorify

from .hmmer_server import generate_idmap, server_functional, load_server, shutdown_server
from .hmmer_search import iter_hits, refine_hit
from .hmmer_seqio import iter_fasta_seqs

from .hmmer_search import SCANTYPE_MEM, SCANTYPE_DISK, QUERY_TYPE_SEQ, QUERY_TYPE_HMM, DB_TYPE_SEQ, DB_TYPE_HMM

class HmmerSearcher:

    cpu = None
    usemem = None
    scantype = None
    servermode = None
    no_refine = None

    db = None
    dbtype = None
    querytype = None
    translate = None

    resume = None
    no_file_comments = None

    evalue = score = qcov = Z = maxhits = maxseqlen = tempdir = None
    excluded_taxa = None

    temp_dir = None
    
    ##
    def __init__(self, args):
        self.cpu = args.cpu
        
        self.usemem = args.usemem
        if self.usemem:
            self.scantype = SCANTYPE_MEM
        else:
            self.scantype = SCANTYPE_DISK
            
        self.servermode = args.servermode
        self.no_refine = args.no_refine
        
        self.db = args.db
        self.dbtype = args.dbtype
        self.querytype = args.qtype
        self.translate = args.translate

        self.resume = args.resume
        self.no_file_comments = args.no_file_comments

        self.maxhits = args.maxhits
        self.maxseqlen = args.maxseqlen
        
        self.evalue = args.evalue
        self.score = args.score
        self.qcov = args.qcov
        
        self.Z = args.Z
        
        self.temp_dir = args.temp_dir

        self.excluded_taxa = args.excluded_taxa
        
        return

    ##
    def search(self, in_file, seed_orthologs_file, hmm_hits_file):
        
        host, port, dbpath, idmap = self.setup_hmm_search()

        # Start HMM SCANNING sequences (if requested)                                                                                                                              
        if not pexists(hmm_hits_file):
            self.dump_hmm_matches(in_file, hmm_hits_file, dbpath, port, idmap)

        if not self.no_refine and not pexists(seed_orthologs_file):
            if self.db == 'viruses':
                print('Skipping seed ortholog detection in "viruses" database')
                
            elif self.db in EGGNOG_DATABASES:
                self.refine_matches(in_file, seed_orthologs_file, hmm_hits_file)
                
            else:
                print('refined hits not available for custom hmm databases.')

        shutdown_server()

        return

    ##
    def setup_hmm_search(self):
        host = 'localhost'
        idmap = None


        connecting_to_server = False
        # If searching against a predefined database name
        if self.db in EGGNOG_DATABASES:
            dbpath, port = get_db_info(self.db)
            print(dbpath)
            db_present = [pexists(dbpath + "." + ext)
                          for ext in 'h3f h3i h3m h3p idmap'.split()]

            if False in db_present:
                print(db_present)
                print(colorify('Database %s not present. Use download_eggnog_database.py to fetch it' % (self.db), 'red'))
                raise ValueError('Database not found')

            if not self.no_refine:
                if not pexists(pjoin(get_data_path(), 'OG_fasta')):
                    print(colorify('Database data/OG_fasta/ not present. Use download_eggnog_database.py to fetch it', 'red'))
                    raise ValueError('Database fasta sequences not found')

            if self.scantype == SCANTYPE_MEM:
                idmap_file = dbpath + '.idmap'
                end_port = 53200

        # If searching against a custom hmm database
        elif os.path.isfile(self.db + '.h3f'):

            print(colorify(f"Preparing to query custom database {self.db}", 'green'))
            dbpath = self.db
            
            if self.scantype == SCANTYPE_MEM:
                idmap_file = dbpath + ".idmap"
                if not pexists(idmap_file):
                    if generate_idmap(self.db):
                        idmap_file = self.db + ".idmap"
                        print("idmap succesfully created!", file=sys.stderr)
                    else:
                        raise ValueError("idmap could not be created!")
                port = 53000
                end_port = 53200
            else:
                idmap_file = None
                port = None

        # If searching against a emapper hmm server
        elif ":" in self.db:
            dbname, host, port = map(str.strip, self.db.split(":"))
            self.scantype = SCANTYPE_MEM
            port = int(port)
            if dbname in EGGNOG_DATABASES:
                dbfile, port = get_db_info(dbname)
                self.db = dbname
            else:
                dbfile = dbname

            idmap_file = dbfile + '.idmap'
            if not pexists(idmap_file):
                raise ValueError("idmap file not found: %s" % idmap_file)

            dbpath = host
            if not server_functional(host, port, self.dbtype):
                print(colorify("eggnog-mapper server not found at %s:%s" % (host, port), 'red'))
                exit(1)
            connecting_to_server = True
        else:
            raise ValueError('Invalid database name/server')


        # If memory based searches requested, start server
        if self.scantype == SCANTYPE_MEM and not connecting_to_server:
            master_db, worker_db = None, None
            for try_port in range(port, end_port, 2):
                print(colorify("Loading server at localhost, port %s-%s" %
                               (try_port, try_port + 1), 'lblue'))
                dbpath, master_db, worker_db = load_server(
                    dbpath, try_port, try_port + 1, self.cpu)
                port = try_port
                ready = False
                for _ in range(TIMEOUT_LOAD_SERVER):
                    print("Waiting for server to become ready..."+str(host)+str(try_port))
                    time.sleep(1)
                    if not master_db.is_alive() or not worker_db.is_alive():
                        master_db.terminate()
                        master_db.join()
                        worker_db.terminate()
                        worker_db.join()
                        break
                    elif server_functional(host, port, self.dbtype):
                        ready = True
                        break
                if ready:
                    dbpath = host
                    break
                
        elif self.scantype == SCANTYPE_MEM:
            print(colorify("DB Server already running or not needed!", 'yellow'))
            dbpath = host

        # Preload seqid map to translate hits from hmmpgmd
        if self.scantype == SCANTYPE_MEM:
            print(colorify("Reading idmap %s" % idmap_file, color='lblue'))
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
            print(str(len(idmap)) + " names loaded")

        # If server mode, just listen for connections and exit when interrupted
        if self.servermode:
            while True:
                print(colorify("Server ready listening at %s:%s and using %d CPU cores" % (host, port, self.cpu), 'green'))
                print(colorify("Use `emapper.py -d %s:%s:%s (...)` to search against this server" % (self.db, host, port), 'lblue'))
                time.sleep(10)
            raise emapperException()
        else:
            return host, port, dbpath, idmap
    

    ##
    def dump_hmm_matches(self, fasta_file, hits_file, dbpath, port, idmap):
        hits_header = ("#query_name", "hit", "evalue", "sum_score", "query_length",
                       "hmmfrom", "hmmto", "seqfrom", "seqto", "query_coverage")

        # Cache previous results if resuming is enabled
        VISITED = set()
        if self.resume and pexists(hits_file):
            print(colorify("Resuming previous run. Reading computed output from %s" % hits_file, 'yellow'))
            VISITED = set([line.split('\t')[0].strip()
                           for line in open(hits_file) if not line.startswith('#')])
            print(str(len(VISITED)) + ' queries skipped')
            OUT = open(hits_file, 'a')
        else:
            OUT = open(hits_file, 'w')

        print(colorify("Sequence mapping starts now!", 'green'))
        if not self.no_file_comments:
            print(get_call_info(), file=OUT)
            print('# ' + '\t'.join(hits_header), file=OUT)
            
        total_time = 0
        last_time = time.time()
        start_time = time.time()
        qn = 0 # in case nothing to loop bellow
        for qn, (name, elapsed, hits, querylen, seq) in enumerate(iter_hits(
                                                            fasta_file,
                                                            self.translate,
                                                            self.querytype,
                                                            self.dbtype,
                                                            self.scantype,
                                                            dbpath,
                                                            port,
                                                            evalue_thr=self.evalue,
                                                            score_thr=self.score,
                                                            qcov_thr=self.qcov,
                                                            fixed_Z=self.Z,
                                                            max_hits=self.maxhits,
                                                            skip=VISITED,
                                                            maxseqlen=self.maxseqlen,
                                                            cpus=self.cpu,
                                                            base_tempdir=self.temp_dir)):

            if elapsed == -1:
                # error occurred
                print('\t'.join(
                    [name] + ['ERROR'] * (len(hits_header) - 1)), file=OUT)
            elif not hits:
                print('\t'.join([name] + ['-'] * (len(hits_header) - 1)), file=OUT)
            else:
                for hitindex, (hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore) in enumerate(hits):
                    hitname = hid
                    if idmap:
                        hitname = idmap[hid][0]

                    print('\t'.join(map(str, [name, hitname, heval, hscore,
                                                     int(querylen), int(hmmfrom),
                                                     int(hmmto), int(sqfrom),
                                                     int(sqto),
                                                     float(sqto - sqfrom) / querylen])), file=OUT)
            OUT.flush()

            # monitoring
            total_time += time.time() - last_time
            last_time = time.time()
            if qn and (qn % 25 == 0):
                print(qn + \
                       1, total_time, "%0.2f q/s" % ((float(qn + 1) / total_time)), file=sys.stderr)
                sys.stderr.flush()

        # Writes final stats
        elapsed_time = time.time() - start_time
        if not self.no_file_comments:
            print('# %d queries scanned' % (qn + 1), file=OUT)
            print('# Total time (seconds): '+str(elapsed_time), file=OUT)
            print('# Rate:', "%0.2f q/s" % ((float(qn + 1) / elapsed_time)), file=OUT)
        OUT.close()
        print(colorify(" Processed queries:%s total_time:%s rate:%s" %\
                       (qn+1, elapsed_time, "%0.2f q/s" % ((float(qn+1) / elapsed_time))), 'lblue'))

        return


    ##
    def refine_matches(self, fasta_file, refine_file, hits_file):
        refine_header = map(str.strip, '''#query_name, best_hit_eggNOG_ortholog,
                            best_hit_evalue, best_hit_score'''.split(','))

        print(colorify("Hit refinement starts now", 'green'))
        start_time = time.time()
        print(get_oglevels_file())
        og2level = dict([tuple(map(str.strip, line.split('\t')))
                         for line in gopen(get_oglevels_file())])
        OUT = open(refine_file, "w")

        if not self.no_file_comments:
            print(get_call_info(), file=OUT)
            print('\t'.join(refine_header), file=OUT)

        qn = 0 # in case no hits in loop bellow
        for qn, r in enumerate(self.process_nog_hits_file(hits_file, fasta_file, og2level,
                                                     translate=self.translate,
                                                     cpu=self.cpu,
                                                     excluded_taxa=self.excluded_taxa,
                                                     base_tempdir=self.temp_dir)):
            if qn and (qn % 25 == 0):
                total_time = time.time() - start_time
                print(str(qn + 1)+" "+str(total_time)+" %0.2f q/s (refinement)" % ((float(qn + 1) / total_time)), file=sys.stderr)
                sys.stderr.flush()
            query_name = r[0]
            best_hit_name = r[1]
            if best_hit_name == '-' or best_hit_name == 'ERROR':
                continue
            best_hit_evalue = float(r[2])
            best_hit_score = float(r[3])
            print('\t'.join(map(str, (query_name, best_hit_name,
                                             best_hit_evalue, best_hit_score))), file=OUT)
            #OUT.flush()

        elapsed_time = time.time() - start_time
        if not self.no_file_comments:
            print('# %d queries scanned' % (qn + 1), file=OUT)
            print('# Total time (seconds): '+str(elapsed_time), file=OUT)
            print('# Rate: '+"%0.2f q/s" % ((float(qn + 1) / elapsed_time)), file=OUT)
        OUT.close()
        print(colorify(" Processed queries:%s total_time:%s rate:%s" %\
                       (qn+1, elapsed_time, "%0.2f q/s" % ((float(qn+1) / elapsed_time))), 'lblue'))


    ##
    def process_nog_hits_file(self, hits_file, query_fasta, og2level, skip_queries=None,
                              translate=False, cpu=1, excluded_taxa=None, base_tempdir=None):

        sequences = {name: seq for name, seq in iter_fasta_seqs(
            query_fasta, translate=translate)}
        cmds = []
        visited_queries = set()

        if skip_queries:
            visited_queries.update(skip_queries)

        tempdir = mkdtemp(prefix='emappertmp_phmmer_', dir=base_tempdir)

        for line in gopen(hits_file):
            if line.startswith('#'):
                continue

            fields = list(map(str.strip, line.split('\t')))
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
            for r in pool.imap(refine_hit, cmds):
                yield r
            pool.terminate()

        shutil.rmtree(tempdir)
    
## END
