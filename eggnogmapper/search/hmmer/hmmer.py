##
## CPCantalapiedra 2020

import sys
import time
from tempfile import mkdtemp
import multiprocessing
import shutil

from os.path import join as pjoin, isdir as pisdir, isfile as pisfile

from ...common import \
    get_oglevels_file, get_OG_fasta_path, cleanup_og_name, \
    gopen, get_hmmer_databases, get_pfam_clans_file

from ...emapperException import EmapperException

from ...utils import colorify

from ..hits_io import parse_seeds

from .hmmer_server import shutdown_server_by_pid, create_servers, check_servers
from .hmmer_search import iter_hits, refine_hit
from .hmmer_seqio import iter_fasta_seqs

from .hmmer_search import SCANTYPE_MEM, SCANTYPE_DISK
from .hmmer_search_hmmpgmd import QUERY_TYPE_HMM, DB_TYPE_HMM
from .hmmer_setup import setup_hmm_search, SETUP_TYPE_EGGNOG, SETUP_TYPE_CUSTOM, SETUP_TYPE_REMOTE
from .hmmer_idmap import load_idmap_idx
from .hmmer_overlaps import process_overlaps, \
    CLEAN_OVERLAPS_ALL, CLEAN_OVERLAPS_CLANS, CLEAN_OVERLAPS_HMMSEARCH_ALL, CLEAN_OVERLAPS_HMMSEARCH_CLANS

class HmmerSearcher:

    name = "hmmer"
    
    call_info = None
    cpu = None
    usemem = None
    port = None
    end_port = None
    num_servers = None
    num_workers = None
    cpus_per_worker = None
    scantype = None
    setup_type = None
    servers = None

    servers_list = None
    db = None
    dbtype = None
    qtype = None
    translate = None
    trans_table = None

    resume = None
    no_file_comments = None

    evalue = score = qcov = Z = maxhits = report_no_hits = maxseqlen = cut_ga = None
    clean_overlaps = None
    excluded_taxa = None

    hmmcmd_temp_dir = None
    phmmer_temp_dir = None

    # Results
    queries = hits = no_hits = None
    hits_dict = None
    
    ##
    def __init__(self, args):
        self.call_info = args.call_info
        
        self.cpu = args.cpu
        
        self.usemem = args.usemem
        self.port = args.port
        self.end_port = args.end_port
        self.num_servers = args.num_servers
        self.num_workers = args.num_workers
        self.cpus_per_worker = args.cpus_per_worker
        
        if self.usemem or ":" in args.db or args.servers_list is not None:
            self.scantype = SCANTYPE_MEM
        else:
            self.scantype = SCANTYPE_DISK
        
        self.db = args.db
        self.servers_list = args.servers_list

        if args.qtype == QUERY_TYPE_HMM and args.dbtype == DB_TYPE_HMM:
            raise EmapperException("ERROR: HMM to HMM search is not supported.")
        
        self.dbtype = args.dbtype
        self.qtype = args.qtype
        
        self.translate = args.translate
        self.trans_table = args.trans_table

        self.resume = args.resume
        self.no_file_comments = args.no_file_comments

        self.maxhits = args.maxhits
        self.report_no_hits = args.report_no_hits
        self.maxseqlen = args.maxseqlen
        self.cut_ga = args.cut_ga
        self.clean_overlaps = args.clean_overlaps
        
        self.evalue = args.evalue
        self.score = args.score
        
        self.qcov = args.qcov
        
        self.Z = args.Z

        self.hmmcmd_temp_dir = mkdtemp(prefix='emappertmp_hmmcmd_', dir=args.temp_dir)
        self.phmmer_temp_dir = mkdtemp(prefix='emappertmp_phmmer_', dir=args.temp_dir)

        self.excluded_taxa = args.excluded_taxa
        
        return

    ##
    def get_call_info(self):
        return self.call_info
    
    ##
    def clear(self):
        if self.hmmcmd_temp_dir is not None and pisdir(self.hmmcmd_temp_dir):
            try:
                shutil.rmtree(self.hmmcmd_temp_dir)
            except OSError as err:
                print(f"Warning: OS error while removing {self.hmmcmd_temp_dir}", file = sys.stderr)
                print(f"OS error: {err}", file = sys.stderr)
                
        if self.phmmer_temp_dir is not None and pisdir(self.phmmer_temp_dir):
            try:
                shutil.rmtree(self.phmmer_temp_dir)
            except OSError as err:
                print(f"Warning: OS error while removing {self.phmmer_temp_dir}", file = sys.stderr)
                print(f"OS error: {err}", file = sys.stderr)

        # Shutdown server, If a temp local server was set up
        if (self.setup_type == SETUP_TYPE_EGGNOG or self.setup_type == SETUP_TYPE_CUSTOM) and self.scantype == SCANTYPE_MEM:
            for dbpath, port, master_pid, workers_pids in self.servers:
                shutdown_server_by_pid(master_pid, workers_pids)
            self.servers = []
                
                
        return
    
    ##
    def search_hmm_matches(self, in_file, hmm_hits_file, silent = False):
        
        # Prepare HMM database and/or server
        dbname, dbpath, host, port, end_port, idmap_file, self.setup_type = setup_hmm_search(self.db, self.scantype, self.dbtype, self.qtype,
                                                                                             self.port, self.end_port, self.servers_list, silent)
        
        if (self.setup_type == SETUP_TYPE_EGGNOG or self.setup_type == SETUP_TYPE_CUSTOM) and self.scantype == SCANTYPE_MEM:
            dbpath, host, port, self.servers = create_servers(self.dbtype, dbpath, host, port, end_port,
                                                              self.num_servers, self.num_workers, self.cpus_per_worker,
                                                              silent)
            
        elif self.setup_type == SETUP_TYPE_REMOTE and self.scantype == SCANTYPE_MEM:
            dbpath, host, port, self.servers = check_servers(self.dbtype, self.qtype, dbpath, host, port, self.servers_list)

            
        # Search for HMM hits (OG)
        if self.servers is None:
            hosts = [(dbpath, port)]
        else:
            hosts = [(dbpath, port) for dbpath, port, master_pid, workers_pids in self.servers] # I cannot use master_db and workers for later multiprocessing
            
        self.dump_hmm_matches(in_file, hmm_hits_file, dbpath, port, hosts, idmap_file, silent)

        # Shutdown server, If a temp local server was set up
        if (self.setup_type == SETUP_TYPE_EGGNOG or self.setup_type == SETUP_TYPE_CUSTOM) and self.scantype == SCANTYPE_MEM:
            for dbpath, port, master_pid, workers_pids in self.servers:
                shutdown_server_by_pid(master_pid, workers_pids)
            self.servers = []
            
        return
        
    ##
    def search(self, in_file, seed_orthologs_file, hmm_hits_file):

        hits = None
        
        print(f"hmmer.py:search DB: {self.db}")
        
        # Prepare HMM database and/or server
        dbname, dbpath, host, port, end_port, idmap_file, self.setup_type = setup_hmm_search(self.db, self.scantype, self.dbtype, self.qtype,
                                                                                             self.port, self.end_port, self.servers_list)

        print(f"hmmer.py:search DB: {self.db}, name {dbname}, path {dbpath}, host {host}, port {port}, endport {end_port}, idmap {idmap_file}")
        
        if (self.setup_type == SETUP_TYPE_EGGNOG or self.setup_type == SETUP_TYPE_CUSTOM) and self.scantype == SCANTYPE_MEM:
            dbpath, host, port, self.servers = create_servers(self.dbtype, dbpath, host, port, end_port,
                                                              self.num_servers, self.num_workers, self.cpus_per_worker)

        elif self.setup_type == SETUP_TYPE_REMOTE and self.scantype == SCANTYPE_MEM:
            dbpath, host, port, self.servers = check_servers(self.dbtype, self.qtype, dbpath, host, port, self.servers_list)
            
            
        # Search for HMM hits (OG)
        if self.servers is None:
            hosts = [(dbpath, port)]
        else:
            hosts = [(dbpath, port) for dbpath, port, master_pid, workers_pids in self.servers] # I cannot use master_db and workers for later multiprocessing
            
        self.dump_hmm_matches(in_file, hmm_hits_file, dbpath, port, hosts, idmap_file)
        
        # Search for seed orthologs within the HMM hits
        if dbname == 'viruses':
            print('Skipping seed ortholog detection in "viruses" database')

        elif dbname in get_hmmer_databases():
            hits = self.refine_matches(dbname, in_file, seed_orthologs_file, hmm_hits_file)

        else:
            print(f'Could not find {dbname} among eggnog databases. Skipping seed ortholog detection.')

        # Shutdown server, If a temp local server was set up
        if (self.setup_type == SETUP_TYPE_EGGNOG or self.setup_type == SETUP_TYPE_CUSTOM) and self.scantype == SCANTYPE_MEM:
            for dbpath, port, master_pid, workers_pids in self.servers:
                shutdown_server_by_pid(master_pid, workers_pids)
            self.servers = []
            
        return hits
    
    ##
    def dump_hmm_matches(self, in_file, hits_file, dbpath, port, servers, idmap_file, silent = False):
        hits_header = ("query_name", "hit", "evalue", "sum_score", "query_length",
                       "hmmfrom", "hmmto", "seqfrom", "seqto", "query_coverage")

        CLANS_FILE = get_pfam_clans_file()
        
        # Cache previous results if resuming is enabled
        VISITED = set()
        if self.resume and pisfile(hits_file):
            print(colorify("Resuming previous run. Reading computed output from %s" % hits_file, 'yellow'))
            VISITED = set([line.split('\t')[0].strip()
                           for line in open(hits_file) if not line.startswith('#')])
            print(str(len(VISITED)) + ' queries skipped')
            OUT = open(hits_file, 'a')
        else:
            OUT = open(hits_file, 'w')

        if not self.no_file_comments:
            print(self.get_call_info(), file=OUT)
            if self.resume == False or not pisfile(hits_file):
                print('# ' + '\t'.join(hits_header), file=OUT)
        
        total_time = 0
        last_time = time.time()
        start_time = time.time()

        # Loading the DB identifiers will also taken into account for total_time
        idmap_idx = None
        if idmap_file:
            idmap_idx = load_idmap_idx(idmap_file)

        if silent == False:
            print(colorify("Sequence mapping starts now!", 'green'))

        if self.clean_overlaps is not None and self.clean_overlaps in [CLEAN_OVERLAPS_HMMSEARCH_ALL, CLEAN_OVERLAPS_HMMSEARCH_CLANS]:
            namedhits = []
            
        qn = -1 # in case nothing to loop bellow
        for name, elapsed, hits, querylen, seq in iter_hits(in_file,
                                                            self.translate,
                                                            self.qtype,
                                                            self.dbtype,
                                                            self.scantype,
                                                            dbpath,
                                                            port,
                                                            servers,
                                                            evalue_thr=self.evalue,
                                                            score_thr=self.score,
                                                            qcov_thr=self.qcov,
                                                            fixed_Z=self.Z,
                                                            max_hits=self.maxhits,
                                                            skip=VISITED,
                                                            maxseqlen=self.maxseqlen,
                                                            cut_ga=self.cut_ga,
                                                            cpus=self.cpu,
                                                            base_tempdir=self.hmmcmd_temp_dir,
                                                            silent=silent,
                                                            trans_table=self.trans_table):

            if elapsed == -1:
                # error occurred. hits should contain a single element with the error msg. e.g. hits = ["ERROR_MSG"]
                print('\t'.join([name] + hits * (len(hits_header) - 1)), file=sys.stderr)
                print('\t'.join([name] + ['-'] * (len(hits_header) - 1)), file=OUT)
            elif not hits and self.report_no_hits == True:
                print('\t'.join([name] + ['-'] * (len(hits_header) - 1)), file=OUT)
            else:
                
                if self.clean_overlaps is not None and self.clean_overlaps in [CLEAN_OVERLAPS_ALL, CLEAN_OVERLAPS_CLANS]:
                    hits = process_overlaps(hits, self.clean_overlaps, CLANS_FILE, idmap_idx)

                elif self.clean_overlaps is not None and self.clean_overlaps in [CLEAN_OVERLAPS_HMMSEARCH_ALL, CLEAN_OVERLAPS_HMMSEARCH_CLANS]:
                    namedhits.append((name, querylen, hits))
                    
                # output
                if self.clean_overlaps is None or self.clean_overlaps not in [CLEAN_OVERLAPS_HMMSEARCH_ALL, CLEAN_OVERLAPS_HMMSEARCH_CLANS]:
                    self.output_hits(name, querylen, hits, OUT, idmap_idx)
                    
            OUT.flush()

            qn += 1

            # monitoring
            total_time += time.time() - last_time
            last_time = time.time()
            if qn and (qn % 25 == 0):
                if silent == False:
                    print(qn + 1, total_time, "%0.2f q/s" % ((float(qn + 1) / total_time)), file=sys.stderr)
                    sys.stderr.flush()

        if self.clean_overlaps is not None and self.clean_overlaps in [CLEAN_OVERLAPS_HMMSEARCH_ALL, CLEAN_OVERLAPS_HMMSEARCH_CLANS]:
            if silent == False:
                sys.stderr.write("Postprocessing overlapping hits...\n")
            namedhits = process_overlaps(namedhits, self.clean_overlaps, CLANS_FILE, idmap_idx)
            for (name, querylen, hits) in namedhits:
                self.output_hits(name, querylen, hits, OUT, idmap_idx)

        # Writes final stats
        elapsed_time = time.time() - start_time
        if not self.no_file_comments:
            print('## %d queries scanned' % (qn + 1), file=OUT)
            print('## Total time (seconds): '+str(elapsed_time), file=OUT)
            print('## Rate:', "%0.2f q/s" % ((float(qn + 1) / elapsed_time)), file=OUT)
        OUT.close()
        if silent == False:
            print(colorify(" Processed queries:%s total_time:%s rate:%s" %\
                           (qn+1, elapsed_time, "%0.2f q/s" % ((float(qn+1) / elapsed_time))), 'lblue'))

        return

    
    ##
    def output_hits(self, name, querylen, hits, OUT, idmap_idx = None):
        for hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore in hits:
            hitname = hid
            if idmap_idx:
                hitname = idmap_idx[hid][0]

            print('\t'.join(map(str, [name, hitname, f'{heval:.1e}', f'{hscore:.1f}', int(querylen),
                                      int(hmmfrom), int(hmmto), int(sqfrom), int(sqto),
                                      float(sqto - sqfrom) / querylen])), file=OUT)

        return


    ##
    def refine_matches(self, dbname, in_file, refine_file, hits_file):
        refine_header = map(str.strip, '''query_name, best_hit_eggNOG_ortholog,
                            best_hit_evalue, best_hit_score'''.split(','))

        print(colorify("Hit refinement starts now", 'green'))
        start_time = time.time()

        # Cache previous results if resuming is enabled
        last_resumed_query = None
        if self.resume == True:
            if pisfile(hits_file):
                if pisfile(refine_file):
                    hits_parser = parse_seeds(refine_file)
                    for hit in hits_parser:
                        yield hit
                        last_resumed_query = hit[0]
            else:
                raise EmapperException(f"Couldn't find hits file {hits_file} to resume.")
                
            OUT = open(refine_file, 'a')
        else:
            OUT = open(refine_file, 'w')

        if not self.no_file_comments:
            print(self.get_call_info(), file=OUT)
            if self.resume == False:
                print('# ' + '\t'.join(refine_header), file=OUT)

        qn = -1 # in case no hits in loop bellow
        sequences = {name: seq for name, seq in iter_fasta_seqs(in_file, translate=self.translate, trans_table=self.trans_table)}
        self.queries = set(sequences.keys())
        for qn, r in enumerate(self.process_nog_hits_file(dbname, hits_file, sequences,
                                                          last_resumed_query, cpu=self.cpu,
                                                          excluded_taxa=self.excluded_taxa)):
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
            
            yield [query_name, best_hit_name, best_hit_evalue, best_hit_score]
            #OUT.flush()

        elapsed_time = time.time() - start_time
        if not self.no_file_comments:
            print('## %d queries scanned' % (qn + 1), file=OUT)
            print('## Total time (seconds): '+str(elapsed_time), file=OUT)
            print('## Rate: '+"%0.2f q/s" % ((float(qn + 1) / elapsed_time)), file=OUT)
        OUT.close()
        print(colorify(" Processed queries:%s total_time:%s rate:%s" %\
                       (qn+1, elapsed_time, "%0.2f q/s" % ((float(qn+1) / elapsed_time))), 'lblue'))
        return


    ##
    def process_nog_hits_file(self, dbname, hits_file, sequences, last_resumed_query,
                              cpu=1, excluded_taxa=None):

        cmds = []
        visited_queries = set()

        # semaphore to start processing new hits
        last_resumed_query_found = False if last_resumed_query is not None else True
        
        for line in gopen(hits_file):
            if line.startswith('#'):
                continue
            
            fields = list(map(str.strip, line.split('\t')))
            seqname = fields[0]

            if last_resumed_query is not None:
                if seqname == last_resumed_query:
                    last_resumed_query_found = True
                    continue
                else:
                    if last_resumed_query_found == False:
                        continue
                    else:
                        last_resumed_query = None # start parsing new queries
                

            # if there is no hit (check stderr for possible errors of each specific query)
            if fields[1] == '-':
                continue

            if seqname in visited_queries:
                continue

            hitname = cleanup_og_name(fields[1])
            # level = og2level[hitname]

            seq = sequences[seqname]
            visited_queries.add(seqname)
            target_fasta = get_OG_fasta_path(dbname, hitname)
            # print(f"hmmer.py:process_nog_hits_file Target FASTA: {target_fasta}")
            # target_fasta = pjoin(get_fasta_path(), level, "%s.fa" % hitname)
            
            cmds.append([seqname, seq, target_fasta, excluded_taxa, self.phmmer_temp_dir])

        if cmds is not None and len(cmds) > 0:
            # multiprocessing.set_start_method("spawn")
            pool = multiprocessing.Pool(cpu)
            for r in pool.imap(refine_hit, cmds):
                yield r
            pool.close()
            pool.terminate()

        return
    
## END
