##
## CPCantalapiedra 2020

import sys
import time
from tempfile import mkdtemp
import multiprocessing
import shutil

from os.path import exists as pexists
from os.path import join as pjoin

from ..common import EGGNOG_DATABASES, get_oglevels_file, get_fasta_path, get_call_info, cleanup_og_name, gopen
from ..utils import colorify

from .hmmer_server import shutdown_server_by_pid, create_servers, check_servers
from .hmmer_search import iter_hits, refine_hit
from .hmmer_seqio import iter_fasta_seqs

from .hmmer_search import SCANTYPE_MEM, SCANTYPE_DISK, QUERY_TYPE_SEQ
from .hmmer_setup import setup_hmm_search, SETUP_TYPE_EGGNOG, SETUP_TYPE_CUSTOM, SETUP_TYPE_REMOTE
from .hmmer_idmap import load_idmap_idx

class HmmerSearcher:

    cpu = None
    usemem = None
    num_servers = None
    num_workers = None
    cpus_per_worker = None
    scantype = None

    db = None
    dbtype = None
    qtype = None
    translate = None

    resume = None
    no_file_comments = None

    evalue = score = qcov = Z = maxhits = report_no_hits = maxseqlen = cut_ga = None
    clean_overlaps = None
    excluded_taxa = None

    temp_dir = None
    
    ##
    def __init__(self, args):
        self.cpu = args.cpu
        
        self.usemem = args.usemem
        self.num_servers = args.num_servers
        self.num_workers = args.num_workers
        self.cpus_per_worker = args.cpus_per_worker
        
        if self.usemem or ":" in args.db or args.servers_list is not None:
            self.scantype = SCANTYPE_MEM
        else:
            self.scantype = SCANTYPE_DISK
        
        self.db = args.db
        self.servers_list = args.servers_list
        self.dbtype = args.dbtype
        self.qtype = args.qtype
        self.translate = args.translate

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
        
        self.temp_dir = args.temp_dir

        self.excluded_taxa = args.excluded_taxa
        
        return
    
        
    ##
    def search_hmm_matches(self, in_file, hmm_hits_file):
        
        annot = None
        
        # Prepare HMM database and/or server
        dbname, dbpath, host, port, end_port, idmap_file, setup_type = setup_hmm_search(self.db, self.scantype, self.dbtype, self.qtype, self.servers_list)

        servers = None
        if (setup_type == SETUP_TYPE_EGGNOG or setup_type == SETUP_TYPE_CUSTOM) and self.scantype == SCANTYPE_MEM:
            dbpath, host, port, servers = create_servers(self.dbtype, dbpath, host, port, end_port,
                                                         self.num_servers, self.num_workers, self.cpus_per_worker)
            
        elif setup_type == SETUP_TYPE_REMOTE and self.scantype == SCANTYPE_MEM:
            dbpath, host, port, servers = check_servers(self.dbtype, self.qtype, dbpath, host, port, self.servers_list)

            
        # Search for HMM hits (OG)
        # if not pexists(hmm_hits_file): This avoids resuming the previous run
        if servers is None:
            hosts = [(dbpath, port)]
        else:
            hosts = [(dbpath, port) for dbpath, port, master_pid, workers_pids in servers] # I cannot use master_db and workers for later multiprocessing
            
        self.dump_hmm_matches(in_file, hmm_hits_file, dbpath, port, hosts, idmap_file)

        # Shutdown server, If a temp local server was set up
        if (setup_type == SETUP_TYPE_EGGNOG or setup_type == SETUP_TYPE_CUSTOM) and self.scantype == SCANTYPE_MEM:
            for dbpath, port, master_pid, workers_pids in servers:
                shutdown_server_by_pid(master_pid, workers_pids)
            
        return
        
    ##
    def search(self, in_file, seed_orthologs_file, hmm_hits_file):

        annot = None
        
        # Prepare HMM database and/or server
        dbname, dbpath, host, port, end_port, idmap_file, setup_type = setup_hmm_search(self.db, self.scantype, self.dbtype, self.qtype, self.servers_list)

        servers = None
        if (setup_type == SETUP_TYPE_EGGNOG or setup_type == SETUP_TYPE_CUSTOM) and self.scantype == SCANTYPE_MEM:
            dbpath, host, port, servers = create_servers(self.dbtype, dbpath, host, port, end_port,
                                                         self.num_servers, self.num_workers, self.cpus_per_worker)

        elif setup_type == SETUP_TYPE_REMOTE and self.scantype == SCANTYPE_MEM:
            dbpath, host, port, servers = check_servers(self.dbtype, self.qtype, dbpath, host, port, self.servers_list)
            
            
        # Search for HMM hits (OG)
        # if not pexists(hmm_hits_file): This avoids resuming the previous run
        if servers is None:
            hosts = [(dbpath, port)]
        else:
            hosts = [(dbpath, port) for dbpath, port, master_pid, workers_pids in servers] # I cannot use master_db and workers for later multiprocessing
            
        self.dump_hmm_matches(in_file, hmm_hits_file, dbpath, port, hosts, idmap_file)
        
        # Search for seed orthologs within the HMM hits
        if dbname == 'viruses':
            print('Skipping seed ortholog detection in "viruses" database')
            annot = False

        elif dbname in EGGNOG_DATABASES:
            if not pexists(seed_orthologs_file):
                self.refine_matches(in_file, seed_orthologs_file, hmm_hits_file)

        else:
            print(f'Could not find {dbname} among eggnog databases.')
            annot = False

        # Shutdown server, If a temp local server was set up
        if (setup_type == SETUP_TYPE_EGGNOG or setup_type == SETUP_TYPE_CUSTOM) and self.scantype == SCANTYPE_MEM:
            for dbpath, port, master_pid, workers_pids in servers:
                shutdown_server_by_pid(master_pid, workers_pids)
            
        return annot



    ##
    def dump_hmm_matches(self, in_file, hits_file, dbpath, port, servers, idmap_file):
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

        if not self.no_file_comments:
            print(get_call_info(), file=OUT)
            print('# ' + '\t'.join(hits_header), file=OUT)
        
        total_time = 0
        last_time = time.time()
        start_time = time.time()

        # Loading the DB identifiers will also taken into account for total_time
        idmap_idx = None
        if idmap_file:
            idmap_idx = load_idmap_idx(idmap_file)

        print(colorify("Sequence mapping starts now!", 'green'))
        
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
                                                            base_tempdir=self.temp_dir):

            if elapsed == -1:
                # error occurred. hits should contain a single element with the error msg. e.g. hits = ["ERROR_MSG"]
                print('\t'.join([name] + hits * (len(hits_header) - 1)), file=OUT)
            elif not hits and self.report_no_hits == True:
                print('\t'.join([name] + ['-'] * (len(hits_header) - 1)), file=OUT)
            else:

                if self.clean_overlaps == True:
                    clean_doms = self.process_overlaps(hits)
                    self.output_hits(name, querylen, clean_doms, OUT, idmap_idx)
                    
                else:
                    self.output_hits(name, querylen, hits, OUT, idmap_idx)
                    
            OUT.flush()

            qn += 1

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
    def process_overlaps(self, hits):
        # for hit in hits:
        #     print(hit)
        # sorted_hits = sorted(hits, key=lambda x: float(x[1]))
        # for hit in sorted_hits:
        #     print(hit)
        # sys.exit(1)
        clean_doms = []
        total_range = set()
        
        for hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore in hits:
            hmmfrom, hmmto, sqfrom, sqto = map(int, [hmmfrom, hmmto, sqfrom, sqto])
            new_span = set(range(sqfrom, sqto+1))
            
            total_overlap = new_span & total_range
            if len(total_overlap) > 0:
                best = True
                tmp_clean_doms = []
                tmp_overlapping = []

                for phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore in clean_doms:
                    prev_span = set(range(psqfrom, psqto+1))
                    overlap = new_span & prev_span
                    if len(overlap) > 0 and best == True:
                        if heval > pheval:
                            best = False
                        tmp_overlapping.append([phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore])
                    else:
                        tmp_clean_doms.append([phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore])
                    
                if best == True:
                    tmp_clean_doms.append([hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore])
                else:
                    tmp_clean_doms.extend(tmp_overlapping)

                # update clean_doms and total_range
                clean_doms = tmp_clean_doms
                for phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore in clean_doms:
                    clean_span = set(range(psqfrom, psqto+1))
                    total_range.update(clean_span)
            else:
                clean_doms.append([hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore])
                total_range.update(new_span)
        
        return clean_doms
    
        
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
    def refine_matches(self, in_file, refine_file, hits_file):
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
        for qn, r in enumerate(self.process_nog_hits_file(hits_file, in_file, og2level,
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
            target_fasta = pjoin(get_fasta_path(), level, "%s.fa" % hitname)
            
            cmds.append([seqname, seq, target_fasta, excluded_taxa, tempdir])

        if cmds:
            pool = multiprocessing.Pool(cpu)
            for r in pool.imap(refine_hit, cmds):
                yield r
            pool.terminate()

        shutil.rmtree(tempdir)
    
## END
