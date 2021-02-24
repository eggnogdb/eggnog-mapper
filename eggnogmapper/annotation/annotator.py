##
## CPCantalapiedra 2019

import sys
import time
import multiprocessing
    
from ..emapperException import EmapperException
from ..common import get_call_info, TAX_SCOPE_AUTO, TAX_SCOPE_AUTO_BROAD
from ..utils import colorify
from ..search.hmmer.hmmer_seqio import iter_fasta_seqs

from .db_sqlite import get_eggnog_db
from .ncbitaxa.ncbiquery import get_ncbi
from .pfam.pfam_modes import run_pfam_mode, PFAM_TRANSFER_NARROWEST_OG, PFAM_TRANSFER_SEED_ORTHOLOG, PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO

from .annotator_worker import annotate_hit_line_mem, annotate_hit_line_ondisk
from . import output

ANNOTATIONS_HEADER = output.ANNOTATIONS_HEADER

PFAM_COL = -1 # position of PFAMs annotations in list of annotations

##
class Annotator:

    hits = hits_dict = None
    annotations = annotations_dict = None
    orthologs = None
    
    annot = report_orthologs = None
    
    dbmem = None

    no_file_comments = cpu = None

    # options for pfam hmmpgmd searches
    num_servers = num_workers = cpus_per_worker = port = end_port = None

    seed_ortholog_score = seed_ortholog_evalue = None
    tax_scope_mode = tax_scope_id = target_taxa = target_orthologs = excluded_taxa = None
    TAXONOMIC_RESOLUTION = None
    go_evidence = go_excluded = None
    pfam_realign = pfam_transfer = queries_fasta = translate = temp_dir = None
    md5 = None
        
    ##
    def __init__(self, args, annot, report_orthologs):

        self.annot = annot
        self.report_orthologs = report_orthologs

        self.dbmem = args.dbmem
        
        self.no_file_comments = args.no_file_comments
        self.cpu = args.cpu
        self.num_servers = args.num_servers
        self.num_workers = args.num_workers
        self.cpus_per_worker = args.cpus_per_worker
        self.port = args.port
        self.end_port = args.end_port
        self.seed_ortholog_score = args.seed_ortholog_score
        self.seed_ortholog_evalue = args.seed_ortholog_evalue

        self.tax_scope_mode = args.tax_scope_mode
        if self.tax_scope_mode == "auto":
            self.TAXONOMIC_RESOLUTION = TAX_SCOPE_AUTO
        elif self.tax_scope_mode == "auto_broad":
            self.TAXONOMIC_RESOLUTION = TAX_SCOPE_AUTO_BROAD
        else:
            # self.TAXONOMIC_RESOLUTION = None
            pass
        
        self.tax_scope_id = args.tax_scope_id
                
        self.target_taxa = args.target_taxa
        self.target_orthologs = args.target_orthologs
        self.excluded_taxa = args.excluded_taxa
            
        self.go_evidence = args.go_evidence
        self.go_excluded = args.go_excluded
        
        self.pfam_realign = args.pfam_realign
        self.pfam_transfer = args.pfam_transfer
        
        self.queries_fasta = args.input
        self.translate = args.translate
        self.temp_dir = args.temp_dir
        
        self.md5 = args.md5
        
        return

    #
    def get_hits(self):
        return self.hits

    def get_hits_dict(self):
        hits_dict = None
        if self.hits_dict is not None:
            hits_dict = self.hits_dict
        elif self.hits is not None:
            hits_dict = {hit[0]:hit for hit in self.hits}
            self.hits_dict = hits_dict
        # else: None
        return hits_dict
    
    #
    def get_annotations(self):
        return self.annotations
    
    def get_annotations_dict(self):
        annotations_dict = None
        if self.annotations_dict is not None:
            annotations_dict = self.annotations_dict
        elif self.annotations is not None:
            annotations_dict = {annot[0]:annot for annot in self.annotations}
            self.annotations_dict = annotations_dict
        # else: None
        return annotations_dict

    #
    def get_orthologs(self):
        return self.orthologs
    
    ##
    def annotate(self, hits_gen_func, store_hits, annot_file, orthologs_file, pfam_file):

        if self.report_orthologs == True or self.annot == True:            
            ##
            # md5 hashes
            if self.md5 == True:
                print(colorify("Creating md5 hashes of input sequences", 'green'))
                md5_queries = md5_seqs(self.queries_fasta)
            else:
                md5_queries = None

            md5_field = (self.md5 == True and md5_queries is not None)

            ##
            # Annotations
            print(colorify("Functional annotation of refined hits starts now", 'green'))

            #
            # Prepare output files and print headers and call info
            ORTHOLOGS_OUT = None
            if self.report_orthologs == True:
                ORTHOLOGS_OUT = open(orthologs_file, "w")            
                output.output_orthologs_header(ORTHOLOGS_OUT, self.no_file_comments)

            ANNOTATIONS_OUT = None
            if self.annot == True:
                ANNOTATIONS_OUT = open(annot_file, "w")
                output.output_annotations_header(ANNOTATIONS_OUT, self.no_file_comments, md5_field)

            # closures to generate output

            output_orthologs_f = None
            if self.report_orthologs == True:
                ncbi = get_ncbi(usemem = True)
                output_orthologs_f = output.output_orthologs_closure(ORTHOLOGS_OUT, ncbi)

            output_annotations_f = None
            if self.annot == True:
                output_annotations_f = output.output_annotations_closure(ANNOTATIONS_OUT, md5_field, md5_queries)

            ##
            # Obtain annotations
            qn, elapsed_time = self._annotate(hits_gen_func, store_hits, pfam_file, output_orthologs_f, output_annotations_f)

            ##
            # Output footer and close files
            if self.report_orthologs == True:
                output.output_orthologs_footer(ORTHOLOGS_OUT, self.no_file_comments, qn, elapsed_time)
                ORTHOLOGS_OUT.close()

            if self.annot == True:
                output.output_annotations_footer(ANNOTATIONS_OUT, self.no_file_comments, qn, elapsed_time)
                ANNOTATIONS_OUT.close()
            
        return

    ##
    def _annotate(self, hits_gen_func, store_hits, pfam_file, output_ortho_f, output_annot_f):
        
        start_time = time.time()
        
        if self.dbmem == True:
            all_orthologs, all_annotations, qn = self._annotate_dbmem(hits_gen_func, store_hits, pfam_file)
        else:
            all_orthologs, all_annotations, qn = self._annotate_ondisk(hits_gen_func, store_hits, pfam_file)

        elapsed_time = time.time() - start_time
        print(colorify(f" Processed queries:{qn} total_time:{elapsed_time} rate:{(float(qn) / elapsed_time):.2f} q/s", 'lblue'))
            
        ##
        # PFAMs annotation
        if self.annot == True and self.pfam_realign in [PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO] and all_annotations is not None and len(all_annotations) > 0:
            all_annotations = run_pfam_mode(self.pfam_realign, all_annotations, self.queries_fasta, self.translate,
                                            self.cpu, self.num_servers, self.num_workers, self.cpus_per_worker, self.port, self.end_port,
                                            self.temp_dir, pfam_file)
            
            elapsed_time = time.time() - start_time
            print(colorify(f" Processed queries:{qn} total_time:{elapsed_time} rate:{(float(qn) / elapsed_time):.2f} q/s", 'lblue'))

        ##
        # Output rows
        if self.report_orthologs == True:
            output_ortho_f(all_orthologs)
            
        if self.annot == True:
            output_annot_f(all_annotations)

        self.orthologs = all_orthologs
        self.annotations = all_annotations
            
        return qn, elapsed_time


    ##
    def _annotate_dbmem(self, hits_gen_func, store_hits, pfam_file):
        all_orthologs = {}
        all_annotations = []

        ##
        # Load sqlite DBs into memory
        
        start_time = time.time() # do not take into account time to load the db into memory
        eggnog_db = get_eggnog_db(usemem = True)
        ncbi = get_ncbi(usemem = True)
        total_time = time.time() - start_time
        print(colorify(f"Time to load the DB into memory: {total_time}", "lblue"), file=sys.stderr)
        sys.stderr.flush()
        
        ##
        # Annotate hits
        
        start_time = time.time() # do not take into account time to load the db into memory
        
        qn = 0
        try:
            for result in map(annotate_hit_line_mem, self.iter_hit_lines(hits_gen_func, store_hits)):
                qn += 1
                if qn and (qn % 100 == 0):
                    total_time = time.time() - start_time
                    print(f"{qn} {total_time} {(float(qn) / total_time):.2f} q/s (func. annotation)", file=sys.stderr)
                    sys.stderr.flush()

                if result:
                    self._process_annot_result(result, all_orthologs, all_annotations)

        except EmapperException:
            raise
        except Exception as e:
            # import traceback
            # traceback.print_exc()
            raise EmapperException(f"Error: annotation went wrong for query number {qn}. "+str(e))
        finally:
            eggnog_db.close()
            ncbi.close()

        elapsed_time = time.time() - start_time
        print(colorify(f" All queries processed. Time to perform queries:{elapsed_time} rate:{(float(qn) / elapsed_time):.2f} q/s", 'lblue'))
                    
        return all_orthologs, all_annotations, qn

    
    ##
    def _annotate_ondisk(self, hits_gen_func, store_hits, pfam_file):

        all_orthologs = {}
        all_annotations = []
        
        multiprocessing.set_start_method("spawn")
        
        pool = multiprocessing.Pool(self.cpu)        

        start_time = time.time() # do not take into account time to load the pool of processes
        
        qn = 0
        try:
            for result in pool.imap(annotate_hit_line_ondisk, self.iter_hit_lines(hits_gen_func, store_hits)):
                qn += 1
                if qn and (qn % 100 == 0):
                    total_time = time.time() - start_time
                    print(f"{qn} {total_time} {(float(qn) / total_time):.2f} q/s (func. annotation)", file=sys.stderr)
                    sys.stderr.flush()

                if result:
                    self._process_annot_result(result, all_orthologs, all_annotations)      

        except EmapperException:
            raise
        except Exception as e:
            # import traceback
            # traceback.print_exc()
            raise EmapperException(f"Error: annotation went wrong for query number {qn}. "+str(e))
        finally:
            pool.close()
            pool.terminate() # it should remove the global eggnog_db variables also

        elapsed_time = time.time() - start_time
        print(colorify(f" All queries processed. Time to perform queries:{elapsed_time} rate:{(float(qn) / elapsed_time):.2f} q/s", 'lblue'))
            
        return all_orthologs, all_annotations, qn

    
    ##
    def _process_annot_result(self, result, all_orthologs, all_annotations):
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations, 
         narr_og_name, narr_og_cat, narr_og_desc,
         best_og_name, best_og_cat, best_og_desc,                     
         match_nogs_names, all_orthologies, annot_orthologs) = result

        if self.report_orthologs == True and all_orthologies is not None and annot_orthologs is not None:
            # filter co-orthologs to keep only target_orthologs: "all", "one2one", ...
            if query_name in all_orthologs:
                query_orthologs = all_orthologs[query_name]
            else:
                query_orthologs = {}
                all_orthologs[query_name] = query_orthologs

            for target in all_orthologies:
                if target in query_orthologs:
                    query_orthologs[target].update(all_orthologies[target])
                else:
                    query_orthologs[target] = set(all_orthologies[target])
            if "annot_orthologs" in query_orthologs:
                query_orthologs["annot_orthologs"].update(annot_orthologs)
            else:
                query_orthologs["annot_orthologs"] = set(annot_orthologs)

        if self.annot == True and annotations is not None:
            # prepare annotations for printing
            annot_columns = [query_name, best_hit_name, str(best_hit_evalue), str(best_hit_score),
                             ",".join(match_nogs_names), 
                             narr_og_name, narr_og_cat.replace('\n', ''), narr_og_desc.replace('\n', ' ')]

            # # If tax_scope_mode == narrowest there is no need to output best_og colums
            # if self.tax_scope_id is not None or self.tax_scope_mode != "narrowest":
            annot_columns.extend([best_og_name, best_og_cat.replace('\n', ''), best_og_desc.replace('\n', ' ')])
            
            for h in ANNOTATIONS_HEADER:
                if h in annotations:
                    annot_columns.append(','.join(sorted(annotations[h])))
                else:
                    annot_columns.append('-')

            all_annotations.append(annot_columns)

        return
    

    ##
    def parse_hits(self, filename):

        def _parse_hits():
            for line in open(filename, 'r'):
                if line.startswith('#') or not line.strip():
                    continue

                line = list(map(str.strip, line.split('\t')))
                # query, target, evalue, score
                if len(line) == 4: # short hits
                    hit = [line[0], line[1], float(line[2]), float(line[3])]
                elif len(line) == 11:
                    hit = [line[0], line[1], float(line[2]), float(line[3]),
                           int(line[4]), int(line[5]), int(line[6]), int(line[7]),
                           float(line[8]), float(line[9]), float(line[10])]

                yield hit
                            
        return _parse_hits
    
    ##
    def iter_hit_lines(self, hits_gen_func, store_hits = True):

        if store_hits: self.hits = []
        
        for hit in hits_gen_func():
            
            if store_hits: self.hits.append(hit)
            
            yield_tuple = (hit, self.annot, self.seed_ortholog_score, self.seed_ortholog_evalue,
                           self.tax_scope_mode, self.tax_scope_id, self.TAXONOMIC_RESOLUTION,
                           self.target_taxa, self.target_orthologs, self.excluded_taxa,
                           self.go_evidence, self.go_excluded, self.pfam_transfer)
            
            yield yield_tuple
            
        return


def md5_seqs(fasta_file):
    from hashlib import md5
    md5_queries = {}
        
    for name, seq in iter_fasta_seqs(fasta_file):
        md5_seq = md5(seq.encode('utf-8')).hexdigest()
        md5_queries[name] = md5_seq
            
    return md5_queries

## END
