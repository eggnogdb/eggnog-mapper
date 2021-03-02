##
## CPCantalapiedra 2019

import sys
import time
import multiprocessing
from collections import defaultdict

from ..emapperException import EmapperException
from ..common import get_call_info
from ..utils import colorify
from ..search.hmmer.hmmer_seqio import iter_fasta_seqs

from .db_sqlite import get_eggnog_db
from .ncbitaxa.ncbiquery import get_ncbi
from .pfam.pfam_modes import run_pfam_mode, PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO

from .annotator_worker import annotate_hit_line_mem, annotate_hit_line_ondisk
from . import output

ANNOTATIONS_HEADER = output.ANNOTATIONS_HEADER

PFAM_COL = -1 # position of PFAMs annotations in list of annotations

##
class Annotator:

    # hits = hits_dict = None
    # annotations = annotations_dict = None
    # orthologs = None
    
    annot = report_orthologs = None
    
    dbmem = None

    no_file_comments = cpu = None

    # options for pfam hmmpgmd searches
    num_servers = num_workers = cpus_per_worker = port = end_port = None

    seed_ortholog_score = seed_ortholog_evalue = None
    tax_scope_mode = tax_scope_ids = target_taxa = target_orthologs = excluded_taxa = None
    
    go_evidence = go_excluded = None
    pfam_realign = queries_fasta = translate = temp_dir = None
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
        self.tax_scope_ids = args.tax_scope_ids
                
        self.target_taxa = args.target_taxa
        self.target_orthologs = args.target_orthologs
        self.excluded_taxa = args.excluded_taxa
            
        self.go_evidence = args.go_evidence
        self.go_excluded = args.go_excluded
        
        self.pfam_realign = args.pfam_realign
        
        self.queries_fasta = args.input
        self.translate = args.translate
        self.temp_dir = args.temp_dir
        
        self.md5 = args.md5
        
        return

    # #
    # def get_hits(self):
    #     return self.hits

    # def get_hits_dict(self):
    #     hits_dict = None
    #     if self.hits_dict is not None:
    #         hits_dict = self.hits_dict
    #     elif self.hits is not None:
    #         hits_dict = {hit[0]:hit for hit in self.hits}
    #         self.hits_dict = hits_dict
    #     # else: None
    #     return hits_dict
    
    # #
    # def get_annotations(self):
    #     return self.annotations
    
    # def get_annotations_dict(self):
    #     annotations_dict = None
    #     if self.annotations_dict is not None:
    #         annotations_dict = self.annotations_dict
    #     elif self.annotations is not None:
    #         annotations_dict = {annot[0]:annot for annot in self.annotations}
    #         self.annotations_dict = annotations_dict
    #     # else: None
    #     return annotations_dict

    # #
    # def get_orthologs(self):
    #     return self.orthologs
    
    ##
    def annotate(self, hits_gen_func, annot_file, orthologs_file, pfam_file):

        ncbi = None
        try:
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
                # Prepare taxa restrictions
                # target_taxa are used to restrict the species from which to retrieve co-ortholog proteins
                # the opposite is excluded_taxa
                # In both cases, we need to normalize the list of taxa, if any
                
                if self.target_taxa is not None:
                    ncbi = get_ncbi(usemem = False)
                    self.target_taxa = normalize_target_taxa(self.target_taxa, ncbi)
                    
                if self.excluded_taxa is not None:
                    ncbi = get_ncbi(usemem = False)
                    self.excluded_taxa = normalize_target_taxa(self.excluded_taxa, ncbi)

                if ncbi is not None: ncbi.close() # close it, so that for orthologs we can load it into memory
                
                ##
                # Annotations

                # Obtain annotations
                annots_generator = self._annotate(hits_gen_func)

                ##
                # PFAM realign
                # Note that this needs all the annotations at once,
                # and therefore breaks the generators pipeline
                if (self.annot == True and
                    self.pfam_realign in [PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO] and
                    annots_generator is not None):
                    
                    annots_generator = run_pfam_mode(self.pfam_realign, annots_generator,
                                                     self.queries_fasta, self.translate,
                                                     self.cpu, self.num_servers,
                                                     self.num_workers, self.cpus_per_worker,
                                                     self.port, self.end_port,
                                                     self.temp_dir, pfam_file)
                
                ##
                # Output
                
                if self.annot == True:
                    annots_generator = output.output_annotations(annots_generator,
                                                                 annot_file,
                                                                 self.no_file_comments,
                                                                 md5_field,
                                                                 md5_queries)
                    
                if self.report_orthologs == True:
                    annots_generator = output.output_orthologs(annots_generator,
                                                               orthologs_file,
                                                               self.no_file_comments)

        finally:
            if ncbi is not None: ncbi.close()
            
        return annots_generator


    def _annotate(self, hits_gen_func):
        if self.dbmem == True:
            annots_generator = self._annotate_dbmem(hits_gen_func)
        else:
            annots_generator = self._annotate_ondisk(hits_gen_func)
        return annots_generator


    ##
    def _annotate_dbmem(self, hits_gen_func):
        try:
            ##
            # Load sqlite DBs into memory
            start_time = time.time()
            eggnog_db = get_eggnog_db(usemem = True)
            total_time = time.time() - start_time
            print(colorify(f"Time to load the DB into memory: {total_time}", "lblue"), file=sys.stderr)
            sys.stderr.flush()

            ##
            # Annotate hits
            for result in map(annotate_hit_line_mem, self.iter_hit_lines(hits_gen_func)):
                yield result

        except EmapperException:
            raise
        except Exception as e:
            # import traceback
            # traceback.print_exc()
            raise EmapperException(f"Error: annotation went wrong. "+str(e))
        finally:
            eggnog_db.close()
                    
        return

    
    ##
    def _annotate_ondisk(self, hits_gen_func):
        
        pool = multiprocessing.Pool(self.cpu)
        chunk_size = 1
        # "my recommendation is targeting 10 ms chunk processing time"
        # https://stackoverflow.com/a/43817408/2361653
        # As our tasks take no less than 0.1 secs, a large chunk_size makes no sense at all
        # Note that this makes q/s an approximation until all tasks have been finished
        
        try:
            for result in pool.imap(annotate_hit_line_ondisk, self.iter_hit_lines(hits_gen_func), chunk_size):
                yield result

        except EmapperException:
            raise
        except Exception as e:
            import traceback
            traceback.print_exc()
            raise EmapperException(f"Error: annotation failed. "+str(e))
        finally:
            pool.close()
            pool.terminate() # it should remove the global eggnog_db variables also
            
        return
    

    ##
    def parse_hits(self, filename):
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
        return
    
    ##
    def iter_hit_lines(self, hits_gen_func):
        
        for hit in hits_gen_func:
            
            yield_tuple = (hit, self.annot, self.seed_ortholog_score, self.seed_ortholog_evalue,
                           self.tax_scope_mode, self.tax_scope_ids,
                           self.target_taxa, self.target_orthologs, self.excluded_taxa,
                           self.go_evidence, self.go_excluded)
            
            yield yield_tuple
            
        return

##
def normalize_target_taxa(target_taxa, ncbi):
    """
    Receives a list of taxa IDs and/or taxa names and returns a set of expanded taxids numbers
    """
    expanded_taxa = set()
    
    for taxon in target_taxa:
        taxid = ""
        try:
            taxid = int(taxon)
        except ValueError:
            taxid = ncbi.get_name_translator([taxon])[taxon][0]
        else:
            taxon = ncbi.get_taxid_translator([taxid])[taxid]

        if taxid is not None:
            species = ncbi.get_descendant_taxa(taxid, intermediate_nodes = True)
            for sp in species:
                expanded_taxa.add(sp)

    return expanded_taxa


def md5_seqs(fasta_file):
    from hashlib import md5
    md5_queries = {}
        
    for name, seq in iter_fasta_seqs(fasta_file):
        md5_seq = md5(seq.encode('utf-8')).hexdigest()
        md5_queries[name] = md5_seq
            
    return md5_queries

## END
