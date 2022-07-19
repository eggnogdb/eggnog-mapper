##
## CPCantalapiedra 2019

from os.path import isfile as pisfile
import sys
import time
import multiprocessing
from collections import defaultdict, Counter

from ..emapperException import EmapperException
from ..common import get_call_info, get_data_path, ITYPE_PROTS, ITYPE_CDS
from ..utils import colorify
from ..search.hmmer.hmmer_seqio import iter_fasta_seqs

from .annotator_worker_novel_fams import annotate_hit_line
from . import output_novel_fams

ANNOTATIONS_HEADER = output_novel_fams.ANNOTATIONS_HEADER

PFAM_COL = -1 # position of PFAMs annotations in list of annotations

##
class AnnotatorNovelFams:
    
    annot = report_orthologs = None
    
    dbmem = None

    no_file_comments = cpu = None

    # options for pfam hmmpgmd searches
    num_servers = num_workers = cpus_per_worker = port = end_port = None

    seed_ortholog_score = seed_ortholog_evalue = None
    tax_scope_mode = tax_scope_ids = target_taxa = target_orthologs = excluded_taxa = None
    
    go_evidence = go_excluded = None
    pfam_realign = trans_table = temp_dir = None
    md5 = None

    resume = None
        
    ##
    def __init__(self, args, annot, excel, report_orthologs):

        self.annot = annot
        self.report_orthologs = report_orthologs
        self.excel = excel

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
        
        self.trans_table = args.trans_table
        self.itype = args.itype
        
        self.temp_dir = args.temp_dir
        
        self.md5 = args.md5

        self.resume = args.resume
        
        return

    
    ##
    def annotate(self, hits_gen_func, annot_file, excel_file, orthologs_file, pfam_file, queries_file):

        annots_generator = None
        # ncbi = None
        
        try:
            if self.report_orthologs == True or self.annot == True:            
                ##
                # md5 hashes
                if self.md5 == True:
                    print(colorify("Creating md5 hashes of input sequences", 'green'))
                    
                    if self.itype == ITYPE_PROTS:
                        translate = False
                    elif self.itype == ITYPE_CDS:
                        translate = True
                        
                    md5_queries = md5_seqs(queries_file, translate, self.trans_table)
                else:
                    md5_queries = None

                md5_field = (self.md5 == True and md5_queries is not None)

                ##
                # target_taxa and excluded_taxa not available yet for novel fams
                
                ##
                # Annotations
                
                # If resume, create generator of previous annotations
                annots_parser = None
                if self.resume == True:
                    annots_parser = parse_annotations(self.annot, annot_file)
                
                # Obtain annotations
                annots_generator = self._annotate(hits_gen_func, annots_parser)

                ##
                # PFAM realign
                # no pfam mode for novel fams yet
                
                ##
                # Output
                
                if self.annot == True:
                    annots_generator = output_novel_fams.output_annotations(annots_generator,
                                                                            annot_file,
                                                                            self.resume,
                                                                            self.no_file_comments,
                                                                            md5_field,
                                                                            md5_queries)

                    # excel output not available yet for novel_fams
                # no orthologs report for novel fams

                # unpack the annotations removing the "exists" or "skip"
                # boolean used when --resume
                
                annots_generator = unpack_annotations(annots_generator)

        finally:
            pass
            # if ncbi is not None: ncbi.close()
            
        return annots_generator

    
    ##
    def _annotate(self, hits_gen_func, annots_parser):
        
        pool = multiprocessing.Pool(self.cpu)
        chunk_size = 1
        # "my recommendation is targeting 10 ms chunk processing time"
        # https://stackoverflow.com/a/43817408/2361653
        # As our tasks take no less than 0.1 secs, a large chunk_size makes no sense at all
        # Note that this makes q/s an approximation until all tasks have been finished
        
        try:
            for result in pool.imap(annotate_hit_line,
                                    self.iter_hit_lines(hits_gen_func, annots_parser),
                                    chunk_size):
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
    def iter_hit_lines(self, hits_gen_func, annots_parser):

        curr_annot = None
        if annots_parser is not None:
            try:
                curr_annot = next(annots_parser)
            except StopIteration:
                curr_annot = None
        
        for hit in hits_gen_func:

            annotation = None
            
            if curr_annot is not None:
                if hit[0] == curr_annot[0]:
                    annotation = curr_annot
                    try:
                        curr_annot = next(annots_parser)
                    except StopIteration:
                        curr_annot = None
            
            yield_tuple = (hit, self.annot, self.seed_ortholog_score, self.seed_ortholog_evalue,
                           self.tax_scope_mode, self.tax_scope_ids,
                           self.target_taxa, self.target_orthologs, self.excluded_taxa,
                           self.go_evidence, self.go_excluded, get_data_path(), annotation)
            
            yield yield_tuple
            
        return

##
def unpack_annotations(annots):
    for (hit, annotation), skip in annots:
        yield hit, annotation

##
def parse_annotations(annot, annot_file):
    
    if annot == True:
        with open(annot_file, 'r') as annot_f:
            for line in annot_f:
                if line.startswith("#"): continue
                
                hit, annotation = parse_annotation_line(line)
                
                # this assumes the annotated hit it is also present
                # in orthologs_file
                yield annotation

    # no report of orthologs for novel fams
                
    else:
        pass # no annotations then
    
    return

def parse_annotation_line(line):
    hit = None
    annotation = None

    data = list(map(str.strip, line.split("\t")))

    query_name = data[0]
    best_hit_name = data[1]
    best_hit_evalue = float(data[2])
    best_hit_score = float(data[3])
    hit = [query_name, best_hit_name, best_hit_evalue, best_hit_score]
    novel_fam = data[4]
    
    annotation = (query_name, best_hit_name, best_hit_evalue, best_hit_score, novel_fam)
    
    return hit, annotation


##
def md5_seqs(fasta_file, translate, trans_table):
    from hashlib import md5
    md5_queries = {}
        
    for name, seq in iter_fasta_seqs(fasta_file, translate=translate, trans_table=trans_table):
        md5_seq = md5(seq.encode('utf-8')).hexdigest()
        md5_queries[name] = md5_seq
            
    return md5_queries

## END
