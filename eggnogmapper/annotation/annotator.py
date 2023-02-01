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

from .db_sqlite import get_eggnog_db
from .ncbitaxa.ncbiquery import get_ncbi
from .pfam.pfam_modes import run_pfam_mode, PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO

from .annotator_worker import annotate_hit_line_mem, annotate_hit_line_ondisk
from . import output

ANNOTATIONS_HEADER = output.ANNOTATIONS_HEADER

PFAM_COL = -1 # position of PFAMs annotations in list of annotations

##
class Annotator:
    
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
        ncbi = None
        
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

                # If resume, create generator of previous annotations
                annots_parser = None
                if self.resume == True:
                    annots_parser = parse_annotations(self.annot, annot_file,
                                                      self.report_orthologs, orthologs_file)
                
                # Obtain annotations
                annots_generator = self._annotate(hits_gen_func, annots_parser)

                ##
                # PFAM realign
                # Note that this needs all the annotations at once,
                # and therefore breaks the generators pipeline
                if (self.annot == True and
                    self.pfam_realign in [PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO] and
                    annots_generator is not None):

                    if self.itype == ITYPE_PROTS:
                        translate = False
                    elif self.itype == ITYPE_CDS:
                        translate = True
                        
                    annots_generator = run_pfam_mode(self.pfam_realign, annots_generator,
                                                     queries_file, self.resume,
                                                     translate, self.trans_table,
                                                     self.cpu, self.num_servers,
                                                     self.num_workers, self.cpus_per_worker,
                                                     self.port, self.end_port,
                                                     self.temp_dir, pfam_file)
                
                ##
                # Output
                
                if self.annot == True:
                    annots_generator = output.output_annotations(annots_generator,
                                                                 annot_file,
                                                                 self.resume,
                                                                 self.no_file_comments,
                                                                 md5_field,
                                                                 md5_queries)

                    if self.excel == True:
                        annots_generator = output.output_excel(annots_generator,
                                                               excel_file,
                                                               self.resume,
                                                               self.no_file_comments,
                                                               md5_field,
                                                               md5_queries)
                        
                    
                if self.report_orthologs == True:
                    annots_generator = output.output_orthologs(annots_generator,
                                                               orthologs_file,
                                                               self.resume,
                                                               self.no_file_comments)

                # unpack the annotations removing the "exists" or "skip"
                # boolean used when --resume
                
                annots_generator = unpack_annotations(annots_generator)

        finally:
            if ncbi is not None: ncbi.close()
            
        return annots_generator


    def _annotate(self, hits_gen_func, annots_parser):
        print(colorify(f"Functional annotation of hits...", "lgreen"), file=sys.stderr)
        if self.dbmem == True:
            annots_generator = self._annotate_dbmem(hits_gen_func, annots_parser)
        else:
            annots_generator = self._annotate_ondisk(hits_gen_func, annots_parser)
        return annots_generator


    ##
    def _annotate_dbmem(self, hits_gen_func, annots_parser):
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
            for result in map(annotate_hit_line_mem,
                              self.iter_hit_lines(hits_gen_func, annots_parser)):
                yield result

        except EmapperException:
            raise
        except Exception as e:
            import traceback
            traceback.print_exc()
            raise EmapperException(f"Error: annotation went wrong. "+str(e))
        finally:
            eggnog_db.close()
                    
        return

    
    ##
    def _annotate_ondisk(self, hits_gen_func, annots_parser):
        
        pool = multiprocessing.Pool(self.cpu)
        chunk_size = 1
        # "my recommendation is targeting 10 ms chunk processing time"
        # https://stackoverflow.com/a/43817408/2361653
        # As our tasks take no less than 0.1 secs, a large chunk_size makes no sense at all
        # Note that this makes q/s an approximation until all tasks have been finished
        
        try:
            for result in pool.imap(annotate_hit_line_ondisk,
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
def parse_annotations(annot, annot_file, report_orthologs, orthologs_file):
    
    if annot == True:
        with open(annot_file, 'r') as annot_f:
            for line in annot_f:
                if line.startswith("#"): continue
                
                hit, annotation = parse_annotation_line(line)
                
                # this assumes the annotated hit it is also present
                # in orthologs_file
                yield annotation
            
    elif report_orthologs == True:
        
        prev_query = None
        with open(orthologs_file, 'r') as orth_f:
            for line in orth_f:
                if line.startswith("#"): continue

                # just yield that query already exists in file
                query_name = line[0]
                if prev_query is not None and query_name != prev_query:
                    annotation = (prev_query, None, None, None,
                                  None, None, None, None, None)
                    yield annotation

                prev_query = query_name
                
            # last query
            if prev_query is not None:
                annotation = (prev_query, None, None, None,
                              None, None, None, None, None)
                yield annotation
                
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
    
    annotations = defaultdict(Counter)
    
    for i, field in enumerate(data[8:]):
        if i < len(ANNOTATIONS_HEADER):
            field_name = ANNOTATIONS_HEADER[i]
            annotations[field_name] = field.split(",")
        
    og_cat_desc = ("-", data[5], data[6])

    max_annot_lvl = data[7]
    
    match_nog_names = data[4].split(",")
    
    all_orthologies = None
    annot_orthologs = None
    
    annotation = (query_name, best_hit_name, best_hit_evalue, best_hit_score,
                  annotations,
                  og_cat_desc, max_annot_lvl, match_nog_names,
                  all_orthologies, annot_orthologs)
    
    return hit, annotation


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


def md5_seqs(fasta_file, translate, trans_table):
    from hashlib import md5
    md5_queries = {}
        
    for name, seq in iter_fasta_seqs(fasta_file, translate=translate, trans_table=trans_table):
        md5_seq = md5(seq.encode('utf-8')).hexdigest()
        md5_queries[name] = md5_seq
            
    return md5_queries

## END
