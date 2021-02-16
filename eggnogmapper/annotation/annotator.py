##
## CPCantalapiedra 2019

from collections import Counter
import sys
import time
import multiprocessing
    
from ..emapperException import EmapperException
from ..common import get_call_info, TAX_SCOPE_AUTO, TAX_SCOPE_AUTO_BROAD
from ..utils import colorify
from ..vars import LEVEL_NAMES, LEVEL_DEPTH
from ..search.hmmer.hmmer_seqio import iter_fasta_seqs

from . import annota
from . import db_sqlite
from . import orthologs as ortho
from . import output
from .pfam.pfam_modes import run_pfam_mode, PFAM_TRANSFER_NARROWEST_OG, PFAM_TRANSFER_SEED_ORTHOLOG, PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO
from .ncbitaxa.ncbiquery import NCBITaxa

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


    ##
    def annotate(self, seed_orthologs_file, annot_file, orthologs_file, pfam_file):

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
            output_orthologs_f = output.output_orthologs_closure(ORTHOLOGS_OUT, NCBITaxa())
            output_annotations_f = output.output_annotations_closure(ANNOTATIONS_OUT, md5_field, md5_queries)

            ##
            # Obtain annotations
            qn, elapsed_time = self._annotate(seed_orthologs_file, pfam_file, output_orthologs_f, output_annotations_f)
            db_sqlite.close()

            ##
            # Output footer and close files
            if self.report_orthologs == True:
                output.output_orthologs_footer(ORTHOLOGS_OUT, self.no_file_comments, qn, elapsed_time)
                ORTHOLOGS_OUT.close()

            if self.annot == True:
                output.output_annotations_footer(ANNOTATIONS_OUT, self.no_file_comments, qn, elapsed_time)
                ANNOTATIONS_OUT.close()

            print(colorify("Functional annotation of refined hits starts now", 'green'))
            
        return

    ##
    def _annotate(self, seed_orthologs_file, pfam_file, output_ortho_f, output_annot_f):
        
        start_time = time.time()
        
        if self.dbmem == True:
            all_orthologs, all_annotations, qn = self._annotate_dbmem(seed_orthologs_file, pfam_file)
        else:
            all_orthologs, all_annotations, qn = self._annotate_ondisk(seed_orthologs_file, pfam_file)

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
            
        return qn, elapsed_time


    ##
    def _annotate_dbmem(self, seed_orthologs_file, pfam_file):
        all_orthologs = {}
        all_annotations = []

        start_time = time.time() # do not take into account time to load the db into memory
        db_sqlite.connect(usemem = True)
        total_time = time.time() - start_time
        print(colorify(f"Time to load the DB into memory: {total_time}", "lblue"), file=sys.stderr)
        sys.stderr.flush()        

        start_time = time.time() # do not take into account time to load the db into memory
        
        qn = 0
        try:
            for result in map(annotate_hit_line, self.iter_hit_lines(seed_orthologs_file)):
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
            db_sqlite.close()

        elapsed_time = time.time() - start_time
        print(colorify(f" All queries processed. Time to perform queries:{elapsed_time} rate:{(float(qn) / elapsed_time):.2f} q/s", 'lblue'))
                    
        return all_orthologs, all_annotations, qn

    
    ##
    def _annotate_ondisk(self, seed_orthologs_file, pfam_file):

        all_orthologs = {}
        all_annotations = []
        
        # multiprocessing.set_start_method("spawn")
        pool = multiprocessing.Pool(self.cpu)

        start_time = time.time() # do not take into account time to load the pool of processes
        
        qn = 0
        try:
            for result in pool.imap(annotate_hit_line_process, self.iter_hit_lines(seed_orthologs_file)):
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
            pool.terminate()

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
    def iter_hit_lines(self, filename):
        
        for line in open(filename, 'r'):
            if line.startswith('#') or not line.strip():
                continue
            
            yield_tuple = (line, self.annot, self.seed_ortholog_score, self.seed_ortholog_evalue,
                           self.tax_scope_mode, self.tax_scope_id, self.TAXONOMIC_RESOLUTION,
                           self.target_taxa, self.target_orthologs, self.excluded_taxa,
                           self.go_evidence, self.go_excluded, self.pfam_transfer)
            
            yield yield_tuple
            
        return


##
def annotate_hit_line_process(arguments):
    # should connect also if no previous connection
    # exists in this Pool process (worker)
    db_sqlite.connect()

    return annotate_hit_line(arguments)

# annotate_hit_line is outside the class because must be pickable
##
def annotate_hit_line(arguments):
    
    line, annot, seed_ortholog_score, seed_ortholog_evalue, \
        tax_scope_mode, tax_scope_id, tax_resolution, \
        target_taxa, target_orthologs, excluded_taxa, \
        go_evidence, go_excluded, \
        pfam_transfer = arguments
    
    try:
        if not line.strip() or line.startswith('#'):
            return None
        
        ##
        # Split fields of search results
        r = list(map(str.strip, line.split('\t')))

        query_name = r[0]
        best_hit_name = r[1]
        best_hit_evalue = float(r[2])
        best_hit_score = float(r[3])
        
        ##
        # Filter by empty hit, error, evalue and/or score
        if filter_out(best_hit_name, best_hit_evalue, best_hit_score, seed_ortholog_evalue, seed_ortholog_score):
            return None
        
        ##
        # Retrieve OGs (orthologs groups) the hit belongs to
        match_nogs = get_member_ogs(best_hit_name)
        if not match_nogs:
            return None
            
        ##
        # Obtain names of OGs, narrowest OG, and the best OG according to tax_scope
        match_nogs_names, narr_og_id, narr_og_level, best_og_id, best_og_level = parse_nogs(match_nogs, tax_scope_mode, tax_scope_id, tax_resolution)
            
        annot_levels = set()
        annot_levels.add(best_og_level)
        
        if best_og_id is None:
            annotations = None
            all_orthologies = None
            orthologs = None
            best_og_name = "-"
            best_og_cat = "-"
            best_og_desc = "-"
            narr_og_name = "-"
            narr_og_cat = "-"
            narr_og_desc = "-"
            
        else:
            best_og_name = f"{best_og_id}@{best_og_level}|{LEVEL_NAMES.get(best_og_level, best_og_level)}"
            best_og_cat, best_og_desc = get_og_description(best_og_id, best_og_level)
            
            narr_og_name = f"{narr_og_id}@{narr_og_level}|{LEVEL_NAMES.get(narr_og_level, narr_og_level)}"
            narr_og_cat, narr_og_desc = get_og_description(narr_og_id, narr_og_level)

            ##
            # Normalize target_taxa if any
            if target_taxa is not None:
                target_taxa = normalize_target_taxa(target_taxa)
            else:
                target_taxa = None

            if excluded_taxa is not None:
                excluded_taxa = normalize_target_taxa(excluded_taxa)
            else:
                excluded_taxa = None

            ##
            # Retrieve co-orthologs of seed ortholog
            # annot_levels are used to restrict the speciation events retrieved
            # target_taxa are used to restrict the species from which to retrieve co-ortholog proteins
            try:
                all_orthologies, best_OG = ortho.get_member_orthologs(best_hit_name, annot_levels, match_nogs_names)
                if best_OG is not None:
                    best_og_name = best_OG
                    best_og_id = best_OG.split("|")[0].split("@")[0]
                    best_og_level = best_OG.split("|")[0].split("@")[1]
                    if best_og_id == "seed_ortholog":
                        best_og_cat = "-"
                        best_og_desc = "-"
                    else:
                        best_og_cat, best_og_desc = get_og_description(best_og_id, best_og_level)

            except Exception as e:
                # import traceback
                # traceback.print_exc()
                raise e
            else:
                # filter co-orthologs to keep only target_orthologs: "all", "one2one", ...
                orthologs = _filter_orthologs(all_orthologies, target_orthologs, target_taxa, excluded_taxa)

            ##
            # Retrieve annotations of co-orthologs
            if annot == True and orthologs is not None and len(orthologs) > 0:

                annotations = annota.summarize_annotations(orthologs,
                                                           annotations_fields = ANNOTATIONS_HEADER,
                                                           target_go_ev = go_evidence,
                                                           excluded_go_ev = go_excluded)

                if pfam_transfer == PFAM_TRANSFER_NARROWEST_OG:
                    if best_og_level == narr_og_level:
                        narr_orthologies = all_orthologies
                    else:
                        narr_annot_levels = set()
                        narr_annot_levels.add(narr_og_level)
                        narr_orthologies, _ = ortho.get_member_orthologs(best_hit_name, narr_annot_levels, match_nogs_names)

                    # filter co-orthologs to keep only target_orthologs: "all", "one2one", ...
                    narr_orthologs = _filter_orthologs(narr_orthologies, target_orthologs, target_taxa, excluded_taxa)

                    pfam_annotations = db_sqlite.get_pfam_annotations(','.join(['"%s"' % n for n in narr_orthologs]))
                    if pfam_annotations is not None and len(pfam_annotations) > 0:
                        annotations["PFAMs"] = Counter()
                        for pfam_annotation in pfam_annotations:
                            annotations["PFAMs"].update([str(x).strip() for x in pfam_annotation[0].split(",")])
                    else:
                        annotations["PFAMs"] = Counter()

                elif pfam_transfer == PFAM_TRANSFER_SEED_ORTHOLOG:
                    pfam_annotations = db_sqlite.get_pfam_annotations('"'+best_hit_name+'"')
                    if pfam_annotations is not None and len(pfam_annotations) > 0:
                        pfam_annotations = Counter(list(pfam_annotations[0][0].split(",")))
                        annotations["PFAMs"] = pfam_annotations
                    else:
                        annotations["PFAMs"] = Counter()                    
                else: # pfam_transfer == PFAM_TRANSFER_BEST_OG
                    pass

            else:
                annotations = {}

    except Exception as e:
        raise EmapperException(f"Error: annotation went wrong for line \"{line.strip()}\". "+str(e))

    # WARNING: do NOT close db connections, because it becomes super slow
    # finally:
    #     db_sqlite.close()
    
    return (query_name, best_hit_name, best_hit_evalue, best_hit_score,
            annotations,
            narr_og_name, narr_og_cat, narr_og_desc,
            best_og_name, best_og_cat, best_og_desc,
            match_nogs_names, all_orthologies, orthologs)


def _filter_orthologs(all_orthologies, target_orthologs, target_taxa, excluded_taxa):
    orthologs = sorted(all_orthologies[target_orthologs])
    if excluded_taxa is not None:
        orthologs = [o for o in orthologs if int(o.split(".")[0]) not in excluded_taxa]
    if target_taxa is not None:
        orthologs = [o for o in orthologs if int(o.split(".")[0]) in target_taxa]
    return orthologs
            
##
def parse_nogs(match_nogs, tax_scope_mode, tax_scope_id, tax_resolution):        
    match_nogs_names = []
    best_og_id = None
    best_og_level = None
    best_og_depth = None
    narr_og_id = None
    narr_og_level = None
    narr_og_depth = None

    lvl_depths = set(LEVEL_DEPTH.keys())
    
    tax_scope_id_2 = tax_scope_id
    
    for nog in sorted(match_nogs, key=lambda x: LEVEL_DEPTH[x.split("@")[1]]):
        nog_id = nog.split("@")[0]
        nog_tax_id = nog.split("@")[1]

        nog_name = f"{nog}|{LEVEL_NAMES.get(nog_tax_id, nog_tax_id)}"
        match_nogs_names.append(nog_name)

        nog_depth = LEVEL_DEPTH[nog_tax_id]

        # Obtain narrowest OG
        if narr_og_depth is None or nog_depth >= narr_og_depth:
            narr_og_id = nog_id
            narr_og_level = nog_tax_id
            narr_og_depth = nog_depth
            if tax_scope_id is None and tax_scope_mode == "narrowest":
                best_og_id = narr_og_id
                best_og_level = narr_og_level
                best_og_depth = narr_og_depth

        # Obtain best OG based on tax scope
        if tax_scope_id is None:
            if tax_scope_mode in {"auto", "auto_broad"}:
                for filter_tax_id in tax_resolution:
                    if filter_tax_id == nog_tax_id:
                        best_og_id = nog_id
                        best_og_level = nog_tax_id
                        best_og_depth = nog_depth
                        break

            elif tax_scope_mode == "narrowest":
                pass # Already processed above

            else:
                raise EmapperException(f"Unrecognized tax scope mode {tax_scope_mode}")    

        else: # tax_scope_id is not None:
            for i, filter_tax_id in enumerate(tax_scope_id_2):
                if filter_tax_id == nog_tax_id:
                    best_og_id = nog_id
                    best_og_level = nog_tax_id
                    best_og_depth = nog_depth
                    # only leave IDs of equal or more priority (left-most in the list)
                    tax_scope_id_2 = tax_scope_id_2[:i+1]
                    break

    if best_og_id is None:
        if tax_scope_mode == "none":
            pass
        elif tax_scope_mode == "narrowest":
            best_og_id = narr_og_id
            best_og_level = narr_og_level
            best_og_depth = narr_og_depth
        elif tax_scope_mode in {"auto", "auto_broad"}:
            for nog in sorted(match_nogs, key=lambda x: LEVEL_DEPTH[x.split("@")[1]]):
                nog_tax_id = nog.split("@")[1]
                for filter_tax_id in tax_resolution:
                    if filter_tax_id == nog_tax_id:
                        best_og_id = nog_id
                        best_og_level = nog_tax_id
                        best_og_depth = nog_depth
                        break                
        else:
            raise EmapperException(f"Error. Unrecognized tax scope mode {tax_scope_mode}.")
    
    # print(match_nogs_names)
    # print(f"Best OG: {best_og_id}-{best_og_level}")

    return match_nogs_names, narr_og_id, narr_og_level, best_og_id, best_og_level

        
##
def filter_out(hit_name, hit_evalue, hit_score, threshold_evalue, threshold_score):
    """
    Filter hit if ERROR, by score or by evalue
    """
    if hit_name == '-' or hit_name == 'ERROR':
        return True
    
    if hit_score < threshold_score or hit_evalue > threshold_evalue:
        return True
    
    return False

##
def normalize_target_taxa(target_taxa):
    """
    Receives a list of taxa IDs and/or taxa names and returns a set of expanded taxids numbers
    """
    ncbi = NCBITaxa()
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


def get_member_ogs(name):
    match = db_sqlite.get_member_ogs(name)
    ogs = None
    if match:
        ogs = [str(x).strip() for x in match[0].split(',')]
    return ogs


def get_og_description(og, level):
    best = ['-', '-', '-']
    
    for og, nm, desc, cat in db_sqlite.get_ogs_description(og, level):
        desc = desc.strip()
        if desc and desc != 'N/A' and desc != 'NA':
            best = [nm, cat, desc]
            break
    
    return best[1], best[2]


def md5_seqs(fasta_file):
    from hashlib import md5
    md5_queries = {}
        
    for name, seq in iter_fasta_seqs(fasta_file):
        md5_seq = md5(seq.encode('utf-8')).hexdigest()
        md5_queries[name] = md5_seq
            
    return md5_queries

## END
