##
## CPCantalapiedra 2019

import sys
import time
import multiprocessing

from ..emapperException import EmapperException
from ..common import get_call_info, TAXONOMIC_RESOLUTION
from ..utils import colorify
from ..vars import LEVEL_PARENTS, LEVEL_NAMES, LEVEL_DEPTH

# from ..orthologs.orthology import normalize_target_taxa
from . import annota
from . import db_sqlite
from . import orthologs as ortho

HIT_HEADER = ["#query_name",
              "seed_eggNOG_ortholog",
              "seed_ortholog_evalue",
              "seed_ortholog_score",
              "eggNOG OGs",
              "narr_og_name",
              "narr_og_cat",
              "narr_og_desc"]

BEST_OG_HEADER = ["best_og_name",
                  "best_og_cat",
                  "best_og_desc"]

ANNOTATIONS_HEADER = ['Preferred_name',
                      'GOs',
                      'EC',
                      'KEGG_ko',
                      'KEGG_Pathway',
                      'KEGG_Module',
                      'KEGG_Reaction',
                      'KEGG_rclass',
                      'BRITE',
                      'KEGG_TC',
                      'CAZy',
                      'BiGG_Reaction']

##
def get_annotator(args, annot, report_orthologs):
    annotator = None

    annotator = Annotator(args, annot, report_orthologs)
    
    return annotator

##
class Annotator:

    annot = report_orthologs = None

    no_file_comments = cpu = None

    seed_ortholog_score = seed_ortholog_evalue = None
    tax_scope_mode = tax_scope_id = target_taxa = target_orthologs = excluded_taxa = None
    go_evidence = go_excluded = None
    
    ##
    def __init__(self, args, annot, report_orthologs):

        self.annot = annot
        self.report_orthologs = report_orthologs
        
        self.no_file_comments = args.no_file_comments
        self.cpu = args.cpu
        self.seed_ortholog_score = args.seed_ortholog_score
        self.seed_ortholog_evalue = args.seed_ortholog_evalue

        self.tax_scope_mode, self.tax_scope_id = self.__parse_tax_scope(args.tax_scope)
                
        self.target_taxa = args.target_taxa
        self.target_orthologs = args.target_orthologs
        self.excluded_taxa = args.excluded_taxa
        self.go_evidence = args.go_evidence
        self.go_excluded = args.go_excluded
        
        return


    ##
    # Parses tax_scope command line argument
    # to define tax_scope_mode and tax_scope_id (one or more tax IDs)
    def __parse_tax_scope(self, tax_scope):
        tax_scope_mode = None
        tax_scope_id = None

        tax_scope_fields = tax_scope.strip().split(",")
        tax_scope_mode = tax_scope_fields[0]

        # Auto
        if tax_scope_mode == "auto":
            # if len(tax_scope_fields) > 1:
            #     tax_scope_mode = "auto_custom"
            #     tax_scope_id = tax_scope_fields[1:]
            # else:
            tax_scope_id = None

        # Narrowest
        elif tax_scope_mode == "narrowest":
            tax_scope_id = None

        # Tax IDs
        else:
            # Only the specified tax ID
            if len(tax_scope_fields) == 1:
                tax_scope_mode = "none"
                tax_scope_id = [tax_scope_fields[0]]

            # Tax ID lis, with or without mode for those not found in the list
            elif len(tax_scope_fields) > 1:
                last_pos = tax_scope_fields[-1]
                if last_pos in ["narrowest", "none"]:
                    tax_scope_mode = last_pos
                    tax_scope_id = tax_scope_fields[:-1]
                else:
                    tax_scope_mode = "none"
                    tax_scope_id = tax_scope_fields
            else:
                raise EmapperException(f"Error: unrecognized tax scope format {tax_scope}.")

        return tax_scope_mode, tax_scope_id

    
    ##
    def annotate(self, seed_orthologs_file, annot_file, orthologs_file):
        
        print(colorify("Functional annotation of refined hits starts now", 'green'))
        
        all_orthologs, all_annotations, qn, elapsed_time = self._annotate(seed_orthologs_file)

        # Output orthologs
        if self.report_orthologs:
            ORTHOLOGS = open(orthologs_file, "w")
            for (query_name, orthologs) in all_orthologs:
                print('\t'.join(map(str, (query_name, ','.join(orthologs)))), file=ORTHOLOGS)
            ORTHOLOGS.close()

        # Output annotations
        if self.annot:
            OUT = open(annot_file, "w")

            if not self.no_file_comments:
                print(get_call_info(), file=OUT)
                print('\t'.join(HIT_HEADER), end="\t", file=OUT)
                if self.tax_scope_id is not None or self.tax_scope_mode != "narrowest":
                    print('\t'.join(BEST_OG_HEADER), end="\t", file=OUT)                    
                print('\t'.join(ANNOTATIONS_HEADER), file=OUT)

            for annot_columns in all_annotations:
                print('\t'.join(annot_columns), file=OUT)

            if not self.no_file_comments:
                print('# %d queries scanned' % (qn), file=OUT)
                print('# Total time (seconds):', elapsed_time, file=OUT)
                print('# Rate:', "%0.2f q/s" % ((float(qn) / elapsed_time)), file=OUT)

            OUT.close()

        return

    ##
    def _annotate(self, seed_orthologs_file):

        all_orthologs = []
        all_annotations = []
        
        start_time = time.time()
        
        pool = multiprocessing.Pool(self.cpu)

        qn = 0
        try:
            for result in pool.imap(annotate_hit_line, self.iter_hit_lines(seed_orthologs_file)):
                qn += 1
                if qn and (qn % 500 == 0):
                    total_time = time.time() - start_time
                    print(f"{pq} {total_time} {(float(qn) / total_time):.2f} q/s (func. annotation)", file=sys.stderr)
                    # print(qn, total_time, "%0.2f q/s (func. annotation)" % ((float(qn) / total_time)), file=sys.stderr)
                    sys.stderr.flush()

                if result:
                    (query_name, best_hit_name, best_hit_evalue, best_hit_score,
                     annotations, 
                     narr_og_name, narr_og_cat, narr_og_desc,
                     best_og_name, best_og_cat, best_og_desc,                     
                     match_nogs_names, orthologs) = result

                    if self.report_orthologs:
                        all_orthologs.append((query_name, orthologs))

                    if self.annot:
                        # prepare annotations for printing
                        annot_columns = [query_name, best_hit_name, str(best_hit_evalue), str(best_hit_score),
                                         ",".join(match_nogs_names), 
                                         narr_og_name, narr_og_cat.replace('\n', ''), narr_og_desc.replace('\n', ' ')]

                        if self.tax_scope_id is not None or self.tax_scope_mode != "narrowest":
                            annot_columns.extend([best_og_name, best_og_cat.replace('\n', ''), best_og_desc.replace('\n', ' ')])

                        for h in ANNOTATIONS_HEADER:
                            if h in annotations:
                                annot_columns.append(','.join(sorted(annotations[h])))
                            else:
                                annot_columns.append('-')

                        all_annotations.append(annot_columns)

        except EmapperException:
            raise
        except Exception as e:
            raise EmapperException(f"Error: annotation went wrong for query number {qn}. "+str(e))
        finally:
            pool.terminate()

        elapsed_time = time.time() - start_time

        print(colorify(f" Processed queries:{qn} total_time:{elapsed_time} rate:{(float(qn) / elapsed_time):.2f} q/s", 'lblue'))
        
        return all_orthologs, all_annotations, qn, elapsed_time
    
    ##
    def iter_hit_lines(self, filename):
        
        for line in open(filename, 'r'):
            if line.startswith('#') or not line.strip():
                continue
            
            yield_tuple = (line, self.seed_ortholog_score, self.seed_ortholog_evalue,
                           self.tax_scope_mode, self.tax_scope_id, self.target_taxa, self.target_orthologs, self.excluded_taxa,
                           self.go_evidence, self.go_excluded)
            
            yield yield_tuple
            
        return

    
# annotate_hit_line is outside the class because must be pickable
##
def annotate_hit_line(arguments):

    # should connect also if no previous connection
    # exists in this Pool process (worker)
    db_sqlite.connect()

    line, seed_ortholog_score, seed_ortholog_evalue, tax_scope_mode, tax_scope_id, target_taxa, target_orthologs, excluded_taxa, go_evidence, go_excluded = arguments
    
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
        # Obtain names of OGs and the best OG according to tax_scope
        match_nogs_names, narr_og_id, narr_og_level, best_og_id, best_og_level = parse_nogs(match_nogs, tax_scope_mode, tax_scope_id)
        
        annot_levels = set()
        annot_levels.add(best_og_level)
                
        # swallowest_level = f"{best_og_id}@{best_og_level}|{LEVEL_NAMES.get(best_og_level, best_og_level)}"
        if best_og_id is None:
            best_og_name = "-"
            best_og_cat = "-"
            best_og_desc = "-"
        else:
            best_og_name = f"{best_og_id}@{best_og_level}|{LEVEL_NAMES.get(best_og_level, best_og_level)}"
            best_og_cat, best_og_desc = get_og_description(best_og_id)

        narr_og_name = f"{narr_og_id}@{narr_og_level}|{LEVEL_NAMES.get(narr_og_level, narr_og_level)}"
        narr_og_cat, narr_og_desc = get_og_description(narr_og_id)

        ##
        # Normalize target_taxa if any
        if target_taxa != 'all':
            target_taxa = normalize_target_taxa(target_taxa)
        else:
            target_taxa = None

        ##
        # Retrieve co-orthologs of seed ortholog
        # annot_levels are used to restrict the speciation events retrieved
        # target_taxa are used to restrict the species from which to retrieve co-ortholog proteins
        try:
            all_orthologies = ortho.get_member_orthologs(best_hit_name, target_taxa=target_taxa, target_levels=annot_levels)

        except Exception as e:
            # print(str(e))
            orthologs = None
            status = 'Error'
        else:
            # filter co-orthologs to keep only target_orthologs: "all", "one2one", ...
            orthologs = sorted(all_orthologies[target_orthologs])
            if excluded_taxa:
                orthologs = [o for o in orthologs if not o.startswith("%s." % excluded_taxa)]
            status = 'OK'

        ##
        # Retrieve annotations of co-orthologs
        if orthologs:
            annotations = annota.summarize_annotations(orthologs,
                                                       annotations_fields = ANNOTATIONS_HEADER,
                                                       target_go_ev = go_evidence,
                                                       excluded_go_ev = go_excluded)
        else:
            annotations = {}

    except Exception as e:
        raise EmapperException(f"Error: annotation went wrong for line \"{line.strip()}\". "+str(e))
    
    finally:
        db_sqlite.close()
    
    return (query_name, best_hit_name, best_hit_evalue, best_hit_score,
            annotations,
            narr_og_name, narr_og_cat, narr_og_desc,
            best_og_name, best_og_cat, best_og_desc,
            match_nogs_names, orthologs)


##
def parse_nogs(match_nogs, tax_scope_mode, tax_scope_id):        
    match_nogs_names = []
    best_og_id = None
    best_og_level = None
    best_og_depth = None
    narr_og_id = None
    narr_og_level = None
    narr_og_depth = None

    lvl_depths = set(LEVEL_DEPTH.keys())
    
    # tax_scope_id_2 = [x for x in tax_scope_id]
    tax_scope_id_2 = tax_scope_id
    
    for nog in sorted(match_nogs, key=lambda x: LEVEL_DEPTH[x.split("@")[1]]):
        nog_id = nog.split("@")[0]
        nog_level = nog.split("@")[1]

        nog_name = f"{nog}|{LEVEL_NAMES.get(nog_level, nog_level)}"
        match_nogs_names.append(nog_name)

        nog_depth = LEVEL_DEPTH[nog_level]

        # Obtain narrowest OG
        if narr_og_depth is None or nog_depth >= narr_og_depth:
            narr_og_id = nog_id
            narr_og_level = nog_level
            narr_og_depth = nog_depth
            if tax_scope_id is None and tax_scope_mode == "narrowest":
                best_og_id = narr_og_id
                best_og_level = narr_og_level
                best_og_depth = narr_og_depth

        # Obtain best OG based on tax scope
        if tax_scope_id is None:
            if tax_scope_mode == "auto":
                for level in TAXONOMIC_RESOLUTION:
                    if level == nog_level:
                        best_og_id = nog_id
                        best_og_level = nog_level
                        best_og_depth = nog_depth
                        break

            elif tax_scope_mode == "narrowest":
                pass # Already processed above

            else:
                raise EmapperException(f"Unrecognized tax scope mode {tax_scope_mode}")    

        else: # tax_scope_id is not None:
            for i, level in enumerate(tax_scope_id_2):
                if level == nog_level:
                    best_og_id = nog_id
                    best_og_level = nog_level
                    best_og_depth = nog_depth
                    tax_scope_id_2 = tax_scope_id_2[:i+1]
                    break

    if best_og_id is None:
        if tax_scope_mode == "none":
            pass
        elif tax_scope_mode == "narrowest":
            best_og_id = narr_og_id
            best_og_level = narr_og_level
            best_og_depth = narr_og_depth
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
    from ete3 import NCBITaxa
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

        species = ncbi.get_descendant_taxa(taxid, collapse_subspecies=False)
        for sp in species:
            expanded_taxa.add(sp)

    return expanded_taxa


def get_member_ogs(name):
    match = db_sqlite.get_member_ogs(name)
    ogs = None
    if match:
        ogs = [str(x).strip() for x in match[0].split(',')]
    return ogs


def get_og_description(og):
    best = [None, '', '']
    
    for og, nm, desc, cat in db_sqlite.get_ogs_description(og):
        desc = desc.strip()
        if desc and desc != 'N/A' and desc != 'NA':
            best = [nm, cat, desc]
            break
    
    return best[1], best[2]


## END
