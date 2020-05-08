##
## CPCantalapiedra 2019

import sys
import time
import multiprocessing

from ..common import get_version, get_db_version, get_call_info, pexists, TAXONOMIC_RESOLUTION, cleanup_og_name
from ..emapperException import EmapperException
from ..utils import colorify
from ..vars import LEVEL_PARENTS, LEVEL_NAMES, LEVEL_DEPTH

# from ..orthologs.orthology import normalize_target_taxa
from . import annota
from .orthologs import get_member_orthologs

HIT_HEADER = ["#query_name",
              "seed_eggNOG_ortholog",
              "seed_ortholog_evalue",
              "seed_ortholog_score",
              "taxonomic_scope",
              "best_tax_level",
              "COG Functional cat.",
              "eggNOG free text desc.",              
              "eggNOG OGs"]

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
    tax_scope = target_taxa = target_orthologs = excluded_taxa = None
    go_evidence = go_excluded = None

    ##
    def __init__(self, args, annot, report_orthologs):

        self.annot = annot
        self.report_orthologs = report_orthologs
        
        self.no_file_comments = args.no_file_comments
        self.cpu = args.cpu
        self.seed_ortholog_score = args.seed_ortholog_score
        self.seed_ortholog_evalue = args.seed_ortholog_evalue
        self.tax_scope = args.tax_scope
        self.target_taxa = args.target_taxa
        self.target_orthologs = args.target_orthologs
        self.excluded_taxa = args.excluded_taxa
        self.go_evidence = args.go_evidence
        self.go_excluded = args.go_excluded
        
        return

    ##
    def annotate(self, seed_orthologs_file, annot_file):
        
        print(colorify("Functional annotation of refined hits starts now", 'green'))
        
        all_orthologs, all_annotations, qn, elapsed_time = self._annotate(seed_orthologs_file, annot_file)

        # Output orthologs
        if self.report_orthologs:
            ORTHOLOGS = open(annot_file+".orthologs", "w")
            for (query_name, orthologs) in all_orthologs:
                print('\t'.join(map(str, (query_name, ','.join(orthologs)))), file=ORTHOLOGS)
            ORTHOLOGS.close()

        # Output annotations
        if self.annot:
            OUT = open(annot_file, "w")

            if not self.no_file_comments:
                print(get_call_info(), file=OUT)
                print('\t'.join(HIT_HEADER + ANNOTATIONS_HEADER), file=OUT)

            for annot_columns in all_annotations:
                print('\t'.join(annot_columns), file=OUT)

            if not self.no_file_comments:
                print('# %d queries scanned' % (qn), file=OUT)
                print('# Total time (seconds):', elapsed_time, file=OUT)
                print('# Rate:', "%0.2f q/s" % ((float(qn) / elapsed_time)), file=OUT)

            OUT.close()

        return

    ##
    def _annotate(self, seed_orthologs_file, annot_file):

        all_orthologs = []
        all_annotations = []
        
        start_time = time.time()
        
        pool = multiprocessing.Pool(self.cpu)

        qn = 0
        for result in pool.imap(annotate_hit_line, self.iter_hit_lines(seed_orthologs_file)):
            qn += 1
            if qn and (qn % 500 == 0):
                total_time = time.time() - start_time
                print(f"{pq} {total_time} {(float(qn) / total_time):.2f} q/s (func. annotation)", file=sys.stderr)
                # print(qn, total_time, "%0.2f q/s (func. annotation)" % ((float(qn) / total_time)), file=sys.stderr)
                sys.stderr.flush()

            if result:
                (query_name, best_hit_name, best_hit_evalue, best_hit_score,
                 annotations, annot_level_max, swallowest_level,
                 og_cat, og_desc, match_nogs_names, orthologs) = result

                if self.report_orthologs:
                    all_orthologs.append((query_name, orthologs))
                    
                if self.annot:
                    # prepare annotations for printing
                    annot_columns = [query_name, best_hit_name, str(best_hit_evalue), str(best_hit_score),
                                     annot_level_max, swallowest_level, og_cat.replace('\n', ''), og_desc.replace('\n', ' '),
                                     ",".join(match_nogs_names)]
                
                    for h in ANNOTATIONS_HEADER:
                        if h in annotations:
                            annot_columns.append(','.join(sorted(annotations[h])))
                        else:
                            annot_columns.append('')

                    all_annotations.append(annot_columns)

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
                   self.tax_scope, self.target_taxa, self.target_orthologs, self.excluded_taxa,
                   self.go_evidence, self.go_excluded)
            
            yield yield_tuple
            
        return

# annotate_hit_line is outside the class because must be pickable
##
def annotate_hit_line(arguments):
    try:
        return _annotate_hit_line(arguments)
    except:
        import traceback
        traceback.print_exc(file=sys.stdout)
        raise

    return

##
def _annotate_hit_line(arguments):

    # should connect also if no previous connection
    # exists in this Pool process (worker)
    annota.connect()

    line, seed_ortholog_score, seed_ortholog_evalue, tax_scope, target_taxa, target_orthologs, excluded_taxa, go_evidence, go_excluded = arguments

    if not line.strip() or line.startswith('#'):
        return None
    
    r = list(map(str.strip, line.split('\t')))

    query_name = r[0]
    
    best_hit_name = r[1]
    if best_hit_name == '-' or best_hit_name == 'ERROR':
        return None

    best_hit_evalue = float(r[2])
    best_hit_score = float(r[3])
    if best_hit_score < seed_ortholog_score or best_hit_evalue > seed_ortholog_evalue:
        return None

    match_nogs = annota.get_member_ogs(best_hit_name)
    if not match_nogs:
        return None

    match_levels = set()
    swallowest_og = None
    swallowest_level = None
    lvl_depths = set(LEVEL_DEPTH.keys())
    
    for nog in match_nogs:
        nog_lvls = LEVEL_PARENTS[nog.split("@")[1]]
        match_levels.update(nog_lvls)

        # detect swallowest OG
        nog_lvl = sorted(set(nog_lvls) & set(lvl_depths), key=lambda x: LEVEL_DEPTH[x], reverse=True)[0]
        nog_depth = LEVEL_DEPTH[nog_lvl]
        if swallowest_level is None or nog_depth > swallowest_depth:
            swallowest_depth = nog_depth
            swallowest_level = nog_lvl
            swallowest_og = nog.split("@")[0]

    swallowest_level = LEVEL_NAMES.get(swallowest_level, swallowest_level)

    og_cat, og_desc = annota.get_deepest_og_description(swallowest_og)

    match_nogs_names = [nog+"|"+LEVEL_NAMES.get(nog.split("@")[1], nog.split("@")[1]) for nog in
                        sorted(match_nogs, key=lambda x: LEVEL_DEPTH[x.split("@")[1]])]
    
    annot_levels = set()
    if tax_scope == "auto":
        for level in TAXONOMIC_RESOLUTION:
            if level in match_levels:
                annot_levels.add(level)
                annot_level_max = LEVEL_NAMES.get(level, level)
                break
    else:
        annot_levels.add(tax_scope)
        annot_level_max = LEVEL_NAMES.get(tax_scope, tax_scope)
    
    if target_taxa != 'all':
        target_taxa = normalize_target_taxa(target_taxa)
    else:
        target_taxa = None

    try:
        all_orthologies = get_member_orthologs(best_hit_name, target_taxa=target_taxa, target_levels=annot_levels)
        
    except Exception as e:
        # print(str(e))
        orthologs = None
        status = 'Error'
    else:
        orthologs = sorted(all_orthologies[target_orthologs])
        
        # print("annotator: orthologs")
        # print(orthologs)
        if excluded_taxa:
            orthologs = [o for o in orthologs if not o.startswith("%s." % excluded_taxa)]
        status = 'OK'

    if orthologs:
        annotations = annota.summarize_annotations(orthologs,
                                                   target_go_ev=go_evidence,
                                                   excluded_go_ev=go_excluded)
    else:
        annotations = {}

    annota.close()
    
    return (query_name, best_hit_name, best_hit_evalue, best_hit_score,
            annotations, annot_level_max, swallowest_level,
            og_cat, og_desc, match_nogs_names, orthologs)


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

## END
