##
## CPCantalapiedra 2019

import sys
import time
import multiprocessing

from ..common import get_version, get_db_version, get_call_info, pexists, TAXONOMIC_RESOLUTION, cleanup_og_name
from ..emapperException import EmapperException
from ..utils import colorify
from ..vars import LEVEL_PARENTS, LEVEL_NAMES, LEVEL_DEPTH

from ..orthologs.orthology import normalize_target_taxa
from . import annota

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
def get_annotator(args):
    annotator = None

    annotator = Annotator(args)
    
    return annotator

##
class Annotator:

    no_file_comments = report_orthologs = None

    cpu = None

    seed_ortholog_score = seed_ortholog_evalue = None
    tax_scope = target_taxa = target_orthologs = excluded_taxa = None
    go_evidence = go_excluded = None

    ##
    def __init__(self, args):

        self.no_file_comments = args.no_file_comments
        self.report_orthologs = args.report_orthologs
        self.cpu = args.cpu
        self.seed_ortholog_score = args.seed_ortholog_score
        self.seed_ortholog_evalue = args.seed_ortholog_evalue
        self.tax_scope = args.tax_scope
        self.target_taxa = args.target_taxa
        self.target_orthologs = args.target_orthologs
        self.excluded_taxa = args.excluded_taxa
        self.go_evidence = args.go_evidence
        self.go_excluded = args.go_excluded
        
        annota.connect()
        
        return

    ##
    def annotate(self, seed_orthologs_file, annot_file, hmm_hits_file):

        start_time = time.time()
        
        seq2bestOG = {}
        seq2annotOG = {}
        if pexists(hmm_hits_file):
            seq2bestOG = self.get_seq_hmm_matches(hmm_hits_file)
            seq2annotOG = annota.get_ogs_annotations(set([v[0] for v in seq2bestOG.values()]))

        print(colorify("Functional annotation of refined hits starts now", 'green'))

        OUT = open(annot_file, "w")

        if self.report_orthologs:
            ORTHOLOGS = open(annot_file+".orthologs", "w")

        if not self.no_file_comments:
            print(get_call_info(), file=OUT)
            print('\t'.join(HIT_HEADER + ANNOTATIONS_HEADER), file=OUT)
            
        qn = 0
        pool = multiprocessing.Pool(self.cpu)
        # annotate_hit_line is outside the class because must be pickable
        for result in pool.imap(annotate_hit_line, self.iter_hit_lines(seed_orthologs_file)):
            qn += 1
            if qn and (qn % 500 == 0):
                total_time = time.time() - start_time
                print(qn, total_time, "%0.2f q/s (func. annotation)" % (
                    (float(qn) / total_time)), file=sys.stderr)
                sys.stderr.flush()

            if result:
                (query_name, best_hit_name, best_hit_evalue, best_hit_score,
                 annotations, annot_level_max, swallowest_og, swallowest_level, match_nogs, orthologs) = result

                if self.report_orthologs:
                    print('\t'.join(map(str, (query_name, ','.join(orthologs)))), file=ORTHOLOGS)

                # prepare annotations for printing
                annot_columns = [query_name,
                                 best_hit_name,
                                 str(best_hit_evalue),
                                 str(best_hit_score)]

                annot_columns.append(annot_level_max)
                
                annot_columns.append(swallowest_level)
                
                og_cat, og_desc = annota.get_deeper_og_description(swallowest_og)
                
                annot_columns.extend([og_cat.replace('\n', ''),
                                      og_desc.replace('\n', ' ')])

                match_nogs_names = [nog+"|"+LEVEL_NAMES.get(nog.split("@")[1], nog.split("@")[1]) for nog in
                                    sorted(match_nogs, key=lambda x: LEVEL_DEPTH[x.split("@")[1]])]
                
                annot_columns.append(",".join(match_nogs_names))
                
                for h in ANNOTATIONS_HEADER:
                    if h in annotations:
                        annot_columns.append(','.join(sorted(annotations[h])))
                    else:
                        annot_columns.append('')
                        
                print('\t'.join(annot_columns), file=OUT)

            #OUT.flush()

        pool.terminate()

        elapsed_time = time.time() - start_time
        if not self.no_file_comments:
            print('# %d queries scanned' % (qn), file=OUT)
            print('# Total time (seconds):', elapsed_time, file=OUT)
            print('# Rate:', "%0.2f q/s" % ((float(qn) / elapsed_time)), file=OUT)
        OUT.close()

        if self.report_orthologs:
            ORTHOLOGS.close()

        print(colorify(" Processed queries:%s total_time:%s rate:%s" %\
                       (qn, elapsed_time, "%0.2f q/s" % ((float(qn) / elapsed_time))), 'lblue'))
        
        return

    ##
    def get_seq_hmm_matches(self, hits_file):
        # annota.connect()
        print(colorify("Reading HMM matches", 'green'))
        seq2oginfo = {}
        start_time = time.time()
        hitnames = set()
        if pexists(hits_file):
            for line in open(hits_file):
                
                if not line.strip() or line.startswith('#'):
                    continue
                
                (query, hit, evalue, sum_score, query_length, hmmfrom, hmmto,
                 seqfrom, seqto, q_coverage) = map(str.strip, line.split('\t'))

                if query not in seq2oginfo and hit not in ['ERROR', '-']:
                    hitname = cleanup_og_name(hit)
                    seq2oginfo[query] = [hitname, evalue, sum_score, query_length,
                                         hmmfrom, hmmto, seqfrom, seqto,
                                         q_coverage]
        
        return seq2oginfo
    
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
        all_orthologies = annota.get_member_orthologs(best_hit_name, target_taxa=target_taxa, target_levels=annot_levels)
        
        # print("annotator: all_orthologies")
        # print(all_orthologies)
    except Exception:
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

    return (query_name, best_hit_name, best_hit_evalue, best_hit_score,
            annotations, annot_level_max, swallowest_og, swallowest_level, match_nogs, orthologs)

## END
