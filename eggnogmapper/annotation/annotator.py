##
## CPCantalapiedra 2019

import sys
import time
import multiprocessing

from ..common import get_version, get_db_version, get_call_info, pexists, TAXONOMIC_RESOLUTION
from ..emapperException import EmapperException
from ..utils import colorify
from ..vars import LEVEL_PARENTS, LEVEL_NAMES, LEVEL_DEPTH

from . import annota

HIT_HEADER = ["#query_name",
              "seed_eggNOG_ortholog",
              "seed_ortholog_evalue",
              "seed_ortholog_score",
              "best_tax_level", ]

ANNOTATIONS_HEADER = ['Preferred_name', 'GOs', 'EC',
                      'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass',
                      'BRITE', 'KEGG_TC', 'CAZy', 'BiGG_Reaction']

HIT_OG_HEADER = ["taxonomic scope", "eggNOG OGs", "best eggNOG OG",
                 "COG Functional cat.", "eggNOG free text desc."]

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
    #     return self.annotate_hits_file(seed_orthologs_file, annot_file, hmm_hits_file)

    # ##
    # def annotate_hits_file(self, seed_orthologs_file, annot_file, hmm_hits_file):

        start_time = time.time()
        
        seq2bestOG = {}
        seq2annotOG = {}
        if pexists(hmm_hits_file):
            seq2bestOG = get_seq_hmm_matches(hmm_hits_file)
            seq2annotOG = annota.get_ogs_annotations(set([v[0] for v in seq2bestOG.values()]))

        print(colorify("Functional annotation of refined hits starts now", 'green'))

        OUT = open(annot_file, "w")

        if self.report_orthologs:
            ORTHOLOGS = open(annot_file+".orthologs", "w")

        if not self.no_file_comments:
            print(get_call_info(), file=OUT)
            print('\t'.join(HIT_HEADER + ANNOTATIONS_HEADER + HIT_OG_HEADER), file=OUT)
            
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
                 annotations, annot_level_max, swallowest_level, match_nogs, orthologs) = result
                if query_name in seq2bestOG:
                    (hitname, evalue, score, qlength, hmmfrom, hmmto, seqfrom,
                     seqto, q_coverage) = seq2bestOG[query_name]
                    bestOG = '%s|%s|%s' %(hitname, evalue, score)
                    og_cat, og_desc = seq2annotOG.get(hitname, ['', ''])
                else:
                    bestOG = 'NA|NA|NA'
                    og_cat, og_desc = annota.get_best_og_description(match_nogs)

                if self.report_orthologs:
                    print('\t'.join(map(str, (query_name, ','.join(orthologs)))), file=ORTHOLOGS)

                # prepare annotations for printing
                annot_columns = [query_name,
                                 best_hit_name,
                                 str(best_hit_evalue),
                                 str(best_hit_score),
                                 LEVEL_NAMES[swallowest_level]]

                for h in ANNOTATIONS_HEADER:
                    if h in annotations:
                        annot_columns.append(','.join(sorted(annotations[h])))
                    else:
                        annot_columns.append('')

                annot_columns.extend([annot_level_max,
                                        ','.join(match_nogs),
                                        bestOG,
                                        og_cat.replace('\n', ''),
                                        og_desc.replace('\n', ' ')])

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
    for nog in match_nogs:
        match_levels.update(LEVEL_PARENTS[nog.split("@")[1]])

    swallowest_level = sorted(match_levels & set(LEVEL_DEPTH.keys()),
                              key=lambda x: LEVEL_DEPTH[x], reverse=True)[0]

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
        target_taxa = orthology.normalize_target_taxa(target_taxa)
    else:
        target_taxa = None

    try:
        all_orthologies = annota.get_member_orthologs(best_hit_name, target_taxa=target_taxa, target_levels=annot_levels)
    except Exception:
        orthologs = None
        status = 'Error'
    else:
        orthologs = sorted(all_orthologies[target_orthologs])
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
            annotations, annot_level_max, swallowest_level, match_nogs, orthologs)

## END
