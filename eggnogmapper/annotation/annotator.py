##
## CPCantalapiedra 2019

import os, sys
import time
import multiprocessing
import subprocess
from tempfile import NamedTemporaryFile
    
from ..emapperException import EmapperException
from ..common import get_call_info, TAXONOMIC_RESOLUTION, get_pfam_db, HMMFETCH
from ..utils import colorify
from ..vars import LEVEL_PARENTS, LEVEL_NAMES, LEVEL_DEPTH

# from ..orthologs.orthology import normalize_target_taxa
from . import annota
from . import db_sqlite
from . import orthologs as ortho
from .pfam import PfamAligner, get_pfam_args, parse_pfam_file, group_queries_pfams, get_hmmsearch_args, parse_hmmsearch_file

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
                      'BiGG_Reaction',
                      'PFAMs']

PFAM_COL = -1 # position of PFAMs annotations in list of annotations

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
    pfam = queries_fasta = translate = temp_dir = None
    
    ##
    def __init__(self, args, annot, report_orthologs):

        self.annot = annot
        self.report_orthologs = report_orthologs
        
        self.no_file_comments = args.no_file_comments
        self.cpu = args.cpu
        self.seed_ortholog_score = args.seed_ortholog_score
        self.seed_ortholog_evalue = args.seed_ortholog_evalue

        self.tax_scope_mode = args.tax_scope_mode
        self.tax_scope_id = args.tax_scope_id
                
        self.target_taxa = args.target_taxa
        self.target_orthologs = args.target_orthologs
        self.excluded_taxa = args.excluded_taxa
        self.go_evidence = args.go_evidence
        self.go_excluded = args.go_excluded
        
        self.pfam = args.pfam
        self.queries_fasta = args.input
        self.translate = args.translate
        self.temp_dir = args.temp_dir
        
        return

    
    ##
    def annotate(self, seed_orthologs_file, annot_file, orthologs_file, pfam_file):
        
        print(colorify("Functional annotation of refined hits starts now", 'green'))
        
        all_orthologs, all_annotations, qn, elapsed_time = self._annotate(seed_orthologs_file, pfam_file)

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
                
                # If tax_scope_mode == narrowest there is no need to output best_og colums
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
    def _annotate(self, seed_orthologs_file, pfam_file):

        all_orthologs = []
        all_annotations = []
        
        start_time = time.time()

        # multiprocessing.set_start_method("spawn")
        pool = multiprocessing.Pool(self.cpu)

        qn = 0
        try:
            for result in pool.imap(annotate_hit_line, self.iter_hit_lines(seed_orthologs_file)):
                qn += 1
                if qn and (qn % 500 == 0):
                    total_time = time.time() - start_time
                    print(f"{qn} {total_time} {(float(qn) / total_time):.2f} q/s (func. annotation)", file=sys.stderr)
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

                        # If tax_scope_mode == narrowest there is no need to output best_og colums
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

        ##
        # PFAMs annotation
        
        aligned_pfams = None
        if self.pfam == 'denovo':
            print(colorify("de novo search of PFAM domains", 'green'))
            # filter fasta file to have only annotated queries
            queries = [annot_columns[0] for annot_columns in all_annotations]
            fasta_file = filter_fasta_file(queries, self.queries_fasta)

            # align those queries to whole PFAM to carry out a de novo annotation
            pfam_args, infile = get_pfam_args(self.cpu, fasta_file.name, self.translate, self.temp_dir)
            pfam_aligner = PfamAligner(pfam_args)
            pfam_aligner.align_whole_pfam(infile, pfam_file, silent = True)
            aligned_pfams = parse_pfam_file(pfam_file)

        elif self.pfam == 'align':
            print(colorify("re-aligning PFAM domains from orthologs to queries", 'green'))

            aligned_pfams = self.pfam_align_serial(all_annotations, PFAM_COL, pfam_file)

        elif self.pfam == 'alignp':
            print(colorify("re-aligning PFAM domains from orthologs to queries, in parallel mode", 'green'))

            aligned_pfams = self.pfam_align_parallel(all_annotations, PFAM_COL, pfam_file)    
            
        if aligned_pfams is not None:
            for annot_columns in all_annotations:
                query_name = annot_columns[0]
                if query_name in aligned_pfams:
                    annot_columns[PFAM_COL] = ",".join(sorted(aligned_pfams[query_name]))
                else:
                    annot_columns[PFAM_COL] = "-"

            elapsed_time = time.time() - start_time

            print(colorify(f" Processed queries:{qn} total_time:{elapsed_time} rate:{(float(qn) / elapsed_time):.2f} q/s", 'lblue'))
        
        return all_orthologs, all_annotations, qn, elapsed_time
    
    ##
    def iter_hit_lines(self, filename):
        
        for line in open(filename, 'r'):
            if line.startswith('#') or not line.strip():
                continue
            
            yield_tuple = (line, self.annot, self.seed_ortholog_score, self.seed_ortholog_evalue,
                           self.tax_scope_mode, self.tax_scope_id, self.target_taxa, self.target_orthologs, self.excluded_taxa,
                           self.go_evidence, self.go_excluded)
            
            yield yield_tuple
            
        return

    ##
    def pfam_align_serial(self, annotations, PFAM_COL, pfam_file):
        aligned_pfams = None

        for queries_pfams_group in group_queries_pfams(annotations, PFAM_COL):
            # fetch pfams to a new HMM file
            # fetch queries to a new Fasta file

            arguments = queries_pfams_group, self.queries_fasta, get_pfam_db(), self.temp_dir, self.translate, pfam_file, self.cpu
            alignments = query_pfam_annotate(arguments)
            if alignments is not None:
                if aligned_pfams is None:
                    aligned_pfams = alignments
                else:
                    aligned_pfams.update(alignments)

        return aligned_pfams


    ##
    def wrap_group_queries_pfams(self, annotations, PFAM_COL, pfam_file):
        queries_pfams_groups = group_queries_pfams(annotations, PFAM_COL)

        for queries_pfams_group in group_queries_pfams(annotations, PFAM_COL):
            yield (queries_pfams_group, self.queries_fasta, get_pfam_db(), self.temp_dir, self.translate, pfam_file, 1)
            # cpu = 1 since we are already parallelizing
            
        return

    ##
    def pfam_align_parallel(self, annotations, PFAM_COL, pfam_file):
        aligned_pfams = None

        pool = multiprocessing.Pool(self.cpu)
        
        try:
            
            for alignments in pool.imap(query_pfam_annotate, self.wrap_group_queries_pfams(annotations, PFAM_COL, pfam_file)):
                if alignments is not None:
                    if aligned_pfams is None:
                        aligned_pfams = alignments
                    else:
                        aligned_pfams.update(alignments)                
                
        except EmapperException:
            raise
        except Exception as e:
            raise EmapperException(f"Error: annotation went wrong for pfam alignment in parallel. "+str(e))
        finally:
            pool.terminate()

        return aligned_pfams
    

##
def query_pfam_annotate(arguments):
    queries_pfams_group, queries_fasta, pfam_db, temp_dir, translate, pfam_file, cpu = arguments
    
    aligned_pfams = None

    fasta_file, hmm_file = filter_fasta_hmm_files(queries_pfams_group, queries_fasta, pfam_db)

    if fasta_file is None or hmm_file is None:
        pass
    else:
        # output file for this group
        P = NamedTemporaryFile(mode='w')

        # align queries to the new HMM file
        pfam_args, infile = get_hmmsearch_args(cpu, fasta_file.name, hmm_file.name, translate, temp_dir)
        pfam_aligner = PfamAligner(pfam_args)
        pfam_aligner.align_whole_pfam(infile, P.name, silent = True)

        aligned_pfams = parse_hmmsearch_file(P.name)

        # Append contents of output file for this group into pfam_file,
        # which is the file reporting all the pfam hits together
        with open(pfam_file, 'a') as pfamf:
            pfamf.write(open(P.name, 'r').read())

        P.close()

    if fasta_file is not None:
        fasta_file.close()
    if hmm_file is not None:
        hmm_file.close()

    return aligned_pfams

##
# Create a new temp fasta file including only
# the queries in the list "queries"
def filter_fasta_file(queries, orig_fasta):
    
    Q = NamedTemporaryFile(mode='w')
    with open(orig_fasta, 'r') as origf:
        found = False
        for line in origf:
            if line.startswith(">"):
                orig_query = line.strip()[1:].split(" ")[0]
                if orig_query in queries:
                    found = True
                    print(f">{orig_query}", file=Q)
                else:
                    found = False
            else:
                if found == True:
                    print(f"{line.strip()}", file=Q)
    Q.flush()
    
    return Q

##
def filter_fasta_hmm_files(queries_pfams, orig_fasta, orig_hmm):
    
    new_fasta = new_hmm = None
    
    queries = queries_pfams[0]
    pfams = queries_pfams[1]
    
    # Process fasta queries
    Q = filter_fasta_file(queries, orig_fasta)
    
    # Process pfams
    P = NamedTemporaryFile(mode='w')
    for pfam in pfams:
        print(pfam, file=P)
    P.flush()

    H = None
    if os.stat(P.name).st_size > 0:
        H = NamedTemporaryFile(mode='w')
        cmd = f"{HMMFETCH} -f {orig_hmm} {P.name}"
        cp = subprocess.run(cmd, shell=True, stdout=H, stderr=subprocess.DEVNULL)
        if os.stat(H.name).st_size == 0:
            H = None
    
    return Q, H
    
    
# annotate_hit_line is outside the class because must be pickable
##
def annotate_hit_line(arguments):
    # should connect also if no previous connection
    # exists in this Pool process (worker)
    db_sqlite.connect()
    
    line, annot, seed_ortholog_score, seed_ortholog_evalue, \
        tax_scope_mode, tax_scope_id, \
        target_taxa, target_orthologs, excluded_taxa, \
        go_evidence, go_excluded = arguments
    
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
        match_nogs_names, narr_og_id, narr_og_level, best_og_id, best_og_level = parse_nogs(match_nogs, tax_scope_mode, tax_scope_id)
        
        annot_levels = set()
        annot_levels.add(best_og_level)
        
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
        if annot == True and orthologs is not None and len(orthologs) > 0:
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
            if tax_scope_mode == "auto":
                for filter_tax_id in TAXONOMIC_RESOLUTION:
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
        elif tax_scope_mode == "auto":
            for nog in sorted(match_nogs, key=lambda x: LEVEL_DEPTH[x.split("@")[1]]):
                nog_tax_id = nog.split("@")[1]
                for filter_tax_id in TAXONOMIC_RESOLUTION:
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
