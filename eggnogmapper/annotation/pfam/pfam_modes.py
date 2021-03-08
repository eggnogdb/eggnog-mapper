##
## CPCantalapiedra 2020

from collections import Counter

from ...emapperException import EmapperException
from ...utils import colorify

from .pfam_denovo import pfam_align_denovo
from .pfam_scan import pfam_align_parallel_scan

PFAM_REALIGN_NONE = 'none'
PFAM_REALIGN_REALIGN = 'realign'
PFAM_REALIGN_DENOVO = 'denovo'

def run_pfam_mode(pfam_search_mode, annots_generator, queries_fasta, resume, translate,
                  cpu, num_servers, num_workers, cpus_per_worker, port, end_port,
                  temp_dir, pfam_file):

    ##
    # 1) Align queries to PFAMs
    
    aligned_pfams = None
    all_annotations = None
    
    if pfam_search_mode == PFAM_REALIGN_DENOVO:
        print(colorify("De novo scan of PFAM domains", 'lgreen'))

        all_annotations, queries_pfams = load_all_annotations(annots_generator)

        aligned_pfams = pfam_align_denovo(queries_pfams,
                                          queries_fasta,
                                          resume,
                                          translate,
                                          cpu,
                                          num_servers,
                                          num_workers,
                                          cpus_per_worker,
                                          port,
                                          end_port,
                                          temp_dir,
                                          pfam_file)

    elif pfam_search_mode == PFAM_REALIGN_REALIGN:
        print(colorify("Re-aligning queries to PFAM domains from orthologs", 'lgreen'))

        all_annotations, queries_pfams = load_all_annotations(annots_generator)
        
        aligned_pfams = pfam_align_parallel_scan(queries_pfams,
                                                 queries_fasta,
                                                 resume,
                                                 translate,
                                                 cpu,
                                                 temp_dir,
                                                 pfam_file)

    else:
        raise EmapperException(f"Unrecognized pfam search mode {pfam_search_mode}.")

    ##
    # 2) Add found pfams to annotations output
    
    if aligned_pfams is not None and all_annotations is not None:
        
        for (hit, annotation), exists in all_annotations:

            # if --resume and annotation exists, skip pfam realignment            
            if exists == False:
                (query_name, best_hit_name, best_hit_evalue, best_hit_score,
                 annotations,
                 (narr_og_name, narr_og_cat, narr_og_desc),
                 (best_og_name, best_og_cat, best_og_desc),
                 match_nog_names,
                 all_orthologies, annot_orthologs) = annotation
            
                if query_name in aligned_pfams:
                    annotations["PFAMs"] = Counter(aligned_pfams[query_name])
                else:
                    annotations["PFAMs"] = None
                
            yield ((hit, annotation), exists)

    return

##
def load_all_annotations(annots_generator):

    all_annotations = []
    queries_pfams = []

    for (hit, annotation), exists in annots_generator:
        all_annotations.append(((hit, annotation), exists))

        # if --resume and annotation exists, skip pfam realignment
        if exists == False:
            (query_name, best_hit_name, best_hit_evalue, best_hit_score,
             annotations,
             (narr_og_name, narr_og_cat, narr_og_desc),
             (best_og_name, best_og_cat, best_og_desc),
             match_nog_names,
             all_orthologies, annot_orthologs) = annotation
            
            if "PFAMs" in annotations:
                queries_pfams.append((query_name, list(annotations["PFAMs"])))
        
    return all_annotations, queries_pfams

## END
