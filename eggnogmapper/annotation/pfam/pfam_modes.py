##
## CPCantalapiedra 2020

from ...emapperException import EmapperException
from ...utils import colorify

from .pfam_denovo import pfam_align_denovo, pfam_align_denovo_search
from .pfam_search import pfam_align_serial, pfam_align_parallel
from .pfam_scan import pfam_align_parallel_scan

PFAM_COL = -1 # position of PFAMs annotations in list of annotations

PFAM_TRANSFER_BEST_OG = 'best_og'
PFAM_TRANSFER_NARROWEST_OG = 'narrowest_og'
PFAM_TRANSFER_SEED_ORTHOLOG = 'seed_ortholog'

PFAM_REALIGN_NONE = 'none'
PFAM_REALIGN_REALIGN = 'realign'
PFAM_REALIGN_DENOVO = 'denovo'

def run_pfam_mode(pfam_search_mode, all_annotations, queries_fasta, translate,
                  cpu, num_servers, num_workers, cpus_per_worker, port, end_port,
                  temp_dir, pfam_file):
    aligned_pfams = None
    if pfam_search_mode == PFAM_REALIGN_DENOVO:
        print(colorify("de novo scan of PFAM domains", 'green'))

        aligned_pfams = pfam_align_denovo(all_annotations, queries_fasta, translate,
                                          cpu, num_servers, num_workers, cpus_per_worker, port, end_port,
                                          temp_dir, pfam_file)

    # elif pfam_search_mode == 'denovo_search':

    #     print(colorify("de novo search of PFAM domains", 'green'))

    #     aligned_pfams = pfam_align_denovo_search(all_annotations, queries_fasta, translate, cpu, temp_dir, pfam_file)

    # elif pfam_search_mode == 'align':
    #     print(colorify("re-aligning PFAM domains from orthologs to queries", 'green'))

    #     aligned_pfams = pfam_align_serial(all_annotations, PFAM_COL, queries_fasta, translate, cpu, temp_dir, pfam_file)

    # elif pfam_search_mode == 'alignp':
    #     print(colorify("re-aligning PFAM domains from orthologs to queries, in parallel mode", 'green'))

    #     aligned_pfams = pfam_align_parallel(all_annotations, PFAM_COL, queries_fasta, translate, cpu, temp_dir, pfam_file)

    elif pfam_search_mode == PFAM_REALIGN_REALIGN:
        print(colorify("re-aligning PFAM domains from orthologs to queries, in parallel mode, using hmmscan", 'green'))

        aligned_pfams = pfam_align_parallel_scan(all_annotations, PFAM_COL, queries_fasta, translate, cpu, temp_dir, pfam_file)

    else:
        raise EmapperException(f"Unrecognized pfam search mode {pfam_search_mode}.")

    # Add found pfams to annotations output
    if aligned_pfams is not None:
        for annot_columns in all_annotations:
            query_name = annot_columns[0]
            if query_name in aligned_pfams:

                # if "CG50_07170" in query_name:
                #     print(f"annotator.py:annotate {aligned_pfams[query_name]}")

                annot_columns[PFAM_COL] = ",".join(sorted(aligned_pfams[query_name]))
            else:
                annot_columns[PFAM_COL] = "-"

    return all_annotations

## END
