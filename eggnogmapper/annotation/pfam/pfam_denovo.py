##
## CPCantalapiedra 2020

import os

from .pfam_common import filter_fasta_file
from .pfam import PfamAligner, get_pfam_args, parse_pfam_file, parse_hmmsearch_file

##
def pfam_align_denovo(all_annotations, queries_fasta, translate,
                      cpu, num_servers, num_workers, cpus_per_worker, port, end_port,
                      temp_dir, pfam_file):
    aligned_pfams = None
    
    # filter fasta file to have only annotated queries
    queries = [annot_columns[0] for annot_columns in all_annotations]
    fasta_file = filter_fasta_file(queries, queries_fasta, temp_dir)

    # align those queries to whole PFAM to carry out a de novo annotation
    pfam_args, infile = get_pfam_args(cpu, num_servers, num_workers, cpus_per_worker, port, end_port,
                                      fasta_file.name, translate, temp_dir)
    pfam_aligner = PfamAligner(pfam_args)
    pfam_aligner.align_whole_pfam(infile, pfam_file, silent = True)
    aligned_pfams = parse_pfam_file(pfam_file)

    if fasta_file is not None:
        fasta_file.close()
        if os.path.isfile(f"{fasta_file.name}.map"):
            os.remove(f"{fasta_file.name}.map")
        if os.path.isfile(f"{fasta_file.name}.seqdb"):
            os.remove(f"{fasta_file.name}.seqdb")

    return aligned_pfams
                
##
def pfam_align_denovo_search(all_annotations, queries_fasta, translate, cpu, temp_dir, pfam_file):
    aligned_pfams = None
    
    # filter fasta file to have only annotated queries
    queries = [annot_columns[0] for annot_columns in all_annotations]
    fasta_file = filter_fasta_file(queries, queries_fasta, temp_dir)

    # align those queries to whole PFAM to carry out a de novo annotation
    pfam_args, infile = get_pfam_args(cpu, fasta_file.name, translate, temp_dir, force_seqdb = True)
    pfam_aligner = PfamAligner(pfam_args)
    pfam_aligner.align_whole_pfam(infile, pfam_file, silent = True)
    aligned_pfams = parse_hmmsearch_file(pfam_file)

    if fasta_file is not None:
        fasta_file.close()
        if os.path.isfile(f"{fasta_file.name}.map"):
            os.remove(f"{fasta_file.name}.map")
        if os.path.isfile(f"{fasta_file.name}.seqdb"):
            os.remove(f"{fasta_file.name}.seqdb")

    return aligned_pfams

## END
