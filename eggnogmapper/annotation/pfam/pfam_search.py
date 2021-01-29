##
## CPCantalapiedra 2020

import multiprocessing
from tempfile import NamedTemporaryFile

from ...emapperException import EmapperException
from ...common import get_pfam_db

from .pfam import get_hmmsearch_args, PfamAligner, parse_hmmsearch_file
from .pfam_common import filter_fasta_hmm_files, group_queries_pfams, wrap_group_queries_pfams


##
def pfam_align_serial(annotations, PFAM_COL, queries_fasta, translate, cpu, temp_dir, pfam_file):
    aligned_pfams = None

    for queries_pfams_group in group_queries_pfams(annotations, PFAM_COL):
        # fetch pfams to a new HMM file
        # fetch queries to a new Fasta file

        arguments = queries_pfams_group, queries_fasta, get_pfam_db(), temp_dir, translate, pfam_file, cpu
        alignments = query_pfam_annotate(arguments)
        if alignments is not None:
            if aligned_pfams is None:
                aligned_pfams = alignments
            else:
                for k,v in alignments.items():
                    if k in aligned_pfams:
                        aligned_pfams[k] = aligned_pfams[k] | v
                    else:
                        aligned_pfams[k] = v
                # aligned_pfams.update(alignments)

    return aligned_pfams


##
def pfam_align_parallel(annotations, PFAM_COL, queries_fasta, translate, cpu, temp_dir, pfam_file):
    aligned_pfams = None

    pool = multiprocessing.Pool(cpu)

    try:

        for alignments in pool.imap(query_pfam_annotate,
                                    wrap_group_queries_pfams(annotations, PFAM_COL, queries_fasta, get_pfam_db(), translate, temp_dir, pfam_file)):
            if alignments is not None:
                if aligned_pfams is None:
                    aligned_pfams = alignments
                else:
                    for k,v in alignments.items():
                        if k in aligned_pfams:
                            aligned_pfams[k] = aligned_pfams[k] | v
                        else:
                            aligned_pfams[k] = v
                    # aligned_pfams.update(alignments)                

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

    fasta_file, hmm_file = filter_fasta_hmm_files(queries_pfams_group, queries_fasta, pfam_db, temp_dir)

    if fasta_file is None or hmm_file is None:
        pass
    else:
        # output file for this group
        P = NamedTemporaryFile(mode='w', dir=temp_dir)

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

## END
