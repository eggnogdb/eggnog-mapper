##
## CPCantalapiedra 2020

import os
import multiprocessing
import subprocess
from tempfile import NamedTemporaryFile

from ...emapperException import EmapperException
from ...common import get_pfam_db
from ...common import HMMPRESS

from .pfam_common import filter_fasta_hmm_files, wrap_group_queries_pfams
from .pfam import PfamAligner, get_hmmscan_args, parse_hmmscan_file

##
def pfam_align_parallel_scan(annotations, PFAM_COL, queries_fasta, translate, cpu, temp_dir, pfam_file):
    aligned_pfams = None

    pool = multiprocessing.Pool(cpu)

    try:

        for alignments in pool.imap(query_pfam_annotate_scan,
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

                # if "1105367.SAMN02673274.CG50_08330" in alignments:
                #     print(f"annotator.py:pfam_align_parallel_scan: {alignments['1105367.SAMN02673274.CG50_08330']}")

    except EmapperException:
        raise
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise EmapperException(f"Error: annotation went wrong for pfam alignment in parallel. "+str(e))
    finally:
        pool.terminate()

    return aligned_pfams


##
def query_pfam_annotate_scan(arguments):
    queries_pfams_group, queries_fasta, pfam_db, temp_dir, translate, pfam_file, cpu = arguments
    
    aligned_pfams = None

    fasta_file, hmm_file = filter_fasta_hmm_files(queries_pfams_group, queries_fasta, pfam_db, temp_dir)
        
    if fasta_file is None or hmm_file is None:
        pass
    else:

        # create hmmdb
        cmd = f"{HMMPRESS} {hmm_file.name}"
        cp = subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        # output file for this group
        P = NamedTemporaryFile(mode='w', dir=temp_dir)

        # align queries to the new HMM file
        pfam_args, infile = get_hmmscan_args(cpu, fasta_file.name, hmm_file.name, translate, temp_dir)
        pfam_aligner = PfamAligner(pfam_args)
        pfam_aligner.align_whole_pfam(infile, P.name, silent = True)

        aligned_pfams = parse_hmmscan_file(P.name)

        # if "1105367.SAMN02673274.CG50_07170" in aligned_pfams:
        #     print(f"annotator.py:query_pfam_annotate_scan {aligned_pfams}")

        # Append contents of output file for this group into pfam_file,
        # which is the file reporting all the pfam hits together
        with open(pfam_file, 'a') as pfamf:
            pfamf.write(open(P.name, 'r').read())

        P.close()

    if fasta_file is not None:
        fasta_file.close()
    if hmm_file is not None:
        if os.path.isfile(f"{hmm_file.name}.h3m"):
            os.remove(f"{hmm_file.name}.h3m")
        if os.path.isfile(f"{hmm_file.name}.h3p"):
            os.remove(f"{hmm_file.name}.h3p")
        if os.path.isfile(f"{hmm_file.name}.h3f"):
            os.remove(f"{hmm_file.name}.h3f")
        if os.path.isfile(f"{hmm_file.name}.h3i"):
            os.remove(f"{hmm_file.name}.h3i")
            
        hmm_file.close()

    return aligned_pfams

## END
