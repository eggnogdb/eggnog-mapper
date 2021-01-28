
import os
import subprocess

from tempfile import NamedTemporaryFile, mkdtemp
import uuid
from collections import defaultdict

from ...common import *

from .hmmer_seqio import iter_fasta_seqs
from .hmmer_search_hmmpgmd import iter_seq_hits, iter_hmm_hits, QUERY_TYPE_SEQ, QUERY_TYPE_HMM, DB_TYPE_SEQ, DB_TYPE_HMM


SCANTYPE_MEM = "mem"
SCANTYPE_DISK = "disk"


def iter_hits(source, translate, query_type, dbtype, scantype, host, port, servers = None,
              evalue_thr=None, score_thr=None, max_hits=None, return_seq=False,
              skip=None, maxseqlen=None, fixed_Z=None, qcov_thr=None, cut_ga=False, cpus=1,
              base_tempdir=None, silent=False):

    try:
        max_hits = int(max_hits)
        if max_hits == 0: # unlimited hits
            max_hits = None
    except Exception:
        max_hits = None

    ## On mem searches
    # "hmmscan- and phmmer-like" modes
    if scantype == SCANTYPE_MEM and query_type == QUERY_TYPE_SEQ:
        return iter_seq_hits(source, translate, cpus, servers, dbtype=dbtype, evalue_thr=evalue_thr, score_thr=score_thr,
                             max_hits=max_hits, skip=skip, maxseqlen=maxseqlen, cut_ga=cut_ga, silent=silent)

    # "hmmsearch"-like mode
    elif scantype == SCANTYPE_MEM and query_type == QUERY_TYPE_HMM and dbtype == DB_TYPE_SEQ:
        return iter_hmm_hits(source, cpus, servers, dbtype=dbtype, evalue_thr=evalue_thr, score_thr=score_thr,
                             max_hits=max_hits, skip=skip, maxseqlen=maxseqlen, fixed_Z=fixed_Z, cut_ga=cut_ga, silent=silent)

    ## On disk searches
    # hmmscan mode
    elif scantype == SCANTYPE_DISK and query_type == QUERY_TYPE_SEQ and dbtype == DB_TYPE_HMM:
        return hmmscan(source, translate, host, evalue_thr=evalue_thr, score_thr=score_thr, max_hits=max_hits, 
                       cpus=cpus, maxseqlen=maxseqlen, base_tempdir=base_tempdir, cut_ga=cut_ga, silent=silent)

    # hmmsearch mode
    elif scantype == SCANTYPE_DISK and query_type == QUERY_TYPE_HMM and dbtype == DB_TYPE_SEQ:
        # host is the fasta_file in this case
        # source is the hmm_file in this case
        return hmmsearch(host, translate, source, evalue_thr=evalue_thr, score_thr=score_thr, max_hits=max_hits, 
                         cpus=cpus, maxseqlen=maxseqlen, base_tempdir=base_tempdir, cut_ga=cut_ga, silent=silent)

    # phmmer mode
    elif scantype == SCANTYPE_DISK and query_type == QUERY_TYPE_SEQ and dbtype == DB_TYPE_SEQ:
        return phmmer(source, translate, host, evalue_thr=evalue_thr, score_thr=score_thr, max_hits=max_hits, 
                      cpus=cpus, maxseqlen=maxseqlen, base_tempdir=base_tempdir, cut_ga=cut_ga, silent=silent)
    
    elif query_type == QUERY_TYPE_HMM and dbtype == DB_TYPE_HMM:
        raise Exception("HMM to HMM search is not supported.")        
    
    else:
        raise ValueError('not supported')

    
def safe_cast(v):
    try:
        return float(v)
    except ValueError:
        return v.strip()


##
def hmmscan(fasta_file, translate, hmm_file, cpus=1, evalue_thr=None,
            score_thr=None, max_hits=None, fixed_Z=None, maxseqlen=None, cut_ga=False,
            base_tempdir=None, silent=False):
    
    cmd = HMMSCAN
    return hmmcommand(cmd, fasta_file, translate, hmm_file, cpus, evalue_thr,
                      score_thr, max_hits, fixed_Z, maxseqlen, cut_ga,
                      base_tempdir, silent)


##
def hmmsearch(fasta_file, translate, hmm_file, cpus=1, evalue_thr=None,
              score_thr=None, max_hits=None, fixed_Z=None, maxseqlen=None, cut_ga=False,
              base_tempdir=None, silent=False):
    
    cmd = HMMSEARCH
    return hmmcommand(cmd, fasta_file, translate, hmm_file, cpus, evalue_thr,
                      score_thr, max_hits, fixed_Z, maxseqlen, cut_ga,
                      base_tempdir, silent)


##
def phmmer(fasta_file, translate, fasta_target_file, cpus=1, evalue_thr=None,
           score_thr=None, max_hits=None, fixed_Z=None, maxseqlen=None, cut_ga=False,
           base_tempdir=None, silent=False):
    
    cmd = PHMMER
    return hmmcommand(cmd, fasta_file, translate, fasta_target_file, cpus, evalue_thr,
                      score_thr, max_hits, fixed_Z, maxseqlen, cut_ga,
                      base_tempdir, silent)


##
def hmmcommand(hmmer_cmd, fasta_file, translate, hmm_file, cpus=1, evalue_thr=None,
               score_thr=None, max_hits=None, fixed_Z=None, maxseqlen=None, cut_ga=False,
               base_tempdir=None, silent=False):
    
    if not hmmer_cmd:
        raise ValueError(f'{hmmer_cmd} not found in path')
    
    tempdir = mkdtemp(prefix='emappertmp_hmmcmd_', dir=base_tempdir)

    Q = None
    R = None
    
    ##
    # Translate FASTA nts to aas if needed
    
    if translate or maxseqlen is not None:
        if translate:
            if silent == False:
                print('translating query input file')
        Q = NamedTemporaryFile(mode='w', dir=tempdir)
        has_records = False
        for name, seq in iter_fasta_seqs(fasta_file, translate=translate, silent=silent):
            if maxseqlen is None or len(seq) <= maxseqlen:
                has_records = True
                print(f">{name}\n{seq}", file=Q)
                # Q.write(f">{name}\n{seq}".encode())
        if has_records == False:
            sys.stderr.write(f"No records after maxseqlen filtering for file {fasta_file}.\n")
            shutil.rmtree(tempdir)
            if Q is not None:
                Q.close()
            return
        Q.flush()
        fasta_file = Q.name

    if hmmer_cmd == PHMMER:
        # if cmd in phmmer, hmm_file is actually a fasta file of sequences
        if translate or maxseqlen:
            if translate:
                if silent == False:
                    print('translating target fasta file')
            R = NamedTemporaryFile(mode='w', dir=tempdir)
            for name, seq in iter_fasta_seqs(hmm_file, translate=translate, silent=silent):
                if maxseqlen is None or len(seq) <= maxseqlen:
                    print(f">{name}\n{seq}", file=R)
                    # Q.write(f">{name}\n{seq}".encode())
            R.flush()
            hmm_file = R.name

    # if report_no_hits == True:
    # we need the list of queries to be sure that all queries without hits
    # are reported, which would be consistent with the results from
    # hmmpgmd
    queries_dict = {}
    if hmmer_cmd == HMMSCAN or hmmer_cmd == PHMMER:
        for name, seq in iter_fasta_seqs(fasta_file, translate=False, silent=silent): # seqs were already translated if needed
            queries_dict[name] = len(seq)
    elif hmmer_cmd == HMMSEARCH:
        with open(hmm_file, 'r') as hmm_data:
            for line in hmm_data:
                if line.startswith("NAME"):
                    name = line.split()[-1]

                if line.startswith("LENG"):
                    queries_dict[name] = int(line.split()[-1]) # query length

    else:
        raise Exception(f"cmd {hmmer_cmd} is not supported")

    ##
    # Run command
    
    OUT = NamedTemporaryFile(mode='w+', dir=tempdir)

    if cut_ga:
        cut_ga_p = " --cut_ga"
    else:
        cut_ga_p = ""
        
    cmd = '%s %s --cpu %s -o /dev/null --domtblout %s %s %s' % (
        hmmer_cmd, cut_ga_p, cpus, OUT.name, hmm_file, fasta_file)
    if silent == False:
        print(f'# {cmd}')
    sts = subprocess.call(cmd, shell=True)
    # if sts != 0:
    #     print(f"X {cmd}")
    #     print(f"FASTA: {fasta_file}")
    #     print(os.stat(fasta_file).st_size)
    #     print("Printing fasta file")
    #     with open(fasta_file, 'r') as f:
    #         for line in f:
    #             print(f)
        
    #     print(f"HMMFILE: {hmm_file}")
    # else:
    #     print(f"# {cmd}")        
        

    ##
    # Process output
    
    byquery = defaultdict(list)

    last_query = None
    last_hitname = None
    hit_list = []
    hit_ids = set()
    last_query_len = None
    queries_with_hits = set()
    if sts == 0:
        for line in OUT:
            
            # domain
            if line.startswith('#'):
                continue
            fields = line.split()

            # This is the list of fields from 0 to 21
            # The code was problematic with string values similar to int or floats
            # so I will be not used anymore, but the list is useful to keep track
            # of the available fields
            #
            # (hitname, hacc, tlen, qname, qacc, qlen, evalue, score, bias, didx,
            #  dnum, c_evalue, i_evalue, d_score, d_bias, hmmfrom, hmmto, seqfrom,
            #  seqto, env_from, env_to, acc) = list(map(safe_cast, fields[:22]))
            
            hitname = str(fields[0]).strip()
            qname = str(fields[3]).strip()
            qlen = int(fields[5])
            evalue = float(fields[6])
            score = float(fields[7])
            d_score = float(fields[13])
            hmmfrom = int(fields[15])
            hmmto = int(fields[16])
            seqfrom = int(fields[17])
            seqto = int(fields[18])
            
            # If a new query is being processed,
            # report the results of the previous one
            if last_query and qname != last_query:
                if len(hit_list) > 0:
                    yield last_query, 0, hit_list, last_query_len, None
                    queries_with_hits.add(last_query)

                last_query = qname
                last_hitname = None
                hit_list = []
                hit_ids = set()
                last_query_len = None

            if last_query_len and last_query_len != qlen:
                raise ValueError("Inconsistent qlen when parsing hmmscan output")
            last_query_len = qlen

            # Filter thresholds
            if (cut_ga == True or \
               ((evalue_thr is None or evalue <= evalue_thr) and \
               (score_thr is None or score >= score_thr))) and \
               (max_hits is None or last_hitname == hitname or len(hit_ids) < max_hits):

                # if "CG50_09910" in qname:
                #     print(f"hmmer_search.py: PASS FILTER {hitname}")
                    
                hit_list.append([hitname, evalue, score, hmmfrom,
                                 hmmto, seqfrom, seqto, d_score])
                hit_ids.add(hitname)
                last_hitname = hitname
            
            # if (evalue_thr is None or evalue <= evalue_thr) and \
            #    (score_thr is None or score >= score_thr) and \
            #    (max_hits is None or last_hitname == hitname or len(hit_ids) < max_hits):

            #     hit_list.append([hitname, evalue, score, hmmfrom,
            #                      hmmto, seqfrom, seqto, d_score])
            #     hit_ids.add(hitname)
            #     last_hitname = hitname
                
            last_query = qname

        # Finally, report results of the last processed query
        if last_query and len(hit_list) > 0:
            yield last_query, 0, hit_list, last_query_len, None
            # if report_no_hits == True:
            queries_with_hits.add(last_query)

    OUT.close()
    # if translate:
    if Q is not None:
        Q.close()
    if R is not None:
        R.close()
    shutil.rmtree(tempdir)

    # if report_no_hits == True:
    # report queries without hits
    queries_without_hits = set(queries_dict.keys()) ^ queries_with_hits
    for query in queries_without_hits:
        qlen = queries_dict[query]
        yield query, 0, [], qlen, None

    return


def refine_hit(args):
    seqname, seq, group_fasta, excluded_taxa, tempdir = args
    F = NamedTemporaryFile(delete=True, dir=tempdir, mode='w+')
    F.write(f'>{seqname}\n{seq}')
    F.flush()

    best_hit = get_best_hit(F.name, group_fasta, excluded_taxa, tempdir)
    F.close()

    return [seqname] + best_hit


def get_best_hit(target_seq, target_og, excluded_taxa, tempdir):
    if not PHMMER:
        raise ValueError('phmmer not found in path')

    tempout = pjoin(tempdir, uuid.uuid4().hex)
    cmd = "%s --cpu 1 --incE 0.001 -E 0.001 -o /dev/null --noali --tblout %s %s %s" %\
          (PHMMER, tempout, target_seq, target_og)

    # print cmd
    status = os.system(cmd)
    best_hit_found = False
    if status == 0:
        # take the best hit
        for line in open(tempout):
            if line.startswith('#'):
                continue
            else:
                fields = line.split()
                best_hit_name = fields[0]
                best_hit_evalue = float(fields[4])
                best_hit_score = float(fields[5])
                if not excluded_taxa or not best_hit_name.startswith("%s." % (excluded_taxa)):
                    best_hit_found = True
                    break
        os.remove(tempout)
    else:
        raise ValueError('Error running PHMMER')

    if not best_hit_found:
        best_hit_evalue = '-'
        best_hit_score = '-'
        best_hit_name = '-'

    return [best_hit_name, best_hit_evalue, best_hit_score]

## END
