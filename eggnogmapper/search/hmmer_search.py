#!/usr/bin/env python2

import sys
import os
import socket
import struct
import math
import re
import time
import subprocess
import pickle
import multiprocessing

from tempfile import NamedTemporaryFile, mkdtemp
import uuid
from collections import defaultdict, Counter

from .hmmer_seqio import iter_fasta_seqs
from ..annotation import annota
from ..common import *


SCANTYPE_MEM = "mem"
SCANTYPE_DISK = "disk"

QUERY_TYPE_SEQ = "seq"
QUERY_TYPE_HMM = "hmm"

DB_TYPE_SEQ = "seqdb"
DB_TYPE_HMM = "hmmdb"

B62_IDENTITIES = {'A': 4, 'B': 4, 'C': 9, 'D': 6, 'E': 5, 'F': 6, 'G': 6, 'H': 8,
                  'I': 4, 'K': 5, 'L': 4, 'M': 5, 'N': 6, 'P': 7, 'Q': 5, 'R': 5,
                  'S': 4, 'T': 5, 'V': 4, 'W': 11, 'X': -1, 'Y': 7, 'Z': 4}


def iter_hits(source, translate, query_type, dbtype, scantype, host, port, servers = None,
              evalue_thr=None, score_thr=None, max_hits=None, return_seq=False,
              skip=None, maxseqlen=None, fixed_Z=None, qcov_thr=None, cut_ga=False, cpus=1,
              base_tempdir=None):

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
                             max_hits=max_hits, skip=skip, maxseqlen=maxseqlen, cut_ga=cut_ga)

    # "hmmsearch"-like mode
    elif scantype == SCANTYPE_MEM and query_type == QUERY_TYPE_HMM and dbtype == DB_TYPE_SEQ:
        return iter_hmm_hits(source, cpus, servers, dbtype=dbtype, evalue_thr=evalue_thr, score_thr=score_thr,
                             max_hits=max_hits, skip=skip, fixed_Z=fixed_Z, cut_ga=cut_ga)

    ## On disk searches
    # hmmscan mode
    elif scantype == SCANTYPE_DISK and query_type == QUERY_TYPE_SEQ and dbtype == DB_TYPE_HMM:
        return hmmscan(source, translate, host, evalue_thr=evalue_thr, score_thr=score_thr, max_hits=max_hits, cpus=cpus, maxseqlen=maxseqlen, base_tempdir=base_tempdir)

    # hmmsearch mode
    elif scantype == SCANTYPE_DISK and query_type == QUERY_TYPE_HMM and dbtype == DB_TYPE_SEQ:
        # host is the fasta_file in this case
        # source is the hmm_file in this case
        return hmmsearch(host, translate, source, evalue_thr=evalue_thr, score_thr=score_thr, max_hits=max_hits, cpus=cpus, maxseqlen=maxseqlen, base_tempdir=base_tempdir)

    # phmmer mode
    elif scantype == SCANTYPE_DISK and query_type == QUERY_TYPE_SEQ and dbtype == DB_TYPE_SEQ:
        return phmmer(source, translate, host, evalue_thr=evalue_thr, score_thr=score_thr, max_hits=max_hits, cpus=cpus, maxseqlen=maxseqlen, base_tempdir=base_tempdir)
    
    elif query_type == QUERY_TYPE_HMM and dbtype == DB_TYPE_HMM:
        raise Exception("HMM to HMM search is not supported.")        
    
    else:
        raise ValueError('not supported')

    
def safe_cast(v):
    try:
        return float(v)
    except ValueError:
        return v.strip()


def unpack_hit(bindata, z):
    (name, acc, desc, window_length, sort_key, score, pre_score, sum_score,
     pvalue, pre_pvalue, sum_pvalue, nexpected, nregions, nclustered,
     noverlaps, nenvelopes, ndom, flags, nreported, nincluded, best_domain,
     seqidx, subseq_start, dcl, offset) = struct.unpack("3Q I 4x d 3f 4x 3d f 9I 4Q", bindata)

    evalue = math.exp(pvalue) * z
    return name, evalue, sum_score, ndom


def unpack_stats(bindata):
    (elapsed, user, sys, Z, domZ, Z_setby, domZ_setby, nmodels, nseqs,
     n_past_msv, n_past_bias, n_past_vit, n_past_fwd, nhits, nreported,
     nincluded) = struct.unpack("5d 2I 9q", bindata)

    return elapsed, nhits, Z, domZ


def scan_hits(data, address="127.0.0.1", port=51371, evalue_thr=None,
              score_thr=None, max_hits=None, fixed_Z=None):
    
    # print(data)
    
    s = socket.socket()
    try:
        s.connect((address, port))
    except Exception as e:
        print(address, port, e)
        raise
    s.sendall(data.encode())

    status = s.recv(16)
    st, msg_len = struct.unpack("I 4x Q", status)
    elapsed, nreported = 0, 0
    hits = defaultdict(list)
    if st == 0:
        binresult = b''
        while len(binresult) < msg_len:
            binresult += s.recv(4096)

        elapsed, nreported, Z, domZ = unpack_stats(binresult[0:120])
        if fixed_Z:
            Z = fixed_Z

        hitdata = defaultdict(dict)

        hits_start = 120 # First 120 bits are the stats
        hits_end = hits_start

        for hitid in range(nreported):
            hits_end = hits_start + 152
            name, evalue, score, ndom = unpack_hit(binresult[hits_start:hits_end], Z)
            hitdata[hitid] = {"name": name, "evalue": evalue, "score":score, "ndom": ndom, "doms": []}
            hits_start += 152

        next_start = hits_end
        reported_hits = []
        for hitid in range(nreported):
            hit = hitdata[hitid]
            
            for domid in range(hit["ndom"]):
                dombit = binresult[next_start:next_start + 72]
                dom = struct.unpack("4i 5f 4x d 2i Q 8x", dombit)
                is_reported = dom[10]
                is_included = dom[11]
                hit["doms"].append(dom)
                next_start += 72

            lasthitname = None
            for domid in range(hit["ndom"]):
                alibit = struct.unpack("7Q I 4x 3Q 3I 4x 6Q I 4x Q", binresult[next_start:next_start+168])

                (rfline, mmline, csline, model, mline, aseq, ppline, N,
                hmmname, hmmacc, hmmdesc, hmmfrom, hmmto, M, sqname, sqacc,
                sqdesc,sqfrom, sqto, L, memsize, mem) = alibit
                next_start += 168
                # Skip alignment
                # ....
                next_start += (memsize)
                d = hit["doms"][domid]
                bitscore = d[8]
                ievalue = math.exp(d[9] * Z)
                cevalue = math.exp(d[9] * domZ)
                hitname = hit["name"]
                evalue = hit["evalue"]
                score = hit["score"]
                
                if (evalue_thr is None or evalue <= evalue_thr) and \
                    (score_thr is None or score >= score_thr) and \
                    (max_hits is None or len(reported_hits)+1 <= max_hits or lasthitname == hitname):
                    
                    reported_hits.append((hitname, evalue, score, hmmfrom,
                                          hmmto, sqfrom, sqto, bitscore))

                    lasthitname = hitname

                # IMPORTANT: WE MUST NOT BREAK BECAUSE OF THE BYTES LOGIC
                # READING HMMPGMD RESPONSE
                # elif (max_hits is not None and len(reported_hits)+1 > max_hits and lasthitname != hitname):
                #     break
                    

    else:
        ret = s.recv(4096).decode().strip()
        s.close()        
        raise ValueError('hmmpgmd error: %s' % ret)

    s.close()
    
    return elapsed, reported_hits


def iter_hmm_file(hmmfile, skip):
    hmmer_version = None
    model = ''
    name = 'Unknown'
    leng = None
    with open(hmmfile) as HMMFILE:
        for line in HMMFILE:
            if hmmer_version is None:
                hmmer_version = line
                
            if line.startswith("NAME"):
                name = line.split()[-1]
                model = ''
                leng = None
            if line.startswith("LENG"):
                leng = int(line.split()[-1])
                
            model += line
            if line.strip() == '//':
                if skip and name in skip:
                    continue
                yield hmmer_version, name, leng, model

    return
    
def iter_hmm(hmm):
    hmm_num, hmmer_version, name, leng, model, servers, dbtype, evalue_thr, score_thr, max_hits, fixed_Z, cut_ga = hmm

    num_servers = len(servers)
    num_server = hmm_num % num_servers
    host, port = servers[num_server]
    
    if cut_ga == True:
        cut_ga = " --cut_ga"
    else:
        cut_ga = ""
        
    data = f'@--{dbtype} 1 {cut_ga}\n{hmmer_version}\n{model}'    
    etime, hits = scan_hits(data, host, port,
                            evalue_thr=evalue_thr, score_thr=score_thr,
                            max_hits=max_hits, fixed_Z=fixed_Z)

    return name, etime, hits, leng, None
    
def iter_hmm_hits(hmmfile, cpus, servers, dbtype=DB_TYPE_HMM,
                  evalue_thr=None, score_thr=None,
                  max_hits=None, skip=None, fixed_Z=None, cut_ga=False):
    
    pool = multiprocessing.Pool(cpus)
    # hmms = [[hmmnum, hmmer_version, name, leng, model, servers, dbtype, evalue_thr, score_thr, max_hits, fixed_Z, skip, cut_ga]
    #         for hmmnum, (hmmer_version, name, leng, model) in
    #         enumerate(iter_hmm_file(hmmfile))]
    
    # for r in pool.imap(iter_hmm, hmms):
    for r in pool.imap(iter_hmm, ([hmmnum, hmmer_version, name, leng, model, servers, dbtype, evalue_thr, score_thr, max_hits, fixed_Z, cut_ga]
                                  for hmmnum, (hmmer_version, name, leng, model) in
                                  enumerate(iter_hmm_file(hmmfile, skip)))):
        yield r
    pool.terminate()
    return
    # with open(hmmfile) as HMMFILE:
    #     for line in HMMFILE:

    #         if hmmer_version is None:
    #             hmmer_version = line
                
    #         if line.startswith("NAME"):
    #             name = line.split()[-1]
    #             model = ''
    #             leng = None
    #         if line.startswith("LENG"):
    #             leng = int(line.split()[-1])
                
    #         model += line
    #         if line.strip() == '//':
    #             if skip and name in skip:
    #                 continue
    #             else:                    
    #                 data = f'@--{dbtype} 1\n{hmmer_version}\n{model}'

    #                 etime, hits = scan_hits(data, host, port,
    #                                         evalue_thr=evalue_thr, score_thr=score_thr,
    #                                         max_hits=max_hits, fixed_Z=fixed_Z)

    #                 yield name, etime, hits, leng, None

def iter_seq(seq):
    seqnum, name, seq, servers, dbtype, evalue_thr, score_thr, max_hits, maxseqlen, fixed_Z, skip, cut_ga = seq

    num_servers = len(servers)
    num_server = seqnum % num_servers
    host, port = servers[num_server]
    
    if skip and name in skip:
        return

    if maxseqlen and len(seq) > maxseqlen:
        return name, -1, [], len(seq), None

    if not seq:
        return

    if cut_ga == True:
        cut_ga = " --cut_ga"
    else:
        cut_ga = ""

    seq = re.sub("-.", "", seq)
    data = '@--%s 1%s\n>%s\n%s\n//' % (dbtype, cut_ga, name, seq)
    etime, hits = scan_hits(data, host, port, evalue_thr=evalue_thr,
                            score_thr=score_thr, max_hits=max_hits,
                            fixed_Z=fixed_Z)

    return name, etime, hits, len(seq), None


def iter_seq_hits(src, translate, cpus, servers, dbtype, evalue_thr=None,
                  score_thr=None, max_hits=None, maxseqlen=None, fixed_Z=None,
                  skip=None, cut_ga=False):

    pool = multiprocessing.Pool(cpus)
    # seqs = [[seqnum, name, seq, servers, dbtype, evalue_thr, score_thr, max_hits, maxseqlen, fixed_Z, skip, cut_ga]
    #         for seqnum, (name, seq) in
    #         enumerate(iter_fasta_seqs(src, translate=translate))]
    
    # for r in pool.imap(iter_seq, seqs):
    for r in pool.imap(iter_seq, ([seqnum, name, seq, servers, dbtype, evalue_thr, score_thr, max_hits, maxseqlen, fixed_Z, skip, cut_ga]
                                  for seqnum, (name, seq) in
                                  enumerate(iter_fasta_seqs(src, translate=translate)))):
        yield r
    pool.terminate()
    return


def get_hits(name, record, address="127.0.0.1", port=51371, dbtype=DB_TYPE_HMM, qtype=QUERY_TYPE_SEQ,
             evalue_thr=None, score_thr = None, max_hits=None):

    if qtype == QUERY_TYPE_SEQ:
        seq = re.sub("-.", "", record)
        data = f'@--{dbtype} 1\n>{name}\n{seq}\n//'
    elif qtype == QUERY_TYPE_HMM:
        data = f'@--{dbtype} 1\n{record}\n//'        
    else:
        raise Exception(f"Unrecognized query type {qtype}.")
    
    etime, hits = scan_hits(data, address=address, port=port,
                            evalue_thr=evalue_thr, score_thr=score_thr, max_hits=max_hits)
    
    return name, etime, hits


##
def hmmscan(fasta_file, translate, hmm_file, cpus=1, evalue_thr=None,
            score_thr=None, max_hits=None, fixed_Z=None, maxseqlen=None,
            base_tempdir=None):
    
    cmd = HMMSCAN
    return hmmcommand(cmd, fasta_file, translate, hmm_file, cpus, evalue_thr,
            score_thr, max_hits, fixed_Z, maxseqlen,
            base_tempdir)


##
def hmmsearch(fasta_file, translate, hmm_file, cpus=1, evalue_thr=None,
            score_thr=None, max_hits=None, fixed_Z=None, maxseqlen=None,
            base_tempdir=None):
    
    cmd = HMMSEARCH
    return hmmcommand(cmd, fasta_file, translate, hmm_file, cpus, evalue_thr,
            score_thr, max_hits, fixed_Z, maxseqlen,
            base_tempdir)


##
def phmmer(fasta_file, translate, fasta_target_file, cpus=1, evalue_thr=None,
            score_thr=None, max_hits=None, fixed_Z=None, maxseqlen=None,
            base_tempdir=None):
    
    cmd = PHMMER
    return hmmcommand(cmd, fasta_file, translate, fasta_target_file, cpus, evalue_thr,
            score_thr, max_hits, fixed_Z, maxseqlen,
            base_tempdir)


##
def hmmcommand(hmmer_cmd, fasta_file, translate, hmm_file, cpus=1, evalue_thr=None,
            score_thr=None, max_hits=None, fixed_Z=None, maxseqlen=None,
            base_tempdir=None):
    
    if not hmmer_cmd:
        raise ValueError(f'{hmmer_cmd} not found in path')
    
    tempdir = mkdtemp(prefix='emappertmp_hmmcmd_', dir=base_tempdir)

    ##
    # Translate FASTA nts to aas if needed
    
    if translate or maxseqlen:
        if translate:
            print('translating query input file')
        Q = NamedTemporaryFile(mode='w')
        for name, seq in iter_fasta_seqs(fasta_file, translate=translate):
            if maxseqlen is None or len(seq) <= maxseqlen:
                print(f">{name}\n{seq}", file=Q)
                # Q.write(f">{name}\n{seq}".encode())
        Q.flush()
        fasta_file = Q.name

    if hmmer_cmd == PHMMER:
        # if cmd in phmmer, hmm_file is actually a fasta file of sequences
        if translate or maxseqlen:
            if translate:
                print('translating target fasta file')
            R = NamedTemporaryFile(mode='w')
            for name, seq in iter_fasta_seqs(hmm_file, translate=translate):
                if maxseqlen is None or len(seq) <= maxseqlen:
                    print(f">{name}\n{seq}", file=R)
                    # Q.write(f">{name}\n{seq}".encode())
            R.flush()
            hmm_file = R.name

    # we need the list of queries to be sure that all queries without hits
    # are reported, which would be consistent with the results from
    # hmmpgmd
    queries_dict = {}
    if hmmer_cmd == HMMSCAN or hmmer_cmd == PHMMER:
        for name, seq in iter_fasta_seqs(fasta_file, translate=False): # seqs were already translated if needed
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
    
    OUT = NamedTemporaryFile(dir=tempdir, mode='w+')
    
    cmd = '%s --cpu %s -o /dev/null --domtblout %s %s %s' % (
        hmmer_cmd, cpus, OUT.name, hmm_file, fasta_file)
    # print '#', cmd
    sts = subprocess.call(cmd, shell=True)

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

            (hitname, hacc, tlen, qname, qacc, qlen, evalue, score, bias, didx,
             dnum, c_evalue, i_evalue, d_score, d_bias, hmmfrom, hmmto, seqfrom,
             seqto, env_from, env_to, acc) = list(map(safe_cast, fields[:22]))

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
                raise ValueError(
                    "Inconsistent qlen when parsing hmmscan output")
            last_query_len = qlen

            # Filter thresholds
            if (evalue_thr is None or evalue <= evalue_thr) and \
               (score_thr is None or score >= score_thr) and \
               (max_hits is None or last_hitname == hitname or len(hit_ids) < max_hits):

                hit_list.append([hitname, evalue, score, hmmfrom,
                                 hmmto, seqfrom, seqto, d_score])
                hit_ids.add(hitname)
                last_hitname = hitname
                
            last_query = qname

        # Finally, report results of the last processed query
        if last_query and len(hit_list) > 0:
            yield last_query, 0, hit_list, last_query_len, None
            queries_with_hits.add(last_query)

    OUT.close()
    if translate:
        Q.close()
    shutil.rmtree(tempdir)

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
    cmd = "%s --incE 0.001 -E 0.001 -o /dev/null --noali --tblout %s %s %s" %\
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
