

import sys
import socket
import struct
import math
import re
import multiprocessing

from collections import defaultdict

from ...common import *

from .hmmer_seqio import iter_fasta_seqs

QUERY_TYPE_SEQ = "seq"
QUERY_TYPE_HMM = "hmm"

DB_TYPE_SEQ = "seqdb"
DB_TYPE_HMM = "hmmdb"


def iter_seq_hits(src, translate, cpus, servers, dbtype, evalue_thr=None,
                  score_thr=None, max_hits=None, maxseqlen=None, fixed_Z=None,
                  skip=None, cut_ga=False, silent=False):
        
    pool = multiprocessing.Pool(cpus)
    sys.stderr.write(f"Searching queries with a pool of {cpus} CPUs\n")
    for r in pool.imap(iter_seq, ([seqnum, name, seq, servers, dbtype, evalue_thr, score_thr, max_hits, maxseqlen, fixed_Z, skip, cut_ga]
                                  for seqnum, (name, seq) in
                                  enumerate(iter_fasta_seqs(src, translate=translate, silent=silent)))):
        yield r
    pool.terminate()
    return


def iter_seq(seq):
    seqnum, name, seq, servers, dbtype, evalue_thr, score_thr, max_hits, maxseqlen, fixed_Z, skip, cut_ga = seq

    num_servers = len(servers)
    num_server = seqnum % num_servers
    host, port = servers[num_server]
    # sys.stderr.write(f"Search seq {seqnum} on server {num_server} at {host}:{port}\n")
    
    if skip and name in skip:
        return name, -1, ["SKIPPED"], len(seq), None

    if maxseqlen and len(seq) > maxseqlen:
        return name, -1, ["SEQ_TOO_LARGE "+str(len(seq))], len(seq), None

    if not seq:
        return name, -1, ["NO_SEQ_FOUND"], len(seq), None

    if cut_ga == True:
        cut_ga_p = " --cut_ga"
    else:
        cut_ga_p = ""
        
    seq = re.sub("-.", "", seq)
    data = '@--%s 1%s\n>%s\n%s\n//' % (dbtype, cut_ga_p, name, seq)    
    etime, hits = scan_hits(data, host, port, cut_ga=cut_ga, evalue_thr=evalue_thr,
                            score_thr=score_thr, max_hits=max_hits,
                            fixed_Z=fixed_Z)

    return name, etime, hits, len(seq), None


def iter_hmm_hits(hmmfile, cpus, servers, dbtype=DB_TYPE_HMM,
                  evalue_thr=None, score_thr=None,
                  max_hits=None, skip=None, maxseqlen=None,
                  fixed_Z=None, cut_ga=False, silent=False):
        
    pool = multiprocessing.Pool(cpus)
    for r in pool.imap(iter_hmm, ([hmmnum, hmmer_version, name, leng, model, servers, dbtype, evalue_thr, score_thr,
                                   max_hits, maxseqlen, fixed_Z, skip, cut_ga]
                                  for hmmnum, (hmmer_version, name, leng, model) in
                                  enumerate(iter_hmm_file(hmmfile, skip, silent=silent)))):
        yield r
    pool.terminate()
    return


def iter_hmm_file(hmmfile, skip, silent=False):
    if silent == False:
        sys.stderr.write(f"Parsing hmm file {hmmfile}...\n")
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

    if silent == False:
        sys.stderr.write(f"hmm file {hmmfile} parsing complete.\n")
    return


def iter_hmm(hmm):
    hmm_num, hmmer_version, name, leng, model, servers, dbtype, evalue_thr, score_thr, max_hits, maxseqlen, fixed_Z, skip, cut_ga = hmm

    if skip and name in skip:
        return name, -1, ["SKIPPED"], leng, None

    if maxseqlen and leng > maxseqlen:
        return name, -1, ["SEQ_TOO_LARGE "+str(leng)], leng, None
    
    num_servers = len(servers)
    num_server = hmm_num % num_servers
    host, port = servers[num_server]

    if cut_ga == True:
        cut_ga_p = " --cut_ga"
    else:
        cut_ga_p = ""
        
    data = f'@--{dbtype} 1 {cut_ga_p}\n{hmmer_version}\n{model}'    
    etime, hits = scan_hits(data, host, port,
                            cut_ga=cut_ga, evalue_thr=evalue_thr, score_thr=score_thr,
                            max_hits=max_hits, fixed_Z=fixed_Z)

    return name, etime, hits, leng, None


def get_hits(name, record, address="127.0.0.1", port=51371, dbtype=DB_TYPE_HMM, qtype=QUERY_TYPE_SEQ,
             evalue_thr=None, score_thr = None, max_hits=None):

    if qtype == QUERY_TYPE_SEQ:
        seq = re.sub("-.", "", record)
        data = f'@--{dbtype} 1\n>{name}\n{seq}\n//'
    elif qtype == QUERY_TYPE_HMM:
        data = f'@--{dbtype} 1\n{record}\n//'        
    else:
        raise Exception(f"Unrecognized query type {qtype}.")

    cut_ga = False
    etime, hits = scan_hits(data, address=address, port=port,
                            cut_ga=cut_ga, evalue_thr=evalue_thr, score_thr=score_thr, max_hits=max_hits)
    
    return name, etime, hits


def scan_hits(data, address="127.0.0.1", port=51371, cut_ga=False, evalue_thr=None,
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

        #The first hit section contains all of the sequence matches
        #but not domain mataches.
        for hitid in range(nreported):
            hits_end = hits_start + 152
            name, evalue, score, ndom = unpack_hit(binresult[hits_start:hits_end], Z)
            hitdata[hitid] = {"name": name, "evalue": evalue, "score":score, "ndom": ndom, "doms": []}
            hits_start += 152
        
        next_start = hits_end
        reported_hits = []
        
        #Now, unpack the domain hits for the sequences
        for hitid in range(nreported):
            hit = hitdata[hitid]
            
            for domid in range(hit["ndom"]):
                dombit = binresult[next_start:next_start + 72]
                dom = struct.unpack("4i 5f 4x d 2i Q 8x", dombit)
                # is_reported = dom[10]
                # is_included = dom[11]
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
                is_reported = d[10] == 1
                is_included = d[11] == 1
                
                bitscore = d[8]
                ievalue = math.exp(d[9] * Z)
                cevalue = math.exp(d[9] * domZ)
                hitname = hit["name"]
                evalue = hit["evalue"]
                score = hit["score"]
                
                if (cut_ga == True or \
                    ((evalue_thr is None or evalue <= evalue_thr) and \
                    (score_thr is None or score >= score_thr))) and \
                    (max_hits is None or len(reported_hits)+1 <= max_hits or lasthitname == hitname) and \
                    is_reported and is_included:
                    
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

## END
