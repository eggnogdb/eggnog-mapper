#!/usr/bin/env python

__description__ = 'Reads a fasta file containing protein sequences and searches for significant in a memory based hmmpgmd database.'
__author__ = 'Jaime Huerta Cepas'
__license__ = "GPL v2"

import sys
import socket
import struct
import math
import re
import time 
import os
import subprocess
import cPickle
from tempfile import NamedTemporaryFile
from collections import defaultdict
from hashlib import md5

from Bio import SeqIO

from server import DBDATA
import server

from refine import refine_hit

BASEPATH = os.path.split(os.path.abspath(__file__))[0]
HMMSEARCH = 'hmmsearch'
HMMSCAN = 'hmmscan'

B62_IDENTITIES = {'A': 4,
                  'B': 4,
                  'C': 9,
                  'D': 6,
                  'E': 5,
                  'F': 6,
                  'G': 6,
                  'H': 8,
                  'I': 4,
                  'K': 5,
                  'L': 4,
                  'M': 5,
                  'N': 6,
                  'P': 7,
                  'Q': 5,
                  'R': 5,
                  'S': 4,
                  'T': 5,
                  'V': 4,
                  'W': 11,
                  'X': -1,
                  'Y': 7,
                  'Z': 4}

def unpack_hit(bindata, z):
    (name, acc, desc, window_length, sort_key, score, pre_score, sum_score,
     pvalue, pre_pvalue, sum_pvalue, nexpected, nregions, nclustered, noverlaps,
     nenvelopes, ndom, flags, nreported, nincluded, best_domain, seqidx, subseq_start,
     dcl, offset) = struct.unpack("3Q I 4x d 3f 4x 3d f 9I 4Q", bindata)
    
    # print (name, acc, desc, window_length, sort_key, score, pre_score, sum_score,
    #        pvalue, pre_pvalue, sum_pvalue, nexpected, nregions, nclustered, noverlaps,
    #     nenvelopes, ndom, flags, nreported, nincluded, best_domain, seqidx, subseq_start,
    #        dcl, offset)
    
    evalue = math.exp(pvalue) * z
    return name, evalue, sum_score, ndom

def unpack_stats(bindata):
    (elapsed, user, sys, Z, domZ, Z_setby, domZ_setby, nmodels, nseqs,
     n_past_msv, n_past_bias, n_past_vit, n_past_fwd, nhits, nreported, nincluded) = struct.unpack("5d 2I 9q", bindata)

    return elapsed, nhits, Z, domZ

def unpack_domain():
    pass
    
def unpack_ali():
    pass
    
    
def scan_hits(data, address="127.0.0.1", port=51371, evalue_thr=None, max_hits=None, fixed_Z=None):
    hits = []
    hit_models = set()
    s = socket.socket()
    try:
        s.connect((address, port)) 
    except Exception, e:
        print address, port, e
        raise
    s.sendall(data)

    status = s.recv(16)
    st, msg_len = struct.unpack("I 4x Q", status)
    elapsed, nreported = 0, 0
    if st == 0:
        binresult = ''
        while len(binresult) < msg_len:
            binresult += s.recv(4096)

        elapsed, nreported, Z, domZ = unpack_stats(binresult[0:120])
        if fixed_Z:
            Z = fixed_Z
        
        hits_start = 120
        hits_end = hits_start + (152 * nreported)
        dom_start = hits_end
        
        for hitblock in xrange(hits_start, hits_end, 152):
            name, evalue, score, ndom = unpack_hit(binresult[hitblock: hitblock + 152], Z)
            if ndom:
                dom_end = dom_start + (72 * ndom)
                dombit = binresult[dom_start:dom_end]
                dom = struct.unpack( "4i 5f 4x d 2i Q 8x" * ndom, dombit)
                
                alg_start = dom_end
                dom_start = dom_end
                ndomkeys = 13
                for d in xrange(ndom):                   
                    # Decode domain info
                    off = d * ndomkeys
                    # ienv = dom[off]
                    # jenv = dom[ off + 1 ]
                    iali = dom[ off + 2 ]
                    # jali = dom[ off + 3 ]
                    #ievalue = math.exp(dom[ off + 9 ]) * Z
                    #cevalue = math.exp(dom[ off + 9 ]) * domZ
                    bitscore = dom[ off + 8 ]
                    is_reported = dom[ off + 10 ]
                    is_included = dom[ off + 11 ]
                    

                    # decode the alignment
                    alibit = binresult[alg_start : alg_start + 168]
                    
                    (rfline, mmline, csline, model, mline, aseq, ppline, N, hmmname, hmmacc,
                     hmmdesc, hmmfrom, hmmto, M, sqname, sqacc, sqdesc,
                     sqfrom, sqto, L, memsize, mem) = struct.unpack( "7Q I 4x 3Q 3I 4x 6Q I 4x Q", alibit)
                    # next domain start pos
                    alg_start += 168 + memsize
                    dom_start = alg_start
                        
                    if evalue_thr is None or evalue <= evalue_thr:
                        hit_models.add(name)
                        hits.append((name, evalue, score, hmmfrom, hmmto, sqfrom, sqto, bitscore))

            if max_hits and len(hit_models) == max_hits:
                break
    else:
        s.close()
        raise ValueError('hmmpgmd error: %s' %data[:50])

    s.close()
    return  elapsed, hits

def iter_hmm_hits(hmmfile, host, port, dbtype="hmmdb",
                  evalue_thr=None, max_hits=None, skip=None, maxseqlen=None, fixed_Z=None): 
    HMMFILE = open(hmmfile)    
    with open(hmmfile) as HMMFILE:
        while HMMFILE.tell() != os.fstat(HMMFILE.fileno()).st_size:
            model = ''
            name = 'Unknown'
            leng = None
            for line in HMMFILE:
                if line.startswith("NAME"):
                    name = line.split()[-1]
                if line.startswith("LENG"):
                    hmm_leng = int(line.split()[-1])
                model += line
                if line.strip() == '//':
                    break
                    
            if skip and name in skip:
                continue
                
            data = '@--%s 1\n%s' %(dbtype, model)
            etime, hits = scan_hits(data, host, port, evalue_thr=evalue_thr, max_hits=max_hits, fixed_X=fixed_Z)
            yield name, etime, hits, hmm_leng, None
    
def iter_seq_hits(src, host, port, dbtype,
                  evalue_thr=None, max_hits=None, skip=None, maxseqlen=None, fixed_Z=None):    
    for seqnum, record in enumerate(SeqIO.parse(src, "fasta")):        
        name = record.id
        if skip and name in skip:
            continue
            
        if maxseqlen and len(record.seq) > maxseqlen:           
            yield name, -1, [], len(record.seq), None
            continue

        if not record.seq:
            continue

        seq = str(record.seq)
        seq = re.sub("-.", "", seq)
        data = '@--%s 1\n>%s\n%s\n//' %(dbtype, name, seq)
        etime, hits = scan_hits(data, host, port, evalue_thr=evalue_thr, max_hits=max_hits, fixed_Z=fixed_Z)

        #max_score = sum([B62_IDENTITIES.get(nt, 0) for nt in seq])        
        yield name, etime, hits, len(seq), None

def iter_hits(source, query_type, dbtype, scantype, host, port,
              evalue_thr=None, max_hits=None, return_seq=False, skip=None, maxseqlen=None, fixed_Z=None, qcov_thr=None, fixex_Z=None):
    try:
        max_hits = int(max_hits)
    except Exception:
        max_hits = None
    
    if scantype == 'mem' and query_type == "seq":
        return iter_seq_hits(source, host, port, dbtype=dbtype, evalue_thr=evalue_thr, max_hits=max_hits)
    elif scantype == 'mem' and query_type == "hmm" and dbtype == "seqdb":
        return iter_hmm_hits(src, )
    elif scantype == 'disk' and query_type == "seq":
        return hmmscan(source, host)
    else:
        raise ValueError('not supported')         
            
def get_hits(name, seq, address="127.0.0.1", port=51371, dbtype='hmmdb', evalue_thr=None, max_hits=None):    
    seq = re.sub("-.", "", seq)    
    data = '@--%s 1\n>%s\n%s\n//' %(dbtype, name, seq)    
    etime, hits = scan_hits(data, address=address, port=port, evalue_thr=evalue_thr, max_hits=max_hits)
    print etime
    return name, etime, hits
            
def server_up(host, port):
    import socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex((host, port))
    sock.close()
    if result == 0:
        return True
    else: 
        return False

def server_functional(host, port, dbtype):
    try:
        get_hits("test", "TESTSEQ", host, port, dbtype)
    except Exception, e:
        print 'Server not ready', e
        return False
    return True
    
def safe_cast(v):
    try:
        return float(v)
    except ValueError:
        return v.strip()
        
def hmmscan(query_file, database_path, ncpus=10, fixed_Z=None):
    #F = NamedTemporaryFile()
    #F.write(fasta)
    #F.flush()
    OUT = NamedTemporaryFile()
    cmd = '%s --cpu %s -o /dev/null --domtblout %s %s %s' %(HMMSCAN, ncpus, OUT.name, database_path, query_file)
    
    #print cmd
    sts = subprocess.call(cmd, shell=True)
    byquery = defaultdict(list)

    last_query = None
    hit_list = []
    last_query_len = None
    if sts == 0:
        for line in OUT:
            # TBLOUT
            #['#', '---', 'full', 'sequence', '----', '---', 'best', '1', 'domain', '----', '---', 'domain', 'number', 'estimation', '----']
            #['#', 'target', 'name', 'accession', 'query', 'name', 'accession', 'E-value', 'score', 'bias', 'E-value', 'score', 'bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description', 'of', 'target']
            #['#-------------------', '----------', '--------------------', '----------', '---------', '------', '-----', '---------', '------', '-----', '---', '---', '---', '---', '---', '---', '---', '---', '---------------------']
            #['delNOG20504', '-', '553220', '-', '1.3e-116', '382.9', '6.2', '3.4e-116', '381.6', '6.2', '1.6', '1', '1', '0', '1', '1', '1', '1', '-']
            #fields = line.split() # output is not tab delimited! Should I trust this split?
            #hit, _, query, _ , evalue, score, bias, devalue, dscore, dbias = fields[0:10]
            
            #DOMTBLOUT
            #                                                                             --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
            # target name        accession   tlen query name            accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
            #------------------- ---------- -----  -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
            #Pkinase              PF00069.22   264 1000565.METUNv1_02451 -            858   4.5e-53  180.2   0.0   1   1   2.4e-56   6.6e-53  179.6   0.0     1   253   580   830   580   838 0.89 Protein kinase domain
            
            
            if line.startswith('#'):
                continue
            fields = line.split()
            (hitname, hacc, tlen, qname, qacc, qlen, evalue, score, bias, didx, dnum, c_evalue,
             i_evalue, d_score, d_bias, hmmfrom, hmmto, seqfrom, seqto, env_from, env_to, acc) = map(safe_cast, fields[:22])
            if last_query and qname != last_query:
                yield last_query, -1, hit_list, last_query_len, None
                hit_list = []
                last_query = qname
                last_query_len = None
            hit_list.append([hitname, evalue, score, hmmfrom, hmmto, seqfrom, seqto, d_score])
            last_query = qname
            if last_query_len and last_query_len != qlen:
                raise ValuerError("Inconsistent qlen when parsing hmmscan output")
            last_query_len = qlen

        yield last_query, 0, hit_list, last_query_len, None
            
    OUT.close()


def hmmsearch(query_hmm, target_db, ncpus=10):
    OUT = NamedTemporaryFile()
    cmd = '%s --cpu %s -o /dev/null -Z 1000000 --tblout %s %s %s' %(HMMSEARCH, ncpus, OUT.name, query_hmm, target_db)

    sts = subprocess.call(cmd, shell=True)
    byquery = defaultdict(list)
    if sts == 0:
        for line in OUT:
            #['#', '---', 'full', 'sequence', '----', '---', 'best', '1', 'domain', '----', '---', 'domain', 'number', 'estimation', '----']
            #['#', 'target', 'name', 'accession', 'query', 'name', 'accession', 'E-value', 'score', 'bias', 'E-value', 'score', 'bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description', 'of', 'target']
            #['#-------------------', '----------', '--------------------', '----------', '---------', '------', '-----', '---------', '------', '-----', '---', '---', '---', '---', '---', '---', '---', '---', '---------------------']
            #['delNOG20504', '-', '553220', '-', '1.3e-116', '382.9', '6.2', '3.4e-116', '381.6', '6.2', '1.6', '1', '1', '0', '1', '1', '1', '1', '-']
            if line.startswith('#'): continue
            fields = line.split() # output is not tab delimited! Should I trust this split?
            hit, _, query, _ , evalue, score, bias, devalue, dscore, dbias = fields[0:10]
            evalue, score, bias, devalue, dscore, dbias = map(float, [evalue, score, bias, devalue, dscore, dbias])
            byquery[query].append([query, evalue, score])
            
    OUT.close()
    return byquery

def get_nogname(name):
    if len(name) == 5:
        return "ENOG41%s" %name
    return name

def load_nog_lineages():
    import cPickle
    if os.path.exists('NOG_hierarchy.pkl'):
        nog2lineage = cPickle.load(open('NOG_hierarchy.pkl'))
    else:
        nog2lineage = {}
        for line in open('../build_db/NOG_hierarchy.tsv'):
            fields = line.strip().split('\t')
            nog2lineage[fields[0].split('@')[0]] = map(lambda x: tuple(x.split('@')), fields[2].split(','))
        cPickle.dump(nog2lineage, open('NOG_hierarchy.pkl', 'wb'), protocol=2)
    return nog2lineage

def generate_idmap(dbpath):
    cmd = """hmmstat %s |grep -v '#'|awk '{print $1" "$2}' > %s""" %(dbpath, dbpath+'.idmap')
    print 'Generating idmap in', dbpath+'.idmap'
    return os.system(cmd) == 0
        
def main(args):
    master_db, worker_db = None, None
    
    if args.db in DBDATA:
        db = DBDATA
        port = DBDATA[args.db]["client_port"]
        host = "localhost"
        idmap = cPickle.load(open(DBDATA[args.db]['idmap'], 'rb'))
        scantype = 'mem'
    elif os.path.isfile(args.db):
        if args.usemem:
            if not args.idmap:
                if generate_idmap(args.db):
                    args.idmap = args.db+".idmap"
                    print >>sys.stderr, "idmap succesfully created!"
                else:
                    print >>sys.stderr, "idmap could not be created!"
                    idmap = None
            else:
                idmap = None
                
            for try_port in range(52000, 53000, 2):
                print >>sys.stderr, "Loading server at localhost, port", try_port, try_port+1
                dbpath, master_db, worker_db = server.load_server(args.db, try_port, try_port+1, args.cpu)
                port = try_port
                host = 'localhost'
                ready = False
                while 1:
                    print >>sys.stderr, "Waiting for server to become ready..."
                    time.sleep(1)
                    if not server.check_pid(master_db.pid) or not server.check_pid(worker_db.pid):
                        break
                    elif server_functional(host, port, args.dbtype):
                        ready = True
                        break                   
                if ready:
                    break
            scantype = 'mem'
        else:
            idmap = None
            port = None
            host = args.db
            scantype = 'disk'
    else:
        db = None
        scantype = 'mem'
        if ":" in args.db:
            host, port = args.db.split(":")
            host = host.strip()            
            port = int(port)
        else:
            host = 'localhost'
            port = int(args.db)
        idmap = None

    if port and not server_up(host, port):
        print >>sys.stderr, "hmmpgmd Server not found at %s:%s" %(host, port)
        exit(-1)
        
    if args.idmap:
        print >>sys.stderr, "Reading idmap"
        idmap = {}
        for _lnum, _line in enumerate(open(args.idmap)):
            if not _line.strip():
                continue                
            try:
                _seqid, _seqname = map(str, _line.strip().split(' '))
            except ValueError:
                if _lnum == 0:
                    continue # idmap generate by esl_reformat has an info line at begining
                else:
                    raise                   
            _seqid = int(_seqid)
            idmap[_seqid] = [_seqname]
        print >>sys.stderr, len(idmap), "names loaded" 
    
    VISITED = set()
    if args.output:
        if args.resume:
            print "Resuming previous run. Reading computed output from", args.output
            VISITED = set([line.split('\t')[0].strip() for line in open(args.output) if not line.startswith('#')])
            print len(VISITED), 'processed queries skipped'
            OUT = open(args.output, 'a')
        else:
            OUT = open(args.output, 'w')
    else:
        OUT = sys.stdout


    HEADER = ['query', 'hit', 'e-value', 'sum_score', 'query_length', 'hmmfrom', 'hmmto', 'seqfrom', 'seqto', 'q_coverage']
    print >>OUT, '# ' + time.ctime()
    print >>OUT, '# ' + ' '.join(sys.argv)
    print >>OUT, '# ' + '\t'.join(HEADER)
        
    print >>sys.stderr, "Analysis starts now."
    total_time = 0
    last_time = time.time()
    for qn, (name, elapsed, hits, querylen, seq) in enumerate(iter_hits(args.fastafile[0],
                                                                            args.qtype,
                                                                            args.dbtype,
                                                                            scantype,
                                                                            host,
                                                                            port,
                                                                            evalue_thr=args.evalue,
                                                                            qcov_thr=args.qcov,
                                                                            fixed_Z=args.Z,                                                                                
                                                                            max_hits=args.maxhits,
                                                                            skip=VISITED,
                                                                            maxseqlen=args.maxseqlen)):

        if elapsed == -1:
            # error occured 
            print >>OUT, '\t'.join([name] + ['ERROR'] * len(HEADER))
        elif not hits:            
            print >>OUT, '\t'.join([name] + ['-'] * len(HEADER))
        else:
            for hitindex, (hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore) in enumerate(hits):
                hitname = hid
                level = "-"
                if idmap:                    
                    hitname = idmap[hid][0]                    
                print >>OUT, '\t'.join(map(str, [name, hitname, heval, hscore, querylen, hmmfrom, hmmto, sqfrom, sqto, (sqto-sqfrom)/querylen]))                                                        
        OUT.flush()

        # monitoring
        total_time += time.time() - last_time
        last_time = time.time()
        if qn and (qn % 25 == 0):
            print >>sys.stderr, qn, total_time, "%0.2f q/s" %((float(qn)/total_time))
            sys.stderr.flush()

    # finish
    print >>sys.stderr, qn, total_time, "%0.2f q/s" %((float(qn)/total_time))
    sys.stderr.flush()
    print >>OUT, '# %d queries scanned' %(qn + 1)
    print >>OUT, '# Total time (seconds):', total_time
    if args.output:
        OUT.close()

    if master_db:
        master_db.terminate()
    if worker_db:
        worker_db.terminate()
        
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    # server
    parser.add_argument('--db', required=True, dest='db', help='specify the target database for sequence searches. Choose among: euk,bact,arch, host:port, or local hmmpressed database')
    parser.add_argument('--usemem', action="store_true", help='If a local hmmpressed databased is provided as target, this flag allows to store the database in memory prior to all computations. Database is unloaded when finished.')
    parser.add_argument('--cpu', type=int, default=1)
    
    parser.add_argument('--dbtype', dest="dbtype", choices=["hmmdb", "seqdb"], default="hmmdb")
    parser.add_argument('--qtype',  choices=["hmm", "seq"], default="seq")
    parser.add_argument('--idmap', dest='idmap', type=str)
    
    parser.add_argument('--evalue', dest='evalue', default=0.001, type=float, help="e-value threshold")
    parser.add_argument('--maxhits', dest='maxhits', type=int, help="max number of hits to report")
    parser.add_argument('--maxseqlen', type=int, help="exclude query sequences larger than `maxseqlen`")
    parser.add_argument('--qcov', type=float, help="min query coverage (from 0 to 1)")
    parser.add_argument('--Z', dest='Z', type=float, help='NOT IMPLEMENTED YET')
    
    parser.add_argument('--output', type=str, help="output file")    
    parser.add_argument('--resume', action="store_true", help="Resumes a previous execution skipping reported hits in the output file.")
    
    parser.add_argument('fastafile', metavar="fastafile", nargs=1, help='query file')

    #parser.add_argument('--refine', action="store_true", dest='refine', help="Refine hits searching best protein within the matching group")
    #parser.add_argument('--cache', type=str)
    args = parser.parse_args()
    main(args)

