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

BASEPATH = os.path.split(os.path.abspath(__file__))[0]
HMMSEARCH = 'hmmsearch'
HMMSCAN = 'hmmscan'

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
    
    
def scan_hits(data, address="127.0.0.1", port=51371, evalue_thr=None, max_hits=None):
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

def iter_hits(msf, msfformat='fasta', address="127.0.0.1", port=51371, dbtype='hmmdb', evalue_thr=None, max_hits=None, return_seq=False, skip=None, maxseqlen=None, cache=None):
    if cache:
        seqnum2md5 = {}
        for seqnum, record in enumerate(SeqIO.parse(msf, msfformat)):
            seqnum2md5[seqnum] = md5(str(record.seq)).hexdigest()
        
    try:
        max_hits = int(max_hits)
    except Exception:
        max_hits = None
        
    for seqnum, record in enumerate(SeqIO.parse(msf, msfformat)):                
        name = record.id
        if skip and name in skip:
            continue
        if maxseqlen and len(record.seq) > maxseqlen:           
            yield name, -1, [], None
            continue

        if not record.seq:
            continue

        if cache and seqnum2md5[seqnum] in cached_hits:
            name, evalue, score, hmmfrom, hmmto, sqfrom, sqto, bitscore = cached_hits[seqnum2md5[seqnum]]
            if evalue_thr is None or evalue <= evalue_thr:
                for h in cached_hits:
                    hit_models.add(name)
                    hits.append((name, evalue, score, hmmfrom, hmmto, sqfrom, sqto, bitscore))
                    if max_hits and len(hit_models) == max_hits:
                        break
        else:
            seq = str(record.seq)
            seq = re.sub("-.", "", seq)
            data = '@--%s 1\n>%s\n%s\n//' %(dbtype, name, seq)
            etime, hits = scan_hits(data, address=address, port=port, evalue_thr=evalue_thr, max_hits=max_hits)
            
        if return_seq: 
            yield name, etime, hits, seq
        else:
            yield name, etime, hits, None

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

def hmmscan(fasta, database_path, ncpus=10):
    F = NamedTemporaryFile()
    F.write(fasta)
    F.flush()
    OUT = NamedTemporaryFile()
    cmd = '%s --cpu %s -o /dev/null -Z 10000000 --tblout %s %s %s' %(HMMSCAN, ncpus, OUT.name, database_path, F.name)
    #print cmd
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
            byquery[query].append([hit, evalue, score])
            
    OUT.close()
    F.close()
    return byquery

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
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='host', default='127.0.0.1')
    parser.add_argument('-p', dest='port', default=51371, type=int)
    parser.add_argument('--db', dest='db', required=True, choices=DBDATA.keys(), help='specify the target database for sequence searches')
    
    parser.add_argument('--leveldb', dest='level', help='specify a specific taxonomic level database')
    
    parser.add_argument('--evalue', dest='evalue', default=0.001, type=float, help="e-value threshold")
    parser.add_argument('--maxhits', dest='maxhits', type=int, help="max number of hits to report")
    parser.add_argument('--output', type=str, help="output file")
    parser.add_argument('--maxseqlen', type=int, help="exclude query sequences larger than `maxseqlen`")
    parser.add_argument('--resume', action="store_true", help="Resumes a previous execution skipping reported hits in the output file.")
    
    parser.add_argument('fastafile', metavar="fastafile", nargs=1, help='query file')

    parser.add_argument('--refine', action="store_true", dest='refine', help="Refine hits using EggNOG hierarchical group lineages (experimental)")
    parser.add_argument('--refine_method', type=int, default=2, help="(experimental)")
    parser.add_argument('--cache', type=str)
    
    args = parser.parse_args()
    VISITED = set()
    if args.output:
        if args.resume:
            print "Resuming previous run. Reading computed output from", args.output
            VISITED = set([line.split('\t')[0].strip() for line in open(args.output) if not line.startswith('#')])
            print len(VISITED), 'processed queries'
            OUT = open(args.output, 'a')
        else:
            OUT = open(args.output, 'w')
    else:
        OUT = sys.stdout

    if args.refine:        
        from pymongo import MongoClient
        mongoCli = MongoClient()
        mongoDB = mongoCli.eggnog4_1
        db_nog_lineages = mongoDB.nog_lineages

    args.port = DBDATA[args.db]['client_port']
    idmap = cPickle.load(open(DBDATA[args.db]['idmap'], 'rb'))
    
    if not server_up(args.host, args.port):
        print >>sys.stderr, "hmmpgmd Server not found at %s:%s" %(args.host, args.port)
        exit(-1)
        

    print >>OUT, '# ' + time.ctime()
    print >>OUT, '# ' + ' '.join(sys.argv)
    print >>OUT, '# ' + '\t'.join(['query', 'hit', 'e-value', 'sum_score', 'hmmfrom', 'hmmto', 'seqfrom', 'seqto', 'domain_score'])
        
    total_time = 0
    print >>sys.stderr, "Analysis starts now" 
    for qn, (name, elapsed, hits, seq) in enumerate(iter_hits(args.fastafile[0], address=args.host, port=args.port, dbtype='hmmdb',
                                                         evalue_thr=args.evalue, max_hits=args.maxhits, return_seq=args.refine, skip=VISITED, maxseqlen=args.maxseqlen)):

        if elapsed >= 0:
            total_time += elapsed
                    
        if args.refine:
            t1 = time.time()
            if args.refine_method == 1:
                all_hits = []
                refine_databases = set()
                for h in hits:
                    hitname = h[0]
                    if idmap: 
                        hitname = idmap[h[0]][0]
                    nog_database_path = os.path.join('/kappa/data/eggnog41/build_db/hmm_db/', '%s@NOG'%hitname, '%s@NOG.hmm'%hitname)
                    refine_databases.add(nog_database_path)
                print refine_databases
                for dbpath in refine_databases:
                    if os.path.exists(dbpath+'.h3f'):
                        for query, refinehits in hmmscan('>%s\n%s' %(name, seq), dbpath).iteritems():
                            for subhit in refinehits:
                                all_hits.append([subhit[2], subhit[1], subhit[0]])
            else:
                all_hits = []
                hmm_models = set()
                hit_names = [idmap[h[0]][0] for h in hits]
                for match in db_nog_lineages.find({'g': {'$in':hit_names} }, {'f':1}):
                    hmm_models.update(match['f'])
                hmm_models.discard('')

                F = NamedTemporaryFile()
                F.write('>%s\n%s' %(name, seq))
                F.flush()
                dbpath = F.name
                for i, hmm_file in enumerate(hmm_models):
                    level = hmm_file.split('.', 1)[0]
                    query_hmm = os.path.join('/kappa/data/eggnog41/build_db/final_files/%s_hmm/%s' %(level, hmm_file))
                    print '\r%05d/%05d' %(i, len(hmm_models)),
                    sys.stdout.flush()
                    for query, refinehits in hmmsearch(query_hmm, dbpath).iteritems():
                        for subhit in refinehits:
                            all_hits.append([subhit[2], subhit[1], subhit[0]])
                F.close()
            print time.time() - t1
            # dump hits
            all_hits.sort(reverse=True) # high scores first
            maxhits = args.maxhits if args.maxhits else len(all_hits)
            for h in all_hits[:maxhits]:
                print >>OUT, '\t'+ '\t'.join(map(str, [name, h[2], h[1], h[0]]))
            print 
        else: 
            if elapsed == -1:
                # error occured 
                print >>OUT, '\t'.join([name, 'ERROR', 'ERROR', 'ERROR', 'ERROR', 'ERROR', 'ERROR', 'ERROR', 'ERROR'])
            elif not hits:
                print >>OUT, '\t'.join([name, '-', '-', '-', '-', '-', '-', '-', '-',])
            else:
                for hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore in hits:
                    hitname = hid
                    if idmap: 
                        hitname = idmap[hid][0]
                    print >>OUT, '\t'.join(map(str, [name, hitname, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore]))
        
        OUT.flush()
        if qn and (qn % 25 == 0):
            print >>sys.stderr, qn, total_time, "%0.2f q/s" %((float(qn)/total_time))
            sys.stderr.flush()
            
    print >>sys.stderr, qn, total_time, "%0.2f q/s" %((float(qn)/total_time))
    sys.stderr.flush()
    print >>OUT, '# %d queries scanned' %(qn + 1)
    print >>OUT, '# Total time (seconds):', total_time
    if args.output:
        OUT.close()


