import sys
import tempfile
import uuid
import os
import time
from collections import Counter
import cPickle
from pymongo import MongoClient
mongoCli = MongoClient()
mongoDB = mongoCli.eggnog4_1
db_speciation = mongoDB.sp_events
db_members = mongoDB.members

PHMMER_BIN = "/g/bork1/huerta/_soft/hmmer-3.1b2/src/phmmer"

def get_best_hit(target_seq, target_og):    
    tempout = str(uuid.uuid4())
    cmd = "%s --incE 0.001 -E 0.001 -o /dev/null --noali --tblout %s %s %s" %(
        PHMMER_BIN, tempout, target_seq, target_og)
    status = os.system(cmd)
    best_hit = None
    if status == 0:
        for line in open(tempout):
            if line.startswith('#'):
                continue
            else:
                best_hit = line.split()
                break
        os.remove(tempout)
    else:
        raise ValueError('Error running')

    if best_hit:
        best_hit_name = best_hit[2]
        best_hit_evalue = best_hit[4]
        best_hit_score = best_hit[5]
        orthologs = sorted(get_grainned_orthologs_by_member([best_hit_name]))
        pname = Counter(get_preferred_names_dict(orthologs).values())
        name_ranking = sorted(pname.items(), key=lambda x:x[1], reverse=True)
        if name_ranking[0][1] > 2:
            best_name = name_ranking[0][0]
        else:
            best_name = '-'
    else:
        best_hit_evalue = '-'
        best_hit_score = '-'
        best_hit_name = '-'
        best_name = '-'
        orthologs = []
        
    return [best_hit_name, best_hit_evalue, best_hit_score, best_name, orthologs]

def get_preferred_names_dict(names):
    query = {'$or': [{"n":n.split('.', 1)[1], "t":int(n.split('.', 1)[0])} for n in names]} 
    return dict([("%s.%s" %(e['t'], e['n']), e['p']) for e in db_members.find(query, {'n':1, 't':1, 'p':1})])
    
def refine_hit(args):
    seqname, seq, group_fasta = args
    F = tempfile.NamedTemporaryFile(dir="./", delete=True)
    F.write('>%s\n%s' %(seqname, seq))
    F.flush()    
    best_hit = get_best_hit(F.name, group_fasta)
    F.close()

    return [seqname] + best_hit

def get_grainned_orthologs_by_member(target_members, target_taxa=None, target_nogs=None):
    t1 = time.time()
    query = {}
    
    target_seqids = set()
    query_taxa = set()
    for m in target_members: 
        taxid, seqid = m.split(".", 1)
        query_taxa.add(taxid)
        target_seqids.add(seqid)

    if len(target_members) >1 :
        query = {'m': {"$in":target_members}}
    else:
        query = {'m': target_members[0]}

    if target_taxa:
        target_taxa = set(map(str, target_taxa))
    else:
        target_taxa = set()
    
    # all_taxa = query_taxa | target_taxa
        
    orthology = {}
    target_members = set(target_members)
    for event in db_speciation.find(query, {'z':1, 'm':1, 'n':1, 'l':1}):
        if target_nogs and event['n'] not in target_nogs:
            continue
        
        all_seqs = [m.split(".", 1) for m in event['m']]
        by_sp1, by_sp2 = {}, {}
        for _sp, _side in [(by_sp1, all_seqs[:event['z']]),
                           (by_sp2, all_seqs[event['z']:])]:
            for t, s in _side:
                if not target_taxa or t in target_taxa or t in query_taxa:
                    mid = "%s.%s" %(t, s)
                    _sp.setdefault(t, set()).add(mid)

        # merge by side1 coorthologs
        targets = target_taxa or by_sp2.keys()
        for sp1, co1 in by_sp1.iteritems():
            if target_members & co1:
                key1 = (sp1, tuple(sorted((co1))))
                for sp2 in targets:
                    if sp2 not in by_sp2: continue
                    co2 = by_sp2[sp2]
                    key2 = (sp2, tuple(sorted(co2)))
                    orthology.setdefault(key1, set()).add(key2)
                    
        # merge by side1 coorthologs
        targets = target_taxa or by_sp1.keys()
        for sp1, co1 in by_sp2.iteritems():
            if target_members & co1:
                key1 = (sp1, tuple(sorted((co1))))
                for sp2 in targets:
                    if sp2 not in by_sp1: continue
                    co2 = by_sp1[sp2]
                    key2 = (sp2, tuple(sorted(co2)))
                    orthology.setdefault(key1, set()).add(key2)
                    
    #as_list = [((k[0], get_sp(k[0]), k[1]), [(t2, get_sp(t2), co2) for t2,co2 in v])
    #            for k,v in orthology.iteritems()]

    all_orthologs = set()
    for k, v in orthology.iteritems():
        all_orthologs.update(k[1])
        for t2, co2 in v:
            all_orthologs.update(co2)
    
    #print time.time() - t1
    return all_orthologs

def get_sp(a):
    return a

def process_hits_file(hits_file, query_fasta):
    from ete3 import SeqGroup
    from multiprocessing import Pool
    FASTA_PATH = 'OG_fasta/'    
    print "Loading OG data..."
    og2level = cPickle.load(open("og2level.pkl"))
    
    aln = SeqGroup(query_fasta)
    cmds = []
    visited_queries = set()
    for line in open(hits_file):
        if line.startswith('#'):
            continue

        fields = line.split('\t')
        seq = aln.get_seq(fields[0].strip())
        seqname = fields[0]
        hitname = fields[1]
        if seqname in visited_queries:
            continue
        visited_queries.add(seqname)
        level = og2level.get(hitname, 'unknown')            
        target_fasta = os.path.join(FASTA_PATH, level, "%s.fa" %hitname)
        cmds.append([seqname, seq, target_fasta])
    print len(cmds)
    pool = Pool(10)
    result = []
    process = pool.map_async(refine_hit, cmds, callback=result.append)
    process.wait()
    #result.append( refine_hit(cmds[0]))
    
    return result
                    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--hits', dest='hitsfile', type=str, help="hits file", required=True)
    parser.add_argument('-f', dest='fastafile', type=str, help="fasta file", required=True)
    parser.add_argument('-o', dest='output', type=str, help="output")

    args = parser.parse_args()
    
    results = process_hits_file(args.hitsfile, args.fastafile)
    if args.output:
        OUT = open(args.output, "w")
    else:
        OUT = sys.stdout
    print >>OUT, '\t'.join("#query_seq, best_hit_eggNOG_ortholog, best_hit_evalue, best_hit_score, predicted_name, strict_orthologs".split(','))
    for r in results[0]:
        print >>OUT, '\t'.join(map(str, (r[0], r[1], r[2], r[3], r[4], ','.join(r[5]))))
        
    if args.output:
        OUT.close()
