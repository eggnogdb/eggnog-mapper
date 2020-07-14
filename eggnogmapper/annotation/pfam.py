##
## CPCantalapiedra 2020

import os

from ..search.hmmer import HmmerSearcher
from ..search.hmmer_search import SCANTYPE_MEM, SCANTYPE_DISK, QUERY_TYPE_SEQ, DB_TYPE_HMM, DB_TYPE_SEQ, QUERY_TYPE_HMM
from ..common import get_pfam_db

##
def get_pfam_args(cpu, fasta_file):

    query_number = len([1 for line in open(fasta_file) if line.startswith(">")])

    print(f"pfam.py:get_pfam_args {cpu}, {query_number}")

    if query_number < 100:
        usemem = False
        num_servers = cpu
        num_workers = 1
        cpus_per_worker = 1
        scan_type = SCANTYPE_DISK
        db = get_pfam_db()
        infile = fasta_file
        dbtype = DB_TYPE_HMM
        qtype = QUERY_TYPE_SEQ
        
    elif query_number >= 100 and query_number < 15000:
        usemem = True
        num_servers = cpu
        num_workers = 1
        cpus_per_worker = 1
        scan_type = SCANTYPE_MEM
        db = get_pfam_db()
        infile = fasta_file
        dbtype = DB_TYPE_HMM
        qtype = QUERY_TYPE_SEQ
        
    elif query_number >= 15000:
        if mapfile(fasta_file):
            usemem = True
            num_servers = cpu
            num_workers = 1
            cpus_per_worker = 1
            scan_type = SCANTYPE_MEM
            db = fasta_file
            infile = get_pfam_db()
            dbtype = DB_TYPE_SEQ
            qtype = QUERY_TYPE_HMM            
        else:
            usemem = True
            num_servers = cpu
            num_workers = 1
            cpus_per_worker = 1
            scan_type = SCANTYPE_MEM
            db = get_pfam_db()
            infile = fasta_file
            dbtype = DB_TYPE_HMM
            qtype = QUERY_TYPE_SEQ
    
    return usemem, num_servers, num_workers, cpus_per_worker, scan_type, db, infile, dbtype, qtype

def mapfile(fasta_file):
    exists = False

    exists = os.path.exists(fasta_file+".map") and os.path.exists(fasta_file+".seqdb")
    
    return exists


def parse_pfam_file(pfam_file):
    pfams = {}

    with open(pfam_file, 'r') as pfamf:
        for line in pfamf:
            if line.startswith("#"): continue
            query, pfam, evalue, score, qlen, hmmfrom, hmmto, seqfrom, seqto, qcov = map(str.strip, line.split("\t"))
            if query in pfams:
                pfams[query].add(pfam)
            else:
                pfams[query] = {pfam}

    return pfams
    
##
class PfamAligner:

    args = None
    
    def __init__(self, args):
        self.args = args
        return

    
    ##
    def align_whole_pfam(self, queries_fasta, pfam_file):

        # hmmscan
        s = HmmerSearcher(self.args)
        s.search_hmm_matches(queries_fasta, pfam_file)

        return

    
    ##
    def align_query_pfams(self, query_name, queries_fasta, pfams_list):
        aligned_pfams = None # [pfam|score|evalue|start|end, ...]

        # create a file with query sequence only
        # create a hmm database with only pfams_list
        # hmmscan

        return aligned_pfams

## END
