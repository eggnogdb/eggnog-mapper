##
## CPCantalapiedra 2020

import os

from argparse import Namespace

from ..search.hmmer import HmmerSearcher
from ..search.hmmer_search import SCANTYPE_MEM, SCANTYPE_DISK, QUERY_TYPE_SEQ, DB_TYPE_HMM, DB_TYPE_SEQ, QUERY_TYPE_HMM
from ..common import get_pfam_db, get_call_info


##
def get_hmmsearch_args(cpu, fasta_file, hmm_file, translate, temp_dir):
    usemem = False
    num_servers = cpu
    num_workers = 1
    cpus_per_worker = 1
    scan_type = SCANTYPE_DISK
    db = fasta_file
    infile = hmm_file
    dbtype = DB_TYPE_SEQ
    qtype = QUERY_TYPE_HMM

    pfam_args = Namespace(call_info = get_call_info(),
                          cpu = cpu,
                          usemem = usemem,
                          num_servers = num_servers,
                          num_workers = num_workers,
                          cpus_per_worker = cpus_per_worker,
                          scan_type = scan_type,
                          db = db,
                          servers_list = None,
                          dbtype = dbtype,
                          qtype = qtype,
                          translate = translate,
                          resume = False,
                          no_file_comments = False,
                          maxhits = 0, # unlimited
                          report_no_hits = False,
                          maxseqlen = 5000,
                          cut_ga = True,
                          clean_overlaps = "clans",
                          evalue = 1E-10,
                          score = None,
                          qcov = 0,
                          Z = 40000000,
                          temp_dir = temp_dir,
                          excluded_taxa = None)
    
    return pfam_args, infile
    
##
def get_pfam_args(cpu, fasta_file, translate, temp_dir):

    query_number = len([1 for line in open(fasta_file) if line.startswith(">")])
    
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
        
    else: #elif query_number >= 15000:
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
    
    pfam_args = Namespace(call_info = get_call_info(),
                          cpu = cpu,
                          usemem = usemem,
                          num_servers = num_servers,
                          num_workers = num_workers,
                          cpus_per_worker = cpus_per_worker,
                          scan_type = scan_type,
                          db = db,
                          servers_list = None,
                          dbtype = dbtype,
                          qtype = qtype,
                          translate = translate,
                          resume = False,
                          no_file_comments = False,
                          maxhits = 0, # unlimited
                          report_no_hits = False,
                          maxseqlen = 5000,
                          cut_ga = True,
                          clean_overlaps = "clans",
                          evalue = 1E-10,
                          score = None,
                          qcov = 0,
                          Z = 40000000,
                          temp_dir = temp_dir,
                          excluded_taxa = None)
    
    # return usemem, num_servers, num_workers, cpus_per_worker, scan_type, db, infile, dbtype, qtype    
    return pfam_args, infile


def mapfile(fasta_file):
    exists = False

    exists = os.path.exists(fasta_file+".map") and os.path.exists(fasta_file+".seqdb")
    
    return exists


def group_queries_pfams(all_annotations, PFAM_COL):
    queries_pfams_groups = []

    # Extract list of queries for each pfam
    pfams_queries = {}
    for annot_columns in all_annotations:
        query_name = annot_columns[0]
        pfams = annot_columns[PFAM_COL]
        
        if pfams == "-": continue
        
        for pfam in pfams.split(","):
            if pfam in pfams_queries:
                pfams_queries[pfam].add(query_name)
            else:
                pfams_queries[pfam] = {query_name}

    # Re-group Pfams with the same list of queries
    queries_pfams_keys = {}
    for pfam, queries in pfams_queries.items():
        pq_key = ",".join(queries)
        if pq_key in queries_pfams_keys:
            queries_pfams_keys[pq_key]["pfams"].add(pfam)
        else:
            queries_pfams_keys[pq_key] = {"pfams":{pfam}, "queries":queries}

    # Return tuples of pfams-queries
    return [(x["queries"], x["pfams"]) for x in queries_pfams_keys.values()]


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


def parse_hmmsearch_file(pfam_file):
    pfams = {}

    with open(pfam_file, 'r') as pfamf:
        for line in pfamf:
            if line.startswith("#"): continue
            pfam, query, evalue, score, qlen, hmmfrom, hmmto, seqfrom, seqto, qcov = map(str.strip, line.split("\t"))
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
    def align_whole_pfam(self, infile, pfam_file, silent = False):

        # hmmscan
        s = HmmerSearcher(self.args)
        s.search_hmm_matches(infile, pfam_file, silent)

        return

## END
