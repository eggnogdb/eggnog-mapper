##
## CPCantalapiedra 2020

import os
import subprocess
from argparse import Namespace

from ...common import get_pfam_db, get_call_info, ESL_REFORMAT

from ...search.hmmer.hmmer import HmmerSearcher
from ...search.hmmer.hmmer_search import SCANTYPE_MEM, SCANTYPE_DISK, QUERY_TYPE_SEQ, DB_TYPE_HMM, DB_TYPE_SEQ, QUERY_TYPE_HMM
from ...search.hmmer.hmmer_setup import DEFAULT_PORT, DEFAULT_END_PORT

##
def get_hmmscan_args(cpu, fasta_file, hmm_file, translate, temp_dir):
    usemem = False
    scan_type = SCANTYPE_DISK
    db = hmm_file
    infile = fasta_file
    dbtype = DB_TYPE_HMM
    qtype = QUERY_TYPE_SEQ

    # port, end_port, num_servers, num_workers, cpus_per_worker are not really used
    # for SCANTYPE_DISK
    pfam_args = Namespace(call_info = get_call_info(),
                          cpu = cpu,
                          usemem = usemem,
                          port = -1,
                          end_port = -1,
                          num_servers = -1,
                          num_workers = -1,
                          cpus_per_worker = -1,
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
                          maxseqlen = None, #5000,
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
def get_hmmsearch_args(cpu, fasta_file, hmm_file, translate, temp_dir):
    usemem = False
    scan_type = SCANTYPE_DISK
    db = fasta_file
    infile = hmm_file
    dbtype = DB_TYPE_SEQ
    qtype = QUERY_TYPE_HMM

    # port, end_port, num_servers, num_workers, cpus_per_worker are not really used
    # for SCANTYPE_DISK
    pfam_args = Namespace(call_info = get_call_info(),
                          cpu = cpu,
                          usemem = usemem,
                          port = -1,
                          end_port = -1,
                          num_servers = -1,
                          num_workers = -1,
                          cpus_per_worker = -1,
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
                          maxseqlen = None, #5000,
                          cut_ga = True,
                          clean_overlaps = "hmmsearch_clans",
                          evalue = 1E-10,
                          score = None,
                          qcov = 0,
                          Z = 40000000,
                          temp_dir = temp_dir,
                          excluded_taxa = None)
    
    return pfam_args, infile
    
##
def get_pfam_args(cpu, num_servers, num_workers, cpus_per_worker, port, end_port, fasta_file, translate, temp_dir, force_seqdb = False):

    query_number = len([1 for line in open(fasta_file) if line.startswith(">")])

    # port, end_port, num_servers, num_workers, cpus_per_worker are not really used
    # for SCANTYPE_DISK, only for SCANTYPE_MEM
    port = port
    end_port = end_port
    num_servers = num_servers
    num_workers = num_workers
    cpus_per_worker = cpus_per_worker
    
    # If query number < 100, use hmmscan
    if query_number < 100:
        usemem = False
        scan_type = SCANTYPE_DISK
        db = get_pfam_db()
        infile = fasta_file
        dbtype = DB_TYPE_HMM
        qtype = QUERY_TYPE_SEQ
        clean_overlaps = "clans"

    # if query number between 100 and 15000, use hmmpgmd (hmmsearch)
    elif query_number >= 100 and query_number < 15000:

        # if not mapfile(fasta_file):
        #     create_fasta_hmmpgmd_db(fasta_file)
        #     print(f"CREATED FASTA FILE DB {fasta_file}")
        # usemem = True
        # scan_type = SCANTYPE_MEM
        # db = fasta_file
        # infile = get_pfam_db()
        # dbtype = DB_TYPE_SEQ
        # qtype = QUERY_TYPE_HMM
        # clean_overlaps = "hmmsearch_clans"

        usemem = True
        scan_type = SCANTYPE_MEM
        db = get_pfam_db()
        infile = fasta_file
        dbtype = DB_TYPE_HMM
        qtype = QUERY_TYPE_SEQ
        clean_overlaps = "clans"
        
                
    else: #elif query_number >= 15000:
        if mapfile(fasta_file):
            usemem = True
            scan_type = SCANTYPE_MEM
            db = fasta_file
            infile = get_pfam_db()
            dbtype = DB_TYPE_SEQ
            qtype = QUERY_TYPE_HMM
            clean_overlaps = "hmmsearch_clans"
            
        else:
            if force_seqdb == True and create_fasta_hmmpgmd_db(fasta_file):
                usemem = True
                scan_type = SCANTYPE_MEM
                db = fasta_file
                infile = get_pfam_db()
                dbtype = DB_TYPE_SEQ
                qtype = QUERY_TYPE_HMM
                clean_overlaps = "hmmsearch_clans"

                print("CREATED ESL_REFORMAT DB. USING HMMSEARCH")
                
            else:
                usemem = True
                scan_type = SCANTYPE_MEM
                db = get_pfam_db()
                infile = fasta_file
                dbtype = DB_TYPE_HMM
                qtype = QUERY_TYPE_SEQ
                clean_overlaps = "clans"
    
    pfam_args = Namespace(call_info = get_call_info(),
                          cpu = cpu,
                          usemem = usemem,
                          port = port,
                          end_port = end_port,
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
                          maxseqlen = None, #5000,
                          cut_ga = True,
                          clean_overlaps = clean_overlaps,
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


def create_fasta_hmmpgmd_db(fasta_file):
    cmd = f"{ESL_REFORMAT} hmmpgmd {fasta_file} > {fasta_file}.seqdb"
    cp = subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # the previous command should create also a f"{fasta_file}.map" file
    
    return mapfile(fasta_file)



def parse_pfam_file(pfam_file):
    pfams = {}

    with open(pfam_file, 'r') as pfamf:
        for line in pfamf:
            if line.startswith("#"): continue
            query, pfam, evalue, score, qlen, hmmfrom, hmmto, seqfrom, seqto, qcov = map(str.strip, line.split("\t"))

            # if "CG50_07170" in query:
            #     print(f"pfam.py:parse_pfam_file: {pfam}")
                
            if query in pfams:
                pfams[query].add(pfam)
            else:
                pfams[query] = {pfam}
                
            # if "CG50_07170" in query:
            #     print(f"pfam.py:parse_pfam_file: {pfams}")

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

def parse_hmmscan_file(pfam_file):
    pfams = {}

    with open(pfam_file, 'r') as pfamf:
        for line in pfamf:
            if line.startswith("#"): continue
            query, pfam, evalue, score, qlen, hmmfrom, hmmto, seqfrom, seqto, qcov = map(str.strip, line.split("\t"))
                
            if query in pfams:
                pfams[query].add(pfam)
            else:
                pfams[query] = {pfam}

            # if "CG50_08330" in query:
            #     print(f"pfam.py:parse_hmm_scan_file: {pfams}")

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
