##
## CPCantalapiedra 2020

import os
import subprocess
from tempfile import NamedTemporaryFile

from ...common import HMMFETCH

##
# Create a new temp fasta file including only
# the queries in the list "queries"
def filter_fasta_file(queries, orig_fasta, temp_dir):
    
    Q = NamedTemporaryFile(mode='w', dir=temp_dir)
    with open(orig_fasta, 'r') as origf:
        found = False
        for line in origf:
            if line.startswith(">"):
                orig_query = line.strip()[1:].split(" ")[0]
                if orig_query in queries:
                    found = True
                    print(f">{orig_query}", file=Q)
                else:
                    found = False
            else:
                if found == True:
                    print(f"{line.strip()}", file=Q)
    Q.flush()
    return Q


##
def filter_fasta_hmm_files(queries_pfams, orig_fasta, orig_hmm, temp_dir):
    
    new_fasta = new_hmm = None
    
    queries = queries_pfams[0]
    pfams = queries_pfams[1]
    
    # Process fasta queries
    Q = filter_fasta_file(queries, orig_fasta, temp_dir)
    
    # Process pfams
    P = NamedTemporaryFile(mode='w', dir=temp_dir)
    for pfam in pfams:
        print(pfam, file=P)
    P.flush()

    H = None
    if os.stat(P.name).st_size > 0:
        H = NamedTemporaryFile(mode='w', dir=temp_dir)
        cmd = f"{HMMFETCH} -f {orig_hmm} {P.name}"
        cp = subprocess.run(cmd, shell=True, stdout=H, stderr=subprocess.DEVNULL)
        if os.stat(H.name).st_size == 0:
            H = None
            
    P.close()
    
    return Q, H

def group_queries_pfams(all_annotations, PFAM_COL):
    queries_pfams_groups = []

    # Extract list of pfams for each query
    queries_pfams = {}
    for annot_columns in all_annotations:
        query_name = annot_columns[0]
        pfams = annot_columns[PFAM_COL]

        if pfams == "-": continue

        for pfam in sorted(pfams.split(",")):
            if query_name in queries_pfams:
                queries_pfams[query_name].add(pfam)
            else:
                queries_pfams[query_name] = {pfam}
                
    # print(f"pfam.py:group:queries_pfams {queries_pfams['1105367.SAMN02673274.CG50_07170']}")

    # Re-group Pfams with the same list of queries
    queries_pfams_keys = {}
    for query, pfams in queries_pfams.items():
        pq_key = ",".join(pfams)
        if pq_key in queries_pfams_keys:
            queries_pfams_keys[pq_key]["queries"].add(query)
        else:
            queries_pfams_keys[pq_key] = {"queries":{query}, "pfams":pfams}
            
        # if "CG50_09910" in query:
        #     print(f"pfam.py:group_queries_pfams {queries_pfams_keys[pq_key]}")

    # Return tuples of pfams-queries
    return [(x["queries"], x["pfams"]) for x in queries_pfams_keys.values()]

##
# This is deprecated, since given that one query
# can be processed separatedly for different pfams
# we cannot process overlaps afterwards.
# Thus, this function was changed for the one above
# def group_queries_pfams(all_annotations, PFAM_COL):
#     queries_pfams_groups = []

#     # Extract list of queries for each pfam
#     pfams_queries = {}
#     for annot_columns in all_annotations:
#         query_name = annot_columns[0]
#         pfams = annot_columns[PFAM_COL]
        
#         if pfams == "-": continue
        
#         for pfam in pfams.split(","):
#             if pfam in pfams_queries:
#                 pfams_queries[pfam].add(query_name)
#             else:
#                 pfams_queries[pfam] = {query_name}

#     # Re-group Pfams with the same list of queries
#     queries_pfams_keys = {}
#     for pfam, queries in pfams_queries.items():
#         pq_key = ",".join(queries)
#         if pq_key in queries_pfams_keys:
#             queries_pfams_keys[pq_key]["pfams"].add(pfam)
#         else:
#             queries_pfams_keys[pq_key] = {"pfams":{pfam}, "queries":queries}

#     # Return tuples of pfams-queries
#     return [(x["queries"], x["pfams"]) for x in queries_pfams_keys.values()]

##
def wrap_group_queries_pfams(annotations, PFAM_COL, queries_fasta, pfam_db, translate, temp_dir, pfam_file):
    # queries_pfams_groups = group_queries_pfams(annotations, PFAM_COL)

    for queries_pfams_group in group_queries_pfams(annotations, PFAM_COL):
        yield (queries_pfams_group, queries_fasta, pfam_db, temp_dir, translate, pfam_file, 1)
        # cpu = 1 since we are already parallelizing

    return

## END
