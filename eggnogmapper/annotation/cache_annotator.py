##
# CPCantalapiedra 2020

from os.path import exists as pexists
from os.path import join as pjoin

from ..common import get_call_info
from ..emapperException import EmapperException
from ..utils import colorify
from ..search.hmmer.hmmer_seqio import iter_fasta_seqs


##
def md5_seqs_dict(fasta_file):
    from hashlib import md5
    md5_queries = {}

    for name, seq in iter_fasta_seqs(fasta_file):
        md5_seq = md5(seq.encode('utf-8')).hexdigest()
        md5_queries[md5_seq] = {"name":name, "seq":seq, "found":0}

    return md5_queries

##
class CacheAnnotator:

    no_file_comments = cpu = None
    
    queries_fasta = None
    temp_dir = None
    
    ##
    def __init__(self, args):
        
        self.no_file_comments = args.no_file_comments
        self.cpu = args.cpu
        
        self.queries_fasta = args.input
        self.temp_dir = args.temp_dir
        
        return

    ##
    def annotate(self, cache_file, annot_file, no_annot_file):
        
        print(colorify("Functional annotation from cached files starts now", 'green'))

        md5_queries = md5_seqs_dict(self.queries_fasta)
        
        if pexists(cache_file):
            OUT = open(annot_file, "w")
            if not self.no_file_comments:
                print(get_call_info(), file=OUT)

            with open(cache_file, 'r') as cached_annot:
                for line in cached_annot:
                    if line.startswith("#query"):
                        data = line.strip().split("\t")
                        cached_md5 = data[-1]
                        if cached_md5 != "md5":
                            print(colorify("WARNING: last column of cached annotations file should contain the md5 field. The last column name in header is not called 'md5'", 'red'))
                        print(line, file=OUT)

                    if line.startswith("#"): continue
                    data = line.strip().split("\t")
                    cached_md5 = data[-1]
                    if cached_md5 in md5_queries:
                        md5_queries[cached_md5]["found"] = 1
                        print(line, file=OUT)
            OUT.close()
        else:
            print(colorify(f"Skipping cached annotation. {cached_annot_fn} file not found", 'orange'))

        OUT = open(no_annot_file, "w")        
        for name, seq in [(md5_queries[md5]["name"], md5_queries[md5]["seq"]) for md5 in md5_queries if md5_queries[md5]["found"] == 0]:
            print(f">{name}", file=OUT)
            print(seq, file=OUT)
        OUT.close()
            
        return
    
## END
