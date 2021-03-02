##
## CPCantalapiedra 2020

from ..search.hmmer.hmmer_seqio import iter_fasta_seqs

##
# Create a fasta file from search hits from blastx-based gene prediction

def create_prots_file(infile, hits, outfile):
    
    seqs_dict = {name:seq for name, seq in iter_fasta_seqs(infile)}

    with open(outfile, 'w') as OUT:
        for hit in hits:
            query = hit[0]
            query_no_suffix = query[:query.rfind("_")]
            qstart = hit[4]
            qend = hit[5]

            if query_no_suffix in seqs_dict:
                seq = seqs_dict[query_no_suffix]
                print(f">{query}\n{seq[qstart-1:qend]}", file=OUT)
                
            yield hit
            
    return
        
        

## END
