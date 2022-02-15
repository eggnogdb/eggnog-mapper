##
## CPCantalapiedra 2020

from Bio.Seq import Seq

from ..search.hmmer.hmmer_seqio import iter_fasta_seqs

##
# Create a fasta file from search hits from blastx-based gene prediction

def create_prots_file(infile, hits, outfile, translate, table):
    
    seqs_dict = {name:seq for name, seq in iter_fasta_seqs(infile)}

    with open(outfile, 'w') as OUT:
        for hit in hits:
            query = hit[0]
            query_no_suffix = query[:query.rfind("_")]

            if query_no_suffix in seqs_dict:
                seq = seqs_dict[query_no_suffix]
                qstart = hit[4]
                qend = hit[5]

                # reverse strand
                if qend < qstart:
                    qstart = hit[5]
                    qend = hit[4]
                    orf = Seq(seq[qstart-1:qend]).reverse_complement()
                
                else:
                    orf = Seq(seq[qstart-1:qend])

                if translate == True:
                    table = table if table is not None else 1
                    orf = orf.translate(table = table)
                    
                print(f">{query}\n{orf}", file=OUT)
                
            yield hit
            
    return
        
        

## END
