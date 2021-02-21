##
## CPCantalapiedra 2020

from ..search.hmmer.hmmer_seqio import iter_fasta_seqs

##
def create_prots_file(infile, hits, outfile):
    hits_dict = {}
    for hit in hits:
        query = hit[0]
        query_no_suffix = query[:query.rfind("_")]
        qstart = hit[4]
        qend = hit[5]
        if query_no_suffix in hits_dict:
            hits_dict[query_no_suffix].append((query, qstart, qend))
        else:
            hits_dict[query_no_suffix] = [(query, qstart, qend)]

    with open(outfile, 'w') as OUT:
        suffix = 0
        for name, seq in iter_fasta_seqs(infile):
            if name in hits_dict:
                for query, qstart, qend in hits_dict[name]:
                    print(f">{query}\n{seq[qstart-1:qend]}", file=OUT)
                    suffix += 1

    return

## END
