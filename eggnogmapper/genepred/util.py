##
## CPCantalapiedra 2020

from ..search.hmmer.hmmer_seqio import iter_fasta_seqs

##
def create_prots_file(infile, hits, outfile):
    hits_dict = {}
    for hit in hits:
        query = hit[0]
        qstart = hit[4]
        qend = hit[5]
        if query in hits_dict:
            hits_dict[query].append((qstart, qend))
        else:
            hits_dict[query] = [(qstart, qend)]

    with open(outfile, 'w') as OUT:
        suffix = 0
        for name, seq in iter_fasta_seqs(infile):
            if name in hits_dict:
                for qstart, qend in hits_dict[name]:
                    print(f">{name}_{suffix}\n{seq[qstart-1:qend]}", file=OUT)
                    suffix += 1

    return

##
def create_gff_file(infile, searcher_name, version, hits, outfile):
    hits_dict = {}
    with open(outfile, 'w') as OUT:

        print("##gff-version 3", file=OUT)
        print(f"# {version}", file=OUT)

        for hit in sorted(hits, key=lambda x: (x[0],x[4],x[5],x[3])):
            query = hit[0]
            target = hit[1]
            evalue = hit[2]
            score = hit[3]
            qstart = hit[4]
            qend = hit[5]
            sstart = hit[6]
            send = hit[7]
            scov = hit[9]
            if qstart <= qend:
                strand = "+"
            else:
                strand = "-"
                qend = hit[4]
                qstart = hit[5]

            frame = "-" # we cannot know the frame as we align against proteins

            if query in hits_dict:
                hits_dict[query] += 1
            else:
                hits_dict[query] = 0
            suffix = hits_dict[query]

            print(f"{query}\teggNOG-mapper\tCDS\t{qstart}\t{qend}\t{score}\t{strand}\t{frame}\t"
                  f"ID={query}_{suffix};score={score};evalue={evalue};eggNOG_target={target};target_cov={scov};"
                  f"sstart={sstart};send={send};searcher={searcher_name}",
                  file=OUT)

    return

## END
