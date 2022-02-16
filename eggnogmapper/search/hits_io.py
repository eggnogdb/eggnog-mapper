##
## CPCantalapiedra 2019

import time

from ..common import get_call_info

##
# Generator of hits from filename
def parse_seeds(filename):
    for line in open(filename, 'r'):
        if line.startswith('#') or not line.strip():
            continue

        line = list(map(str.strip, line.split('\t')))

        # short hits
        # query, target, evalue, score
        if len(line) == 4: # short hits
            
            hit = [line[0], line[1], float(line[2]), float(line[3])]

        # full hits
        # query, target, evalue, score,
        # qstart, qend, sstart, send
        # pident, qcov, scov
        elif len(line) == 11:
            
            hit = [line[0], line[1], float(line[2]), float(line[3]),
                   int(line[4]), int(line[5]), int(line[6]), int(line[7]),
                   float(line[8]), float(line[9]), float(line[10])]

        yield hit
    return

##
# Receives an iterable of hits to output
# and also returns a generator object of hits
def output_seeds(cmds, hits, out_file, no_file_comments, outfmt_short, change_seeds_coords = False):
    start_time = time.time()

    with open(out_file, 'w') as OUT:

        # comments
        if not no_file_comments:
            print(get_call_info(), file=OUT)
            if cmds is not None:
                for cmd in cmds:
                    print('##'+cmd, file=OUT)

        # header
        if outfmt_short == True:
            print('#'+"\t".join("qseqid sseqid evalue bitscore".split(" ")), file=OUT)
        else:
            print('#'+"\t".join(("qseqid sseqid evalue bitscore qstart qend "
                                 "sstart send pident qcov scov").split(" ")), file=OUT)
            
            
        qn = 0
        for hit in hits:
            # change seeds coordinates relative to the ORF, not to the contig (to use them for the .seed_orthologs file)
            # this is done mainly for blastx hits, for which we want to register the ORF-relative coordinates in the seeds
            # but use the contig-relative coordinates for the GFF
            if change_seeds_coords == True:
                orig_hit = hit
                hit = change_seed_coords(orig_hit)
            else:
                orig_hit = hit

            print('\t'.join(map(str, hit)), file=OUT)

            # always yield the hit
            yield orig_hit
            qn += 1
            
        elapsed_time = time.time() - start_time
        if not no_file_comments:
            print('## %d queries scanned' % (qn), file=OUT)
            print('## Total time (seconds):', elapsed_time, file=OUT)
            print('## Rate:', "%0.2f q/s" % ((float(qn) / elapsed_time)), file=OUT)
    return

def change_seed_coords(hit):
    [query, target, evalue, score, qstart, qend, sstart, send, pident, qcov, scov] = hit
    if qstart <= qend:
        qend = qend - (qstart - 1)
        qstart = 1
    else:
        qstart = qstart - (qend - 1)
        qend = 1
    return [query, target, evalue, score, qstart, qend, sstart, send, pident, qcov, scov]

## END
