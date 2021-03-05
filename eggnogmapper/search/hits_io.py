##
## CPCantalapiedra 2019

import time

from ..common import get_call_info

##
# Generator of hits from filename
def parse_hits(filename):
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
def output_hits(cmds, hits, out_file, resume, no_file_comments, outfmt_short):
    start_time = time.time()
    
    if resume == True:
        file_mode = 'a'
    else:
        file_mode = 'w'

    with open(out_file, file_mode) as OUT:

        # comments
        if not no_file_comments:
            print(get_call_info(), file=OUT)
            if cmds is not None:
                for cmd in cmds:
                    print('##'+cmd, file=OUT)

        # header (only first time, not for further resume)
        if file_mode == 'w':
            if outfmt_short == True:
                print('#'+"\t".join("qseqid sseqid evalue bitscore".split(" ")), file=OUT)
            else:
                print('#'+"\t".join(("qseqid sseqid evalue bitscore qstart qend "
                                     "sstart send pident qcov scov").split(" ")), file=OUT)
            
            
        qn = 0
        # rows
        # (hit, skip): hits are wrapped in a tuple with a boolean flag
        # to be output or not
        for hit, skip in hits:
            # only print the hit if not already present and --resume
            if skip == False:
                print('\t'.join(map(str, hit)), file=OUT)

            # always yield the hit    
            yield hit
            qn += 1
            
        elapsed_time = time.time() - start_time
        if not no_file_comments:
            print('## %d queries scanned' % (qn), file=OUT)
            print('## Total time (seconds):', elapsed_time, file=OUT)
            print('## Rate:', "%0.2f q/s" % ((float(qn) / elapsed_time)), file=OUT)
    return

## END
