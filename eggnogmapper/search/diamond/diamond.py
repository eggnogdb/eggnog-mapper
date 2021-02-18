##
## CPCantalapiedra 2019

from os.path import join as pjoin
import shutil
import subprocess
from tempfile import mkdtemp, mkstemp
import uuid

from ...emapperException import EmapperException
from ...common import DIAMOND, get_eggnog_dmnd_db, get_call_info, ITYPE_CDS, ITYPE_PROTS, ITYPE_GENOME, ITYPE_META
from ...utils import colorify, translate_cds_to_prots

from ..hmmer.hmmer_seqio import iter_fasta_seqs

SENSMODE_FAST = "fast"
SENSMODE_MID_SENSITIVE = "mid-sensitive"
SENSMODE_SENSITIVE = "sensitive"
SENSMODE_MORE_SENSITIVE = "more-sensitive"
SENSMODE_VERY_SENSITIVE = "very-sensitive"
SENSMODE_ULTRA_SENSITIVE = "ultra-sensitive"
# sens modes in diamond 0.9.24
# SENSMODES = [SENSMODE_FAST, SENSMODE_SENSITIVE, SENSMODE_MORE_SENSITIVE]
# sens modes in diamond 2.0.4
SENSMODES = [SENSMODE_FAST, SENSMODE_MID_SENSITIVE, SENSMODE_SENSITIVE, SENSMODE_MORE_SENSITIVE, SENSMODE_VERY_SENSITIVE, SENSMODE_ULTRA_SENSITIVE]

OVERLAP_TOL_FRACTION = 1/3

def create_diamond_db(dbprefix, in_fasta):
    cmd = (
        f'{DIAMOND} makedb --in {in_fasta} --db {dbprefix}'
    )

    print(colorify('  '+cmd, 'yellow'))
    try:
        completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
    except subprocess.CalledProcessError as cpe:
        raise EmapperException("Error running diamond: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])
        
    return

class DiamondSearcher:

    name = "diamond"
    
    # Command
    cpu = tool = dmnd_db = temp_dir = no_file_comments = None
    matrix = gapopen = gapextend = None
    block_size = index_chunks = None

    # Filters
    pident_thr = score_thr = evalue_thr = query_cov = subject_cov = None

    # Output format from diamond
    outfmt_short = False

    in_file = None
    itype = None
    translate = None
    query_gencode = None

    # Results
    queries = hits = no_hits = None

    ##
    def __init__(self, args):
        
        self.itype = args.itype
        self.translate = args.translate
        self.query_gencode = args.trans_table

        self.dmnd_db = args.dmnd_db if args.dmnd_db else get_eggnog_dmnd_db()

        self.cpu = args.cpu

        self.sensmode = args.sensmode
        
        self.query_cov = args.query_cover
        self.subject_cov = args.subject_cover

        self.matrix = args.matrix
        self.gapopen = args.gapopen
        self.gapextend = args.gapextend
        self.block_size = args.dmnd_block_size
        self.index_chunks = args.dmnd_index_chunks

        self.pident_thr = args.pident
        self.evalue_thr = args.dmnd_evalue
        self.score_thr = args.dmnd_score
        # self.excluded_taxa = args.excluded_taxa if args.excluded_taxa else None

        self.outfmt_short = args.outfmt_short
        
        self.temp_dir = args.temp_dir
        self.no_file_comments = args.no_file_comments
        
        return

    ##
    def search(self, in_file, seed_orthologs_file, hits_file = None):
        # DiamondSearcher does not use the "hits_file"
        # but we need this wrapper to define the search interface for Emapper
        return self._search(in_file, seed_orthologs_file)

    def _search(self, in_file, seed_orthologs_file):
        if not DIAMOND:
            raise EmapperException("%s command not found in path" % (DIAMOND))

        self.in_file = in_file
        
        tempdir = mkdtemp(prefix='emappertmp_dmdn_', dir=self.temp_dir)
        try:
            output_file = pjoin(tempdir, uuid.uuid4().hex)
            cmd = self.run_diamond(in_file, tempdir, output_file)
            self.hits = self.parse_diamond(output_file)            
            self.output_diamond(cmd, self.hits, seed_orthologs_file)

        except Exception as e:
            raise e
        finally:
            shutil.rmtree(tempdir)
            
        return

    ##
    def get_hits(self):
        return self.hits

    ##
    def get_no_hits(self):
        if self.hits is not None:
            hit_queries = set([x[0] for x in self.hits])
            self.queries = set({name for name, seq in iter_fasta_seqs(self.in_file)})
            self.no_hits = set(self.queries).difference(hit_queries)
            
        return self.no_hits

    ##
    def run_diamond(self, fasta_file, tempdir, output_file):
        ##
        # search type
        if self.itype == ITYPE_CDS and self.translate == True:
            tool = 'blastp'
            handle, query_file = mkstemp(dir = tempdir, text = True)
            translate_cds_to_prots(fasta_file, query_file)
        elif self.itype == ITYPE_CDS or self.itype == ITYPE_GENOME or self.itype == ITYPE_META:
            tool = 'blastx'
            query_file = fasta_file
        elif self.itype == ITYPE_PROTS:
            tool = 'blastp'
            query_file = fasta_file
        else:
            raise EmapperException(f"Unrecognized --itype {self.itype}.")

        ##
        #prepare command
        cmd = (
            f'{DIAMOND} {tool} -d {self.dmnd_db} -q {query_file} '
            f'--threads {self.cpu} -o {output_file} '
        )
        
        if self.sensmode != SENSMODE_FAST: cmd += f' --{self.sensmode}'

        if self.evalue_thr is not None: cmd += f' -e {self.evalue_thr}'
        if self.score_thr is not None: cmd += f' --min-score {self.score_thr}'
        if self.pident_thr is not None: cmd += f' --id {self.pident_thr}'
        if self.query_cov is not None: cmd += f' --query-cover {self.query_cov}'
        if self.subject_cov is not None: cmd += f' --subject-cover {self.subject_cov}'

        if self.query_gencode: cmd += f' --query-gencode {self.query_gencode}'
        if self.matrix: cmd += f' --matrix {self.matrix}'
        if self.gapopen: cmd += f' --gapopen {self.gapopen}'
        if self.gapextend: cmd += f' --gapextend {self.gapextend}'
        if self.block_size: cmd += f' --block-size {self.block_size}'
        if self.index_chunks: cmd += f' -c {self.index_chunks}'

        if self.itype == ITYPE_CDS or self.itype == ITYPE_PROTS:
            cmd += " --top 3 "
        else: # self.itype == ITYPE_GENOME or self.itype == ITYPE_META: i.e. gene prediction
            cmd += " --max-target-seqs 0 --max-hsps 0 "

        ##
        # output format
        OUTFMT_SHORT = " --outfmt 6 qseqid sseqid evalue bitscore"
        OUTFMT_LONG = " --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp"
        if self.itype == ITYPE_GENOME or self.itype == ITYPE_META: # i.e. gene prediction
            cmd += OUTFMT_LONG
        else:
            if self.outfmt_short == True:
                cmd += OUTFMT_SHORT
            else:
                cmd += OUTFMT_LONG

        # NOTE about short output format:
        # diamond should run faster if no pident, qcov, scov values are used either as filter or to be output
        # This is because it needs to compute them, whereas using only evalue and score does not need to recompute.
        # Therefore, the fastest way to obtain diamond alignments is using the OUTFMT_SHORT format and
        # not using --id, --query-cover, --subject-cover thresholds. Of course, does not always fit our needs.

        ##
        # run command
        print(colorify('  '+cmd, 'yellow'))
        try:
            completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
        except subprocess.CalledProcessError as cpe:
            raise EmapperException("Error running diamond: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])
        
        return cmd


    ##
    def parse_diamond(self, raw_dmnd_file):
        if self.itype == ITYPE_CDS or self.itype == ITYPE_PROTS:
            return self._parse_diamond(raw_dmnd_file)
        else: #self.itype == ITYPE_GENOME or self.itype == ITYPE_META:
            return self._parse_genepred(raw_dmnd_file)
        
    ##
    def _parse_diamond(self, raw_dmnd_file):
        hits = []

        visited = set()
        with open(raw_dmnd_file, 'r') as raw_f:
            for line in raw_f:
                if not line.strip() or line.startswith('#'):
                    continue

                fields = list(map(str.strip, line.split('\t')))
                # fields are defined in run_diamond
                # OUTFMT_SHORT = " --outfmt 6 qseqid sseqid evalue bitscore"
                # OUTFMT_LONG = " --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp"
                
                query = fields[0]

                # only one result per query
                if query in visited:
                    continue
                visited.add(query)
                
                hit = fields[1]

                if self.outfmt_short == True:
                    evalue = float(fields[2])
                    score = float(fields[3])
                    hits.append([query, hit, evalue, score])
                else:
                    pident = float(fields[2])
                    qstart = int(fields[6])
                    qend = int(fields[7])
                    sstart = int(fields[8])
                    send = int(fields[9])
                    evalue = float(fields[10])
                    score = float(fields[11])
                    qcov = float(fields[12])
                    scov = float(fields[13])
                    hits.append([query, hit, evalue, score, qstart, qend, sstart, send, qcov, scov])
            
        return hits

    ##
    def _parse_genepred(self, raw_dmnd_file):
        hits = []
        curr_query_hits = []
        prev_query = None
        
        visited = set()
        with open(raw_dmnd_file, 'r') as raw_f:
            for line in raw_f:
                if not line.strip() or line.startswith('#'):
                    continue

                fields = list(map(str.strip, line.split('\t')))
                # fields are defined in run_diamond
                # OUTFMT_LONG = " --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp"

                hit = fields[1]
                pident = float(fields[2])
                evalue = float(fields[10])
                score = float(fields[11])
                
                query = fields[0]
                qstart = int(fields[6])
                qend = int(fields[7])
                sstart = int(fields[8])
                send = int(fields[9])
                qcov = float(fields[12])
                scov = float(fields[13])
                
                hit = [query, hit, evalue, score, qstart, qend, sstart, send, qcov, scov]

                if query == prev_query:
                    if not hit_does_overlap(hit, curr_query_hits):
                        hits.append(hit)
                        curr_query_hits.append(hit)
                else:
                    hits.append(hit)
                    curr_query_hits = [hit]
                    
                prev_query = query
        
        return hits

    ##
    def output_diamond(self, cmd, hits, out_file):
        if self.itype == ITYPE_CDS or self.itype == ITYPE_PROTS:
            return self._output_diamond(cmd, hits, out_file)
        else: #self.itype == ITYPE_GENOME or self.itype == ITYPE_META:
            return self._output_genepred(cmd, hits, out_file)
        
    ##
    def _output_diamond(self, cmd, hits, out_file):
        with open(out_file, 'w') as OUT:

            # comments
            if not self.no_file_comments:
                print(get_call_info(), file=OUT)
                print('#'+cmd, file=OUT)

            # header
            if self.outfmt_short == True:
                print('#'+"\t".join("qseqid sseqid evalue bitscore".split(" ")), file=OUT)
            else:
                print('#'+"\t".join("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp".split(" ")), file=OUT)

            # rows
            for line in hits:
                print('\t'.join(map(str, line)), file=OUT)
                
        return

    ##
    def _output_genepred(self, cmd, hits, out_file):
        queries_suffixes = {}
        with open(out_file, 'w') as OUT:

            # comments
            if not self.no_file_comments:
                print(get_call_info(), file=OUT)
                print('#'+cmd, file=OUT)

            # header
            print('#'+"\t".join("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp".split(" ")), file=OUT)

            # rows
            for line in hits:
                query = line[0]
                if query in queries_suffixes:
                    queries_suffixes[query] += 1
                    suffix = queries_suffixes[query]
                else:
                    suffix = 0
                    queries_suffixes[query] = suffix
                    
                print('\t'.join(map(str, [f"{query}_{suffix}"] + line[1:])), file=OUT)
                
        return

def hit_does_overlap(hit, hits):
    does_overlap = False

    hitstart = hit[4]
    hitend = hit[5]
    if hitstart > hitend:
        hitend = hit[4]
        hitstart = hit[5]

    for o in hits:
        ostart = o[4]
        oend = o[5]
        if ostart > oend:
            oend = o[4]
            ostart = o[5]

        overlap = get_overlap(hitstart, hitend, ostart, oend)

        if overlap is not None and overlap > 0:
            does_overlap = True
            break

    return does_overlap


def get_overlap(hitstart, hitend, ostart, oend, allow_diff_frame = False):
    overlap = None

    # if different frame and not allow different frame to compute overlap
    # return overlap None
    # If allow different frame is True, overlap will be computed
    if abs(hitstart - ostart) % 3 != 0 and allow_diff_frame == False:
        overlap = None
    else:        
        # no overlap
        if hitend <= ostart:
            overlap = hitend - ostart

        # no overlap
        elif hitstart >= oend:
            overlap = oend - hitstart

        # envelopes
        elif (hitstart >= ostart and hitend <= oend) or (ostart >= hitstart and oend <= hitend):
            overlap_start = max(hitstart, ostart)
            overlap_end = min(hitend, oend)
            overlap = overlap_end - (overlap_start - 1)

        # overlap, no envelope
        else:
            hittol = (hitend - (hitstart - 1)) * OVERLAP_TOL_FRACTION
            otol = (oend - (ostart - 1)) * OVERLAP_TOL_FRACTION
            # the tolerance to apply to each end
            # depends on which sequence overhangs on that specific end
            if hitstart < ostart:
                tol1 = hittol
                tol2 = otol
            else:
                tol1 = otol
                tol2 = hittol

            hang_left = abs(hitstart - ostart)
            hang_right = abs(hitend - oend)

            if hang_left > tol1 and hang_right > tol2:
                overlap = -1 # consider as no overlapping
            else:
                overlap_start = max(hitstart, ostart)
                overlap_end = min(hitend, oend)
                overlap = overlap_end - (overlap_start - 1)
            
    return overlap

## END
