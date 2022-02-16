##
## CPCantalapiedra 2019

import os
from os.path import isdir as pisdir, isfile as pisfile
import shutil
import subprocess
from sys import stderr as sys_stderr
from tempfile import mkdtemp, mkstemp

from ...emapperException import EmapperException
from ...common import DIAMOND, get_eggnog_dmnd_db, ITYPE_CDS, ITYPE_PROTS, ITYPE_GENOME, ITYPE_META
from ...utils import colorify, translate_cds_to_prots

from ..hmmer.hmmer_seqio import iter_fasta_seqs

from ..hits_io import output_seeds

SENSMODE_FAST = "fast"
SENSMODE_DEFAULT = "default"
SENSMODE_MID_SENSITIVE = "mid-sensitive"
SENSMODE_SENSITIVE = "sensitive"
SENSMODE_MORE_SENSITIVE = "more-sensitive"
SENSMODE_VERY_SENSITIVE = "very-sensitive"
SENSMODE_ULTRA_SENSITIVE = "ultra-sensitive"
# sens modes in diamond 0.9.24
# SENSMODES = [SENSMODE_FAST, SENSMODE_SENSITIVE, SENSMODE_MORE_SENSITIVE]
# sens modes in diamond 2.0.4
SENSMODES = [SENSMODE_DEFAULT, SENSMODE_FAST, SENSMODE_MID_SENSITIVE, SENSMODE_SENSITIVE, SENSMODE_MORE_SENSITIVE, SENSMODE_VERY_SENSITIVE, SENSMODE_ULTRA_SENSITIVE]

# Diamond --iterate flag will be controlled with --dmnd_iterate, with next options
DMND_ITERATE_YES = "yes"
DMND_ITERATE_NO = "no"
DMND_ITERATE_DEFAULT = DMND_ITERATE_YES

# Diamond --algo option, which can be changed to increase performance when searching small query sets
DMND_ALGO_AUTO = "auto"
DMND_ALGO_0 = "0"
DMND_ALGO_1 = "1"
DMND_ALGO_CTG = "ctg"
DMND_ALGO_DEFAULT = DMND_ALGO_AUTO
        
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
    matrix = frameshift = gapopen = gapextend = None
    block_size = index_chunks = None

    # Filters
    pident_thr = score_thr = evalue_thr = query_cov = subject_cov = None

    # Output format from diamond
    outfmt_short = False

    in_file = None
    itype = None
    translate = None
    query_gencode = None

    allow_overlaps = None
    overlap_tol = None
    
    resume = None

    ##
    def __init__(self, args):
        
        self.itype = args.itype
        self.translate = args.translate
        self.query_gencode = args.trans_table

        self.allow_overlaps = args.allow_overlaps
        self.overlap_tol = args.overlap_tol

        self.dmnd_db = args.dmnd_db if args.dmnd_db else get_eggnog_dmnd_db()

        self.cpu = args.cpu

        self.sensmode = args.sensmode
        self.iterate = args.dmnd_iterate
        self.ignore_warnings = args.dmnd_ignore_warnings
        self.algo = args.dmnd_algo
        
        self.query_cov = args.query_cover
        self.subject_cov = args.subject_cover

        self.matrix = args.matrix
        self.frameshift = args.dmnd_frameshift
        self.gapopen = args.gapopen
        self.gapextend = args.gapextend
        self.block_size = args.dmnd_block_size
        self.index_chunks = args.dmnd_index_chunks

        self.pident_thr = args.pident
        self.evalue_thr = args.dmnd_evalue
        self.score_thr = args.dmnd_score
        # self.excluded_taxa = args.excluded_taxa if args.excluded_taxa else None

        self.outfmt_short = args.outfmt_short
        
        self.temp_dir = mkdtemp(prefix='emappertmp_dmdn_', dir=args.temp_dir)
        self.no_file_comments = args.no_file_comments

        self.resume = args.resume

        self.gff_ID_field = args.decorate_gff_ID_field
        
        return

    ##
    def clear(self):
        if self.temp_dir is not None and pisdir(self.temp_dir):
            try:
                shutil.rmtree(self.temp_dir)
            except OSError as err:
                print(f"Warning: OS error while removing {self.temp_dir}", file = sys_stderr)
                print(f"OS error: {err}", file = sys_stderr)
        return
    
    ##
    def search(self, in_file, seed_orthologs_file, hits_file):
        hits_generator = None
        
        if not DIAMOND:
            raise EmapperException("%s command not found in path" % (DIAMOND))

        self.in_file = in_file
        
        try:
            cmds = None
            
            # 1) either resume from previous hits or run diamond to generate the hits
            if self.resume == True:
                if pisfile(hits_file):
                    pass
                else:
                    raise EmapperException(f"Couldn't find hits file {hits_file} to resume.")
            else:
                cmds = self.run_diamond(in_file, hits_file)

            # 2) parse search hits to seeds orthologs
            if self.itype == ITYPE_CDS or self.itype == ITYPE_PROTS:
                hits_generator = self._parse_diamond(hits_file)
                
            else: #self.itype == ITYPE_GENOME or self.itype == ITYPE_META:
                # parse_genepred (without coordinate change)
                hits_generator = self._parse_genepred(hits_file)

                
            # 3) output seeds
            if self.itype == ITYPE_CDS or self.itype == ITYPE_PROTS:
                change_seeds_coords = False
            else: #self.itype == ITYPE_GENOME or self.itype == ITYPE_META:
                # change seeds coordinates relative to the ORF, not to the contig (to use them for the .seed_orthologs file)
                change_seeds_coords = True
                
            hits_generator = output_seeds(cmds, hits_generator,
                                          seed_orthologs_file, 
                                          self.no_file_comments, self.outfmt_short,
                                          change_seeds_coords)

        except Exception as e:
            raise e
            
        return hits_generator

    ##
    def run_diamond(self, fasta_file, output_file):
        cmds = []

        handle = None
        ##
        # search type
        if self.itype == ITYPE_CDS and self.translate == True:
            tool = 'blastp'
            handle, query_file = mkstemp(dir = self.temp_dir, text = True)
            translate_cds_to_prots(fasta_file, query_file, self.query_gencode)
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
            f'--threads {self.cpu} -o {output_file} --tmpdir {self.temp_dir}'
        )
        
        if self.sensmode != SENSMODE_DEFAULT: cmd += f' --{self.sensmode}'

        if self.iterate is not None and self.iterate == DMND_ITERATE_YES:
            cmd += f' --iterate'

        if self.ignore_warnings is not None and self.ignore_warnings == True:
            cmd += f' --ignore-warnings'

        if self.algo is not None and self.algo != DMND_ALGO_AUTO:
            cmd += f' --algo {self.algo}'

        if self.evalue_thr is not None: cmd += f' -e {self.evalue_thr}'
        if self.score_thr is not None: cmd += f' --min-score {self.score_thr}'
        if self.pident_thr is not None: cmd += f' --id {self.pident_thr}'
        if self.query_cov is not None: cmd += f' --query-cover {self.query_cov}'
        if self.subject_cov is not None: cmd += f' --subject-cover {self.subject_cov}'

        if self.query_gencode: cmd += f' --query-gencode {self.query_gencode}'
        if self.matrix: cmd += f' --matrix {self.matrix}'
        if self.frameshift is not None: cmd += f' --frameshift {self.frameshift}'
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
            cmds.append(cmd)
        except subprocess.CalledProcessError as cpe:
            raise EmapperException("Error running diamond: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])
        finally:
            if handle is not None:
                os.close(handle)
        
        return cmds
        
    ##
    def _parse_diamond(self, raw_dmnd_file):        

        prev_query = None
        # parse hits
        with open(raw_dmnd_file, 'r') as raw_f:
            for line in raw_f:
                if not line.strip() or line.startswith('#'):
                    continue

                fields = list(map(str.strip, line.split('\t')))
                # fields are defined in run_diamond
                # OUTFMT_SHORT = " --outfmt 6 qseqid sseqid evalue bitscore"
                # OUTFMT_LONG = " --outfmt 6 qseqid sseqid pident length mismatch
                # gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp"
                
                query = fields[0]

                # only one result per query
                if prev_query is not None and query == prev_query:
                    continue
                else:
                    prev_query = query
                
                target = fields[1]

                if self.outfmt_short == True:
                    evalue = float(fields[2])
                    score = float(fields[3])
                    hit = [query, target, evalue, score]
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
                    hit = [query, target, evalue, score, qstart, qend, sstart, send, pident, qcov, scov]

                yield hit # hit
        return

    ##
    def _parse_genepred(self, raw_dmnd_file):
        
        curr_query_hits = []
        prev_query = None
        queries_suffixes = {}
        
        with open(raw_dmnd_file, 'r') as raw_f:
            for line in raw_f:
                if not line.strip() or line.startswith('#'):
                    continue

                fields = list(map(str.strip, line.split('\t')))
                # fields are defined in run_diamond
                # OUTFMT_LONG = " --outfmt 6 qseqid sseqid pident length mismatch
                # gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp"
                
                query = fields[0]
                target = fields[1]
                pident = float(fields[2])
                evalue = float(fields[10])
                score = float(fields[11])
                qstart = int(fields[6])
                qend = int(fields[7])
                sstart = int(fields[8])
                send = int(fields[9])
                qcov = float(fields[12])
                scov = float(fields[13])
                
                hit = [query, target, evalue, score, qstart, qend, sstart, send, pident, qcov, scov]

                if query == prev_query:
                    if self.allow_overlaps == ALLOW_OVERLAPS_ALL:
                        yield [f"{hit[0]}_{suffix}"]+hit[1:] # hit
                        
                    else:
                        if not hit_does_overlap(hit, curr_query_hits, self.allow_overlaps, self.overlap_tol):
                            if query in queries_suffixes:
                                queries_suffixes[query] += 1
                                suffix = queries_suffixes[query]
                            else:
                                suffix = 0
                                queries_suffixes[query] = suffix

                            yield [f"{hit[0]}_{suffix}"]+hit[1:] # hit
                            curr_query_hits.append(hit)
                        
                else:
                    if query in queries_suffixes:
                        queries_suffixes[query] += 1
                        suffix = queries_suffixes[query]
                    else:
                        suffix = 0
                        queries_suffixes[query] = suffix

                    yield [f"{hit[0]}_{suffix}"]+hit[1:] # hit
                    curr_query_hits = [hit]
                    
                prev_query = query
        return

#
ALLOW_OVERLAPS_NONE = "none"
ALLOW_OVERLAPS_OPPOSITE_STRAND = "strand"
ALLOW_OVERLAPS_DIFF_FRAME = "diff_frame"
ALLOW_OVERLAPS_ALL = "all"

def hit_does_overlap(hit, hits, allow_overlaps, overlap_tol):
    does_overlap = False

    hitstart = hit[4]
    hitend = hit[5]
    hit_strand = "+"
    if hitstart > hitend:
        hitend = hit[4]
        hitstart = hit[5]
        hit_strand = "-"

    for o in hits:
        ostart = o[4]
        oend = o[5]
        o_strand = "+"
        if ostart > oend:
            oend = o[4]
            ostart = o[5]
            o_strand = "-"

        same_strand = (hit_strand == o_strand)

        if allow_overlaps == ALLOW_OVERLAPS_OPPOSITE_STRAND and not same_strand:
            continue
        
        same_frame = (abs(hitstart - ostart) % 3 == 0)
        
        if allow_overlaps == ALLOW_OVERLAPS_DIFF_FRAME and (not same_strand or not same_frame):
            continue

        overlap = get_overlap(hitstart, hitend, ostart, oend, overlap_tol)

        if overlap is not None and overlap > 0:
            does_overlap = True
            break

    return does_overlap


def get_overlap(hitstart, hitend, ostart, oend, overlap_tol):
    overlap = None
    
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
        hittol = (hitend - (hitstart - 1)) * overlap_tol
        otol = (oend - (ostart - 1)) * overlap_tol

        overlap_start = max(hitstart, ostart)
        overlap_end = min(hitend, oend)
        overlap = overlap_end - (overlap_start - 1)

        if overlap < min(hittol, otol):
            overlap = -1
            
    return overlap

## END
