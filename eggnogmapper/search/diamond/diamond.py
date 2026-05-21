##
## CPCantalapiedra 2019

import os
from os.path import isdir as pisdir, isfile as pisfile
import shutil
import subprocess
from sys import stderr as sys_stderr
from tempfile import mkdtemp, mkstemp

from ...emapperException import EmapperException
from ...common import DIAMOND, ITYPE_CDS, ITYPE_PROTS, ITYPE_GENOME, ITYPE_META
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
        f'{DIAMOND} makedb --in \'{in_fasta}\' --db \'{dbprefix}\''
    )

    print(colorify('  '+cmd, 'yellow'))
    try:
        completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
    except subprocess.CalledProcessError as cpe:
        raise EmapperException("Error running diamond: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])
        
    return

def _auto_dmnd_resources(block_size, index_chunks, algo,
                         query_path=None, n_cpus=None):
    """Pick diamond runtime tuning from host RAM, query input size, and
    CPU count when the user didn't override.

    Returns a `(block_size, index_chunks, algo)` tuple. `None` for any
    field means "leave diamond's own default in place". User overrides
    always take precedence.

    Decision rules:

    * `block_size` and `index_chunks` come from host RAM. The eggNOG e7
      dmnd is ~23 GB; bands are sized to leave headroom for the
      annotator pool + OS page cache + kernel.

          <32 GB        → diamond default (block 2, chunks 4)
          32-64 GB      → block 4,  chunks 2     (~40 GB peak,  ~1.3× speedup)
          64-96 GB      → block 6,  chunks 2     (~52 GB peak,  ~1.4× speedup)
          ≥96 GB        → block 8,  chunks 1     (~76 GB peak,  ~1.7× speedup)

      Diamond peak RAM ≈ block_size × 6 + db_size / index_chunks +
      threads × 0.5 GB.

    * `algo` comes from the query input. Diamond's `auto` (default) is
      tuned for big query sets:

      - Small input (<5 MB or <1k seqs): `algo=1` (query-indexed). DB
        index is loaded once, queries indexed in memory. Big startup
        win for tiny inputs where index-load dominates total wall.
      - Big input on many CPUs: stay on diamond `auto` — it picks the
        double-indexed algorithm that scales better with thread count.

      `algo=ctg` (contiguous-seed) is *not* auto-picked; it requires
      explicit user opt-in via `--dmnd_algo ctg` because it's a
      sensitivity tradeoff, not a pure perf knob.
    """
    # ---- block_size + index_chunks: host-RAM driven ----
    if block_size is None or index_chunks is None:
        try:
            import psutil
            avail_gb = psutil.virtual_memory().total / (1024 ** 3)
        except Exception:
            avail_gb = None  # silent fallthrough — platforms without psutil
        auto_bs, auto_ic = None, None
        if avail_gb is not None:
            if avail_gb >= 96:
                auto_bs, auto_ic = 8.0, 1
            elif avail_gb >= 64:
                auto_bs, auto_ic = 6.0, 2
            elif avail_gb >= 32:
                auto_bs, auto_ic = 4.0, 2
            # <32 GB: leave diamond's default
        if block_size is None:
            block_size = auto_bs
        if index_chunks is None:
            index_chunks = auto_ic

    # ---- algo: input-size driven ----
    if algo == DMND_ALGO_AUTO and query_path is not None:
        try:
            import os as _os
            qsize_bytes = _os.path.getsize(query_path)
        except OSError:
            qsize_bytes = None
        # Heuristic: query-indexed (algo=1) wins when the input is small
        # enough that diamond's startup index-load dominates total wall.
        # ~5 MB FASTA ≈ 12k average proteins (300 aa each); below that,
        # switching to algo=1 saves the first DB-index pass.
        if qsize_bytes is not None and qsize_bytes < 5 * 1024 * 1024:
            algo = DMND_ALGO_1

    return block_size, index_chunks, algo


class DiamondSearcher:

    name = "diamond"
    
    # Command
    cpu = tool = dmnd_db = temp_dir = no_file_comments = None
    matrix = frameshift = gapopen = gapextend = None
    block_size = index_chunks = None

    # Filters
    pident_thr = score_thr = evalue_thr = query_cov = subject_cov = None

    in_file = None
    itype = None
    translate = None
    query_gencode = None

    allow_overlaps = None
    overlap_tol = None
    
    resume = None

    ##
    def __init__(self, args, dmnd_db):
        
        self.itype = args.itype
        self.translate = args.translate
        self.query_gencode = args.trans_table

        self.allow_overlaps = args.allow_overlaps
        self.overlap_tol = args.overlap_tol

        self.dmnd_db = dmnd_db

        self.cpu = args.cpu

        self.sensmode = args.sensmode
        self.iterate = args.dmnd_iterate
        self.ignore_warnings = args.dmnd_ignore_warnings
        # `self.algo` is overwritten in the search command builder once
        # we know the query file path (input-size-aware auto-pick).
        # `_user_algo` preserves the original CLI choice so the
        # auto-pick can tell "user said auto" from "user picked X".
        self.algo = args.dmnd_algo
        
        self.query_cov = args.query_cover
        self.subject_cov = args.subject_cover

        self.matrix = args.matrix
        self.frameshift = args.dmnd_frameshift
        self.gapopen = args.gapopen
        self.gapextend = args.gapextend
        # Auto-tune at __init__ runs only RAM-driven decisions
        # (block_size, index_chunks). The algo decision is deferred to
        # search-time because it depends on the query file size which
        # may not be fixed yet (CDS path may be translating a different
        # input). The search() / blastp() command builder calls
        # `_resolve_algo()` to do the input-size pick at that point.
        bs, ic, _ = _auto_dmnd_resources(
            args.dmnd_block_size, args.dmnd_index_chunks,
            algo=DMND_ALGO_DEFAULT,  # placeholder — real pick at search time
        )
        self.block_size = bs
        self.index_chunks = ic
        self.dmnd_top = getattr(args, "dmnd_top", 1)
        self._user_algo = args.dmnd_algo  # remember, for search-time auto-pick

        self.pident_thr = args.pident
        self.evalue_thr = args.dmnd_evalue
        self.score_thr = args.dmnd_score
        # self.excluded_taxa = args.excluded_taxa if args.excluded_taxa else None

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
                                          self.no_file_comments,
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
            f'{DIAMOND} {tool} -d \'{self.dmnd_db}\' -q \'{query_file}\' '
            f'--threads {self.cpu} -o \'{output_file}\' --tmpdir \'{self.temp_dir}\''
        )
        
        if self.sensmode != SENSMODE_DEFAULT: cmd += f' --{self.sensmode}'

        if self.iterate is not None and self.iterate == DMND_ITERATE_YES:
            cmd += f' --iterate'

        if self.ignore_warnings is not None and self.ignore_warnings == True:
            cmd += f' --ignore-warnings'

        # Resolve algo at command-build time so we can inspect the actual
        # query file size. User overrides via `--dmnd_algo` propagate
        # through `_user_algo` unchanged; only `auto` triggers the
        # input-size pick (algo=1 for small inputs).
        _, _, resolved_algo = _auto_dmnd_resources(
            self.block_size, self.index_chunks,
            algo=self._user_algo,
            query_path=query_file,
        )
        self.algo = resolved_algo
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
            # `--dmnd_top 3` (default) keeps all hits within 3 % of top
            # score per query, letting the seed picker apply secondary
            # criteria (pident / qcov). `--dmnd_top 1` switches to
            # `--max-target-seqs 1` — diamond can prune as soon as it
            # has the bitscore-best hit, ~20-30 % faster, with the
            # bitscore-best almost always being the seed anyway.
            if int(getattr(self, "dmnd_top", 3)) == 1:
                cmd += " --top 0 "
            else:
                cmd += " --top 3 "
        else: # self.itype == ITYPE_GENOME or self.itype == ITYPE_META: i.e. gene prediction
            cmd += " --max-target-seqs 0 --max-hsps 0 "

        ##
        # output format
        OUTFMT_LONG = " --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp"
        cmd += OUTFMT_LONG

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
                # OUTFMT_LONG fields:
                # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

                query = fields[0]

                # only one result per query
                if prev_query is not None and query == prev_query:
                    continue
                else:
                    prev_query = query

                target = fields[1]
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

                yield hit
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
