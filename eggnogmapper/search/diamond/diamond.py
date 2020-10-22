##
## CPCantalapiedra 2019

from os.path import join as pjoin
import shutil
import subprocess
from tempfile import mkdtemp
import uuid

from ...emapperException import EmapperException
from ...common import DIAMOND, get_eggnog_dmnd_db, get_call_info
from ...utils import colorify

from ..hmmer.hmmer_seqio import iter_fasta_seqs

SENSMODE_FAST = "fast"
SENSMODE_MID_SENSITIVE = "mid-sensitive"
SENSMODE_SENSITIVE = "sensitive"
SENSMODE_MORE_SENSITIVE = "more-sensitive"
SENSMODE_VERY_SENSITIVE = "very-sensitive"
SENSMODE_ULTRA_SENSITIVE = "ultra-sensitive"
SENSMODES = [SENSMODE_FAST, SENSMODE_SENSITIVE, SENSMODE_MORE_SENSITIVE]
# SENSMODES = [SENSMODE_FAST, SENSMODE_MID_SENSITIVE, SENSMODE_SENSITIVE, SENSMODE_MORE_SENSITIVE, SENSMODE_VERY_SENSITIVE, SENSMODE_ULTRA_SENSITIVE]

class DiamondSearcher:

    # Command
    cpu = tool = dmnd_db = temp_dir = no_file_comments = None
    matrix = gapopen = gapextend = None

    # Filters
    score_thr = evalue_thr = query_cov = subject_cov = excluded_taxa = None

    in_file = None

    # Results
    queries = hits = no_hits = None

    ##
    def __init__(self, args):
        
        if args.translate:
            self.tool = 'blastx'
        else:
            self.tool = 'blastp'

        self.dmnd_db = args.dmnd_db if args.dmnd_db else get_eggnog_dmnd_db()

        self.cpu = args.cpu

        self.sensmode = args.sensmode
        
        self.query_cov = args.query_cover
        self.subject_cov = args.subject_cover

        self.matrix = args.matrix
        self.gapopen = args.gapopen
        self.gapextend = args.gapextend

        self.evalue_thr = args.dmnd_evalue
        self.score_thr = args.dmnd_score
        self.excluded_taxa = args.excluded_taxa if args.excluded_taxa else None
        
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
            cmd = self.run_diamond(in_file, output_file)
            self.hits = self.parse_diamond(output_file)            
            self.output_diamond(cmd, self.hits, seed_orthologs_file)

        except Exception as e:
            raise e
        finally:
            shutil.rmtree(tempdir)
            
        return

    ##
    def get_hits(self):
        if self.hits is not None:
            hit_queries = set([x[0] for x in self.hits])
            # no need to translate, we only need seq identifiers
            self.queries = set({name for name, seq in iter_fasta_seqs(in_file)})
            # sequences = {name: seq for name, seq in iter_fasta_seqs(in_file, translate=self.traslate)}
            # self.queries = set(sequences.keys())
            self.no_hits = set(self.queries).difference(hit_queries)
            
        return self.hits, self.no_hits

    ##
    def run_diamond(self, fasta_file, output_file):

        if self.sensmode == "fast":
            cmd = (
                f'{DIAMOND} {self.tool} -d {self.dmnd_db} -q {fasta_file} '
                f'--{self.sensmode} --threads {self.cpu} -e {self.evalue_thr} -o {output_file} '
                f'--query-cover {self.query_cov} --subject-cover {self.subject_cov}'
            )
        else:
            cmd = (
                f'{DIAMOND} {self.tool} -d {self.dmnd_db} -q {fasta_file} '
                f'--{self.sensmode} --threads {self.cpu} -e {self.evalue_thr} -o {output_file} '
                f'--query-cover {self.query_cov} --subject-cover {self.subject_cov}'
            )

        if self.matrix: cmd += ' --matrix {self.matrix}'
        if self.gapopen: cmd += ' --gapopen {self.gapopen}'
        if self.gapextend: cmd += ' --gapextend {self.gapextend}'
        
        if self.excluded_taxa: cmd += " --max-target-seqs 25 "
        else: cmd += " --top 3 "

        print(colorify('  '+cmd, 'yellow'))
        try:
            completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
        except subprocess.CalledProcessError as cpe:
            raise EmapperException("Error running diamond: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])

        return cmd

    ##
    def parse_diamond(self, raw_dmnd_file):
        hits = []

        visited = set()
        with open(raw_dmnd_file, 'r') as raw_f:
            for line in raw_f:
                if not line.strip() or line.startswith('#'):
                    continue

                fields = list(map(str.strip, line.split('\t')))
                query = fields[0]
                hit = fields[1]
                evalue = float(fields[10])
                score = float(fields[11])

                if query in visited:
                    continue

                if evalue > self.evalue_thr or score < self.score_thr:
                    continue

                if self.excluded_taxa and hit.startswith("%s." % self.excluded_taxa):
                    continue

                visited.add(query)

                hits.append([query, hit, evalue, score])
            
        return hits

    ##
    def output_diamond(self, cmd, hits, out_file):
        with open(out_file, 'w') as OUT:
        
            if not self.no_file_comments:
                print(get_call_info(), file=OUT)
                print('#'+cmd, file=OUT)

            for line in hits:
                print('\t'.join(map(str, line)), file=OUT)
                
        return

## END
