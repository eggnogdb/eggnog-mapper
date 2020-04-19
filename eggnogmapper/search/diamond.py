##
## CPCantalapiedra 2019

from os.path import join as pjoin
import shutil
import subprocess
from tempfile import mkdtemp
import uuid

from ..common import DIAMOND, get_eggnog_dmnd_db, get_call_info
from ..utils import colorify

class DiamondSearcher:

    # Command
    cpu = tool = dmnd_db = dmnd_opts = temp_dir = no_file_comments = None

    # Filters
    score_thr = evalue_thr = query_cov = subject_cov = excluded_taxa = None

    ##
    def __init__(self, args):

        self.cpu = args.cpu
        
        if args.translate:
            self.tool = 'blastx'
        else:
            self.tool = 'blastp'

        if args.dmnd_db:
            self.dmnd_db = args.dmnd_db
        else:
            self.dmnd_db = get_eggnog_dmnd_db()
        
        self.dmnd_opts = ''
        if args.matrix is not None:
            self.dmnd_opts += ' --matrix %s' % args.matrix
        if args.gapopen is not None:
            self.dmnd_opts += ' --gapopen %d' % args.gapopen
        if args.gapextend is not None:
            self.dmnd_opts += ' --gapextend %d' % args.gapextend

        self.temp_dir = args.temp_dir
        self.no_file_comments = args.no_file_comments
        
        self.score_thr = args.seed_ortholog_score
        self.evalue_thr = args.seed_ortholog_evalue
        self.query_cov = args.query_cover
        self.subject_cov = args.subject_cover
        self.excluded_taxa = args.excluded_taxa if args.excluded_taxa else None
        
        return

    ##
    def search(self, in_file, seed_orthologs_file):
        if not DIAMOND:
            raise EmapperException("%s command not found in path" % (DIAMOND))

        tempdir = mkdtemp(prefix='emappertmp_dmdn_', dir=self.temp_dir)
        try:
            output_file = pjoin(tempdir, uuid.uuid4().hex)
            cmd = self.run_diamond(in_file, output_file)
            parsed = self.parse_diamond(output_file)
            
            self.output_diamond(cmd, parsed, seed_orthologs_file)

        except subprocess.CalledProcessError as e:
            raise e
        finally:
            shutil.rmtree(tempdir)
            
        return

    ##
    def run_diamond(self, fasta_file, output_file):
               
        cmd = (
            f'{DIAMOND} {self.tool} -d {self.dmnd_db} -q {fasta_file} '
            f'--more-sensitive --threads {self.cpu} -e {self.evalue_thr} -o {output_file} '
            f'--query-cover {self.query_cov} --subject-cover {self.subject_cov} {self.dmnd_opts}'
        )
        
        if self.excluded_taxa:
            cmd += " --max-target-seqs 25 "
        else:
            cmd += " --top 3 "

        print(colorify('  '+cmd, 'yellow'))

        with open(output_file+'.stdout', 'w') as STDOUT:
            subprocess.check_call(cmd, shell=True, stdout=STDOUT)

        return cmd

    ##
    def parse_diamond(self, raw_dmnd_file):
        parsed = []

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

                parsed.append([query, hit, evalue, score])
            
        return parsed

    ##
    def output_diamond(self, cmd, parsed, out_file):
        with open(out_file, 'w') as OUT:
        
            if not self.no_file_comments:
                print(get_call_info(), file=OUT)
                print('#'+cmd, file=OUT)

            for line in parsed:
                print('\t'.join(map(str, line)), file=OUT)
                
        return

## END
