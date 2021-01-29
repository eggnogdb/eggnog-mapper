##
## CPCantalapiedra 2020

from os.path import join as pjoin
import shutil
import subprocess
from tempfile import mkdtemp

from ..emapperException import EmapperException
from ..common import PRODIGAL, ITYPE_GENOME, ITYPE_META
from ..utils import colorify

# This class handles prediction of genes
# using Prodigal v2
class ProdigalPredictor:

    temp_dir = None
    pmode = None
    cpu = None

    outdir = None # dir with prodigal out files
    outgff = outprots = outcds = outorfs = None # prodigal out files
    
    
    def __init__(self, args):

        if args.itype == ITYPE_GENOME:
            self.pmode = "single" # or self.pmode = ""
        elif args.itype == ITYPE_META:
            self.pmode = "meta"
        else:
            raise EmapperException(f"Unsupported input type {args.itype} for ProdigalPredictor")
        self.cpu = args.cpu
        
        self.temp_dir = args.temp_dir
        
        return

    def predict(self, in_file):
        if not PRODIGAL:
            raise EmapperException("%s command not found in path" % (PRODIGAL))

        self.outdir = mkdtemp(prefix='emappertmp_prod_', dir=self.temp_dir)
        try:
            cmd = self.run_prodigal(in_file, self.outdir)

        except Exception as e:
            raise e
        # finally:
        #     shutil.rmtree(tempdir)
        return

    def clear(self):
        shutil.rmtree(self.outdir)
        return

    def run_prodigal(self, in_file, outdir):
        self.outfile = pjoin(outdir, "output.gff")
        self.outprots = pjoin(outdir, "output.faa")
        self.outcds = pjoin(outdir, "output.fna")
        self.outorfs = pjoin(outdir, "output.orfs")
        cmd = (
            f'{PRODIGAL} -i {in_file} -p {self.pmode} '
            f'-o {self.outfile} -f gff '
            f'-a {self.outprots} -d {self.outcds} '
            f'-s {self.outorfs}'
        )

        print(colorify('  '+cmd, 'yellow'))
        try:
            completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
        except subprocess.CalledProcessError as cpe:
            raise EmapperException("Error running prodigal: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])

        return cmd

## END
