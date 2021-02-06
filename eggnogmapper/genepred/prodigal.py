##
## CPCantalapiedra 2020

from os.path import isfile
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

    trans_table = None
    training_genome = training_file = None
    
    outdir = None # dir with prodigal out files
    outgff = outprots = outcds = outorfs = None # prodigal out files

    PMODE_SINGLE = "single"
    PMODE_META = "meta"
    
    
    def __init__(self, args):

        if args.itype == ITYPE_GENOME:
            self.pmode = self.PMODE_SINGLE # or self.pmode = ""
        elif args.itype == ITYPE_META:
            self.pmode = self.PMODE_META
        else:
            raise EmapperException(f"Unsupported input type {args.itype} for ProdigalPredictor")
        self.cpu = args.cpu

        self.trans_table = args.trans_table
        self.training_genome = args.training_genome
        self.training_file = args.training_file
        
        self.temp_dir = args.temp_dir
        
        return

    def predict(self, in_file):
        if not PRODIGAL:
            raise EmapperException("%s command not found in path" % (PRODIGAL))

        self.outdir = mkdtemp(prefix='emappertmp_prod_', dir=self.temp_dir)
        try:
            # Training: run only if the training file does NOT exist
            if self.training_genome is not None and self.training_file is not None:
                if isfile(self.training_file):
                    print(colorify(f'Warning: --training_file {self.training_file} already exists. '
                                   f'Training will be skipped, and prediction will be run using the existing training file.', 'red'))                                    
                else:
                    cmd = self.run_training(self.training_genome, self.training_file, self.outdir)

            # Gene prediction
            cmd = self.run_prodigal(in_file, self.outdir)

        except Exception as e:
            raise e
        # finally:
        #     shutil.rmtree(tempdir)
        return

    def clear(self):
        shutil.rmtree(self.outdir)
        return

    def run_training(self, in_file, training_file, outdir):
        cmd = (
            f'{PRODIGAL} -i {in_file} -t {training_file}'
        )

        if self.trans_table is not None:
            cmd += f' -g {self.trans_table}'

        print(colorify('  '+cmd, 'yellow'))
        try:
            completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
        except subprocess.CalledProcessError as cpe:
            raise EmapperException("Error running prodigal: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])

        return cmd
    
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

        if self.trans_table is not None:
            if self.pmode == self.PMODE_META:
                print(colorify(f'Warning: --trans_table (-g Prodigal option) '
                               f'is ignored by Prodigal when using -p {self.PMODE_META}', 'red'))                
            cmd += f' -g {self.trans_table}'

        if self.training_file is not None and isfile(self.training_file):
            if self.pmode == self.PMODE_META:
                print(colorify(f'Warning: Ignoring --training_file, because Prodigal does not allow training for -p {self.PMODE_META} ', 'red'))                
            else:
                cmd += f' -t {self.training_file}'

        print(colorify('  '+cmd, 'yellow'))
        try:
            completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
        except subprocess.CalledProcessError as cpe:
            raise EmapperException("Error running prodigal: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])

        return cmd

## END
