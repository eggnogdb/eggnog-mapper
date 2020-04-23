##
## JHCepas
## CPCantalapiedra 2020

import unittest
import os, shutil

from .common import run, check_seed_orthologs, check_annotations, check_hmm_hits

# General eggnog-mapper settings
HMM_HITS_SUFFIX = '.emapper.hmm_hits'
SEED_ORTHOLOGS_SUFFIX = '.emapper.seed_orthologs'
ANNOTATIONS_SUFFIX = '.emapper.annotations'

class Test(unittest.TestCase):

    # Tests
    def test_diamond(self):
        '''
        Tests the whole emapper command, run with -m diamond
        '''

        ##
        # Setup test
        
        in_file = "tests/fixtures/test_queries.fa"
        data_dir = "tests/fixtures"
        outdir = "tests/integration/out"
        outprefix = "test"

        # Observed and expected files
        obs_seed_orthologs = os.path.join(outdir, outprefix+SEED_ORTHOLOGS_SUFFIX)
        obs_annotations = os.path.join(outdir, outprefix+ANNOTATIONS_SUFFIX)
        
        exp_seed_orthologs = os.path.join(data_dir, 'test_output.emapper.seed_orthologs')
        exp_annotations = os.path.join(data_dir, 'test_output.emapper.annotations')

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = f'-m diamond -i {in_file} --data_dir {data_dir} --output_dir {outdir} -o {outprefix}'

        st, out, err = run(cmd)
        if st != 0:
            # print(out)
            # print(err)
            print(out.decode("utf-8"))
            print(err.decode("utf-8"))
        assert st == 0 # check exit status is ok

        ##
        # Check test
        
        # Check alignment phase: detection of seed orthologs
        check_seed_orthologs(obs_seed_orthologs, exp_seed_orthologs)
        
        # Check annotation phase
        check_annotations(obs_annotations, exp_annotations)

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return

    def test_hmmer(self):
        '''
        Tests the whole emapper command, run with -m hmmer
        '''
        # ./emapper.py -m hmmer -i tests/fixtures/test_queries.fa --data_dir tests/fixtures -d bact -o bact --output_dir tmp_borrar
        
        ##
        # Setup test
        
        in_file = "tests/fixtures/test_queries.fa"
        data_dir = "tests/fixtures"
        outdir = "tests/integration/out"
        outprefix = "bact"
        database = "bact"

        # Observed and expected files
        obs_hmm_hits = os.path.join(outdir, outprefix+HMM_HITS_SUFFIX)        
        obs_seed_orthologs = os.path.join(outdir, outprefix+SEED_ORTHOLOGS_SUFFIX)
        obs_annotations = os.path.join(outdir, outprefix+ANNOTATIONS_SUFFIX)

        exp_files_dir = "tests/fixtures/hmmer_expected_output/"
        exp_hmm_hits = os.path.join(exp_files_dir, "bact.emapper.hmm_hits")
        exp_seed_orthologs = os.path.join(exp_files_dir, "bact.emapper.seed_orthologs")
        exp_annotations = os.path.join(exp_files_dir, "bact.emapper.annotations")

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = f'-m hmmer -i {in_file} --data_dir {data_dir} -d {database} --output_dir {outdir} -o {outprefix}'

        st, out, err = run(cmd)
        if st != 0:
            # print(out)
            # print(err)
            print(out.decode("utf-8"))
            print(err.decode("utf-8"))
        assert st == 0 # check exit status is ok

        ##
        # Check test

        # Check HMM hits from alignment phase
        check_hmm_hits(obs_hmm_hits, exp_hmm_hits)
        
        # Check seed orthologs from alignment phase
        check_seed_orthologs(obs_seed_orthologs, exp_seed_orthologs)
        
        # Check annotation phase
        check_annotations(obs_annotations, exp_annotations)

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return
    
if __name__ == '__main__':
    unittest.main()

## END
