##
## JHCepas
## CPCantalapiedra 2020

import unittest
import os, shutil

from .common import run, check_seed_orthologs, check_annotations, check_hmm_hits, check_hmm_query_hit

# General eggnog-mapper settings
HMM_HITS_SUFFIX = '.emapper.hmm_hits'
SEED_ORTHOLOGS_SUFFIX = '.emapper.seed_orthologs'
ANNOTATIONS_SUFFIX = '.emapper.annotations'

class Test(unittest.TestCase):

    # Tests for "hmmscan" mode
    def test_hmm_mapper(self):
        '''
        Tests hmm_mapper to search ("hmmscan" mode) on a custom db, on disk
        '''
        # ./hmm_mapper.py -i tests/fixtures/test_queries.fa -d tests/fixtures/hmmer_custom_dbs/bact.hmm -o dummy --output_dir tmp_borrar
        
        ##
        # Setup test
        
        in_file = "tests/fixtures/test_queries.fa"
        outdir = "tests/integration/out"
        outprefix = "test"
        database = "tests/fixtures/hmmer_custom_dbs/bact.hmm"

        # Observed and expected files
        obs_hmm_hits = os.path.join(outdir, outprefix+HMM_HITS_SUFFIX)        

        exp_files_dir = "tests/fixtures/hmmer_expected_output/"
        exp_hmm_hits = os.path.join(exp_files_dir, "bact.emapper.hmm_hits")

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = f'./hmm_mapper.py -i {in_file} -d {database} --output_dir {outdir} -o {outprefix}'

        # print(f"\t{cmd}")
        
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

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return

    
    def test_usemem(self):
        '''
        Tests hmm_mapper to search ("hmmscan" mode) on a custom db, on mem (--usemem)
        '''
        # ./hmm_mapper.py --usemem -i tests/fixtures/test_queries.fa -d tests/fixtures/hmmer_custom_dbs/bact.hmm -o dummy --output_dir tmp_borrar
        
        ##
        # Setup test
        
        in_file = "tests/fixtures/test_queries.fa"
        outdir = "tests/integration/out"
        outprefix = "test"
        database = "tests/fixtures/hmmer_custom_dbs/bact.hmm"

        # Observed and expected files
        obs_hmm_hits = os.path.join(outdir, outprefix+HMM_HITS_SUFFIX)        

        exp_files_dir = "tests/fixtures/hmmer_expected_output/"
        exp_hmm_hits = os.path.join(exp_files_dir, "bact.emapper.hmm_hits")

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = f'./hmm_mapper.py --usemem -i {in_file} -d {database} --output_dir {outdir} -o {outprefix}'

        # print(f"\t{cmd}")
        
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
        # We can not use check_hmm_hits, because the decimal precision changes
        # between using server or not
        # Therefore, we will compare only the first 2 columns: query and hit
        # check_hmm_hits(obs_hmm_hits, exp_hmm_hits)
        check_hmm_query_hit(obs_hmm_hits, exp_hmm_hits)

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return



    # Tests for "hmmsearch" mode
    def test_hmmsearch(self):
        '''
        Tests hmm_mapper to search ("hmmsearch" mode) on a custom db, on disk
        '''
        # ./hmm_mapper.py --qtype hmm -i tests/fixtures/hmmer_custom_dbs/bact.hmm \
        # --dbtype seqdb -d tests/fixtures/test_queries.fa \
        # -o dummy --output_dir tmp_borrar
        
        ##
        # Setup test
        
        database = "tests/fixtures/test_queries.fa"
        outdir = "tests/integration/out"
        outprefix = "test"
        in_file = "tests/fixtures/hmmer_custom_dbs/bact.hmm"

        # Observed and expected files
        obs_hmm_hits = os.path.join(outdir, outprefix+HMM_HITS_SUFFIX)        

        exp_files_dir = "tests/fixtures/hmmer_expected_output/"
        exp_hmm_hits = os.path.join(exp_files_dir, "bact.hmmsearch.hmm_hits")

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = f'./hmm_mapper.py --qtype hmm -i {in_file} --dbtype seqdb -d {database} --output_dir {outdir} -o {outprefix}'

        # print(f"\t{cmd}")
        
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

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return

    
    def test_hmmsearch_usemem(self):
        '''
        Tests hmm_mapper to search ("hmmsearch" mode) on a custom db, on mem (--usemem)
        '''
        # ./hmm_mapper.py --usemem \
        #         --qtype hmm -i tests/fixtures/hmmer_custom_dbs/bact.hmm \
        #         --dbtype seqdb -d tests/fixtures/test_queries.fa \
        #         -o dummy --output_dir tmp_borrar
        
        ##
        # Setup test
        
        database = "tests/fixtures/test_queries.fa"
        outdir = "tests/integration/out"
        outprefix = "test"
        in_file = "tests/fixtures/hmmer_custom_dbs/bact.hmm"

        # Observed and expected files
        obs_hmm_hits = os.path.join(outdir, outprefix+HMM_HITS_SUFFIX)        

        exp_files_dir = "tests/fixtures/hmmer_expected_output/"
        exp_hmm_hits = os.path.join(exp_files_dir, "bact.hmmsearch.hmm_hits")

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = f'./hmm_mapper.py --usemem --qtype hmm -i {in_file} --dbtype seqdb -d {database} --output_dir {outdir} -o {outprefix}'

        # print(f"\t{cmd}")
        
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
        # We can not use check_hmm_hits, because the decimal precision changes
        # between using server or not
        # Therefore, we will compare only the first 2 columns: query and hit
        # check_hmm_hits(obs_hmm_hits, exp_hmm_hits)
        check_hmm_query_hit(obs_hmm_hits, exp_hmm_hits)

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return
    

    # Tests for "phmmer" mode
    def test_phmmer(self):
        '''
        Tests hmm_mapper to search ("phmmer" mode) on a custom db, on disk
        '''
        # ./hmm_mapper.py --qtype seq -i tests/fixtures/test_queries.fa \
        #         --dbtype seqdb -d tests/fixtures/test_queries.fa \
        #         -o dummy --output_dir tmp_borrar
        
        ##
        # Setup test
        
        database = "tests/fixtures/test_queries.fa"
        outdir = "tests/integration/out"
        outprefix = "test"
        in_file = "tests/fixtures/test_queries.fa"

        # Observed and expected files
        obs_hmm_hits = os.path.join(outdir, outprefix+HMM_HITS_SUFFIX)        

        exp_files_dir = "tests/fixtures/hmmer_expected_output/"
        exp_hmm_hits = os.path.join(exp_files_dir, "bact.phmmer.hmm_hits")

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = f'./hmm_mapper.py --qtype seq -i {in_file} --dbtype seqdb -d {database} --output_dir {outdir} -o {outprefix}'

        # print(f"\t{cmd}")
        
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

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return

    
    def test_phmmer_usemem(self):
        '''
        Tests hmm_mapper to search ("phmmer" mode) on a custom db, on mem (--usemem)
        '''
        # ./hmm_mapper.py --usemem \
        #         --qtype seq -i tests/fixtures/test_queries.fa \
        #         --dbtype seqdb -d tests/fixtures/test_queries.fa \
        #         -o dummy --output_dir tmp_borrar
        
        ##
        # Setup test
        
        database = "tests/fixtures/test_queries.fa"
        outdir = "tests/integration/out"
        outprefix = "test"
        in_file = "tests/fixtures/test_queries.fa"

        # Observed and expected files
        obs_hmm_hits = os.path.join(outdir, outprefix+HMM_HITS_SUFFIX)        

        exp_files_dir = "tests/fixtures/hmmer_expected_output/"
        exp_hmm_hits = os.path.join(exp_files_dir, "bact.phmmer.hmm_hits")

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = f'./hmm_mapper.py --usemem --qtype seq -i {in_file} --dbtype seqdb -d {database} --output_dir {outdir} -o {outprefix}'

        # print(f"\t{cmd}")
        
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
        # We can not use check_hmm_hits, because the decimal precision changes
        # between using server or not
        # Therefore, we will compare only the first 2 columns: query and hit
        # check_hmm_hits(obs_hmm_hits, exp_hmm_hits)
        check_hmm_query_hit(obs_hmm_hits, exp_hmm_hits)

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return

    
if __name__ == '__main__':
    unittest.main()

## END
