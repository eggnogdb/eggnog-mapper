##
## CPCantalapiedra 2020

import shutil, os, unittest
from argparse import Namespace

from eggnogmapper.common import DIAMOND
from eggnogmapper.search.diamond import DiamondSearcher

class Test(unittest.TestCase):
    
    def setUp(self):
        os.mkdir("tests/unit/out/")

    def tearDown(self):
        if os.path.exists("tests/unit/out"):
            shutil.rmtree("tests/unit/out")
            
    # Tests
    def test_run_diamond_blastp(self):
        args = Namespace(translate=False,
                         dmnd_db="tests/fixtures/eggnog_proteins.dmnd",
                         cpu=2,
                         query_cover=0,
                         subject_cover=0,
                         matrix=None,
                         gapopen=None,
                         gapextend=None,
                         seed_ortholog_evalue=0.001,
                         seed_ortholog_score=60,
                         excluded_taxa=None,
                         temp_dir=os.getcwd(),
                         no_file_comments=None)
        
        fasta_file = "tests/fixtures/test_queries.fa"
        output_file = "tests/unit/out/test_run_diamond_blastp.seed_orthologs"
                
        searcher = DiamondSearcher(args)
        cmd = searcher.run_diamond(fasta_file, output_file, silent = True)
        
        self.assertIsNotNone(cmd)
        self.assertTrue(cmd.startswith(DIAMOND))
        
        return

    def test_run_diamond_blastx(self):
        args = Namespace(translate=True,
                         dmnd_db="tests/fixtures/eggnog_proteins.dmnd",
                         cpu=2,
                         query_cover=0,
                         subject_cover=0,
                         matrix=None,
                         gapopen=None,
                         gapextend=None,
                         seed_ortholog_evalue=0.001,
                         seed_ortholog_score=60,
                         excluded_taxa=None,
                         temp_dir=os.getcwd(),
                         no_file_comments=None)
        
        fasta_file = "tests/fixtures/test_queries.fna"
        output_file = "tests/unit/out/test_run_diamond_blastx.seed_orthologs"
                
        searcher = DiamondSearcher(args)
        
        cmd = searcher.run_diamond(fasta_file, output_file, silent = True)
        
        self.assertIsNotNone(cmd)
        self.assertTrue(cmd.startswith(DIAMOND))
        
        return
    
if __name__ == '__main__':
    
    unittest.main()

## END
