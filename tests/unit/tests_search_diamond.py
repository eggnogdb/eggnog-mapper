##
## CPCantalapiedra 2020

import os, unittest
from argparse import Namespace

from eggnogmapper.search.diamond import DiamondSearcher

class Test(unittest.TestCase):
    # Tests
    def test_run_diamond_blastp(self):
        args = Namespace(translate=False,
                         dmnd_db=None,
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
        
        print(args)
        searcher = DiamondSearcher(args)
        fasta_file = "tests/fixtures/test_queries.fa"
        output_file = "tests/unit/test_run_diamond_blastp.seed_orthologs"
        searcher.run_diamond(fasta_file, output_file)
        return
    
if __name__ == '__main__':
    unittest.main()
