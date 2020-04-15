##
## JHCepas
## CPCantalapiedra 2020

import unittest
import os, shutil, subprocess

# Settings for this tests
EMAPPER_CMD = './emapper.py'
EXPECTED_SEED_ORTHOLOGS_FN = 'test_output.emapper.seed_orthologs'
EXPECTED_ANNOTATIONS_FN = 'test_output.emapper.annotations'

# General eggnog-mapper settings
SEED_ORTHOLOGS_SUFFIX = '.emapper.seed_orthologs'
ANNOTATIONS_SUFFIX = '.emapper.annotations'

def run(cmd):
    '''
    Runs eggnog-mapper with the arguments specified
    '''
    cmd = EMAPPER_CMD + " " + cmd
    print cmd
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if not out:
        out = ''
    if not err:
        err = ''
    return (process.returncode, out, err)


def _load_orthologs(orthologs_fn):
    '''
    Loads the rows from a seed orthologs output, to a list
    '''
    _rows = []
    with open(orthologs_fn, 'r') as orth_f:
        for line in orth_f:
            if line.startswith("#"): continue
            _rows.append(line)
    return _rows

def _load_annotations(annotations_fn):
    '''
    Loads the rows from an annotations output, to a list
    '''
    _rows = []
    with open(annotations_fn, 'r') as annot_f:
        for line in annot_f:
            if line.startswith("#"): continue
            _rows.append(line)
    return _rows

class Test(unittest.TestCase):

    def _check_seed_orthologs(self, data_dir, outdir, outprefix):
        '''
        Compares the obtained seed orthologs file with the expected one
        '''
        # Check that output file has been created
        out_fn = os.path.join(outdir, outprefix+SEED_ORTHOLOGS_SUFFIX)
        assert os.path.exists(out_fn)

        # Compare expected and observed output
        # Load test output
        out_rows = _load_orthologs(out_fn)
        
        # Load expected output
        exp_fn = os.path.join(data_dir, EXPECTED_SEED_ORTHOLOGS_FN)
        exp_rows = _load_orthologs(exp_fn)

        # compare both files 
        self._basic_rows_comparison(out_rows, exp_rows)
        return

    def _check_annotations(self, data_dir, outdir, outprefix):
        '''
        Compares the obtained annotations file with the expected one
        '''
        # Check that output file has been created
        out_fn = os.path.join(outdir, outprefix+ANNOTATIONS_SUFFIX)
        assert os.path.exists(out_fn)

        # Compare expected and observed output
        # Load test output
        out_rows = _load_annotations(out_fn)
        
        # Load expected output
        exp_fn = os.path.join(data_dir, EXPECTED_ANNOTATIONS_FN)
        exp_rows = _load_annotations(exp_fn)

        # compare both files 
        self._basic_rows_comparison(out_rows, exp_rows)
        return
    
    def _basic_rows_comparison(self, l1, l2):
        '''
        Performs several basic comparisons between 2 lists
        '''
        # check that both lists have the same number of rows
        assert len(l1) == len(l2)
        # check that rows are equal one-by-one
        differences = [i for i, j in zip(l1, l2) if i != j]
        assert len(differences) == 0
        return

    # Tests
    def test_emapper_cmd(self):
        '''
        Tests the whole emapper command, including the expected output files
        '''
        in_file = "tests/fixtures/test_queries.fa"
        data_dir = "tests/fixtures"
        outdir = "tests/out"
        outprefix = "test"
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)
        
        cmd = '-m diamond -i %s --data_dir %s --output_dir %s -o %s' % (in_file, data_dir, outdir, outprefix)

        print "Running eggnog-mapper..."
        st, out, err = run(cmd)
        assert st == 0 # check exit status is ok

        ##
        # Check alignment phase: detection of seed orthologs
        print "Checking seed orthologs..."
        self._check_seed_orthologs(data_dir, outdir, outprefix)

        ##
        # Check annotation phase
        print "Checking annotations..."
        self._check_annotations(data_dir, outdir, outprefix)

        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return

if __name__ == '__main__':
    unittest.main()

## END
