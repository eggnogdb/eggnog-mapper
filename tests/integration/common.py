##
## CPCantalapiedra 2020

import os, subprocess

EMAPPER_CMD = './emapper.py'

def run(cmd):
    '''
    Runs eggnog-mapper with the arguments specified
    '''
    cmd = EMAPPER_CMD + " " + cmd
    
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if not out:
        out = b''
    if not err:
        err = b''
    return (process.returncode, out, err)

##
def _load_commented_file(fn):
    '''
    Loads the rows from a commentd file, skipping comments
    '''
    _rows = []
    with open(fn, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            _rows.append(line)
    return _rows    

def load_hmm_hits(hmm_hits_fn):
    '''
    Loads the rows from a hmm_hits file, to a list
    '''
    return _load_commented_file(hmm_hits_fn)

def load_orthologs(orthologs_fn):
    '''
    Loads the rows from a seed orthologs output, to a list
    '''
    return _load_commented_file(orthologs_fn)


def load_annotations(annotations_fn):
    '''
    Loads the rows from an annotations output, to a list
    '''
    return _load_commented_file(annotations_fn)

##
def check_hmm_hits(obs_out, exp_out):
    '''
    Compares the obtained hmm_hits file with the expected one
    '''
    # Check that output file has been created
    assert os.path.exists(obs_out)

    # Compare expected and observed output
    # Load test output
    obs_rows = load_hmm_hits(obs_out)

    # Load expected output
    exp_rows = load_hmm_hits(exp_out)

    # compare both files 
    _basic_rows_comparison(obs_rows, exp_rows)
    return

def check_seed_orthologs(obs_out, exp_out):
    '''
    Compares the obtained seed orthologs file with the expected one
    '''
    # Check that output file has been created
    assert os.path.exists(obs_out)

    # Compare expected and observed output
    # Load test output
    obs_rows = load_orthologs(obs_out)

    # Load expected output
    exp_rows = load_orthologs(exp_out)

    # compare both files 
    _basic_rows_comparison(obs_rows, exp_rows)
    return

def check_annotations(obs_out, exp_out):
    '''
    Compares the obtained annotations file with the expected one
    '''
    # Check that output file has been created
    assert os.path.exists(obs_out)

    # Compare expected and observed output
    # Load test output
    obs_rows = load_annotations(obs_out)

    # Load expected output
    exp_rows = load_annotations(exp_out)

    # compare both files 
    _basic_rows_comparison(obs_rows, exp_rows)
    return


def _basic_rows_comparison(l1, l2):
    '''
    Performs several basic comparisons between 2 lists
    '''
    # check that both lists have the same number of rows
    assert len(l1) == len(l2)
    # check that rows are equal one-by-one
    differences = [i for i, j in zip(l1, l2) if i != j]
    assert len(differences) == 0
    return

## END
