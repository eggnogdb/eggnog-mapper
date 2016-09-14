import unittest
import subprocess
import os

def run(cmd):
    cmd = './emapper.py ' + cmd
    print cmd
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = process.communicate()
    if not out:
        out = ''
    if not err:
        err = ''
    return (process.returncode, out, err)

class Test(unittest.TestCase):
    def test_executions(self):
        os.system('rm borrame.emapper.*')

        # should run
        st, out, err = run('-i test/polb.fa -d bact:localhost:51500 -o borrame')
        assert st == 0

        # should fail because output files exist
        st, out, err = run('-i test/polb.fa -d bact:localhost:51500 -o borrame')
        assert st == 1
        assert 'Use --resume' in out

        # should run
        st, out, err = run('-i test/polb.fa -d bact:localhost:51500 -o borrame --override')
        assert st == 0

        # should run with no search
        st, out, err = run('-i test/polb.fa -d bact:localhost:51500 -o borrame --resume')
        assert st == 0

        # should run with no search and read prev table
        st, out, err = run('-i test/polb.fa -d bact:localhost:51500 -o borrame --annotate_hits_table borrame.emapper.seed_orthologs --override')
        assert st == 0

        # should run from disk
        os.system('rm borrame.emapper.*')
        st, out, err = run('-i test/polb.fa -d maNOG -o borrame')
        assert st == 0

        # should run from mem
        os.system('rm borrame.emapper.*')
        st, out, err = run('-i test/polb.fa -d maNOG -o borrame --usemem')
        assert st == 0

        # should run from mem
        os.system('rm borrame.emapper.*')
        st, out, err = run('-i test/polb.fa -o borrame -m diamond --cpu 0')
        assert st == 0

if __name__ == '__main__':
    unittest.main()
