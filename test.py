import unittest
import subprocess
import os

def run(cmd):
    cmd = './emapper.py ' + cmd
    print cmd
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
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
        '-i test/polb.fa -d bact:localhost:51500 -o borrame2 --annotate_hits_table borrame.emapper.seed_orthologs'


        # should run from disk
        '-i test/polb.fa -d maNOG -o borrame --override'

        # should
        '-i test/polb.fa -d maNOG -o borrame --override'




if __name__ == '__main__':
    unittest.main()
