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

    def test_annotations(self):
        os.system('rm polb.emapper.* -f')
        st, out, err = run('-i test/polb.fa -d bact -o polb')
        ok = False
        for line in open("polb.emapper.annotations"):
            if line.startswith('362663.ECP_0061'):
                q, hit, evalue, score, pname, gos, ko, bigg, level, groups, bestgroup, cat, desc =  map(str.strip, line.split('\t'))
                assert hit == '362663.ECP_0061'
                assert pname == 'POLB'
                assert set(gos.split(',')) == set("GO:0003674,GO:0003824,GO:0003887,GO:0004518,GO:0004527,GO:0004529,GO:0004536,GO:0005575,GO:0005622,GO:0005623,GO:0005694,GO:0006139,GO:0006259,GO:0006260,GO:0006261,GO:0006281,GO:0006289,GO:0006297,GO:0006301,GO:0006725,GO:0006807,GO:0006950,GO:0006974,GO:0007154,GO:0008150,GO:0008152,GO:0008296,GO:0008408,GO:0009058,GO:0009059,GO:0009432,GO:0009605,GO:0009987,GO:0009991,GO:0016740,GO:0016772,GO:0016779,GO:0016787,GO:0016788,GO:0016796,GO:0016895,GO:0018130,GO:0019438,GO:0019985,GO:0031668,GO:0033554,GO:0034061,GO:0034641,GO:0034645,GO:0034654,GO:0043170,GO:0043226,GO:0043228,GO:0043229,GO:0043232,GO:0044237,GO:0044238,GO:0044249,GO:0044260,GO:0044271,GO:0044424,GO:0044464,GO:0044699,GO:0044763,GO:0045004,GO:0045005,GO:0046483,GO:0050896,GO:0051716,GO:0071496,GO:0071704,GO:0071897,GO:0090304,GO:0090305,GO:1901360,GO:1901362,GO:1901576".split(','))
                assert ko == 'K02336'
                assert cat == 'L'
                assert desc == 'DNA polymerase'
                ok = True
        if not ok:
            ValueError('polB result not obtained')
        else:
            os.system('rm polb.emmaper.* -f')

    
    def test_execution_modes(self):
        os.system('rm borrame.emapper.* -f')

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
        os.system('rm borrame.emapper.* -f')
        st, out, err = run('-i test/polb.fa -d maNOG -o borrame')
        assert st == 0

        # should run from mem
        os.system('rm borrame.emapper.* -f')
        st, out, err = run('-i test/polb.fa -d homNOG -o borrame --usemem')
        assert st == 0

        # should run from mem
        os.system('rm borrame.emapper.* -f')
        st, out, err = run('-i test/polb.fa -o borrame -m diamond --cpu 10')
        assert st == 0

        os.system('rm borrame.emapper.* -f')

if __name__ == '__main__':
    unittest.main()
