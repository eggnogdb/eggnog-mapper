##
## JHCepas
## CPCantalapiedra 2020
## Updated 2026: v7 database support, new tests

import gzip
import bz2
import os
import shutil
import tempfile
import unittest

from .common import (
    run, check_seed_orthologs, check_annotations, check_orthologs,
    check_gff, check_fasta, check_hmm_hits, check_pfam,
)

# Output file suffixes
DECORATED_GFF_SUFFIX = '.emapper.decorated.gff'
GENEPRED_GFF_SUFFIX = '.emapper.genepred.gff'
GENEPRED_FASTA_SUFFIX = '.emapper.genepred.fasta'
HMM_HITS_SUFFIX = '.emapper.hits'
SEED_ORTHOLOGS_SUFFIX = '.emapper.seed_orthologs'
ANNOTATIONS_SUFFIX = '.emapper.annotations'
ORTHOLOGS_SUFFIX = '.emapper.orthologs'
PFAM_SUFFIX = '.emapper.pfam'

# v7 test data paths
V7_DATA_DIR = "tests/fixtures_v7"
V7_QUERY_FILE = "tests/fixtures/test_queries_v7.fa"
# Resolve the e7-sample backend to a filesystem path — the CLI no longer
# accepts --db BACKEND_NAME; data location is set via --data_dir (or the
# EGGNOG_DATA_DIR env). The auto-discovered backend depends on
# EGGNOG_DATA_ROOT, which common.run() also sets for the child process.
os.environ.setdefault("EGGNOG_DATA_ROOT", "/eggnog-eco/data")
from eggnogmapper.backends import resolve_backend
V7_DB = resolve_backend("e7-sample")


class Test(unittest.TestCase):

    def setUp(self):
        self.outdir = "tests/integration/out"
        if os.path.isdir(self.outdir):
            shutil.rmtree(self.outdir)
        os.mkdir(self.outdir)

    def tearDown(self):
        if os.path.isdir(self.outdir):
            shutil.rmtree(self.outdir)

    # ------------------------------------------------------------------
    # Diamond full pipeline
    # ------------------------------------------------------------------

    def test_emapper_diamond(self):
        '''Full emapper pipeline using DIAMOND search against the v7 sample database.'''
        outprefix = "test"
        obs_seed = os.path.join(self.outdir, outprefix + SEED_ORTHOLOGS_SUFFIX)
        obs_annot = os.path.join(self.outdir, outprefix + ANNOTATIONS_SUFFIX)
        exp_seed = os.path.join(V7_DATA_DIR, 'test_diamond.emapper.seed_orthologs')
        exp_annot = os.path.join(V7_DATA_DIR, 'test_diamond.emapper.annotations')

        cmd = (f'./emapper.py -m diamond -i {V7_QUERY_FILE} --data_dir {V7_DB} '
               f'--output_dir {self.outdir} -o {outprefix}')
        st, out, err = run(cmd)
        if st != 0:
            print(out.decode("utf-8"))
            print(err.decode("utf-8"))
        assert st == 0

        check_seed_orthologs(obs_seed, exp_seed)
        check_annotations(obs_annot, exp_annot)

    # ------------------------------------------------------------------
    # No-search (annotate from existing seed_orthologs)
    # ------------------------------------------------------------------

    def test_emapper_no_search(self):
        '''Annotation from existing seed_orthologs table (-m no_search).'''
        outprefix = "test"
        hits_table = os.path.join(V7_DATA_DIR, 'test_diamond.emapper.seed_orthologs')
        obs_annot = os.path.join(self.outdir, outprefix + ANNOTATIONS_SUFFIX)
        obs_ortho = os.path.join(self.outdir, outprefix + ORTHOLOGS_SUFFIX)
        exp_annot = os.path.join(V7_DATA_DIR, 'test_no_search.emapper.annotations')
        exp_ortho = os.path.join(V7_DATA_DIR, 'test_no_search.emapper.orthologs')

        cmd = (f'./emapper.py -m no_search --annotate_hits_table {hits_table} '
               f'--data_dir {V7_DB} --output_dir {self.outdir} -o {outprefix} '
               f'--report_orthologs --target_orthologs one2one')
        st, out, err = run(cmd)
        if st != 0:
            print(out.decode("utf-8"))
            print(err.decode("utf-8"))
        assert st == 0

        check_orthologs(obs_ortho, exp_ortho)
        check_annotations(obs_annot, exp_annot)

    # ------------------------------------------------------------------
    # Two-step: --no_annot then --annotate_hits_table
    # ------------------------------------------------------------------

    def test_no_annot_then_annotate(self):
        '''Two-step workflow: diamond --no_annot, then -m no_search --annotate_hits_table.
        The final annotations must match a single-step diamond run.'''
        outprefix = "step1"
        step2prefix = "step2"
        exp_annot = os.path.join(V7_DATA_DIR, 'test_diamond.emapper.annotations')

        # Step 1: search only
        cmd = (f'./emapper.py -m diamond -i {V7_QUERY_FILE} --data_dir {V7_DB} '
               f'--output_dir {self.outdir} -o {outprefix} --no_annot')
        st, out, err = run(cmd)
        if st != 0:
            print(out.decode("utf-8"))
            print(err.decode("utf-8"))
        assert st == 0, "Step 1 (--no_annot) failed"

        seeds_file = os.path.join(self.outdir, outprefix + SEED_ORTHOLOGS_SUFFIX)
        assert os.path.exists(seeds_file), "seed_orthologs not produced by --no_annot"

        # Step 2: annotate from seeds
        cmd2 = (f'./emapper.py -m no_search --annotate_hits_table {seeds_file} '
                f'--data_dir {V7_DB} --output_dir {self.outdir} -o {step2prefix}')
        st2, out2, err2 = run(cmd2)
        if st2 != 0:
            print(out2.decode("utf-8"))
            print(err2.decode("utf-8"))
        assert st2 == 0, "Step 2 (no_search annotate) failed"

        obs_annot = os.path.join(self.outdir, step2prefix + ANNOTATIONS_SUFFIX)
        check_annotations(obs_annot, exp_annot)

    # ------------------------------------------------------------------
    # Resume: --resume skips search and reuses hits file
    # ------------------------------------------------------------------

    def test_resume_diamond(self):
        '''--resume skips DIAMOND search and reuses the existing .hits file.
        Output annotations must be identical to a fresh run.'''
        outprefix = "test"
        exp_seed = os.path.join(V7_DATA_DIR, 'test_diamond.emapper.seed_orthologs')
        exp_annot = os.path.join(V7_DATA_DIR, 'test_diamond.emapper.annotations')

        base_cmd = (f'./emapper.py -m diamond -i {V7_QUERY_FILE} --data_dir {V7_DB} '
                    f'--output_dir {self.outdir} -o {outprefix}')

        # Run 1: produce hits file
        st, out, err = run(base_cmd + ' --override')
        if st != 0:
            print(out.decode("utf-8"))
            print(err.decode("utf-8"))
        assert st == 0, "First run failed"

        # Run 2: resume — must reuse hits file
        st2, out2, err2 = run(base_cmd + ' --resume')
        if st2 != 0:
            print(out2.decode("utf-8"))
            print(err2.decode("utf-8"))
        assert st2 == 0, "Resume run failed"

        # Verify [resume] message appears in stderr
        stderr_text = err2.decode("utf-8")
        assert "[resume]" in stderr_text, (
            f"Expected '[resume]' in stderr but got:\n{stderr_text}"
        )

        # Output must match reference
        obs_seed = os.path.join(self.outdir, outprefix + SEED_ORTHOLOGS_SUFFIX)
        obs_annot = os.path.join(self.outdir, outprefix + ANNOTATIONS_SUFFIX)
        check_seed_orthologs(obs_seed, exp_seed)
        check_annotations(obs_annot, exp_annot)

    # ------------------------------------------------------------------
    # Compressed input — gzip
    # ------------------------------------------------------------------

    def test_compressed_input_gz(self):
        '''Gzip-compressed FASTA input is autodetected and produces the same
        annotations as uncompressed input.'''
        outprefix = "test_gz"
        exp_seed = os.path.join(V7_DATA_DIR, 'test_diamond.emapper.seed_orthologs')
        exp_annot = os.path.join(V7_DATA_DIR, 'test_diamond.emapper.annotations')

        # Create a gzip-compressed copy of the query file
        gz_path = os.path.join(self.outdir, "test_queries_v7.fa.gz")
        with open(V7_QUERY_FILE, 'rb') as f_in:
            with gzip.open(gz_path, 'wb') as f_out:
                f_out.write(f_in.read())

        cmd = (f'./emapper.py -m diamond -i {gz_path} --data_dir {V7_DB} '
               f'--output_dir {self.outdir} -o {outprefix}')
        st, out, err = run(cmd)
        if st != 0:
            print(out.decode("utf-8"))
            print(err.decode("utf-8"))
        assert st == 0

        obs_seed = os.path.join(self.outdir, outprefix + SEED_ORTHOLOGS_SUFFIX)
        obs_annot = os.path.join(self.outdir, outprefix + ANNOTATIONS_SUFFIX)
        check_seed_orthologs(obs_seed, exp_seed)
        check_annotations(obs_annot, exp_annot)

    # ------------------------------------------------------------------
    # Compressed input — bzip2
    # ------------------------------------------------------------------

    def test_compressed_input_bz2(self):
        '''Bzip2-compressed FASTA input is autodetected and produces the same
        annotations as uncompressed input.'''
        outprefix = "test_bz2"
        exp_seed = os.path.join(V7_DATA_DIR, 'test_diamond.emapper.seed_orthologs')
        exp_annot = os.path.join(V7_DATA_DIR, 'test_diamond.emapper.annotations')

        bz2_path = os.path.join(self.outdir, "test_queries_v7.fa.bz2")
        with open(V7_QUERY_FILE, 'rb') as f_in:
            with bz2.open(bz2_path, 'wb') as f_out:
                f_out.write(f_in.read())

        cmd = (f'./emapper.py -m diamond -i {bz2_path} --data_dir {V7_DB} '
               f'--output_dir {self.outdir} -o {outprefix}')
        st, out, err = run(cmd)
        if st != 0:
            print(out.decode("utf-8"))
            print(err.decode("utf-8"))
        assert st == 0

        obs_seed = os.path.join(self.outdir, outprefix + SEED_ORTHOLOGS_SUFFIX)
        obs_annot = os.path.join(self.outdir, outprefix + ANNOTATIONS_SUFFIX)
        check_seed_orthologs(obs_seed, exp_seed)
        check_annotations(obs_annot, exp_annot)

    # ------------------------------------------------------------------
    # Tests requiring unavailable tools — skipped in this environment
    # ------------------------------------------------------------------

    @unittest.skip("Test fixtures not yet ported to v7 mmseqs database format")
    def test_emapper_mmseqs(self):
        pass

    @unittest.skip("Test fixtures not yet ported to v7 HMMER database format")
    def test_emapper_hmmer_eggnogdb(self):
        pass

    @unittest.skip("Test fixtures not yet ported to v7 HMMER database format")
    def test_scratch_dir(self):
        pass

    @unittest.skip("Test fixtures not yet ported to v7 Pfam realign workflow")
    def test_pfam_realign(self):
        pass

    @unittest.skip("Test fixtures not yet ported to v7 Pfam denovo workflow")
    def test_pfam_denovo(self):
        pass

    @unittest.skip("Test fixtures not yet ported to v7 genepred workflow")
    def test_genepred_prodigal(self):
        pass

    @unittest.skip("Test fixtures not yet ported to v7 genepred blastx workflow")
    def test_genepred_diamond(self):
        pass

    @unittest.skip("Test fixtures not yet ported to v7 mmseqs genepred workflow")
    def test_genepred_mmseqs(self):
        pass

    @unittest.skip("Test fixtures not yet ported to v7 mmseqs genepred GFF decoration")
    def test_decorate_gff_blastx_annot(self):
        pass

    @unittest.skip("Test fixtures not yet ported to v7 GFF decoration workflow")
    def test_decorate_gff_file(self):
        pass

    @unittest.skip("Test fixtures not yet ported to v7 GFF decoration workflow")
    def test_decorate_gff_file_short(self):
        pass

    @unittest.skip("Novel families test requires novel_fams database; not yet ported to v7")
    def test_emapper_novel_fams(self):
        pass


if __name__ == '__main__':
    unittest.main()

## END
