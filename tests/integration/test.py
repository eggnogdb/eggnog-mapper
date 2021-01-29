##
## JHCepas
## CPCantalapiedra 2020

import unittest
import os, shutil

from .common import run, check_gff, check_fasta, check_seed_orthologs, check_annotations, check_hmm_hits, check_orthologs, check_pfam

# General eggnog-mapper settings
GENEPRED_GFF_SUFFIX = '.emapper.genepred.gff'
GENEPRED_FASTA_SUFFIX = '.emapper.genepred.fasta'
HMM_HITS_SUFFIX = '.emapper.hmm_hits'
SEED_ORTHOLOGS_SUFFIX = '.emapper.seed_orthologs'
ANNOTATIONS_SUFFIX = '.emapper.annotations'
ORTHOLOGS_SUFFIX = '.emapper.orthologs'
PFAM_SUFFIX = '.emapper.pfam'

class Test(unittest.TestCase):

    # Tests
    def test_emapper_diamond(self):
        '''
        Tests the whole emapper (-m diamond) command
        '''
        # ./emapper.py -m diamond -i tests/fixtures/test_queries.fa --data_dir tests/fixtures --output_dir tmp -o test

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

        cmd = f'./emapper.py -m diamond -i {in_file} --data_dir {data_dir} --output_dir {outdir} -o {outprefix}'

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

    def test_emapper_mmseqs(self):
        '''
        Tests the whole emapper (-m mmseqs) command
        '''
        # ./emapper.py -m mmseqs -i tests/fixtures/test_queries.fa --data_dir tests/fixtures --output_dir tmp -o test
        
        ##
        # Setup test
        
        in_file = "tests/fixtures/test_queries.fa"
        data_dir = "tests/fixtures"
        outdir = "tests/integration/out"
        outprefix = "test"

        # Observed and expected files
        obs_seed_orthologs = os.path.join(outdir, outprefix+SEED_ORTHOLOGS_SUFFIX)
        obs_annotations = os.path.join(outdir, outprefix+ANNOTATIONS_SUFFIX)
        
        exp_seed_orthologs = os.path.join(data_dir, 'test_mmseqs.emapper.seed_orthologs')
        exp_annotations = os.path.join(data_dir, 'test_mmseqs.emapper.annotations')

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = f'./emapper.py -m mmseqs -i {in_file} --data_dir {data_dir} --output_dir {outdir} -o {outprefix}'

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
    
    
    def test_emapper_no_search(self):
        '''
        Tests annotation (-m no_search) of an existing hits table, and reports orthologs (--report_orthologs, --target_orthologs, --target_taxa and --excluded_taxa)
        '''
        ('./emapper.py -m no_search --annotate_hits_table tests/fixtures/test_output.emapper.seed_orthologs '
         '--data_dir tests/fixtures --output_dir tmp -o tmp --report_orthologs --target_orthologs one2one --target_taxa 72274,1123487 --excluded_taxa 205918,1395571')
        
        ##
        # Setup test
        
        in_file = "tests/fixtures/test_output.emapper.seed_orthologs"
        data_dir = "tests/fixtures"
        outdir = "tests/integration/out"
        outprefix = "test"

        # Observed and expected files
        obs_annotations = os.path.join(outdir, outprefix+ANNOTATIONS_SUFFIX)
        obs_orthologs = os.path.join(outdir, outprefix+ORTHOLOGS_SUFFIX)
        
        exp_annotations = os.path.join(data_dir, 'test_no_search.emapper.annotations')
        exp_orthologs = os.path.join(data_dir, 'test_no_search.emapper.orthologs')

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = (f'./emapper.py -m no_search --annotate_hits_table {in_file} --data_dir {data_dir} --output_dir {outdir} -o {outprefix} --report_orthologs '
               f'--target_orthologs one2one --target_taxa 72274,1123487 --excluded_taxa 205918,1395571')

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

        # Check orthologs
        check_orthologs(obs_orthologs, exp_orthologs)
        
        # Check annotation phase
        check_annotations(obs_annotations, exp_annotations)

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return

    
    def test_emapper_hmmer_eggnogdb(self):
        '''
        Tests the whole emapper (-m hmmer) command against a eggNOG DB
        '''
        # ./emapper.py -m hmmer -i tests/fixtures/test_queries.fa --data_dir tests/fixtures -d bact --output_dir tmp -o bact
        
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

        cmd = f'./emapper.py -m hmmer -i {in_file} --data_dir {data_dir} -d {database} --output_dir {outdir} -o {outprefix}'

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

    
    def test_scratch_dir(self):
        '''
        Tests the use of scratch_dir
        '''
        # ./emapper.py -m hmmer -i tests/fixtures/test_queries.fa --data_dir tests/fixtures -d bact -o bact --output_dir tmp_borrar --scratch_dir tmp_scratch
        
        ##
        # Setup test
        
        in_file = "tests/fixtures/test_queries.fa"
        data_dir = "tests/fixtures"
        outdir = "tests/integration/out"
        scratchdir = "tests/integration/scratch"
        outprefix = "bact"
        database = "bact"

        # Observed and expected files
        obs_hmm_hits = os.path.join(outdir, outprefix+HMM_HITS_SUFFIX)        
        obs_seed_orthologs = os.path.join(outdir, outprefix+SEED_ORTHOLOGS_SUFFIX)
        obs_annotations = os.path.join(outdir, outprefix+ANNOTATIONS_SUFFIX)

        scratch_hmm_hits = os.path.join(scratchdir, outprefix+HMM_HITS_SUFFIX)        
        scratch_seed_orthologs = os.path.join(scratchdir, outprefix+SEED_ORTHOLOGS_SUFFIX)
        scratch_annotations = os.path.join(scratchdir, outprefix+ANNOTATIONS_SUFFIX)

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

        if os.path.isdir(scratchdir):
            shutil.rmtree(scratchdir)
        os.mkdir(scratchdir)

        cmd = f'./emapper.py -m hmmer -i {in_file} --data_dir {data_dir} -d {database} --output_dir {outdir} -o {outprefix} --scratch_dir {scratchdir}'

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
        check_hmm_hits(scratch_hmm_hits, exp_hmm_hits)
        
        # Check seed orthologs from alignment phase
        check_seed_orthologs(obs_seed_orthologs, exp_seed_orthologs)
        check_seed_orthologs(scratch_seed_orthologs, exp_seed_orthologs)
        
        # Check annotation phase
        check_annotations(obs_annotations, exp_annotations)
        check_annotations(scratch_annotations, exp_annotations)

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)

        if os.path.isdir(scratchdir):
            shutil.rmtree(scratchdir)
        
        return

    
    def test_pfam_transfer_narrowest_og(self):
        '''
        Tests --pfam_transfer narrowest_og
        '''
        # ./emapper.py -m diamond -i tests/fixtures/test_pfam_groups.fa --data_dir tests/fixtures --output_dir tmp -o test --pfam_transfer narrowest_og
        
        ##
        # Setup test
        
        in_file = "tests/fixtures/test_pfam_groups.fa"
        data_dir = "tests/fixtures"
        outdir = "tests/integration/out"
        outprefix = "test"

        # Observed and expected files
        obs_seed_orthologs = os.path.join(outdir, outprefix+SEED_ORTHOLOGS_SUFFIX)
        obs_annotations = os.path.join(outdir, outprefix+ANNOTATIONS_SUFFIX)
        
        exp_seed_orthologs = os.path.join(data_dir, 'pfam_transfer_narrowest_og_output', 'test.emapper.seed_orthologs')
        exp_annotations = os.path.join(data_dir, 'pfam_transfer_narrowest_og_output', 'test.emapper.annotations')

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = f'./emapper.py -m diamond -i {in_file} --data_dir {data_dir} --output_dir {outdir} -o {outprefix} --pfam_transfer narrowest_og'

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


    def test_pfam_align_seed(self):
        '''
        Tests --pfam_transfer seed_ortholog --pfam_realign realign
        '''
        # ./emapper.py -m diamond -i tests/fixtures/test_pfam_groups.fa --data_dir tests/fixtures --output_dir tmp -o test --pfam_transfer seed_ortholog --pfam_realign realign
        
        ##
        # Setup test
        
        in_file = "tests/fixtures/test_pfam_groups.fa"
        data_dir = "tests/fixtures"
        outdir = "tests/integration/out"
        outprefix = "test"

        # Observed and expected files
        obs_seed_orthologs = os.path.join(outdir, outprefix+SEED_ORTHOLOGS_SUFFIX)
        obs_annotations = os.path.join(outdir, outprefix+ANNOTATIONS_SUFFIX)
        obs_pfam = os.path.join(outdir, outprefix+PFAM_SUFFIX)
        
        exp_seed_orthologs = os.path.join(data_dir, 'pfam_align_seed_output', 'test.emapper.seed_orthologs')
        exp_annotations = os.path.join(data_dir, 'pfam_align_seed_output', 'test.emapper.annotations')
        exp_pfam = os.path.join(data_dir, 'pfam_align_seed_output', 'test.emapper.pfam')

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = f'./emapper.py -m diamond -i {in_file} --data_dir {data_dir} --output_dir {outdir} -o {outprefix} --pfam_transfer seed_ortholog --pfam_realign realign'

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
        
        # Check alignment phase: detection of seed orthologs
        check_seed_orthologs(obs_seed_orthologs, exp_seed_orthologs)

        # Check PFAM
        check_pfam(obs_pfam, exp_pfam)
        
        # Check annotation phase
        check_annotations(obs_annotations, exp_annotations)

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return
    

    def test_pfam_denovo(self):
        '''
        Tests --pfam_realign denovo
        '''
        # ./emapper.py -m diamond -i tests/fixtures/test_pfam_groups.fa --data_dir tests/fixtures --output_dir tmp -o test --pfam_realign denovo
        
        ##
        # Setup test
        
        in_file = "tests/fixtures/test_pfam_groups.fa"
        data_dir = "tests/fixtures"
        outdir = "tests/integration/out"
        outprefix = "test"

        # Observed and expected files
        obs_seed_orthologs = os.path.join(outdir, outprefix+SEED_ORTHOLOGS_SUFFIX)
        obs_annotations = os.path.join(outdir, outprefix+ANNOTATIONS_SUFFIX)
        obs_pfam = os.path.join(outdir, outprefix+PFAM_SUFFIX)
        
        exp_seed_orthologs = os.path.join(data_dir, 'pfam_denovo_output', 'test.emapper.seed_orthologs')
        exp_annotations = os.path.join(data_dir, 'pfam_denovo_output', 'test.emapper.annotations')
        exp_pfam = os.path.join(data_dir, 'pfam_denovo_output', 'test.emapper.pfam')

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = f'./emapper.py -m diamond -i {in_file} --data_dir {data_dir} --output_dir {outdir} -o {outprefix} --pfam_realign denovo'

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
        
        # Check alignment phase: detection of seed orthologs
        check_seed_orthologs(obs_seed_orthologs, exp_seed_orthologs)

        # Check PFAM
        check_pfam(obs_pfam, exp_pfam)
        
        # Check annotation phase
        check_annotations(obs_annotations, exp_annotations)

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return
    
    def test_genepred_prodigal(self):
        '''
        Test gene prediction with prodigal
        '''
        # ./emapper.py -i tests/fixtures/genepred_contig/contig.fna --itype metagenome --genepred prodigal -m diamond --sensmode sensitive --no_annot
        # --dmnd_db tests/fixtures/genepred_contig/contig.dmnd -o test --output_dir tmp
        
        ##
        # Setup test
        
        in_file = "tests/fixtures/genepred_contig/contig.fna"
        data_dir = "tests/fixtures/genepred_contig/prodigal_out"
        dmnd_db = "tests/fixtures/genepred_contig/contig.dmnd"
        outdir = "tests/integration/out"
        outprefix = "test"

        # Observed and expected files
        obs_genepred_gff = os.path.join(outdir, outprefix+GENEPRED_GFF_SUFFIX)
        obs_genepred_fasta = os.path.join(outdir, outprefix+GENEPRED_FASTA_SUFFIX)
        obs_seed_orthologs = os.path.join(outdir, outprefix+SEED_ORTHOLOGS_SUFFIX)

        exp_genepred_gff = os.path.join(data_dir, 'out.emapper.genepred.gff')
        exp_genepred_fasta = os.path.join(data_dir, 'out.emapper.genepred.fasta')
        exp_seed_orthologs = os.path.join(data_dir, 'out.emapper.seed_orthologs')

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = (f'./emapper.py -i {in_file} --itype metagenome '
               f'--genepred prodigal -m diamond --sensmode sensitive --no_annot '
               f'--dmnd_db {dmnd_db} '
               f'-o {outprefix} --output_dir {outdir}')
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

        # Check GFF from gene prediction
        check_gff(obs_genepred_gff, exp_genepred_gff)
        
        # Check FASTA from gene prediction
        check_fasta(obs_genepred_fasta, exp_genepred_fasta)
        
        # Check alignment phase: detection of seed orthologs
        check_seed_orthologs(obs_seed_orthologs, exp_seed_orthologs)

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return

    def test_genepred_diamond(self):
        '''
        Test gene prediction with diamond
        '''
        # ./emapper.py -i tests/fixtures/genepred_contig/contig.fna --itype metagenome --genepred search -m diamond --sensmode sensitive --no_annot --dmnd_db tests/fixtures/genepred_contig/contig.dmnd -o test --output_dir tmp_borrar
        
        ##
        # Setup test
        
        in_file = "tests/fixtures/genepred_contig/contig.fna"
        data_dir = "tests/fixtures/genepred_contig/diamond_out"
        dmnd_db = "tests/fixtures/genepred_contig/contig.dmnd"
        outdir = "tests/integration/out"
        outprefix = "test"

        # Observed and expected files
        obs_genepred_gff = os.path.join(outdir, outprefix+GENEPRED_GFF_SUFFIX)
        obs_genepred_fasta = os.path.join(outdir, outprefix+GENEPRED_FASTA_SUFFIX)
        obs_seed_orthologs = os.path.join(outdir, outprefix+SEED_ORTHOLOGS_SUFFIX)

        exp_genepred_gff = os.path.join(data_dir, 'out.emapper.genepred.gff')
        exp_genepred_fasta = os.path.join(data_dir, 'out.emapper.genepred.fasta')
        exp_seed_orthologs = os.path.join(data_dir, 'out.emapper.seed_orthologs')

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = (f'./emapper.py -i {in_file} --itype metagenome '
               f'--genepred search -m diamond --sensmode sensitive --no_annot '
               f'--dmnd_db {dmnd_db} '
               f'-o {outprefix} --output_dir {outdir}')
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

        # Check GFF from gene prediction
        check_gff(obs_genepred_gff, exp_genepred_gff)
        
        # Check FASTA from gene prediction
        check_fasta(obs_genepred_fasta, exp_genepred_fasta)
        
        # Check alignment phase: detection of seed orthologs
        check_seed_orthologs(obs_seed_orthologs, exp_seed_orthologs)

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return

    def test_genepred_mmseqs(self):
        '''
        Test gene prediction with mmseqs
        '''

        # ./emapper.py -i tests/fixtures/genepred_contig/contig.fna --itype metagenome --genepred search -m mmseqs --no_annot --mmseqs_db tests/fixtures/genepred_contig/contig.mmseqs/contig.0.hits.mmseqs.db -o test --output_dir tmp_borrar
        
        ##
        # Setup test
        
        in_file = "tests/fixtures/genepred_contig/contig.fna"
        data_dir = "tests/fixtures/genepred_contig/mmseqs_out"
        mmseqs_db = "tests/fixtures/genepred_contig/contig.mmseqs/contig.0.hits.mmseqs.db"
        outdir = "tests/integration/out"
        outprefix = "test"

        # Observed and expected files
        obs_genepred_gff = os.path.join(outdir, outprefix+GENEPRED_GFF_SUFFIX)
        obs_genepred_fasta = os.path.join(outdir, outprefix+GENEPRED_FASTA_SUFFIX)
        obs_seed_orthologs = os.path.join(outdir, outprefix+SEED_ORTHOLOGS_SUFFIX)

        exp_genepred_gff = os.path.join(data_dir, 'out.emapper.genepred.gff')
        exp_genepred_fasta = os.path.join(data_dir, 'out.emapper.genepred.fasta')
        exp_seed_orthologs = os.path.join(data_dir, 'out.emapper.seed_orthologs')

        ##
        # Run test
        
        # Remove (just in case) and recreate the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.mkdir(outdir)

        cmd = (f'./emapper.py -i {in_file} --itype metagenome '
               f'--genepred search -m mmseqs --no_annot '
               f'--mmseqs_db {mmseqs_db} '
               f'-o {outprefix} --output_dir {outdir}')
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

        # Check GFF from gene prediction
        check_gff(obs_genepred_gff, exp_genepred_gff)
        
        # Check FASTA from gene prediction
        check_fasta(obs_genepred_fasta, exp_genepred_fasta)
        
        # Check alignment phase: detection of seed orthologs
        check_seed_orthologs(obs_seed_orthologs, exp_seed_orthologs)

        ##
        # Teardown test
        
        # Remove the output dir
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        
        return
    
if __name__ == '__main__':
    unittest.main()

## END
