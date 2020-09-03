##
## CPCantalapiedra 2020

import errno, os, shutil
from os.path import exists as pexists
from os.path import join as pjoin

from .utils import colorify
from .common import silent_rm, ITYPE_GENOME, ITYPE_META
from .emapperException import EmapperException

from .genepred.prodigal import ProdigalPredictor
from .search.search_modes import get_searcher, SEARCH_MODE_NO_SEARCH
from .annotation.annotator import get_annotator

class Emapper:

    gene_pred_prots_file = hmm_hits_file = seed_orthologs_file = None
    annot_file = orthologs_file = pfam_file = None
    _output_files = None
    
    output_dir = scratch_dir = None
    _current_dir = None

    gene_pred = mode = no_hits_recovery = annot = None
    override = resume = None

    searcher = predictor = None

    ##
    def __init__(self, gene_pred, mode, annot, report_orthologs, prefix, output_dir, scratch_dir, resume, override):

        #
        self.output_dir = output_dir
        self.scratch_dir = scratch_dir
        
        # Output and intermediate files
        self.gene_pred_prots_file = f"{prefix}.emapper.gene_pred.prots.faa"
        self.hmm_hits_file = f"{prefix}.emapper.hmm_hits"
        self.seed_orthologs_file = f"{prefix}.emapper.seed_orthologs"
        self.annot_file = f"{prefix}.emapper.annotations"
        self.orthologs_file = f"{prefix}.emapper.orthologs"
        self.pfam_file = f"{prefix}.emapper.pfam"

        self.gene_pred = gene_pred
        self.mode = mode
        self.no_hits_recovery = "none"
        self.annot = annot
        self.report_orthologs = report_orthologs
        self.resume = resume
        self.override = override

        self._output_files = []
        if gene_pred == True:
            self._output_files.append(self.gene_pred_prots_file)
            
        if mode == SEARCH_MODE_NO_SEARCH:
            self._output_files.extend([self.annot_file, self.pfam_file])
        elif not annot:
            self._output_files.extend([self.hmm_hits_file, self.seed_orthologs_file])
        else:
            self._output_files.extend([self.hmm_hits_file, self.seed_orthologs_file, self.annot_file, self.pfam_file])

        if report_orthologs == True:
            self._output_files.append(self.orthologs_file)

        # force user to decide what to do with existing files
        files_present = set([pexists(pjoin(self.output_dir, fname)) for fname in self._output_files])
        if True in files_present and not resume and not override:
            raise EmapperException("Output files detected in disk. Use --resume or --override to continue")

        if override:
            for outf in self._output_files:
                silent_rm(pjoin(self.output_dir, outf))

        # If using --scratch_dir, change working dir
        # (once finished move them again to output_dir)
        if scratch_dir:
            self._current_dir = scratch_dir
            
            if resume:
                for fname in self._output_files:
                    if pexists(fname):
                        print("   Copying input file %s to scratch dir %s" % (fname, scratch_dir))
                        shutil.copy(fname, scratch_dir)
            
        else:
            self._current_dir = output_dir
            
        return

    ##
    def search(self, args, infile):
        annot = None
        
        s = get_searcher(args, self.mode)
        if s:
            annot = s.search(infile,
                             pjoin(self._current_dir, self.seed_orthologs_file),
                             pjoin(self._current_dir, self.hmm_hits_file))
        return s

    ##
    def annotate(self, args, annotate_hits_table):
        hits_file = None
        if annotate_hits_table:
            if not pexists(annotate_hits_table):
                raise EmapperException("Could not find hits table to annotate: %s" % (annotate_hits_table))
            
            hits_file = annotate_hits_table
        else:
            hits_file = pjoin(self._current_dir, self.seed_orthologs_file)

        if hits_file:
            a = get_annotator(args, self.annot, self.report_orthologs)
            a.annotate(hits_file,
                       pjoin(self._current_dir, self.annot_file),
                       pjoin(self._current_dir, self.orthologs_file),
                       pjoin(self._current_dir, self.pfam_file))
                
        return


    ##
    def gene_prediction(self, args, infile):
        predictor = ProdigalPredictor(args)
        predictor.predict(infile)
        return predictor

    ##
    def recover_no_hits(self):
        final_prots_file = None
        if self.no_hits_recovery == "none":
            pass
        elif self.no_hits_recovery == "alt_orfs":
            if self.gene_pred == False:
                raise EmapperException(f"Hits recovery mode {no_hits_recovery} requires performing the gene prediction step.")
            else:
                hits, no_hits = self.searcher.get_hits()
                alt_orfs = self.predictor.find_alt_orfs(no_hits)
                alt_orfs_fasta = utils.get_fasta(alt_orfs, infile)
                self.searcher.search(args, alt_orfs_fasta)
                hits_2, no_hits_2 = self.searcher.get_hits()
                final_prots = hits + no_hits + hits_2
                final_prots_file = utils.get_fasta(final_prots, infile)

        return final_prots_file
    
    ##
    def run(self, args, infile, annotate_hits_table = None):

        ##
        # Step 0. Gene prediction
        queries_file = None
        if self.gene_pred == True:
            self.predictor = self.gene_prediction(args, infile)
            queries_file = self.predictor.outprots # Use predicted proteins as input for search
            args.translate = False
        else:
            queries_file = infile # Use user input for search
        
        ##
        # Step 1. Sequence search
        if self.mode == SEARCH_MODE_NO_SEARCH:
            # Final file of predicted proteins
            if self.gene_pred == True and queries_file is not None:
                shutil.move(queries_file, pjoin(self._current_dir, self.gene_pred_prots_file))
                
        else: #if self.mode != SEARCH_MODE_NO_SEARCH:
            self.searcher = self.search(args, queries_file)
            # if annot is not None:
            #     self.annot = annot

            ## Step 1.2
            # Recovery of queries without hits
            if self.no_hits_recovery == "none":
                if self.gene_pred == True:
                    # Final file of predicted proteins
                    if queries_file is not None:
                        shutil.move(queries_file, pjoin(self._current_dir, self.gene_pred_prots_file))                
            else:
                final_prots_file = self.recover_no_hits()
                
                # Final file of predicted proteins                
                if final_prots_file is not None:
                    shutil.move(final_prots_file, pjoin(self._current_dir, self.gene_pred_prots_file))

            if self.predictor is not None:
                self.predictor.clear()
            
        ##
        # Step 2. Annotation
        if self.annot == True or self.report_orthologs:
            self.annotate(args, annotate_hits_table)

        ##
        # If running in scratch, move files to real output dir and clean up
        if self.scratch_dir:
            for fname in self._output_files:
                pathname = pjoin(self.scratch_dir, fname)
                if pexists(pathname):
                    print(" Copying result file %s from scratch to %s" % (pathname, self.output_dir))
                    shutil.copy(pathname, self.output_dir)

            print(colorify(f"Data in {self.scratch_dir} will be not removed. Please, clear it manually.", 'red'))

        ##
        # Finalize and exit
        print(colorify('Done', 'green'))
        for fname in self._output_files:
            pathname = pjoin(self.output_dir, fname)            
            colorify('Result files:', 'yellow')
            if pexists(pathname):
                print("   %s" % (pathname))
        
        return

## END
