##
## CPCantalapiedra 2020

import errno, os, shutil
from os.path import exists as pexists
from os.path import join as pjoin

from .utils import colorify
from .common import silent_rm
from .emapperException import EmapperException

from .search.search_modes import get_searcher, SEARCH_MODE_NO_SEARCH
from .annotation.annotator import get_annotator

class Emapper:

    hmm_hits_file = seed_orthologs_file = annot_file = None
    _output_files = None
    
    output_dir = scratch_dir = None
    _current_dir = None

    mode = annot = None
    override = resume = None

    ##
    def __init__(self, mode, annot, report_orthologs, prefix, output_dir, scratch_dir, resume, override):

        #
        self.output_dir = output_dir
        self.scratch_dir = scratch_dir
        
        # Output and intermediate files
        self.hmm_hits_file = f"{prefix}.emapper.hmm_hits"
        self.seed_orthologs_file = f"{prefix}.emapper.seed_orthologs"
        self.annot_file = f"{prefix}.emapper.annotations"
        self.orthologs_file = f"{prefix}.emapper.orthologs"
        self.pfam_file = f"{prefix}.emapper.pfam"

        self.mode = mode
        self.annot = annot
        self.report_orthologs = report_orthologs
        self.resume = resume
        self.override = override
        
        if mode == SEARCH_MODE_NO_SEARCH:
            self._output_files = [self.annot_file, self.pfam_file]
        elif not annot:
            self._output_files = [self.hmm_hits_file, self.seed_orthologs_file]
        else:
            self._output_files = [self.hmm_hits_file, self.seed_orthologs_file, self.annot_file, self.pfam_file]

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
        return annot

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
    def run(self, args, infile, annotate_hits_table = None):
        
        ##
        # Step 1. Sequence search
        if self.mode != SEARCH_MODE_NO_SEARCH:
            annot = self.search(args, infile)
            if annot is not None:
                self.annot = annot

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
