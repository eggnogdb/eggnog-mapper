##
## CPCantalapiedra 2020

import errno, os, shutil

from .utils import colorify
from .common import silent_rm, pexists
from .emapperException import EmapperException

from .search.search_modes import get_searcher, SEARCH_MODE_NO_SEARCH
from .annotation.annotator import get_annotator
from .orthologs import orthology as orthology

class Emapper:

    hmm_hits_file = seed_orthologs_file = annot_file = orthologs_file = None
    _output_files = None
    
    output_dir = scratch_dir = None
    _current_dir = None

    mode = annot = None

    ##
    def __init__(self, mode, annot, prefix, output_dir, scratch_dir, override):

        #
        self.output_dir = output_dir
        self.scratch_dir = scratch_dir
        
        # Output and intermediate files
        
        # self.hmm_hits_file = "%s/%s.emapper.hmm_hits" % (output_dir, prefix)
        # self.seed_orthologs_file = "%s/%s.emapper.seed_orthologs" % (output_dir, prefix)
        # self.annot_file = "%s/%s.emapper.annotations" % (output_dir, prefix)
        # self.orthologs_file = "%s/%s.emapper.predict_orthologs" % (output_dir, prefix)
        self.hmm_hits_file = f"{prefix}.emapper.hmm_hits"
        self.seed_orthologs_file = f"{prefix}.emapper.seed_orthologs"
        self.annot_file = f"{prefix}.emapper.annotations"
        self.orthologs_file = f"{prefix}.emapper.predict_orthologs"

        self.mode = mode
        self.annot = annot
        if mode == SEARCH_MODE_NO_SEARCH:
            self._output_files = [self.annot_file]
        elif not annot:
            self._output_files = [self.hmm_hits_file, self.seed_orthologs_file]
        else:
            self._output_files = [self.hmm_hits_file, self.seed_orthologs_file, self.annot_file]
            
        # self._output_files = [self.hmm_hits_file, self.seed_orthologs_file, self.annot_file]

        # force user to decide what to do with existing files
        files_present = set([pexists(os.path.join(self.output_dir, fname)) for fname in self._output_files])
        if True in files_present and not override:
            raise EmapperException("Output files detected in disk. Use --override to continue")

        if override:
            for outf in self._output_files:
                silent_rm(os.path.join(self.output_dir, outf))

        # If using --scratch_dir, change working dir
        # (once finished move them again to output_dir)
        if scratch_dir:
            self._current_dir = scratch_dir
            # self.hmm_hits_file = "%s/%s.emapper.hmm_hits" % (scratch_dir, prefix)
            # self.seed_orthologs_file = "%s/%s.emapper.seed_orthologs" % (scratch_dir, prefix)
            # self.annot_file = "%s/%s.emapper.annotations" % (scratch_dir, prefix)
            # self.orthologs_file = "%s/%s.emapper.predict_orthologs" % (scratch_dir, prefix)            
        else:
            self._current_dir = output_dir
            
        return


    ##
    def search(self, args, infile):
        s = get_searcher(args, self.mode)
        if s:
            s.search(infile, os.path.join(self._current_dir, self.seed_orthologs_file))
        return

    ##
    def annotate(self, args, annotate_hits_table):
        hits_file = None
        if annotate_hits_table:
            if not os.path.exists(annotate_hits_table):
                raise EmapperException("Could not find hits table to annotate: %s" % (annotate_hits_table))
            
            hits_file = annotate_hits_table
        else:
            hits_file = os.path.join(self._current_dir, self.seed_orthologs_file)

        if hits_file:
            a = get_annotator(args)
            a.annotate(hits_file,
                       os.path.join(self._current_dir, self.annot_file),
                       os.path.join(self._current_dir, self.hmm_hits_file))
                
        return
    
    ##
    def run(self, args, infile, annotate_hits_table = None, predict_ortho = None):

        print("emapper: run")
        ##
        # Step 1. Sequence search
        if self.mode != SEARCH_MODE_NO_SEARCH:
            self.search(args, infile)

        ##
        # Step 2. Annotation
        if self.annot == True:
            self.annotate(args, annotate_hits_table)

        ##
        # Optional step. Orthology prediction
        if predict_ortho:
            orthology.connect()
            dump_orthologs(os.path.join(self._current_dir, self.seed_orthologs_file),
                           os.path.join(self._current_dir, self.orthologs_file, args))

        ##
        # If running in scratch, move files to real output dir and clean up
        if self.scratch_dir:
            for fname in self._output_files:
                pathname = os.path.join(self.scratch_dir, fname)
                if pexists(pathname):
                    print(" Copying result file %s from scratch to %s" % (pathname, self.output_dir))
                    shutil.copy(pathname, self.output_dir) # CPC 2020 should be fname instead of annot_file
                    print("  Cleaning result file %s from scratch dir" %(fname)) # CPC 2020 it says is going to clean but does nothing here

        ##
        # Finalize and exit
        print(colorify('Done', 'green'))
        for fname in self._output_files:
            pathname = os.path.join(self.output_dir, fname)            
            colorify('Result files:', 'yellow')
            if pexists(pathname):
                print("   %s" % (pathname))
        
        return

## END
