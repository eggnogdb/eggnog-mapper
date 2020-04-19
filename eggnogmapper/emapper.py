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

    ##
    def __init__(self, output_dir, scratch_dir, prefix, override):

        #
        self.output_dir = output_dir
        self.scratch_dir = scratch_dir
        
        # Output and intermediate files
        self.hmm_hits_file = "%s/%s.emapper.hmm_hits" % (output_dir, prefix)
        self.seed_orthologs_file = "%s/%s.emapper.seed_orthologs" % (output_dir, prefix)
        self.annot_file = "%s/%s.emapper.annotations" % (output_dir, prefix)
        self.orthologs_file = "%s/%s.emapper.predict_orthologs" % (output_dir, prefix)

        self._output_files = [self.hmm_hits_file, self.seed_orthologs_file, self.annot_file]

        # force user to decide what to do with existing files
        files_present = set([pexists(fname) for fname in self._output_files])
        if True in files_present and not override:
            raise EmapperException("Output files detected in disk. Use --override to continue")

        if override:
            for outf in self._output_files:
                silent_rm(outf)

        # If using --scratch_dir, change working dir
        # (once finished move them again to output_dir)
        if scratch_dir:
            self.hmm_hits_file = "%s/%s.emapper.hmm_hits" % (scratch_dir, prefix)
            self.seed_orthologs_file = "%s/%s.emapper.seed_orthologs" % (scratch_dir, prefix)
            self.annot_file = "%s/%s.emapper.annotations" % (scratch_dir, prefix)
            self.orthologs_file = "%s/%s.emapper.predict_orthologs" % (scratch_dir, prefix)

        return


    ##
    def search(self, args, mode, infile):
        s = get_searcher(args, mode)
        if s:
            s.search(infile, self.seed_orthologs_file)
        return

    ##
    def annotate(self, args, annotate_hits_table):
        hits_file = None
        if annotate_hits_table:
            if not os.path.exists(annotate_hits_table):
                raise EmapperException("Could not find hits table to annotate: %s" % (annotate_hits_table))
                # raise IOError(errno.ENOENT,
                #               os.strerror(errno.ENOENT),
                #               annotate_hits_table)
            hits_file = annotate_hits_table
        else:
            hits_file = self.seed_orthologs_file

        if hits_file:
            a = get_annotator(args)
            a.annotate(hits_file, self.annot_file, self.hmm_hits_file)
                
        return
    
    ##
    def run(self, args, mode, infile, annot = True, annotate_hits_table = None, predict_ortho = None):

        print("emapper: run")
        ##
        # Step 1. Sequence search
        if mode != SEARCH_MODE_NO_SEARCH:
            self.search(args, mode, infile)

        ##
        # Step 2. Annotation
        if annot == True:
            self.annotate(args, annotate_hits_table)

        ##
        # Optional step. Orthology prediction
        if predict_ortho:
            orthology.connect()
            dump_orthologs(self.seed_orthologs_file, self.orthologs_file, args)

        ##
        # If running in scratch, move files to real output dir and clean up
        if args.scratch_dir:
            for fname in self._output_files:
                if pexists(fname):
                    print(" Copying result file %s from scratch to %s" % (fname, args.output_dir))
                    shutil.copy(annot_file, args.output_dir) # CPC 2020 should be fname instead of annot_file
                    print("  Cleaning result file %s from scratch dir" %(fname)) # CPC 2020 it says is going to clean but does nothing here

        ##
        # Finalize and exit
        print(colorify('Done', 'green'))
        for f in self._output_files:
            colorify('Result files:', 'yellow')
            if pexists(f):
                print("   %s" % (f))
        
        return

## END
