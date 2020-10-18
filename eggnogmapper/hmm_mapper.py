##
## CPCantalapiedra 2020

import shutil
from os.path import exists as pexists
from os.path import join as pjoin

from .utils import colorify
from .common import silent_rm
from .emapperException import EmapperException

from .search.hmmer.hmmer import HmmerSearcher

class HmmMapper:

    hmm_hits_file = None
    
    output_dir = scratch_dir = None
    _current_dir = None

    override = resume = None

    ##
    def __init__(self, prefix, output_dir, scratch_dir, resume, override):
        
        #
        self.output_dir = output_dir
        self.scratch_dir = scratch_dir
        
        # Output and intermediate files
        self.hmm_hits_file = f"{prefix}.emapper.hmm_hits"
        
        self.resume = resume
        self.override = override

        output_path = pjoin(self.output_dir, self.hmm_hits_file)

        if pexists(output_path) and not resume and not override:
            raise EmapperException("Output files detected in disk. Use --resume or --override to continue")

        if override:
            silent_rm(output_path)

        # If using --scratch_dir, change working dir
        # (once finished move them again to output_dir)
        if scratch_dir:
            self._current_dir = scratch_dir
            
            if resume:
                if pexists(output_path):
                    print("   Copying input file %s to scratch dir %s" % (output_path, scratch_dir))
                    shutil.copy(output_path, scratch_dir)
                    
        else:
            self._current_dir = output_dir
            
        return


    ##
    def search(self, args, infile):

        args.no_refine = True
        args.excluded_taxa = None
        
        s = HmmerSearcher(args)
        s.search_hmm_matches(infile, pjoin(self._current_dir, self.hmm_hits_file))
        
        return
    

    ##
    def run(self, args, infile):

        print("hmm_mapper: run")
        ##
        # Step 1. Sequence search
        self.search(args, infile)

        ##
        # If running in scratch, move files to real output dir and clean up
        
        if self.scratch_dir:
            scratch_path = pjoin(self.scratch_dir, self.hmm_hits_file)
            if pexists(scratch_path):
                print(" Copying result file %s from scratch to %s" % (scratch_path, self.output_dir))
                shutil.copy(scratch_path, self.output_dir)

            print(colorify(f"Data in {self.scratch_dir} will be not removed. Please, clear it manually.", 'red'))

        ##
        # Finalize and exit
        print(colorify('Done', 'green'))
        colorify('Result files:', 'yellow')
        output_path = pjoin(self.output_dir, self.hmm_hits_file)
        if pexists(output_path):
            print("   %s" % (output_path))
        
        return
    
## END
