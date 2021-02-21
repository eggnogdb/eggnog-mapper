##
## CPCantalapiedra 2020

import errno, os, shutil
from os.path import exists as pexists
from os.path import join as pjoin

from .utils import colorify
from .common import silent_rm, ITYPE_GENOME, ITYPE_META, ITYPE_PROTS
from .emapperException import EmapperException

from .genepred.genepred_modes import GENEPRED_MODE_SEARCH, GENEPRED_MODE_PRODIGAL, get_predictor
from .genepred.util import create_prots_file
from .search.search_modes import get_searcher, SEARCH_MODE_NO_SEARCH, SEARCH_MODE_CACHE
from .annotation.annotators import get_annotator, get_cache_annotator
from .deco.decoration import run_gff_decoration, DECORATE_GFF_NONE

class Emapper:

    genepred_fasta_file = genepred_gff_file = hmm_hits_file = seed_orthologs_file = None
    annot_file = orthologs_file = pfam_file = None
    _output_files = None
    
    output_dir = scratch_dir = None
    _current_dir = None

    itype = genepred = mode = annot = decorate_gff = None
    genepred_is_prodigal = genepred_is_blastx = None
    override = resume = None

    ##
    def __init__(self, itype, genepred, mode, annot, report_orthologs, decorate_gff,
                 prefix, output_dir, scratch_dir, resume, override):

        #
        self.output_dir = output_dir
        self.scratch_dir = scratch_dir
        
        # Output and intermediate files
        self.genepred_fasta_file = f"{prefix}.emapper.genepred.fasta"
        self.genepred_gff_file = f"{prefix}.emapper.genepred.gff"
        self.hmm_hits_file = f"{prefix}.emapper.hmm_hits"
        self.seed_orthologs_file = f"{prefix}.emapper.seed_orthologs"
        self.annot_file = f"{prefix}.emapper.annotations"
        self.no_annot_file = f"{prefix}.emapper.no_annotations.fasta"
        self.orthologs_file = f"{prefix}.emapper.orthologs"
        self.pfam_file = f"{prefix}.emapper.pfam"

        self.genepred = genepred
        self.mode = mode
        
        self.annot = annot
        self.report_orthologs = report_orthologs
        self.decorate_gff = decorate_gff
        self.resume = resume
        self.override = override

        self._output_files = []
        self.itype = itype
        if itype == ITYPE_GENOME or itype == ITYPE_META:
            self._output_files.append(self.genepred_fasta_file)
            
        if itype == ITYPE_GENOME or itype == ITYPE_META or self.decorate_gff != DECORATE_GFF_NONE:
            self._output_files.append(self.genepred_gff_file)

        self.genepred_is_prodigal = (itype == ITYPE_GENOME or itype == ITYPE_META) and genepred == GENEPRED_MODE_PRODIGAL
        self.genepred_is_blastx = (itype == ITYPE_GENOME or itype == ITYPE_META) and genepred == GENEPRED_MODE_SEARCH
            
        if mode == SEARCH_MODE_NO_SEARCH:
            self._output_files.extend([self.annot_file, self.pfam_file])
        elif mode == SEARCH_MODE_CACHE:
            self._output_files.extend([self.annot_file, self.no_annot_file])            
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
    def gene_prediction(self, args, infile):
        queries_file = None
        predictor = get_predictor(args, self.genepred)
        if predictor is not None:
            predictor.predict(infile)
        
        return predictor

    
    ##
    def search(self, args, infile, predictor = None):

        # determine input file: from prediction or infile
        if predictor is None:
            queries_file = infile
        else:
            queries_file = predictor.outprots
            args.translate = False
            args.itype = ITYPE_PROTS

        # search
        searcher = get_searcher(args, self.mode)
        if searcher is not None:
            searcher.search(queries_file,
                            pjoin(self._current_dir, self.seed_orthologs_file),
                            pjoin(self._current_dir, self.hmm_hits_file))

            # If gene prediction from the hits obtained in the search step
            # create a fasta file with the inferred proteins
            if self.genepred_is_blastx == True:
                hits = searcher.get_hits()
                fasta_file = pjoin(self._current_dir, self.genepred_fasta_file)
                create_prots_file(queries_file, hits, fasta_file)
            
        return searcher
    
    
    ##
    def annotate(self, args, searcher, annotate_hits_table, cache_file):
        annotator = None
        if self.annot == True or self.report_orthologs:
            hits_file = None

            if cache_file is not None:
                if not pexists(cache_file):
                    raise EmaperException(f"Could not find cache file: {cache_file}")
                annotator = get_cache_annotator(args)
                if annotator is not None:
                    annotator.annotate(cache_file,
                                       pjoin(self._current_dir, self.annot_file),
                                       pjoin(self._current_dir, self.no_annot_file))
            else:
                annotator = get_annotator(args, self.annot, self.report_orthologs)
                
                annot_in = None # a generator of hits to annotate
                store_hits = False # force or not the annotator to store the hits as annotator.hits
                if annotate_hits_table is not None:
                    if not pexists(annotate_hits_table):
                        raise EmapperException("Could not find the file with the hits table to annotate: %s" % (annotate_hits_table))
                    annot_in = annotator.parse_hits(annotate_hits_table) # function which parses the file and yields hits
                    store_hits = True
                elif searcher is not None:
                    annot_in = annotator.get_hits # function which returns annotator.hits
                    annotator.hits = searcher.get_hits()
                    store_hits = False
                else:
                    raise EmapperException("Could not find hits to annotate.")

                if annot_in is not None and annotator is not None:
                    annotator.annotate(annot_in, store_hits, 
                                       pjoin(self._current_dir, self.annot_file),
                                       pjoin(self._current_dir, self.orthologs_file),
                                       pjoin(self._current_dir, self.pfam_file))
                
        return annotator

    
    ##
    def decorate_gff_f(self, predictor, searcher, annotator):

        gff_outfile = pjoin(self._current_dir, self.genepred_gff_file)
        
        run_gff_decoration(self.decorate_gff, self.genepred_is_prodigal, self.genepred_is_blastx, gff_outfile,
                     predictor, searcher, annotator)
        
        return

    
    ##
    def run(self, args, infile, annotate_hits_table = None, cache_dir = None):

        ##
        # Step 0. Gene prediction
        predictor = self.gene_prediction(args, infile)
        
        ##
        # Step 1. Sequence search
        searcher = self.search(args, infile, predictor)
            
        ##
        # Step 2. Annotation
        annotator = self.annotate(args, searcher, annotate_hits_table, cache_dir)

        ##
        # Decorate GFF
        self.decorate_gff_f(predictor, searcher, annotator)

        ##
        # Clear things
        if predictor is not None:
            shutil.move(predictor.outprots, pjoin(self._current_dir, self.genepred_fasta_file))
            predictor.clear() # removes gene predictor output directory
        
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
