##
## CPCantalapiedra 2020

import errno, os, shutil, psutil
from sys import stderr
import time
from os.path import exists as pexists
from os.path import join as pjoin

from .utils import colorify
from .common import silent_rm, ITYPE_GENOME, ITYPE_META, ITYPE_PROTS, ITYPE_CDS, get_data_path
from .emapperException import EmapperException

from .genepred.genepred_modes import GENEPRED_MODE_SEARCH, GENEPRED_MODE_PRODIGAL, get_predictor
from .genepred.util import create_prots_file
from .search.search_modes import get_searcher, SEARCH_MODE_NO_SEARCH, SEARCH_MODE_CACHE, SEARCH_MODE_NOVEL_FAMS
from .search.hits_io import parse_seeds
from .annotation.annotators import get_annotator, get_annotator_novel_fams, get_cache_annotator
from .deco.decoration import run_gff_decoration, DECORATE_GFF_NONE, create_blastx_hits_gff

class Emapper:

    genepred_fasta_file = genepred_gff_file = searcher_out_file = seed_orthologs_file = None
    annot_file = orthologs_file = pfam_file = None
    _output_files = None
    
    output_dir = scratch_dir = None
    _current_dir = None

    itype = genepred = mode = annot = decorate_gff = None
    genepred_is_prodigal = genepred_is_blastx = None
    override = resume = None

    ##
    def __init__(self, itype, genepred, mode, annot, excel, report_orthologs, decorate_gff,
                 prefix, output_dir, scratch_dir, resume, override):

        #
        self.output_dir = output_dir
        self.scratch_dir = scratch_dir
        
        # Output and intermediate files
        self.genepred_fasta_file = f"{prefix}.emapper.genepred.fasta"
        self.genepred_gff_file = f"{prefix}.emapper.genepred.gff"

        self.search_out_file = f"{prefix}.emapper.hits"
        
        self.seed_orthologs_file = f"{prefix}.emapper.seed_orthologs"
        self.annot_file = f"{prefix}.emapper.annotations"
        self.no_annot_file = f"{prefix}.emapper.no_annotations.fasta"
        self.orthologs_file = f"{prefix}.emapper.orthologs"
        self.pfam_file = f"{prefix}.emapper.pfam"
        self.excel_file = f"{prefix}.emapper.annotations.xlsx"
        self.deco_gff_file = f"{prefix}.emapper.decorated.gff"

        self.genepred = genepred
        self.mode = mode
        
        self.annot = annot
        self.excel = excel
        self.report_orthologs = report_orthologs
        self.decorate_gff = decorate_gff
        self.resume = resume
        self.override = override

        self._output_files = []
        self.itype = itype
        if itype == ITYPE_GENOME or itype == ITYPE_META:
            self._output_files.append(self.genepred_fasta_file)
            self._output_files.append(self.genepred_gff_file)            
            
        if self.decorate_gff != DECORATE_GFF_NONE:
            self._output_files.append(self.deco_gff_file)

        self.genepred_is_prodigal = ((itype == ITYPE_GENOME or itype == ITYPE_META) and
                                     genepred == GENEPRED_MODE_PRODIGAL)
        
        self.genepred_is_blastx = ((itype == ITYPE_GENOME or itype == ITYPE_META) and
                                   genepred == GENEPRED_MODE_SEARCH)
            
        if mode == SEARCH_MODE_NO_SEARCH:
            self._output_files.extend([self.annot_file, self.pfam_file])
        elif mode == SEARCH_MODE_CACHE:
            self._output_files.extend([self.annot_file, self.no_annot_file])            
        elif not annot:
            self._output_files.extend([self.search_out_file, self.seed_orthologs_file])
        else:
            self._output_files.extend([self.search_out_file, self.seed_orthologs_file,
                                       self.annot_file, self.pfam_file])

        if annot == True and excel == True:
            self._output_files.append(self.excel_file)
            
        if report_orthologs == True:
            self._output_files.append(self.orthologs_file)

        # force user to decide what to do with existing files
        files_present = set([pexists(pjoin(self.output_dir, fname)) for fname in self._output_files])
        if True in files_present and not resume and not override:
            raise EmapperException("Output files detected in disk. Use --resume or --override to continue")

        if override == True:
            for outf in self._output_files:
                silent_rm(pjoin(self.output_dir, outf))

        # Some files are not being resumed and will be ovewritten
        if resume == True:
            silent_rm(pjoin(self.output_dir, self.no_annot_file))
            silent_rm(pjoin(self.output_dir, self.excel_file))

        # If using --scratch_dir, change working dir
        # (once finished move them again to output_dir)
        if scratch_dir:
            self._current_dir = scratch_dir
            
            if resume:
                for fname in self._output_files:
                    full_fname = pjoin(self.output_dir, fname)
                    if pexists(full_fname):
                        print("   Copying input file %s to scratch dir %s" % (full_fname, scratch_dir))
                        shutil.copy(full_fname, scratch_dir)
            
        else:
            self._current_dir = output_dir
            
        return

    
    ##
    def gene_prediction(self, args, infile):
        predictor = None

        # --resume skips gene prediction to resume from diamond/mmseqs/hmmer hits directly
        if self.resume == False:
            predictor = get_predictor(args, self.genepred)
            if predictor is not None:
                predictor.predict(infile)
        
        return predictor

    
    ##
    def search(self, args, infile, predictor = None):

        queries_file = None
        
        # determine input file: from prediction or infile
        if predictor is None:
            queries_file = infile
        else:
            queries_file = predictor.outprots
            args.translate = False
            args.itype = ITYPE_PROTS

        # search
        searcher = get_searcher(args, self.mode, get_data_path())
        searcher_name = None
        hits = None
        if searcher is not None:
            searcher_name = searcher.name
            try:
                hits = searcher.search(queries_file,
                                       pjoin(self._current_dir, self.seed_orthologs_file),
                                       pjoin(self._current_dir, self.search_out_file))
            except Exception as e:
                searcher.clear()
                raise(e)

            # If gene prediction from the hits obtained in the search step
            # create GFF and FASTA of predicted CDS/proteins
            if self.genepred_is_blastx == True:
                
                # create GFF of predicted CDS
                print(colorify("Crafting GFF file of CDS ...", "lgreen"))
                gff_outfile = pjoin(self._current_dir, self.genepred_gff_file)
                hits = create_blastx_hits_gff(hits, gff_outfile, searcher_name, args.decorate_gff_ID_field)
                
                # create fasta file of predicted CDS
                print(colorify("Crafting fasta file of CDS ...", "lgreen"))
                fasta_file = pjoin(self._current_dir, self.genepred_fasta_file)
                silent_rm(fasta_file)
                hits = create_prots_file(queries_file, hits, fasta_file, args.translate, args.trans_table)
                queries_file = fasta_file
                if args.translate == True:
                    args.itype = ITYPE_PROTS
                else:
                    args.itype = ITYPE_CDS
                
        return searcher, searcher_name, hits, queries_file
    
    
    ##
    def annotate(self, args, hits, annotate_hits_table, queries_file, cache_file):
        annotated_hits = None
        
        if self.annot == True or self.report_orthologs:

            if cache_file is not None:
                if not pexists(cache_file):
                    raise EmaperException(f"Could not find cache file: {cache_file}")
                
                annotator = get_cache_annotator(args)
                
                if annotator is not None:
                    annotated_hits = annotator.annotate(cache_file,
                                                        pjoin(self._current_dir, self.annot_file),
                                                        pjoin(self._current_dir, self.no_annot_file))
            else:
                
                annot_in = None # a generator of hits to annotate
                
                if annotate_hits_table is not None:
                    if not pexists(annotate_hits_table):
                        raise EmapperException(f"Could not find the file with the hits "
                                               f"table to annotate: {annotate_hits_table}")
                                               
                    # function which parses the file and yields hits
                    annot_in = parse_seeds(annotate_hits_table)
                    
                elif hits is not None:
                    annot_in = hits
                    
                else:
                    raise EmapperException("Could not find hits to annotate.")

                if self.mode == SEARCH_MODE_NOVEL_FAMS:
                    annotator = get_annotator_novel_fams(args, self.annot, self.excel, self.report_orthologs)
                else:
                    annotator = get_annotator(args, self.annot, self.excel, self.report_orthologs)
                
                if annot_in is not None and annotator is not None:
                    annotated_hits = annotator.annotate(annot_in, 
                                                        pjoin(self._current_dir, self.annot_file),
                                                        pjoin(self._current_dir, self.excel_file),
                                                        pjoin(self._current_dir, self.orthologs_file),
                                                        pjoin(self._current_dir, self.pfam_file),
                                                        queries_file)
        else:
            annotated_hits = ((hit, None) for hit in hits) # hits generator without annotations
                
        return annotated_hits

    
    ##
    def decorate_gff_f(self, args, predictor, searcher_name, annotated_hits):

        gff_outfile = pjoin(self._current_dir, self.deco_gff_file)
        if predictor is not None:
            gff_genepred_file = predictor.outgff
            gff_genepred_fasta = predictor.outprots
        else:
            gff_genepred_file = None
            gff_genepred_fasta = None
        
        
        annotated_hits = run_gff_decoration(self.decorate_gff, args.decorate_gff_ID_field,
                                            self.genepred_is_prodigal, self.genepred_is_blastx,
                                            gff_genepred_file, gff_genepred_fasta, gff_outfile,
                                            predictor, searcher_name, annotated_hits)
        
        return annotated_hits


    ##
    def _print_progress(self, n, start_time, mem_monitor):
        total_time = time.time() - start_time
        percen_mem = psutil.virtual_memory().percent
        percen_avail = psutil.virtual_memory().available * 100 / psutil.virtual_memory().total
                
        msg = f"{n} {total_time} {(float(n) / total_time):.2f} q/s "
        if mem_monitor == True:
            msg += f"(% mem usage: {percen_mem:.2f}, % mem avail: {percen_avail:.2f})"
                    
        print(msg, file=stderr)
        stderr.flush()
        
        return total_time
    
    ##
    def run_generator(self, generator, CHUNK_SIZE = 500, mem_monitor = True):
        
        n = 0
        start_time = time.time()
        
        total_time = self._print_progress(n, start_time, mem_monitor)
        
        for item in generator:
            n += 1
            if n and (n % CHUNK_SIZE == 0):
                self._print_progress(n, start_time, mem_monitor)

        total_time = self._print_progress(n, start_time, mem_monitor)
        
        return n, total_time

    ##
    def wrap_up(self, predictor, searcher):
        ##
        # Clear things
        if predictor is not None:
            shutil.move(predictor.outprots, pjoin(self._current_dir, self.genepred_fasta_file))
            shutil.move(predictor.outgff, pjoin(self._current_dir, self.genepred_gff_file))
            predictor.clear() # removes gene predictor output directory
        if searcher is not None:
            searcher.clear()
        
        ##
        # If running in scratch, move files to real output dir and clean up
        if self.scratch_dir:
            for fname in self._output_files:
                full_fname = pjoin(self.scratch_dir, fname)
                if pexists(full_fname):
                    print(" Copying result file %s from scratch to %s" % (full_fname, self.output_dir))
                    shutil.copy(full_fname, self.output_dir)

            print(colorify(f"Data in {self.scratch_dir} will be not removed. Please, clear it manually.", 'red'))

        ##
        # Finalize and exit
        print(colorify('Done', 'green'))
        print(colorify('Result files:', 'yellow'))
        for fname in self._output_files:
            pathname = pjoin(self.output_dir, fname)
            if pexists(pathname):
                print("   %s" % (pathname))
                
        return
    
    ##
    def run(self, args, infile, annotate_hits_table = None, cache_dir = None):

        ##
        # Step 0. Gene prediction
        predictor = self.gene_prediction(args, infile)
        
        ##
        # Step 1. Sequence search
        searcher, searcher_name, hits, queries_file = self.search(args, infile, predictor)
            
        ##
        # Step 2. Annotation
        annotated_hits = self.annotate(args, hits, annotate_hits_table, queries_file, cache_dir)

        ##
        # step 3. Decorate GFF
        annotated_hits = self.decorate_gff_f(args, predictor, searcher_name, annotated_hits)

        ##
        # Run the generators
        n, elapsed_time = self.run_generator(annotated_hits)

        ##
        # Finish
        self.wrap_up(predictor, searcher)
        
        return n, elapsed_time

## END
