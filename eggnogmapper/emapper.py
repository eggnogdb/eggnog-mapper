##
## CPCantalapiedra 2020

import errno, os, shutil, psutil
from sys import stderr
import time
from os.path import exists as pexists
from os.path import join as pjoin

from .utils import colorify
from .common import silent_rm, sort_seeds_file, \
    ITYPE_GENOME, ITYPE_META, ITYPE_PROTS, ITYPE_CDS, get_data_path
from .emapperException import EmapperException

from .genepred.genepred_modes import GENEPRED_MODE_SEARCH, GENEPRED_MODE_PRODIGAL, get_predictor
from .genepred.util import create_prots_file
from .search.search_modes import get_searcher, SEARCH_MODE_NO_SEARCH
from .search.hits_io import parse_seeds
from .annotation.annotator import Annotator
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
        # Phase 7.1c: optional drop log written when --report_dropped is set.
        self.dropped_file = f"{prefix}.emapper.dropped"
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
        elif not annot:
            self._output_files.extend([self.search_out_file, self.seed_orthologs_file])
        else:
            self._output_files.extend([self.search_out_file, self.seed_orthologs_file,
                                       self.annot_file, self.pfam_file])

        if annot == True and excel == True:
            self._output_files.append(self.excel_file)

        if report_orthologs == True:
            self._output_files.append(self.orthologs_file)

        # Phase 7.1c: dropped log is written only when --report_dropped is
        # set. Register unconditionally so --override / --resume handle
        # any leftover from a previous run cleanly (silent_rm is no-op
        # when the file is absent).
        self._output_files.append(self.dropped_file)

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

        predictor = get_predictor(args, self.genepred)
        # --resume skips gene prediction to resume from diamond/mmseqs/hmmer hits directly
        if predictor is not None and self.resume == False:
            predictor.predict(infile)
        elif predictor is not None:
            print(colorify("[resume] Skipping gene prediction — reusing existing files.", "blue"),
                  file=stderr)
        
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
    def annotate(self, args, hits, annotate_hits_table, queries_file):
        annotated_hits = None
        _resumed_ref = [0]

        if self.annot == True or self.report_orthologs:
            annot_in = None  # a generator of hits to annotate

            # --sort_entries: sort a file-based seed input by seed ortholog id
            # so shared seeds form contiguous blocks and each unique seed is
            # annotated once (see AnnotationEngine.annotate_batch dedup). Only
            # applies to file inputs; the in-memory blastx path is left as-is.
            def _sorted_seed_input(path):
                if not getattr(args, "sort_entries", False):
                    return path
                sorted_path = path + ".sorted"
                # On --resume, reuse an existing sorted file instead of
                # re-sorting the whole (potentially huge) input. It is written
                # atomically, so its presence means a previous sort completed;
                # require it to be at least as new as the source so a
                # regenerated seed file invalidates a stale sort.
                if (self.resume
                        and pexists(sorted_path)
                        and os.path.getmtime(sorted_path) >= os.path.getmtime(path)):
                    print(colorify(
                        f"[resume] Reusing sorted seed orthologs: {sorted_path}",
                        "blue"), file=stderr)
                    return sorted_path
                print(colorify(
                    f"[--sort_entries] Sorting {path} by seed ortholog id...",
                    "blue"), file=stderr)
                sort_seeds_file(path, sorted_path,
                                temp_dir=args.temp_dir, parallel=args.cpu)
                return sorted_path

            if annotate_hits_table is not None:
                if not pexists(annotate_hits_table):
                    raise EmapperException(
                        f"Could not find the file with the hits "
                        f"table to annotate: {annotate_hits_table}"
                    )
                annot_in = parse_seeds(_sorted_seed_input(annotate_hits_table))
            elif hits is not None:
                annot_in = hits
            else:
                # Search wrote seed_orthologs eagerly; read from the
                # completed file so annotation is fully decoupled from search.
                seed_file = pjoin(self._current_dir, self.seed_orthologs_file)
                if not pexists(seed_file):
                    raise EmapperException(
                        f"No hits to annotate and seed orthologs file not found: {seed_file}"
                    )
                annot_in = parse_seeds(_sorted_seed_input(seed_file))

            annotator = Annotator(args, self.annot, self.excel, self.report_orthologs)

            if annot_in is not None and annotator is not None:
                annotated_hits, _resumed_ref = annotator.annotate(
                    annot_in,
                    pjoin(self._current_dir, self.annot_file),
                    pjoin(self._current_dir, self.excel_file),
                    pjoin(self._current_dir, self.orthologs_file),
                    pjoin(self._current_dir, self.pfam_file),
                    queries_file,
                    dropped_file=pjoin(self._current_dir, self.dropped_file),
                )
        else:
            # --no_annot: wrap hits as (hit, None) pairs, or yield nothing if
            # no search was run (e.g. -m no_search --no_annot is a no-op).
            annotated_hits = ((hit, None) for hit in hits) if hits is not None else iter([])

        return annotated_hits, _resumed_ref

    
    ##
    def decorate_gff_f(self, args, predictor, searcher_name, annotated_hits):

        gff_outfile = pjoin(self._current_dir, self.deco_gff_file)
        if predictor is not None:
            if self.resume == True:
                gff_genepred_file = pjoin(self._current_dir, self.genepred_gff_file)
                gff_genepred_fasta = pjoin(self._current_dir, self.genepred_fasta_file)
            else:
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
    def _print_progress(self, n, n_new, start_time, proc_start, mem_monitor):
        now = time.time()
        total_time = now - start_time
        percen_mem = psutil.virtual_memory().percent
        percen_avail = psutil.virtual_memory().available * 100 / psutil.virtual_memory().total

        if total_time > 0.005:
            # Measure the rate from the first produced result (proc_start), not
            # from process start, so the one-time setup (DB load, worker fork,
            # --sort_entries sort) is excluded and the number reflects actual
            # annotation throughput instead of being diluted by startup.
            proc_time = (now - proc_start) if proc_start is not None else total_time
            rate = (float(n_new) / proc_time) if proc_time > 0 else 0.0
            n_resumed = n - n_new
            if n_resumed:
                msg = f"{n} ({n_new} new, {n_resumed} resumed) {total_time:.1f}s {rate:.2f} q/s"
            else:
                msg = f"{n} {total_time:.1f}s {rate:.2f} q/s"
            if mem_monitor:
                msg += f" (mem {percen_mem:.1f}% used, {percen_avail:.1f}% avail)"
            print(msg, file=stderr)

        stderr.flush()

        return total_time

    ##
    def run_generator(self, generator, resumed_ref=None, CHUNK_SIZE=500, mem_monitor=True):

        n = 0
        start_time = time.time()
        proc_start = None

        _n_skipped = resumed_ref if resumed_ref is not None else [0]

        self._print_progress(0, 0, start_time, proc_start, mem_monitor)

        for item in generator:
            if proc_start is None:
                # First result: all one-time setup (DB load, worker fork,
                # --sort_entries sort) is done. Anchor the throughput clock here
                # so the reported q/s reflects annotation speed, not startup.
                proc_start = time.time()
                print(colorify(
                    f"[annotation] setup completed in {proc_start - start_time:.1f}s; "
                    "q/s below is measured from here (excludes setup).",
                    "blue"), file=stderr)
            n += 1
            if n % CHUNK_SIZE == 0:
                n_new = n - _n_skipped[0]
                self._print_progress(n, n_new, start_time, proc_start, mem_monitor)

        n_new = n - _n_skipped[0]
        total_time = self._print_progress(n, n_new, start_time, proc_start, mem_monitor)

        return n, total_time

    ##
    def wrap_up(self, predictor, searcher):
        ##
        # Clear things
        if predictor is not None:
            if self.resume == False: # if self.resume == True there aren't new genepred files to move
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
    def run(self, args, infile, annotate_hits_table=None):

        ##
        # Compressed input (gzip/bzip2, detected by magic bytes) is NOT
        # eagerly inflated to a temp file here. Each consumer handles it at
        # the point of use: the FASTA iterators (iter_fasta_seqs/gopen) read
        # gz+bz2 transparently, and the search tools stream gzip natively
        # (decompressing only bzip2). This avoids materialising the whole
        # decompressed input on disk for large gzipped FASTAs.

        ##
        # Step 0. Gene prediction
        predictor = self.gene_prediction(args, infile)

        ##
        # Step 1. Sequence search
        searcher, searcher_name, hits, queries_file = self.search(args, infile, predictor)

        ##
        # Step 2. Annotation
        annotated_hits, _resumed_ref = self.annotate(args, hits, annotate_hits_table, queries_file)

        ##
        # Step 3. Decorate GFF
        annotated_hits = self.decorate_gff_f(args, predictor, searcher_name, annotated_hits)

        ##
        # Run the generators
        n, elapsed_time = self.run_generator(annotated_hits, _resumed_ref)

        ##
        # Finish
        self.wrap_up(predictor, searcher)

        return n, elapsed_time

## END
