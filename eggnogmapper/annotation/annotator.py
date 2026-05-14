##
## CPCantalapiedra 2019

from os.path import isfile as pisfile
import sys
import time
from collections import defaultdict, Counter

from ..emapperException import EmapperException
from ..common import get_call_info, get_data_path, ITYPE_PROTS, ITYPE_CDS
from ..utils import colorify
from ..search.hmmer.hmmer_seqio import iter_fasta_seqs

from .db_sqlite import get_eggnog_db
from .ncbitaxa.ncbiquery import get_ncbi
from .pfam.pfam_modes import run_pfam_mode, PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO

from . import output

ANNOTATIONS_HEADER = output.ANNOTATIONS_HEADER

PFAM_COL = -1 # position of PFAMs annotations in list of annotations

##
class Annotator:
    
    annot = report_orthologs = None

    no_file_comments = cpu = None

    # options for pfam hmmpgmd searches
    num_servers = num_workers = timeout_load_server = cpus_per_worker = port = end_port = None

    seed_ortholog_score = seed_ortholog_evalue = None
    target_taxa = target_orthologs = excluded_taxa = None

    go_evidence = go_excluded = None
    pfam_realign = trans_table = temp_dir = None
    md5 = None

    resume = None

    ##
    def __init__(self, args, annot, excel, report_orthologs):
        """Initialise Annotator from parsed CLI args.

        Args:
            args: ``argparse.Namespace`` produced by ``emapper.py``.
            annot: Whether functional annotation is requested.
            excel: Whether Excel output is requested.
            report_orthologs: Whether the orthologs file is requested.
        """
        self.annot = annot
        self.report_orthologs = report_orthologs
        self.excel = excel

        self.no_file_comments = args.no_file_comments
        self.cpu = args.cpu
        self.num_servers = args.num_servers
        self.num_workers = args.num_workers
        self.timeout_load_server = args.timeout_load_server
        self.cpus_per_worker = args.cpus_per_worker
        self.port = args.port
        self.end_port = args.end_port
        self.seed_ortholog_score = args.seed_ortholog_score
        self.seed_ortholog_evalue = args.seed_ortholog_evalue

        # --tax_scope: "auto-narrow" | "auto-broad" | <clade_name_or_taxid>
        self.tax_scope = args.tax_scope

        # --donor_pool: "closest" | "union"
        self.donor_pool = getattr(args, "donor_pool", "closest") or "closest"

        self.target_taxa = args.target_taxa
        self.target_orthologs = args.target_orthologs
        self.excluded_taxa = args.excluded_taxa

        self.go_evidence = args.go_evidence
        self.go_excluded = args.go_excluded

        self.pfam_realign = args.pfam_realign

        self.trans_table = args.trans_table
        self.itype = args.itype

        self.temp_dir = args.temp_dir

        self.md5 = args.md5

        # Phase 7.1c: opt-in drop log. When True, Annotator opens
        # dropped_file and writes one row per filtered hit.
        self.report_dropped = bool(getattr(args, "report_dropped", False))

        self.resume = args.resume

        return

    def _applied_filters(self) -> dict:
        """Return a dict of the resolved annotation-stage filter values.

        Written to the ``## applied filters:`` block in the
        ``.emapper.annotations`` file.  Defaults are recorded as their
        actual resolved values so a reader does not need to know
        emapper's defaults to interpret the file.

        Returns:
            Mapping of filter name → resolved value.
        """
        return {
            "tax_scope": self.tax_scope,
            "donor_pool": self.donor_pool,
            "target_orthologs": self.target_orthologs,
            "target_taxa": self.target_taxa,
            "excluded_taxa": self.excluded_taxa,
            "seed_ortholog_evalue": self.seed_ortholog_evalue,
            "seed_ortholog_score": self.seed_ortholog_score,
            "go_evidence": self.go_evidence,
            "go_excluded": self.go_excluded,
            "pfam_realign": self.pfam_realign,
        }


    ##
    def annotate(self, hits_gen_func, annot_file, excel_file, orthologs_file,
                 pfam_file, queries_file, dropped_file=None):

        annots_generator = None
        ncbi = None
        # Phase 7.1c: opt-in drop log. Opened once here and threaded into
        # `_annotate_batched`; the file is closed via try/finally so a
        # crash mid-run still leaves a valid (partial) .dropped file.
        self._dropped_handle = None
        self._dropped_writer = None
        if self.report_dropped and dropped_file is not None:
            self._dropped_handle = open(dropped_file, "w")
            self._dropped_handle.write(
                "#query\treason\tseed_ortholog\tevalue\tscore\n"
            )

            def _writer(query, reason, seed, evalue, score):
                self._dropped_handle.write(
                    f"{query}\t{reason}\t{seed}\t{evalue}\t{score}\n"
                )

            self._dropped_writer = _writer
        
        try:
            if self.report_orthologs == True or self.annot == True:            
                ##
                # md5 hashes
                if self.md5 == True:
                    print(colorify("Creating md5 hashes of input sequences", 'green'))
                    
                    if self.itype == ITYPE_PROTS:
                        translate = False
                    elif self.itype == ITYPE_CDS:
                        translate = True
                        
                    md5_queries = md5_seqs(queries_file, translate, self.trans_table)
                else:
                    md5_queries = None

                md5_field = (self.md5 == True and md5_queries is not None)

                ##
                # Prepare taxa restrictions
                # target_taxa are used to restrict the species from which to retrieve co-ortholog proteins
                # the opposite is excluded_taxa
                # In both cases, we need to normalize the list of taxa, if any
                
                if self.target_taxa is not None:
                    ncbi = get_ncbi(usemem = False)
                    self.target_taxa = normalize_target_taxa(self.target_taxa, ncbi)
                    
                if self.excluded_taxa is not None:
                    ncbi = get_ncbi(usemem = False)
                    self.excluded_taxa = normalize_target_taxa(self.excluded_taxa, ncbi)

                if ncbi is not None: ncbi.close() # close it, so that for orthologs we can load it into memory
                
                ##
                # Annotations

                # If resume, create generator of previous annotations
                annots_parser = None
                if self.resume == True:
                    annots_parser = parse_annotations(self.annot, annot_file,
                                                      self.report_orthologs, orthologs_file)
                
                # Obtain annotations
                annots_generator = self._annotate(hits_gen_func, annots_parser)

                ##
                # PFAM realign
                # Note that this needs all the annotations at once,
                # and therefore breaks the generators pipeline
                if (self.annot == True and
                    self.pfam_realign in [PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO] and
                    annots_generator is not None):

                    if self.itype == ITYPE_PROTS:
                        translate = False
                    elif self.itype == ITYPE_CDS:
                        translate = True
                        
                    annots_generator = run_pfam_mode(self.pfam_realign, annots_generator,
                                                     queries_file, self.resume,
                                                     translate, self.trans_table,
                                                     self.cpu, self.num_servers,
                                                     self.num_workers, self.timeout_load_server,
                                                     self.cpus_per_worker,
                                                     self.port, self.end_port,
                                                     self.temp_dir, pfam_file)
                
                ##
                # Output
                
                if self.annot == True:
                    annots_generator = output.output_annotations(annots_generator,
                                                                 annot_file,
                                                                 self.resume,
                                                                 self.no_file_comments,
                                                                 md5_field,
                                                                 md5_queries,
                                                                 applied_filters=self._applied_filters())

                    if self.excel == True:
                        annots_generator = output.output_excel(annots_generator,
                                                               excel_file,
                                                               self.resume,
                                                               self.no_file_comments,
                                                               md5_field,
                                                               md5_queries)
                        
                    
                if self.report_orthologs == True:
                    try:
                        _eggnog_db = get_eggnog_db()
                    except Exception:
                        _eggnog_db = None
                    annots_generator = output.output_orthologs(annots_generator,
                                                               orthologs_file,
                                                               self.resume,
                                                               self.no_file_comments,
                                                               eggnog_db=_eggnog_db)

                # unpack the annotations removing the "exists" or "skip"
                # boolean used when --resume

                annots_generator = unpack_annotations(annots_generator)

        finally:
            if ncbi is not None: ncbi.close()

        # Phase 7.1c: close the .dropped file when the lazy `annots_generator`
        # chain is fully consumed (the consumer drives iteration, not us).
        if self._dropped_handle is not None:
            annots_generator = self._close_dropped_after(annots_generator)

        return annots_generator

    def _close_dropped_after(self, gen):
        """Wrap a generator so the drop-log file is closed only after the
        consumer finishes iterating. Without this the early `finally`
        above (kept for ncbi back-compat) would close the handle while
        downstream batches are still being processed."""
        handle = self._dropped_handle
        try:
            yield from gen
        finally:
            if handle is not None and not handle.closed:
                handle.close()


    def _annotate(self, hits_gen_func, annots_parser):
        """Dispatch annotation through the v7 batch path.

        v5 databases are rejected at AnnotDB.__init__ time
        (eggnogmapper.annotation.db_sqlite); legacy per-hit and dbmem
        paths were removed in Phase 2.
        """
        print(colorify("Functional annotation of hits...", "lgreen"), file=sys.stderr)
        eggnog_db = get_eggnog_db()
        return self._annotate_batched(hits_gen_func, annots_parser, eggnog_db)


    ##
    def _annotate_batched(self, hits_gen_func, annots_parser, eggnog_db):
        """Batch annotation: pre-fetch DB data per batch, CPU work in parallel.

        Only for v7+ integer-encoded databases.
        """
        from .batch_annotate import annotate_batch, get_lineage_cache, _get_engine
        from eggnogmapper.annotator.e7.ceiling import TaxScopeCeilingResolver
        from ..common import get_ncbitaxadb_file

        # Build the per-seed ceiling resolver once for the whole run.
        taxa_db = get_ncbitaxadb_file()
        lineage_cache = get_lineage_cache()
        ceiling_resolver = TaxScopeCeilingResolver.build(
            lineage_cache=lineage_cache,
            mode=self.tax_scope,
            taxa_db_path=taxa_db,
        )

        mode_display = self.tax_scope
        print(
            colorify(
                f"tax_scope ceiling mode: {mode_display!r} "
                f"(donor_pool={self.donor_pool!r})",
                "lblue",
            ),
            file=sys.stderr,
        )

        # Process-pool parallelism for the annotation phase. Each worker
        # inherits the parent's loaded engine (taxid_array,
        # lineage_cache) via fork-COW and reopens its own SQLite
        # connection in the post-fork initializer. Pool is created once
        # for the full mapper run and reused across batches.
        pool = None
        _orig_sigterm = _orig_sigint = None
        n_workers = max(int(getattr(self, "cpu", 1) or 1), 1)
        if n_workers > 1:
            import multiprocessing
            import signal as _signal
            from eggnogmapper.annotator.e7.annotate import (
                _register_worker_engine,
                _worker_init_after_fork,
            )

            engine = _get_engine(eggnog_db, ceiling_resolver, self.donor_pool)
            _register_worker_engine(engine)
            ctx = multiprocessing.get_context("fork")
            # Workers are long-lived (no maxtasksperchild) to avoid
            # os.fork() ENOMEM under memory pressure; cache growth is
            # instead bounded by _OG_CACHE_MAX / _LINEAGE_CACHE_MAX in
            # AnnotationEngine with half-eviction on overflow.
            pool = ctx.Pool(n_workers, initializer=_worker_init_after_fork)
            print(
                colorify(f"Annotation pool: {n_workers} workers (fork)", "lblue"),
                file=sys.stderr,
            )

            # Forward SIGTERM/SIGINT to workers so they exit promptly
            # instead of leaving pool.join() hanging after a
            # job-scheduler kill. Original handlers are restored before
            # re-raising so the process terminates normally.
            def _pool_signal_handler(signum, frame, _pool=pool):
                try:
                    if _orig_sigterm is not None:
                        _signal.signal(_signal.SIGTERM, _orig_sigterm)
                    if _orig_sigint is not None:
                        _signal.signal(_signal.SIGINT, _orig_sigint)
                    _pool.terminate()
                except Exception:
                    pass
                import os as _os
                _os.kill(_os.getpid(), signum)

            try:
                _orig_sigterm = _signal.signal(_signal.SIGTERM, _pool_signal_handler)
                _orig_sigint = _signal.signal(_signal.SIGINT, _pool_signal_handler)
            except (ValueError, OSError):
                # Not in main thread — skip signal forwarding.
                pass

        # Batch large enough to keep all workers continuously fed.
        # AnnotationEngine.annotate_batch uses sub_batch_size=125;
        # 4× headroom ensures imap_unordered pipelines without idle gaps
        # and avoids CPU underutilisation on many-core machines.
        batch_size = max(1000, n_workers * 500)
        batch = []

        try:
            for args_tuple in self.iter_hit_lines(hits_gen_func, annots_parser):
                # Pass through resumed hits immediately
                if args_tuple[-1] is not None:  # already annotated
                    yield ((args_tuple[0], args_tuple[-1]), True)
                    continue

                batch.append(args_tuple)

                if len(batch) >= batch_size:
                    yield from annotate_batch(
                        batch,
                        eggnog_db,
                        annot=self.annot,
                        target_orthologs=self.target_orthologs,
                        target_taxa=self.target_taxa,
                        excluded_taxa=self.excluded_taxa,
                        go_evidence=self.go_evidence,
                        go_excluded=self.go_excluded,
                        seed_ortholog_score=self.seed_ortholog_score,
                        seed_ortholog_evalue=self.seed_ortholog_evalue,
                        ceiling_resolver=ceiling_resolver,
                        donor_pool=self.donor_pool,
                        pool=pool,
                        dropped_writer=self._dropped_writer,
                    )
                    batch = []

            if batch:
                yield from annotate_batch(
                    batch,
                    eggnog_db,
                    annot=self.annot,
                    target_orthologs=self.target_orthologs,
                    target_taxa=self.target_taxa,
                    excluded_taxa=self.excluded_taxa,
                    go_evidence=self.go_evidence,
                    go_excluded=self.go_excluded,
                    seed_ortholog_score=self.seed_ortholog_score,
                    seed_ortholog_evalue=self.seed_ortholog_evalue,
                    ceiling_resolver=ceiling_resolver,
                    donor_pool=self.donor_pool,
                    pool=pool,
                    dropped_writer=self._dropped_writer,
                )

        except EmapperException:
            raise
        except Exception as e:
            import traceback

            traceback.print_exc()
            raise EmapperException(f"Error: batch annotation failed. " + str(e))
        finally:
            import signal as _sig
            import threading as _thr

            # Restore signal handlers before pool shutdown to prevent
            # re-entrant terminate() calls if a signal fires during join.
            for _snum, _orig in (
                (_sig.SIGTERM, _orig_sigterm),
                (_sig.SIGINT, _orig_sigint),
            ):
                if _orig is not None:
                    try:
                        _sig.signal(_snum, _orig)
                    except (ValueError, OSError):
                        pass

            if pool is not None:
                pool.close()
                # Join with a 120 s timeout. If workers are stuck (e.g.
                # blocked on NFS I/O or killed by the OOM killer),
                # force-terminate instead of hanging the parent forever.
                _t = _thr.Thread(target=pool.join, daemon=True)
                _t.start()
                _t.join(timeout=120)
                if _t.is_alive():
                    print(
                        colorify(
                            "WARNING: annotation workers did not exit within"
                            " 120 s; force-terminating.",
                            "red",
                        ),
                        file=sys.stderr,
                    )
                    pool.terminate()
                    pool.join()

        return

    
    ##
    def iter_hit_lines(self, hits_gen_func, annots_parser):
        """Iterate over hits, pairing each with any previously-annotated entry.

        For ``--resume`` runs, ``annots_parser`` provides previously
        computed annotations; matching hits are passed through without
        re-annotating.

        Args:
            hits_gen_func: Iterable of ``[query_name, seed, evalue, score, …]``
                hit tuples from the search phase.
            annots_parser: Generator from ``parse_annotations`` for
                ``--resume`` mode, or ``None``.

        Yields:
            Tuple ``(hit, annot, seed_ortholog_score, seed_ortholog_evalue,
            target_taxa, target_orthologs, excluded_taxa, go_evidence,
            go_excluded, data_path, annotation)`` where ``annotation`` is
            the previously computed annotation (or ``None``).
        """
        curr_annot = None
        if annots_parser is not None:
            try:
                curr_annot = next(annots_parser)
            except StopIteration:
                curr_annot = None

        for hit in hits_gen_func:

            annotation = None

            if curr_annot is not None:
                if hit[0] == curr_annot[0]:
                    annotation = curr_annot
                    try:
                        curr_annot = next(annots_parser)
                    except StopIteration:
                        curr_annot = None

            yield_tuple = (
                hit,
                self.annot,
                self.seed_ortholog_score,
                self.seed_ortholog_evalue,
                self.target_taxa,
                self.target_orthologs,
                self.excluded_taxa,
                self.go_evidence,
                self.go_excluded,
                get_data_path(),
                annotation,
            )

            yield yield_tuple

        return

##
def unpack_annotations(annots):
    for (hit, annotation), skip in annots:
        yield hit, annotation

##
def parse_annotations(annot, annot_file, report_orthologs, orthologs_file):
    
    if annot == True:
        with open(annot_file, 'r') as annot_f:
            for line in annot_f:
                if line.startswith("#"): continue
                
                hit, annotation = parse_annotation_line(line)
                
                # this assumes the annotated hit it is also present
                # in orthologs_file
                yield annotation
            
    elif report_orthologs == True:
        
        prev_query = None
        with open(orthologs_file, 'r') as orth_f:
            for line in orth_f:
                if line.startswith("#"): continue

                # just yield that query already exists in file
                query_name = line[0]
                if prev_query is not None and query_name != prev_query:
                    annotation = (prev_query, None, None, None,
                                  None, None, None, None, None)
                    yield annotation

                prev_query = query_name
                
            # last query
            if prev_query is not None:
                annotation = (prev_query, None, None, None,
                              None, None, None, None, None)
                yield annotation
                
    else:
        pass # no annotations then
    
    return

def parse_annotation_line(line):
    hit = None
    annotation = None

    data = list(map(str.strip, line.split("\t")))

    query_name = data[0]
    best_hit_name = data[1]
    best_hit_evalue = float(data[2])
    best_hit_score = float(data[3])
    hit = [query_name, best_hit_name, best_hit_evalue, best_hit_score]
    
    annotations = defaultdict(Counter)
    
    for i, field in enumerate(data[8:]):
        if i < len(ANNOTATIONS_HEADER):
            field_name = ANNOTATIONS_HEADER[i]
            annotations[field_name] = field.split(",")
        
    og_cat_desc = ("-", data[6], data[7]) # ("-", data[5], data[6])

    max_annot_lvl = data[5] # data[7]
    
    match_nog_names = data[4].split(",")
    
    all_orthologies = None
    annot_orthologs = None

    # Pad to the current 14-element shape so that --resume can replay
    # lines from pre-refactor annotation files through output_annotations_row
    # without hitting the 14-element `else` branch on a 10-element tuple.
    annotations_confidence = None
    tax_ceiling = "-"
    farthest_donor_taxid = "-"
    farthest_donor_lineage = "-"

    annotation = (query_name, best_hit_name, best_hit_evalue, best_hit_score,
                  annotations,
                  og_cat_desc, max_annot_lvl, match_nog_names,
                  all_orthologies, annot_orthologs,
                  annotations_confidence,
                  tax_ceiling,
                  farthest_donor_taxid,
                  farthest_donor_lineage)

    return hit, annotation


##
def normalize_target_taxa(target_taxa, ncbi):
    """
    Receives a list of taxa IDs and/or taxa names and returns a set of expanded taxids numbers
    """
    expanded_taxa = set()
    
    for taxon in target_taxa:
        taxid = ""
        try:
            taxid = int(taxon)
        except ValueError:
            taxid = ncbi.get_name_translator([taxon])[taxon][0]
        else:
            taxon = ncbi.get_taxid_translator([taxid])[taxid]

        if taxid is not None:
            species = ncbi.get_descendant_taxa(taxid, intermediate_nodes = True)
            for sp in species:
                expanded_taxa.add(sp)

    return expanded_taxa


def md5_seqs(fasta_file, translate, trans_table):
    from hashlib import md5
    md5_queries = {}
        
    for name, seq in iter_fasta_seqs(fasta_file, translate=translate, trans_table=trans_table):
        md5_seq = md5(seq.encode('utf-8')).hexdigest()
        md5_queries[name] = md5_seq
            
    return md5_queries

## END
