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
    tax_scope_mode = tax_scope_ids = target_taxa = target_orthologs = excluded_taxa = None
    
    go_evidence = go_excluded = None
    pfam_realign = trans_table = temp_dir = None
    md5 = None

    resume = None
        
    ##
    def __init__(self, args, annot, excel, report_orthologs):

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

        self.tax_scope = args.tax_scope  # Original value (e.g., "auto", "none", taxid)
        self.tax_scope_mode = args.tax_scope_mode
        self.tax_scope_ids = args.tax_scope_ids
                
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

        # Strict scope-OG filter. When True (default), events whose
        # containing OG (sp_events.og_lca) is broader than the seed's
        # tax_scope ceiling are skipped before fetching their orthologs
        # — consistent with the auto tax_scope intent. Set False to
        # restore the legacy permissive behaviour via
        # `--no-scope_strict_og`.
        #
        # `None` is treated as "use default" (True), not as a sentinel
        # for "False". A direct `bool(getattr(..., True))` would silently
        # flip to False when callers (tests, library code) construct an
        # `argparse.Namespace` with `scope_strict_og=None` — that's a
        # silent correctness bug we'd never see in CI.
        _strict = getattr(args, "scope_strict_og", True)
        self.scope_strict_og = True if _strict is None else bool(_strict)

        self.resume = args.resume
        
        return

    def _applied_filters(self):
        """Return a dict of the resolved annotation-stage filter values
        the run is using, for the `## applied filters:` block in
        .emapper.annotations (Phase 7.1a). Defaults are recorded as
        their actual resolved values so a reader doesn't need to know
        emapper's defaults to interpret the file.
        """
        return {
            "tax_scope":              self.tax_scope,
            "tax_scope_mode":         self.tax_scope_mode,
            "target_orthologs":       self.target_orthologs,
            "target_taxa":            self.target_taxa,
            "excluded_taxa":          self.excluded_taxa,
            "seed_ortholog_evalue":   self.seed_ortholog_evalue,
            "seed_ortholog_score":    self.seed_ortholog_score,
            "go_evidence":            self.go_evidence,
            "go_excluded":            self.go_excluded,
            "pfam_realign":           self.pfam_realign,
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
        """Dispatch annotation through the v7 batch path. Legacy per-hit
        and dbmem paths were removed in Phase 2 commit 3."""
        print(colorify("Functional annotation of hits...", "lgreen"), file=sys.stderr)
        eggnog_db = get_eggnog_db()
        if not eggnog_db._int_mode:
            raise EmapperException(
                "This eggnog-mapper version (v3) requires a v7+ "
                "integer-encoded eggnog.db. Rebuild the data dir with "
                "eggnog-builder build-emapper, or downgrade to v2 to use a "
                "legacy database."
            )
        return self._annotate_batched(hits_gen_func, annots_parser, eggnog_db)


    ##
    def _annotate_batched(self, hits_gen_func, annots_parser, eggnog_db):
        """Batch annotation: pre-fetch DB data per batch, CPU work in parallel.
        Only for v7+ integer-encoded databases.
        """
        from .batch_annotate import annotate_batch, get_lineage_cache
        from .tax_scopes.v7_scope import parse_v7_tax_scope

        # v7 databases use lineage-based filtering instead of OG-level filtering.
        # Parse tax_scope for v7 mode.
        v7_tax_scope = None
        v7_tax_scope_auto = False

        tax_scope = self.tax_scope
        if tax_scope and tax_scope.lower() not in ("none",):
            lineage_cache = get_lineage_cache()
            v7_tax_scope, v7_tax_scope_auto = parse_v7_tax_scope(tax_scope, lineage_cache)

            if v7_tax_scope_auto:
                print(colorify("Using auto tax_scope: orthologs filtered by seed ortholog's "
                               "taxonomic domain (Metazoa/Plants/Fungi/Bacteria/Archaea).",
                               "lblue"), file=sys.stderr)
            elif v7_tax_scope:
                print(colorify(f"Using tax_scope lineage filter: {', '.join(sorted(v7_tax_scope))}",
                               "lblue"), file=sys.stderr)

        batch_size = 1000
        batch = []

        # Process-pool parallelism for the annotation phase. Each worker
        # inherits the parent's loaded engine (taxid_array,
        # lineage_cache) via fork-COW and reopens its own SQLite
        # connection in the post-fork initializer. Pool is created once
        # for the full mapper run and reused across batches.
        pool = None
        n_workers = max(int(getattr(self, "cpu", 1) or 1), 1)
        if n_workers > 1:
            import multiprocessing
            from eggnog_annotator.e7.annotate import (
                _register_worker_engine, _worker_init_after_fork,
            )
            from .batch_annotate import _get_engine
            engine = _get_engine(eggnog_db, v7_tax_scope, v7_tax_scope_auto)
            _register_worker_engine(engine)
            ctx = multiprocessing.get_context("fork")
            pool = ctx.Pool(n_workers, initializer=_worker_init_after_fork)
            print(colorify(f"Annotation pool: {n_workers} workers (fork)",
                           "lblue"), file=sys.stderr)

        try:
            for args_tuple in self.iter_hit_lines(hits_gen_func, annots_parser):
                # Pass through resumed hits immediately
                if args_tuple[-1] is not None:  # already annotated
                    yield ((args_tuple[0], args_tuple[-1]), True)
                    continue

                batch.append(args_tuple)

                if len(batch) >= batch_size:
                    yield from annotate_batch(
                        batch, eggnog_db,
                        annot=self.annot,
                        target_orthologs=self.target_orthologs,
                        target_taxa=self.target_taxa,
                        excluded_taxa=self.excluded_taxa,
                        tax_scope_mode=self.tax_scope_mode,
                        tax_scope_ids=None,  # v7 uses lineage filtering, not OG filtering
                        go_evidence=self.go_evidence,
                        go_excluded=self.go_excluded,
                        seed_ortholog_score=self.seed_ortholog_score,
                        seed_ortholog_evalue=self.seed_ortholog_evalue,
                        pool=pool,
                        v7_tax_scope=v7_tax_scope,
                        v7_tax_scope_auto=v7_tax_scope_auto,
                        dropped_writer=self._dropped_writer,
                        scope_strict_og=self.scope_strict_og,
                    )
                    batch = []

            if batch:
                yield from annotate_batch(
                    batch, eggnog_db,
                    annot=self.annot,
                    target_orthologs=self.target_orthologs,
                    target_taxa=self.target_taxa,
                    excluded_taxa=self.excluded_taxa,
                    tax_scope_mode=self.tax_scope_mode,
                    tax_scope_ids=None,  # v7 uses lineage filtering, not OG filtering
                    go_evidence=self.go_evidence,
                    go_excluded=self.go_excluded,
                    seed_ortholog_score=self.seed_ortholog_score,
                    seed_ortholog_evalue=self.seed_ortholog_evalue,
                    pool=pool,
                    v7_tax_scope=v7_tax_scope,
                    v7_tax_scope_auto=v7_tax_scope_auto,
                    dropped_writer=self._dropped_writer,
                    scope_strict_og=self.scope_strict_og,
                )

        except EmapperException:
            raise
        except Exception as e:
            import traceback
            traceback.print_exc()
            raise EmapperException(f"Error: batch annotation failed. "+str(e))
        finally:
            if pool is not None:
                pool.close()
                pool.join()

        return

    
    ##
    def iter_hit_lines(self, hits_gen_func, annots_parser):

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
            
            yield_tuple = (hit, self.annot, self.seed_ortholog_score, self.seed_ortholog_evalue,
                           self.tax_scope_mode, self.tax_scope_ids,
                           self.target_taxa, self.target_orthologs, self.excluded_taxa,
                           self.go_evidence, self.go_excluded, get_data_path(), annotation)
            
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
    
    annotation = (query_name, best_hit_name, best_hit_evalue, best_hit_score,
                  annotations,
                  og_cat_desc, max_annot_lvl, match_nog_names,
                  all_orthologies, annot_orthologs)
    
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
