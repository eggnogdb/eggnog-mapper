#!/usr/bin/env python3

import os, sys, time, traceback
import argparse, multiprocessing

if sys.version_info < (3, 9):
    sys.exit('Sorry, Python < 3.9 is not supported.')
    
# get the path of this script and add it to the "pythonpath"
SCRIPT_PATH = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.insert(0, SCRIPT_PATH)

from eggnogmapper.emapperException import EmapperException
from eggnogmapper.utils import colorify

from eggnogmapper.emapper import Emapper

from eggnogmapper.genepred.genepred_modes import GENEPRED_MODE_SEARCH, GENEPRED_MODE_PRODIGAL

from eggnogmapper.search.search_modes import \
    SEARCH_MODE_NO_SEARCH, SEARCH_MODE_DIAMOND, \
    SEARCH_MODE_MMSEQS2, get_eggnog_dmnd_db

from eggnogmapper.search.diamond.diamond import SENSMODES, SENSMODE_SENSITIVE, \
    ALLOW_OVERLAPS_NONE, ALLOW_OVERLAPS_ALL, ALLOW_OVERLAPS_DIFF_FRAME, ALLOW_OVERLAPS_OPPOSITE_STRAND, \
    DMND_ITERATE_YES, DMND_ITERATE_NO, DMND_ITERATE_DEFAULT, \
    DMND_ALGO_AUTO, DMND_ALGO_0, DMND_ALGO_1, DMND_ALGO_CTG, DMND_ALGO_DEFAULT

from eggnogmapper.search.hmmer.hmmer_setup import DEFAULT_PORT, DEFAULT_END_PORT

from eggnogmapper.annotation.pfam.pfam_modes import PFAM_REALIGN_NONE, PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO

from eggnogmapper.deco.decoration import \
    DECORATE_GFF_NONE, DECORATE_GFF_GENEPRED, DECORATE_GFF_FIELD_DEFAULT

# Ceiling-mode constants for --tax_scope
TAX_SCOPE_AUTO_NARROW = "auto"
TAX_SCOPE_AUTO_BROAD = "auto-broad"

from eggnogmapper.common import existing_file, existing_dir, get_data_path, set_data_path, pexists, \
    get_eggnogdb_file, get_eggnog_mmseqs_db, \
    get_version, get_full_version_info, get_citation, get_call_info, get_db_version, \
    ITYPE_CDS, ITYPE_PROTS, ITYPE_GENOME, ITYPE_META, \
    TIMEOUT_LOAD_SERVER

from eggnogmapper.backends import resolve_backend, DEFAULT_BACKEND


__description__ = ('A program for bulk functional annotation of novel '
                    'sequences using EggNOG database orthology assignments')
__author__ = 'Jaime Huerta Cepas'
__license__ = "GPL v2"

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

def create_arg_parser():
    
    parser = argparse.ArgumentParser(formatter_class=CustomFormatter)

    parser.add_argument('-v', '--version', action='store_true',
                        help="show version and exit.")

    parser.add_argument('--list_taxa', action="store_true",
                        help="List taxa available for --tax_scope, and exit")

    ##
    pg_exec = parser.add_argument_group('Execution Options')
    
    pg_exec.add_argument('--cpu', type=int, default=1, metavar='NUM_CPU',
                        help="Number of CPUs to be used. --cpu 0 to run with all available CPUs.")

    pg_exec.add_argument('--resume', action="store_true",
                          help=("Resumes a previous emapper run, skipping results in existing output files."))

    pg_exec.add_argument('--override', action="store_true",
                        help=(
                            "Overwrites output files if they exist. "
                            "By default, execution is aborted if conflicting files are detected."))
    
    ##
    pg_input = parser.add_argument_group('Input Data Options')

    pg_input.add_argument('-i', dest="input", metavar='FASTA_FILE', type=existing_file,
                          help=f'Input FASTA file containing query sequences (proteins by default; see --itype). '
                          f'Required unless -m {SEARCH_MODE_NO_SEARCH}.')

    pg_input.add_argument('--itype', dest="itype", choices = [ITYPE_CDS, ITYPE_PROTS, ITYPE_GENOME, ITYPE_META],
                          default=ITYPE_PROTS,
                          help=f'Type of data in the input (-i) file.')
    
    pg_input.add_argument('--annotate_hits_table', type=str, metavar='SEED_ORTHOLOGS_FILE',
                          help=f'Annotate TSV formatted table with 4 fields:'
                          f' query, hit, evalue, score. '
                          f' Usually, a .seed_orthologs file from a previous emapper.py run. '
                          f' Requires -m {SEARCH_MODE_NO_SEARCH}.')

    pg_input.add_argument("--data_dir", metavar='DIR', type=existing_dir,
                          help=('Path to eggnog-mapper databases. '
                                'By default, "data/" or the path specified in the '
                                'environment variable EGGNOG_DATA_DIR.')) # DATA_PATH in eggnogmapper.commons

    ##
    pg_genepred = parser.add_argument_group('Gene Prediction Options')
    
    pg_genepred.add_argument('--genepred', dest='genepred', type=str, choices = [GENEPRED_MODE_SEARCH, GENEPRED_MODE_PRODIGAL],
                              default = GENEPRED_MODE_SEARCH,
                              help=(
                                  f'This is applied when --itype {ITYPE_GENOME} or --itype {ITYPE_META}. '
                                  f'{GENEPRED_MODE_SEARCH}: gene prediction is inferred from Diamond/MMseqs2 blastx hits. '
                                  f'{GENEPRED_MODE_PRODIGAL}: gene prediction is performed using Prodigal. '
                              ))

    pg_genepred.add_argument('--trans_table', dest='trans_table', type=str, metavar='TRANS_TABLE_CODE',
                             help=(
                                 f"This option will be used for Prodigal, Diamond or MMseqs2, depending on the mode. "
                                 f"For Diamond searches, this option corresponds to the --query-gencode option. "
                                 f"For MMseqs2 searches, this option corresponds to the --translation-table option. "
                                 f"For Prodigal, this option corresponds to the -g/--trans_table option. "
                                 f"It is also used when --translate, check https://biopython.org/docs/1.75/api/Bio.Seq.html#Bio.Seq.Seq.translate. "
                                 f"Default is the corresponding programs defaults. "
                             ))

    pg_genepred.add_argument('--training_genome', dest='training_genome', type=existing_file, metavar='FILE',
                             help=(
                                 "A genome to run Prodigal in training mode. "
                                 "If this parameter is used, Prodigal will run in two steps: "
                                 "firstly in training mode, and secondly using the training to analize the emapper input data. "
                                 "See Prodigal documentation about Traning mode for more info. "
                                 f"Only used if --genepred {GENEPRED_MODE_PRODIGAL}."
                             ))

    pg_genepred.add_argument('--training_file', dest='training_file', type=str, metavar='FILE',
                             help=(
                                 "A training file from Prodigal training mode. "
                                 "If this parameter is used, Prodigal will run using this training file to analyze the emapper input data. "
                                 f"Only used if --genepred {GENEPRED_MODE_PRODIGAL}."
                             ))

    pg_genepred.add_argument('--allow_overlaps', dest='allow_overlaps', type=str, choices = [ALLOW_OVERLAPS_NONE,
                                                                                             ALLOW_OVERLAPS_OPPOSITE_STRAND,
                                                                                             ALLOW_OVERLAPS_DIFF_FRAME,
                                                                                             ALLOW_OVERLAPS_ALL],
                             default = ALLOW_OVERLAPS_NONE,
                             help = ("When using 'blastx'-based genepred (--genepred search --itype genome/metagenome) "
                                     "this option controls whether overlapping hits are reported or not, "
                                     "or if only those overlapping hits in a different strand or frame are reported. "
                                     "Also, the degree of accepted overlap can be controlled with --overlap_tol."))

    pg_genepred.add_argument('--overlap_tol', dest='overlap_tol', type=float, default=0.0, metavar='FLOAT',
                             help=("This value (0-1) is the proportion such that if (overlap size / hit length) "
                                   "> overlap_tol, hits are considered to overlap. "
                                   "e.g: if overlap_tol is 0.0, any overlap is considered as such. "
                                   "e.g: if overlap_tol is 1.0, one of the hits must overlap entirely to "
                                   "consider that hits do overlap."))
                             
    ##
    pg_search = parser.add_argument_group('Search Options')

    pg_search.add_argument('-m', dest='mode',
                           choices = [SEARCH_MODE_DIAMOND, SEARCH_MODE_MMSEQS2, SEARCH_MODE_NO_SEARCH],
                           default=SEARCH_MODE_DIAMOND,
                           help=(
                               f'{SEARCH_MODE_DIAMOND}: search seed orthologs using diamond (-i is required). '
                               f'{SEARCH_MODE_MMSEQS2}: search seed orthologs using MMseqs2 (-i is required). '
                               f'{SEARCH_MODE_NO_SEARCH}: skip seed orthologs search (--annotate_hits_table is required, unless --no_annot). '
                           ))

    ##
    pg_diamond_mmseqs = parser.add_argument_group('Search filtering common options')

    pg_diamond_mmseqs.add_argument('--pident', dest='pident', type=float, default=None,
                                   help=(
                                       f'Report only alignments above or equal to the given percentage of identity (0-100).'))

    pg_diamond_mmseqs.add_argument('--query_cover', dest='query_cover', type=float, default=None,
                                   help=(
                                       f'Report only alignments above or equal the given percentage of query cover (0-100).'))

    pg_diamond_mmseqs.add_argument('--subject_cover', dest='subject_cover', type=float, default=None,
                                   help=(
                                       f'Report only alignments above or equal the given percentage of subject cover (0-100).'))
    
    pg_diamond_mmseqs.add_argument('--evalue', dest='evalue', type=float, default=0.001,
                                   help='Report only alignments below or equal the e-value threshold.')

    pg_diamond_mmseqs.add_argument('--score', dest='score', type=float, default=None,
                                   help='Report only alignments above or equal the score threshold.')

    ##
    pg_diamond = parser.add_argument_group('Diamond Search Options')

    pg_diamond.add_argument('--dmnd_algo', dest="dmnd_algo", choices = [DMND_ALGO_AUTO, DMND_ALGO_0, DMND_ALGO_1, DMND_ALGO_CTG],
                            default = DMND_ALGO_DEFAULT,
                            help=("Diamond's --algo option, which can be tuned to search small query sets. "
                                  "By default, it is adjusted automatically. "
                                  f"However, the {DMND_ALGO_CTG} option should be activated manually. "
                                  "If you plan to search a small input set of sequences, use --dmnd_algo ctg to make it faster."
                            ))
	
    pg_diamond.add_argument('--dmnd_db', dest="dmnd_db", metavar='DMND_DB_FILE',
		            help="Path to DIAMOND-compatible database")

    pg_diamond.add_argument('--dmnd_sensmode', dest='dmnd_sensmode',
                            choices = SENSMODES,
                            default=SENSMODE_SENSITIVE,
                            help=(
                                "Diamond's sensitivity mode. "
                                "When `--dmnd_iterate yes` (the default), this is the "
                                "*ceiling* sensitivity that diamond escalates to — each "
                                "successive iteration runs at a higher sensitivity, only "
                                "for queries that came back unaligned in the previous "
                                "iteration, until this ceiling is reached. Easy queries "
                                "are caught at the fastest mode; only divergent queries "
                                "pay the higher-sensitivity cost. "
                                "When `--dmnd_iterate no`, this is the single sensitivity "
                                "used for the search. "
                                "emapper's default is "+SENSMODE_SENSITIVE+" "
                                "(diamond's own default is `default`, faster but less "
                                "sensitive — activate with --dmnd_sensmode default)."
                            ))

    pg_diamond.add_argument('--dmnd_iterate', dest='dmnd_iterate', choices = [DMND_ITERATE_YES, DMND_ITERATE_NO],
                            default = DMND_ITERATE_DEFAULT,
                            help=(
                                f"--dmnd_iterate {DMND_ITERATE_YES} --> activates the --iterate option of diamond for iterative searches, "
                                f"from faster, less sensitive modes, up to the *ceiling* sensitivity specified with --dmnd_sensmode "
                                f"(default `{SENSMODE_SENSITIVE}` — easy queries are caught fast, only divergent queries escalate). "
                                f"Available since diamond 2.0.11. --dmnd_iterate {DMND_ITERATE_NO} --> disables the --iterate mode "
                                f"(single pass at --dmnd_sensmode sensitivity)."
                            ))
        
    pg_diamond.add_argument('--dmnd_matrix', dest='dmnd_matrix',
                            choices = ['BLOSUM62', 'BLOSUM90','BLOSUM80','BLOSUM50','BLOSUM45','PAM250','PAM70','PAM30'],
                            default=None, help='Scoring matrix')

    pg_diamond.add_argument('--dmnd_frameshift', dest='dmnd_frameshift', type=int, default=None, 
                            help='Diamond --frameshift/-F option. Not used by default. Recommended by diamond: 15.')
    
    pg_diamond.add_argument('--dmnd_gapopen', dest='dmnd_gapopen', type=int, default=None,
                            help='Gap open penalty')

    pg_diamond.add_argument('--dmnd_gapextend', dest='dmnd_gapextend', type=int, default=None,
                            help='Gap extend penalty')

    pg_diamond.add_argument('--dmnd_block_size', dest='dmnd_block_size', type=float, default=None, metavar='BLOCK_SIZE',
                            help=("Diamond -b/--block-size option. When unset, emapper auto-picks "
                                  "from host RAM: <32 GB→diamond default, 32-48 GB→4, 48-96 GB→6, "
                                  ">=96 GB→8. Larger = faster but uses more RAM. Diamond peak RAM "
                                  "≈ block_size × 6 + db_size / index_chunks + threads × 0.5 GB."))

    pg_diamond.add_argument('--dmnd_index_chunks', dest='dmnd_index_chunks', type=int, default=None, metavar='CHUNKS',
                            help=("Diamond -c/--index-chunks option. When unset, emapper auto-picks "
                                  "from host RAM: <32 GB→diamond default (4), 32-96 GB→2, "
                                  ">=96 GB→1 (full DB resident). Smaller = faster but uses more RAM."))

    pg_diamond.add_argument('--dmnd_top', dest='dmnd_top', type=int, default=1, metavar='N',
                            choices=[1, 3],
                            help=("Diamond hit count per query for protein/CDS searches. "
                                  "Default 1 (=> `--max-target-seqs 1`): emapper only consumes the "
                                  "single bitscore-best hit as the seed, so we let diamond prune as "
                                  "soon as it has it (~20-30 %% faster on typical proteomes). "
                                  "Set to 3 (=> `--top 3`) to keep all hits within 3 %% of top "
                                  "score per query — useful only if you've extended emapper to pick "
                                  "the seed by a non-bitscore criterion (pident / qcov tie-break). "
                                  "No effect for genome/metagenome itype runs."))

    pg_diamond.add_argument('--dmnd_ignore_warnings', action="store_true",
                            help=(
                                "Diamond --ignore-warnings option. "
                                "It avoids Diamond stopping due to warnings (e.g. when a protein contains only ATGC symbols."
                            ))
    
    ##
    pg_mmseqs = parser.add_argument_group('MMseqs2 Search Options')

    pg_mmseqs.add_argument('--mmseqs_db', dest="mmseqs_db", metavar='MMSEQS_DB_FILE',
		           help="Path to MMseqs2-compatible database")

    pg_mmseqs.add_argument('--start_sens', dest='start_sens', default=3, type=float, metavar='START_SENS',
                           help="Starting sensitivity.")

    pg_mmseqs.add_argument('--sens_steps', dest='sens_steps', default=3, type=int, metavar='SENS_STEPS',
                           help="Number of sensitivity steps.")

    pg_mmseqs.add_argument('--final_sens', dest='final_sens', default=7, type=float, metavar='FINAL_SENS',
                           help="Final sensititivy step.")

    pg_mmseqs.add_argument('--mmseqs_sub_mat', dest='mmseqs_sub_mat', default=None, type=str, metavar='SUBS_MATRIX',
                           help="Matrix to be used for --sub-mat MMseqs2 search option. Default=default used by MMseqs2")
    
    ##
    pg_annot = parser.add_argument_group('Annotation Options')
        
    pg_annot.add_argument("--no_annot", action="store_true",
                          help="Skip functional annotation, reporting only hits.")

    pg_annot.add_argument("--report_dropped", action="store_true",
                          help=("Write a TSV log of every query whose hit was dropped before "
                                "or after annotation, with the drop reason. Path: "
                                "<output>.emapper.dropped. Reasons: error_seed (seed name '-' "
                                "or 'ERROR'), evalue_above_threshold, score_below_threshold, "
                                "non_integer_seed_id, no_donor_orthologs (engine returned no "
                                "annotations under the active filters)."))

    pg_annot.add_argument(
        "--donor_pool",
        type=str,
        choices=["closest", "union"],
        default="closest",
        help=(
            "Controls how annotation donors are pooled across cascade tiers. "
            "'closest' (default): the first non-empty tier for each source wins "
            "(closest evolutionary distance). "
            "'union': all tiers are walked and their values unioned; "
            "confidence is the best (lowest) tier seen for each source."
        ),
    )

    pg_annot.add_argument(
        "--lazy_cascade",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Fetch and parse only the ortholog annotations the closest cascade "
            "actually consumes, instead of every ortholog of every seed up front "
            "(byte-identical output, large speedup on big runs). ON by default; "
            "active only with --donor_pool closest (the default) and a no-op with "
            "--donor_pool union. Use --no-lazy_cascade to force the eager path."
        ),
    )

    pg_annot.add_argument('--annot_evalue', dest='annot_evalue', default=0.001, type=float, metavar='MIN_E-VALUE',
                           help='Min E-value expected when searching for seed eggNOG ortholog.'
                           ' Queries not having a significant'
                           ' seed orthologs will not be annotated.')

    pg_annot.add_argument('--annot_score', dest='annot_score', default=None, type=float, metavar='MIN_SCORE',
                           help='Min bit score expected when searching for seed eggNOG ortholog.'
                           ' Queries not having a significant'
                           ' seed orthologs will not be annotated.')
    
    pg_annot.add_argument(
        "--tax_scope",
        type=str,
        default=TAX_SCOPE_AUTO_NARROW,
        help=(
            "Taxonomic ceiling for annotation: only speciation events whose "
            "ev_lca is at or below this ceiling are used for functional "
            "transfer. "
            f"'{TAX_SCOPE_AUTO_NARROW}' (default): per-seed automatic ceiling "
            "— Eukaryota for all eukaryotes, Prokaryota for all prokaryotes. "
            f"'{TAX_SCOPE_AUTO_BROAD}': broad-first variant — tries Metazoa "
            "before Eukaryota for animal seeds. "
            "A clade name (e.g. 'Metazoa') or numeric NCBI taxid: fixed "
            "ceiling for all seeds; hard-fails if the name cannot be resolved."
        ),
    )

    pg_annot.add_argument('--target_orthologs', choices=["one2one", "many2one",
                                                         "one2many","many2many", "all"],
                          default="all",
                          help='defines what type of orthologs (in relation to the seed ortholog) should be used for functional transfer')

    pg_annot.add_argument('--target_taxa', type=str, metavar="LIST_OF_TAX_IDS", 
                          default=None,
                          help=("Only orthologs from the specified comma-separated list of taxa and all its descendants "
                                "will be used for annotation transference. By default, all taxa are used."))

    pg_annot.add_argument('--excluded_taxa', type=str, metavar="LIST_OF_TAX_IDS",
                          default=None, 
                          help=('Orthologs from the specified comma-separated list of taxa and all its descendants will not be '
                                'used for annotation transference. By default, no taxa is excluded.'))

    pg_annot.add_argument("--report_orthologs", action="store_true",
                          help="Output the list of orthologs found for each query to a .orthologs file")
    
    pg_annot.add_argument('--pfam_realign', type=str,
                          choices=(PFAM_REALIGN_NONE, PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO),
                          default=PFAM_REALIGN_NONE,
                          help=(f'Realign queries to PFAM domains. '
                                f'{PFAM_REALIGN_NONE} (default) = no realignment; PFAM domains are '
                                'transferred from the seed orthologs. '
                                f'{PFAM_REALIGN_REALIGN} = realign queries to the PFAM domains transferred '
                                'from orthologs. '
                                f'{PFAM_REALIGN_DENOVO} = realign queries against the whole PFAM database '
                                '(de novo domain discovery). '
                                'The realignment modes run an hmmpgmd server '
                                '(see --num_servers, --num_workers, --port, --end_port).'))

    pg_annot.add_argument('-p', '--port', dest='port', type=int, default=DEFAULT_PORT, metavar='PORT',
                          help=('Port used to setup hmmpgmd server for --pfam_realign modes.'))

    pg_annot.add_argument('--end_port', dest='end_port', type=int, default=DEFAULT_END_PORT, metavar='PORT',
                          help=('Last port to be used to setup hmmpgmd server for --pfam_realign modes.'))

    pg_annot.add_argument('--num_servers', dest='num_servers', type=int, default=1, metavar="NUM_SERVERS",
                          help=("Number of hmmpgmd servers to fire up for --pfam_realign. "
                                "CPUs specified with --cpu will be distributed among servers and workers."))

    pg_annot.add_argument('--num_workers', dest='num_workers', type=int, default=1, metavar="NUM_WORKERS",
                          help=("Number of workers per server (--num_servers) for --pfam_realign. "
                                "CPUs specified with --cpu will be distributed among servers and workers."))

    pg_annot.add_argument('--timeout_load_server', dest='timeout_load_server', type=int, default=TIMEOUT_LOAD_SERVER, metavar="TIMEOUT_LOAD_SERVER",
                          help="Number of attempts to load a server on a specific port. If failed, the next numerical port will be tried.")

    pg_annot.add_argument("--md5", action="store_true",
                          help="Adds the md5 hash of each query as an additional field in annotations output files.")

    pg_annot.add_argument("--annot_batch_size", type=int, default=125, metavar='N',
                          help=("Number of seed orthologs annotated per worker task (sub-batch). Each pool "
                                "task bulk-queries this many seeds at once; the outer batch is sized to give "
                                "every worker ~2 sub-batches. With sorted seeds (the default) the outer batch "
                                "is measured in distinct seeds so all workers stay fed even on highly redundant inputs. "
                                "Default: 125."))

    ##
    pg_out = parser.add_argument_group('Output options')

    pg_out.add_argument('--output', '-o', type=str, metavar='FILE_PREFIX',
                        help="base name for output files")

    pg_out.add_argument("--output_dir", default=os.getcwd(), type=existing_dir, metavar='DIR',
                        help="Where output files should be written")

    pg_out.add_argument("--scratch_dir", metavar='DIR', type=existing_dir,
                        help='Write output files in a temporary scratch dir, move them to the final'
                        ' output dir when finished. Speed up large computations using network file'
                        ' systems.')

    pg_out.add_argument("--temp_dir", default=os.getcwd(), type=existing_dir, metavar='DIR',
                        help="Where temporary files are created. Better if this is a local disk.")

    pg_out.add_argument('--no_file_comments', action="store_true",
                        help="No header lines nor stats are included in the output files")

    pg_out.add_argument('--decorate_gff', type=str, default=DECORATE_GFF_NONE,
                        help=(
                            "Add search hits and/or annotation results to GFF file from gene prediction of a user specified one. "
                            f"{DECORATE_GFF_NONE} = no GFF decoration at all. GFF file from blastx-based gene prediction will be created anyway. "
                            f"{DECORATE_GFF_GENEPRED} = add search hits and/or annotations to GFF file from Prodigal or blastx-based gene prediction. "
                            f"FILE = decorate the specified pre-existing GFF FILE. e.g. --decorage_gff myfile.gff "
                            f"You change the field interpreted as ID of the feature with --decorate_gff_ID_field. "
                        ))

    pg_out.add_argument('--decorate_gff_ID_field', type=str, default=DECORATE_GFF_FIELD_DEFAULT,
                        help="Change the field used in GFF files as ID of the feature.")
    pg_out.add_argument("--excel", action="store_true",
                        help="Output annotations also in .xlsx format.")
        
    return parser


##
def parse_args(parser):
    
    args = parser.parse_args()

    if args.version:
        version = ""
        try:
            version = get_full_version_info()
        except Exception:
            version = get_version()
        print(version)
        sys.exit(0)

    if args.data_dir:
        set_data_path(args.data_dir)
    elif "EGGNOG_DATA_DIR" in os.environ:
        set_data_path(os.environ["EGGNOG_DATA_DIR"])
    else:
        set_data_path(resolve_backend(DEFAULT_BACKEND))

    args.call_info = get_call_info()

    if args.list_taxa:
        print_taxa()
        sys.exit(0)

    if args.cpu == 0:
        args.cpu = multiprocessing.cpu_count()

    if args.resume == True and args.override == True:
        parser.error('Only one of --resume or --override is allowed.')        

    # Gene prediction
    if args.training_genome is not None and args.training_file is None:
        parser.error('"--training_genome requires --training_file"')

    if args.training_genome is None and args.training_file is not None:
        if not os.path.isfile(args.training_file):
            parser.error('"--training_file must point to an existing file, if no --training_genome is provided."')
    
    # Search modes
    if args.mode == SEARCH_MODE_DIAMOND:
        dmnd_db = get_eggnog_dmnd_db(args.dmnd_db, args.mode, get_data_path())
        if not pexists(dmnd_db):
            print(colorify('DIAMOND database %s not present. Use download_eggnog_data.py to fetch it' % dmnd_db, 'red'))
            raise EmapperException()

        if args.input is not None:
            if args.annotate_hits_table is not None:
                print(colorify(f"--annotate_hits_table will be ignored, due to -m {args.mode}", 'blue'))
                args.annotate_hits_table = None
        else:
            # the default -m is diamond, but we will consider -m no_search as default when
            # --annotate_hits_table has been provided and -i has not been provided
            if args.annotate_hits_table is not None:
                print(colorify(f"Assuming -m {SEARCH_MODE_NO_SEARCH}", 'blue'))
                args.mode = SEARCH_MODE_NO_SEARCH
            else:
                parser.error('An input fasta file is required (-i)')

        # Output file required
        if not args.output:
            parser.error('An output project name is required (-o)')
            
    elif args.mode == SEARCH_MODE_MMSEQS2:
        mmseqs_db = args.mmseqs_db if args.mmseqs_db else get_eggnog_mmseqs_db()
        if not pexists(mmseqs_db):
            print(colorify('MMseqs2 database %s not present. Use download_eggnog_data.py to fetch it' % mmseqs_db, 'red'))
            raise EmapperException()

        if not args.input:
            parser.error('An input fasta file is required (-i)')

        # Output file required
        if not args.output:
            parser.error('An output project name is required (-o)')

        if args.annotate_hits_table is not None:
            print(colorify(f"--annotate_hits_table will be ignored, due to -m {SEARCH_MODE_MMSEQS2}", 'blue'))
            args.annotate_hits_table = None
            
    elif args.mode == SEARCH_MODE_NO_SEARCH:
        if args.no_annot == False and not args.annotate_hits_table:
            parser.error(f'No search mode (-m {SEARCH_MODE_NO_SEARCH}) requires a hits table to annotate (--annotate_hits_table FILE.seed_orthologs)')
        if args.md5 == True and args.input is None:
            parser.error(f'--md5 requires an input FASTA file (-i FASTA).')            
            
    else:
        parser.error(f'unrecognized search mode (-m {args.mode})')


    # Search thresholds
    args.dmnd_evalue = args.mmseqs_evalue = args.evalue
    args.dmnd_score = args.mmseqs_score = args.score
    args.qcov = args.query_cover
    
    # Annotation options
    if args.no_annot == False or args.report_orthologs == True:
        if not pexists(get_eggnogdb_file()):
            print(colorify('Annotation database data/eggnog.db not present. Use download_eggnog_data.py to fetch it', 'red'))
            raise EmapperException()

        if args.target_taxa is not None:
            args.target_taxa = args.target_taxa.split(",")
        if args.excluded_taxa is not None:
            args.excluded_taxa = args.excluded_taxa.split(",")
        
    # PFAM annotation options
    
    if args.pfam_realign == PFAM_REALIGN_NONE:
        pass
    elif args.pfam_realign == PFAM_REALIGN_REALIGN or args.pfam_realign == PFAM_REALIGN_DENOVO:
        if not args.input:
            parser.error(f'An input fasta file is required (-i) for --pfam_realign {args.pfam_realign}')
    else:
        raise ValueError(f'Invalid --pfam_realign option {args.pfam_realign}')

    total_workers = args.num_workers * args.num_servers
    if args.cpu < total_workers:
        parser.error(f"Less cpus ({args.cpu}) than total workers ({total_workers}) were specified.")
    if args.cpu % total_workers != 0:
        parser.error(f"Number of cpus ({args.cpu}) must be a multiple of total workers ({total_workers}).")        

    args.cpus_per_worker = int(args.cpu / total_workers)
    
    return args


if __name__ == "__main__":
    __spec__ = None
    try:
        start_time = time.time()
        
        parser = create_arg_parser()
        args = parse_args(parser)
        
        print('# ', get_version())
        print('# emapper.py ', ' '.join(sys.argv[1:]))
            
        emapper = Emapper(args.itype, args.genepred, args.mode, (not args.no_annot),
                          args.excel, args.report_orthologs, args.decorate_gff,
                          args.output, args.output_dir, args.scratch_dir,
                          args.resume, args.override)
        
        n, elapsed_time = emapper.run(args, args.input, args.annotate_hits_table)

        elapsed_time = time.time() - start_time

        addons = [args.mode, args.genepred]
        try:
            db_version = get_db_version()
        except Exception:
            db_version = None
        print(get_citation(addons, db_version=db_version))
        print(f'Total hits processed: {n}')
        print(f'Total time: {elapsed_time:.0f} secs')
        
    except EmapperException as ee:
        print(ee)
        sys.exit(1)
    except Exception as e:
        traceback.print_exc()
        # print(e)
        sys.exit(1)
    else:
        print("FINISHED")
        sys.exit(0)

## END
