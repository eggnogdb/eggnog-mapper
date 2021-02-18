#!/usr/bin/env python3

import os, sys, time, traceback
import argparse, multiprocessing

if sys.version_info < (3,7):
    sys.exit('Sorry, Python < 3.7 is not supported')
    
# get the path of this script and add it to the "pythonpath"
SCRIPT_PATH = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.insert(0, SCRIPT_PATH)

from eggnogmapper.emapperException import EmapperException
from eggnogmapper.emapper import Emapper
from eggnogmapper.genepred.genepred_modes import GENEPRED_MODE_SEARCH, GENEPRED_MODE_PRODIGAL
from eggnogmapper.search.search_modes import SEARCH_MODE_NO_SEARCH, SEARCH_MODE_DIAMOND, SEARCH_MODE_HMMER, SEARCH_MODE_MMSEQS2, SEARCH_MODE_CACHE
from eggnogmapper.search.diamond.diamond import SENSMODES, SENSMODE_SENSITIVE
from eggnogmapper.search.hmmer.hmmer_search import QUERY_TYPE_SEQ, QUERY_TYPE_HMM, DB_TYPE_SEQ, DB_TYPE_HMM
from eggnogmapper.search.hmmer.hmmer_setup import DEFAULT_PORT, DEFAULT_END_PORT
from eggnogmapper.annotation.pfam.pfam_modes import PFAM_TRANSFER_BEST_OG, PFAM_TRANSFER_NARROWEST_OG, PFAM_TRANSFER_SEED_ORTHOLOG, \
    PFAM_REALIGN_NONE, PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO

from eggnogmapper.common import existing_file, existing_dir, set_data_path, pexists, \
    get_eggnogdb_file, get_eggnog_dmnd_db, get_eggnog_mmseqs_db, \
    get_version, get_full_version_info, get_citation, get_call_info, ITYPE_CDS, ITYPE_PROTS, ITYPE_GENOME, ITYPE_META

from eggnogmapper.utils import colorify


__description__ = ('A program for bulk functional annotation of novel '
                    'sequences using EggNOG database orthology assignments')
__author__ = 'Jaime Huerta Cepas'
__license__ = "GPL v2"

def create_arg_parser():
    
    parser = argparse.ArgumentParser()

    parser.add_argument('-v', '--version', action='store_true',
                        help="show version and exit.")

    parser.add_argument('--list_taxa', action="store_true",
                        help="List taxa available for --tax_scope, and exit")

    ##
    pg_exec = parser.add_argument_group('Execution Options')
    
    pg_exec.add_argument('--cpu', type=int, default=1, metavar='NUM_CPU',
                        help="Number of CPUs to be used. --cpu 0 to run with all available CPUs. Default: 2")
    
    ##
    pg_input = parser.add_argument_group('Input Data Options')

    pg_input.add_argument('-i', dest="input", metavar='FASTA_FILE', type=existing_file,
                          help=f'Input FASTA file containing query sequences (proteins by default; see --translate). '
                          f'Required unless -m {SEARCH_MODE_NO_SEARCH} and --annotate_hits_table')

    pg_input.add_argument('--itype', dest="itype", choices = [ITYPE_CDS, ITYPE_PROTS, ITYPE_GENOME, ITYPE_META],
                          default=ITYPE_PROTS,
                          help=f'Input type of the data in the file specified in the -i option')
    
    pg_input.add_argument('--translate', action="store_true",
                          help='Assume input sequences are CDS instead of proteins')

    pg_input.add_argument('--annotate_hits_table', type=str, metavar='SEED_ORTHOLOGS_FILE',
                          help=f'Annotate TSV formatted table with 4 fields:'
                          f' query, hit, evalue, score. Requires -m {SEARCH_MODE_NO_SEARCH}.')

    pg_input.add_argument('-c', '--cache', dest="cache_file", metavar='FILE', type=existing_file,
                          help=f'File containing annotations and md5 hashes of queries, to be used as cache. '
                          f'Required if -m {SEARCH_MODE_CACHE}')
        
    pg_input.add_argument("--data_dir", metavar='DIR', type=existing_dir,
                          help='Path to eggnog-mapper databases.') # DATA_PATH in eggnogmapper.commons

    ##
    pg_genepred = parser.add_argument_group('Gene Prediction Options')
    pg_genepred.add_argument('--genepred', dest='genepred', type=str, choices = [GENEPRED_MODE_SEARCH, GENEPRED_MODE_PRODIGAL],
                              default = GENEPRED_MODE_SEARCH,
                              help=(
                                  f'This applied when --itype {ITYPE_GENOME} or --itype {ITYPE_META}. '
                                  f'{GENEPRED_MODE_SEARCH}: gene prediction is inferred from hits obtained from the search step. '
                                  f'{GENEPRED_MODE_PRODIGAL}: gene prediction is performed from proteins predicted using prodigal. '
                              ))

    pg_genepred.add_argument('--trans_table', dest='trans_table', type=str, metavar='TRANS_TABLE_CODE',
                             help=(
                                 f"This option will be used for Prodigal, Diamond or MMseqs2, depending on the mode. "
                                 f"For Diamond searches, this option corresponds to the --query-gencode option. "
                                 f"For MMseqs2 searches, this option corresponds to the --translation-table option. "
                                 f"For Prodigal, this option corresponds to the -g/--trans_table option. "
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
                             
    ##
    pg_search = parser.add_argument_group('Search Options')

    pg_search.add_argument('-m', dest='mode', 
                           choices = [SEARCH_MODE_DIAMOND, SEARCH_MODE_MMSEQS2, SEARCH_MODE_HMMER, SEARCH_MODE_NO_SEARCH, SEARCH_MODE_CACHE],
                           default=SEARCH_MODE_DIAMOND,
                           help=(
                               f'{SEARCH_MODE_DIAMOND}: search seed orthologs using diamond (-i is required). '
                               f'{SEARCH_MODE_MMSEQS2}: search seed orthologs using MMseqs2 (-i is required). '
                               f'{SEARCH_MODE_HMMER}: search seed orthologs using HMMER. (-i is required). '
                               f'{SEARCH_MODE_NO_SEARCH}: skip seed orthologs search (--annotate_hits_table is required, unless --no_annot). '
                               f'{SEARCH_MODE_CACHE}: skip seed orthologs search and annotate based on cached results (-i and -c are required). '
                               f'Default:{SEARCH_MODE_DIAMOND}'
                           ))

    ##
    pg_diamond_mmseqs = parser.add_argument_group('Search filtering common options')

    pg_diamond_mmseqs.add_argument('--pident', dest='pident', type=float, default=None,
                                   help=(
                                       f'Report only alignments above or equal to the given percentage of identity (0-100). Default=(not used). '
                                       f'No effect if -m {SEARCH_MODE_HMMER}.'))
    
    pg_diamond_mmseqs.add_argument('--query_cover', dest='query_cover', type=float, default=None,
                                   help=(
                                       f'Report only alignments above or equal the given percentage of query cover (0-100). Default=(not used). '
                                       f'No effect if -m {SEARCH_MODE_HMMER}.'))
    
    pg_diamond_mmseqs.add_argument('--subject_cover', dest='subject_cover', type=float, default=None,
                                   help=(
                                       f'Report only alignments above or equal the given percentage of subject cover (0-100). Default=(not used). '
                                       f'No effect if -m {SEARCH_MODE_HMMER}.'))
    
    pg_diamond_mmseqs.add_argument('--evalue', dest='evalue', type=float, default=0.001,
                                   help='Report only alignments below or equal the e-value threshold. Default=0.001')

    pg_diamond_mmseqs.add_argument('--score', dest='score', type=float, default=None,
                                   help='Report only alignments above or equal the score threshold. Default=(not used)')

    ##
    pg_diamond = parser.add_argument_group('Diamond Search Options')
	
    pg_diamond.add_argument('--dmnd_db', dest="dmnd_db", metavar='DMND_DB_FILE',
		            help="Path to DIAMOND-compatible database")

    pg_diamond.add_argument('--sensmode', dest='sensmode', 
                            choices = SENSMODES, 
                            default=SENSMODE_SENSITIVE, help='Sensitive mode')
        
    pg_diamond.add_argument('--matrix', dest='matrix', 
                            choices = ['BLOSUM62', 'BLOSUM90','BLOSUM80','BLOSUM50','BLOSUM45','PAM250','PAM70','PAM30'], 
                            default=None, help='Scoring matrix')

    pg_diamond.add_argument('--gapopen', dest='gapopen', type=int, default=None, 
                            help='Gap open penalty')

    pg_diamond.add_argument('--gapextend', dest='gapextend', type=int, default=None, 
                            help='Gap extend  penalty')

    pg_diamond.add_argument('--block_size', dest='dmnd_block_size', type=float, default=None, metavar='BLOCK_SIZE',
                            help="Diamond -b/--block-size option. Default is the diamond's default.")

    pg_diamond.add_argument('--index_chunks', dest='dmnd_index_chunks', type=int, default=None, metavar='CHUNKS',
                            help="Diamond -c/--index-chunks option. Default is the diamond's default.")

    pg_diamond.add_argument('--outfmt_short', action="store_true",
                            help=(
                                "Diamond output will include only qseqid sseqid evalue and score. "
                                "This could help obtain better performance, if also no --pident, --query_cover or --subject_cover thresholds are used. "
                                "This option is ignored when the diamond search is run in blastx mode for gene prediction (see --genepred)."
                            ))
    
    ##
    pg_mmseqs = parser.add_argument_group('MMseqs2 Search Options')

    pg_mmseqs.add_argument('--mmseqs_db', dest="mmseqs_db", metavar='MMSEQS_DB_FILE',
		           help="Path to MMseqs2-compatible database")

    pg_mmseqs.add_argument('--start_sens', dest='start_sens', default=3, type=float, metavar='START_SENS',
                           help="Starting sensitivity. Default=3")

    pg_mmseqs.add_argument('--sens_steps', dest='sens_steps', default=3, type=int, metavar='SENS_STEPS',
                           help="Number of sensitivity steps. Default=3")

    pg_mmseqs.add_argument('--final_sens', dest='final_sens', default=7, type=float, metavar='FINAL_SENS',
                           help="Final sensititivy step. Default=7")

    pg_mmseqs.add_argument('--mmseqs_sub_mat', dest='mmseqs_sub_mat', default=None, type=str, metavar='SUBS_MATRIX',
                           help="Matrix to be used for --sub-mat MMseqs2 search option. Default=default used by MMseqs2")
    
    ##
    pg_hmmer = parser.add_argument_group('HMMER Search Options')

    pg_hmmer.add_argument('-d', '--database', dest='db', metavar='HMMER_DB_PREFIX',
                          help=('specify the target database for sequence searches. '
                                'Choose among: euk,bact,arch, or a database loaded in a server, db.hmm:host:port (see hmm_server.py)'))

    pg_hmmer.add_argument('--servers_list', dest="servers_list", metavar="FILE",
                          help="A FILE with a list of remote hmmpgmd servers. "
                                "Each row in the file represents a server, in the format 'host:port'. "
                                "If --servers_list is specified, host and port from -d option will be ignored.")
    
    pg_hmmer.add_argument('--qtype',  choices=[QUERY_TYPE_HMM, QUERY_TYPE_SEQ], default=QUERY_TYPE_SEQ,
                          help="Type of input data (-i). "
                          f"Default: {QUERY_TYPE_SEQ}")

    pg_hmmer.add_argument('--dbtype', dest="dbtype",
                          choices=[DB_TYPE_HMM, DB_TYPE_SEQ], default=DB_TYPE_HMM,
                          help="Type of data in DB (-db). "
                          f"Default: {DB_TYPE_HMM}")

    pg_hmmer.add_argument('--usemem', action="store_true",
                          help='''Use this option to allocate the whole database (-d) in memory using hmmpgmd.
                          If --dbtype hmm, the database must be a hmmpress-ed database.
                          If --dbtype seqdb, the database must be a HMMER-format database created with esl-reformat.
                          Database will be unloaded after execution.''')

    pg_hmmer.add_argument('-p', '--port', dest='port', type=int, default=DEFAULT_PORT, metavar='PORT',
                          help=('Port used to setup HMM server, when --usemem. Also used for --pfam_realign modes.'))
    
    pg_hmmer.add_argument('--end_port', dest='end_port', type=int, default=DEFAULT_END_PORT, metavar='PORT',
                          help=('Last port to be used to setup HMM server, when --usemem. Also used for --pfam_realign modes.'))

    pg_hmmer.add_argument('--num_servers', dest='num_servers', type=int, default=1, metavar="NUM_SERVERS",
                          help="When using --usemem, specify the number of servers to fire up."
                          " By default, only 1 server is used. Note that cpus specified with --cpu will be distributed among servers and workers."
                          " Also used for --pfam_realign modes.")
    
    pg_hmmer.add_argument('--num_workers', dest='num_workers', type=int, default=1, metavar="NUM_WORKERS",
                          help="When using --usemem, specify the number of workers per server (--num_servers) to fire up."
                          " By default, cpus specified with --cpu will be distributed among servers and workers."
                          " Also used for --pfam_realign modes.")

    pg_hmmer.add_argument('--hmm_maxhits', dest='maxhits', type=int, default=1, metavar='MAXHITS',
                          help="Max number of hits to report (0 to report all). Default=1.")

    pg_hmmer.add_argument('--report_no_hits', action="store_true",
                          help="Whether queries without hits should be included in the output table.")

    pg_hmmer.add_argument('--hmm_maxseqlen', dest='maxseqlen', type=int, default=5000, metavar='MAXSEQLEN',
                          help="Ignore query sequences larger than `maxseqlen`. Default=5000")

    pg_hmmer.add_argument('--Z', dest='Z', type=float, default=40000000, metavar='DB_SIZE',
                          help='Fixed database size used in phmmer/hmmscan'
                          ' (allows comparing e-values among databases). Default=40,000,000')

    pg_hmmer.add_argument('--cut_ga', action="store_true",
                          help="Adds the --cut_ga to hmmer commands (useful for Pfam mappings, for example). See hmmer documentation.")

    pg_hmmer.add_argument('--clean_overlaps', dest="clean_overlaps", type=str, default=None, metavar="none|all|clans|hmmsearch_all|hmmsearch_clans",
                          help='Removes those hits which overlap, keeping only the one with best evalue. '
                          'Use the "all" and "clans" options when performing a hmmscan type search (i.e. domains are in the database). '
                          'Use the "hmmsearch_all" and "hmmsearch_clans" options when using a hmmsearch type search (i.e. domains are the queries from -i file). '
                          'The "clans" and "hmmsearch_clans" and options will only have effect for hits to/from Pfam.')
    
    ##
    pg_annot = parser.add_argument_group('Annotation Options')
        
    pg_annot.add_argument("--no_annot", action="store_true",
                          help="Skip functional annotation, reporting only hits.")

    pg_annot.add_argument('--dbmem', action="store_true",
                          help='''Use this option to allocate the whole eggnog.db DB in memory.
                          Database will be unloaded after execution.''')
    

    pg_annot.add_argument('--seed_ortholog_evalue', default=0.001, type=float, metavar='MIN_E-VALUE',
                           help='Min E-value expected when searching for seed eggNOG ortholog.'
                           ' Queries not having a significant'
                           ' seed orthologs will not be annotated. Default=0.001')

    pg_annot.add_argument('--seed_ortholog_score', default=60, type=float, metavar='MIN_SCORE',
                           help='Min bit score expected when searching for seed eggNOG ortholog.'
                           ' Queries not having a significant'
                           ' seed orthologs will not be annotated. Default=60')
    
    pg_annot.add_argument("--tax_scope", type=str, default='auto', 
                          help=("Fix the taxonomic scope used for annotation, so only speciation events from a "
                                "particular clade are used for functional transfer. "
                                "By default ('auto'), it is automatically adjusted for every query sequence, with a predefined list of tax IDs. "
                                "'auto_broad' is the same, but using a broader list of tax IDs, aiming for more annotations, yet slower than 'auto'. "
                                "Use 'narrowest' to use the deepest or narrowest taxon among the OGs identified for each hit. "
                                "Use a comma-separated list of tax IDs to choose an OG among the OGs identified for each hit. "
                                "The list of tax IDs can be followed by 'narrowest' or 'none', to specify the behaviour when the tax ID is not found among OGs. "
                                "If only the list of tax IDs is specified, the default behaviour is 'none', "
                                "so that no OG will be used for annotation if not found among the tax IDs. "
                                "If 'narrowest' is specified, if no OG is found among the list of tax IDs, the narrowest OG will be used for annotation. "
                                "An example of list of tax IDs would be 2759,2157,2,1 for euk, arch, bact and root, in that order of preference. "))

    pg_annot.add_argument('--target_orthologs', choices=["one2one", "many2one",
                                                         "one2many","many2many", "all"],
                          default="all",
                          help='defines what type of orthologs (in relation to the seed ortholog) should be used for functional transfer')

    pg_annot.add_argument('--target_taxa', type=str, 
                          default=None,
                          help="Only orthologs from the specified taxa and all its descendants will be used for annotation transference. By default, all taxa are used.")

    pg_annot.add_argument('--excluded_taxa', type=str,
                          default=None, 
                          help='Orthologs from the specified taxa and all its descendants will not be used for annotation transference. By default, no taxa is excluded.')

    pg_annot.add_argument("--report_orthologs", action="store_true",
                          help="Output the list of orthologs found for each query to a .orthologs file")
    
    pg_annot.add_argument('--go_evidence', type=str, choices=('experimental', 'non-electronic', 'all'),
                          default='non-electronic',
                          help='Defines what type of GO terms should be used for annotation. '
                          'experimental = Use only terms inferred from experimental evidence. '
                          'non-electronic = Use only non-electronically curated terms')

    pg_annot.add_argument('--pfam_transfer', type=str, choices=(PFAM_TRANSFER_BEST_OG, PFAM_TRANSFER_NARROWEST_OG, PFAM_TRANSFER_SEED_ORTHOLOG),
                          default=PFAM_TRANSFER_BEST_OG,
                          help='Defines from which orthologs PFAM domains will be transferred. '
                          f'{PFAM_TRANSFER_BEST_OG} = PFAMs will be transferred from orthologs in the best Orthologous Group (OG). '
                          f'{PFAM_TRANSFER_NARROWEST_OG} = PFAMs will be transferred from orthologs in the narrowest OG. '
                          f'{PFAM_TRANSFER_SEED_ORTHOLOG} = PFAMs will be transferred from the seed ortholog. ')

    pg_annot.add_argument('--pfam_realign', type=str, choices=(PFAM_REALIGN_NONE, PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO),
                          default=PFAM_REALIGN_NONE,
                          help='Realign the queries to the PFAM domains. '
                          f'{PFAM_REALIGN_NONE} = no realignment is performed. PFAM annotation will be that transferred as specify in the --pfam_transfer option. '
                          f'{PFAM_REALIGN_REALIGN} = queries will be realigned to the PFAM domains found according to the --pfam_transfer option. '
                          f'{PFAM_REALIGN_DENOVO} = queries will be realigned to the whole PFAM database, ignoring the --pfam_transfer option. '
                          f'Check hmmer options (--num_servers, --num_workers, --port, --end_port) to change how the hmmpgmd server is run. ')

    pg_annot.add_argument("--md5", action="store_true",
                          help="Adds the md5 hash of each query as an additional field in annotations output files.")

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

    pg_out.add_argument('--resume', action="store_true",
                        help="Resumes a previous SEARCH_MODE_HMMER search, skipping reported hits in the output file. "
                        f"Only SEARCH_MODE_HMMER runs (-m {SEARCH_MODE_HMMER}) can be resumed.")
        
    pg_out.add_argument('--override', action="store_true",
                        help="Overwrites output files if they exist.")

    pg_out.add_argument("--temp_dir", default=os.getcwd(), type=existing_dir, metavar='DIR',
                        help="Where temporary files are created. Better if this is a local disk.")

    pg_out.add_argument('--no_file_comments', action="store_true",
                        help="No header lines nor stats are included in the output files")
        
    return parser


##
# Parses tax_scope command line argument
# to define tax_scope_mode and tax_scope_id (one or more tax IDs)
def __parse_tax_scope(tax_scope):
    tax_scope_mode = None
    tax_scope_id = None

    tax_scope_fields = tax_scope.strip().split(",")
    tax_scope_mode = tax_scope_fields[0]

    # Auto
    if tax_scope_mode in {"auto", "auto_broad"}:
        tax_scope_id = None

    # Narrowest
    elif tax_scope_mode == "narrowest":
        tax_scope_id = None

    # Tax IDs
    else:
        # Only the specified tax ID
        if len(tax_scope_fields) == 1:
            tax_scope_mode = "none"
            tax_scope_id = [tax_scope_fields[0]]

        # Tax ID list, with or without mode for those not found in the list
        elif len(tax_scope_fields) > 1:
            last_pos = tax_scope_fields[-1]
            if last_pos in ["narrowest", "auto", "none"]:
                tax_scope_mode = last_pos
                tax_scope_id = tax_scope_fields[:-1]
            else:
                tax_scope_mode = "none"
                tax_scope_id = tax_scope_fields
        else:
            raise EmapperException(f"Error: unrecognized tax scope format {tax_scope}.")

        if tax_scope_id is not None and len(tax_scope_id) > 0:
            tax_scope_id_int = []
            from eggnogmapper.vars import LEVEL_NAMES, LEVEL_DICT
            for tax_id in tax_scope_id:
                if tax_id in LEVEL_NAMES:
                    tax_scope_id_int.append(tax_id)
                elif tax_id in LEVEL_DICT:
                    tax_scope_id_int.append(LEVEL_DICT[tax_id])
                else:
                    raise EmapperException(f"Unrecognized tax ID, tax name or tax_scope mode: '{tax_id}'.")

            tax_scope_id = tax_scope_id_int

    return tax_scope_mode, tax_scope_id


##
def parse_args(parser):
    
    args = parser.parse_args()

    if "EGGNOG_DATA_DIR" in os.environ:
        set_data_path(os.environ["EGGNOG_DATA_DIR"])
    
    if args.data_dir:
        set_data_path(args.data_dir)
        
    if args.version:
        version = ""
        try:
            version = get_full_version_info()
        except Exception:
            version = get_version()
        print(version)
        sys.exit(0)

    args.call_info = get_call_info()

    if args.list_taxa:
        from eggnogmapper.vars import LEVEL_DEPTH, LEVEL_DICT, LEVEL_NAMES, LEVEL_PARENTS
        print("tax_name\ttax_id\tdepth\tparents\tparents_names")
        for tax_name, tax_id in LEVEL_DICT.items():
            depth = LEVEL_DEPTH.get(tax_id, "-")
            parents = LEVEL_PARENTS.get(tax_id, "-")
            parents_names = [LEVEL_NAMES.get(x, "-") for x in parents]
            print(f"{tax_name}\t{tax_id}\t{depth}\t{','.join(parents)}\t{','.join(parents_names)}")
        sys.exit(0)

    if args.cpu == 0:
        args.cpu = multiprocessing.cpu_count()


    # translate
    if args.itype in [ITYPE_GENOME, ITYPE_META, ITYPE_PROTS] and args.translate == True:
        parser.error('"--translate" only can be used with "--itype CDS"')

    # Gene prediction
    if args.training_genome is not None and args.training_file is None:
        parser.error('"--training_genome requires --training_file"')

    if args.training_genome is None and args.training_file is not None:
        if not os.path.isfile(args.training_file):
            parser.error('"--training_file must point to an existing file, if no --training_genome is provided."')
    
    # Search modes
    if args.mode == SEARCH_MODE_DIAMOND:
        dmnd_db = args.dmnd_db if args.dmnd_db else get_eggnog_dmnd_db()
        if not pexists(dmnd_db):
            print(colorify('DIAMOND database %s not present. Use download_eggnog_database.py to fetch it' % dmnd_db, 'red'))
            raise EmapperException()

        if args.input is not None:
            if args.annotate_hits_table is not None:
                print(colorify(f"--annotate_hits_table will be ignored, due to -m {SEARCH_MODE_DIAMOND}", 'blue'))
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

        if args.resume == True:
            print(colorify("Diamond jobs cannot be resumed. --resume will be ignored.", 'blue'))
            args.resume = False

        if args.annotate_hits_table is not None:
            print(colorify(f"--annotate_hits_table will be ignored, due to -m {SEARCH_MODE_DIAMOND}", 'blue'))
            args.annotate_hits_table = None
            
    elif args.mode == SEARCH_MODE_MMSEQS2:
        mmseqs_db = args.mmseqs_db if args.mmseqs_db else get_eggnog_mmseqs_db()
        if not pexists(mmseqs_db):
            print(colorify('MMseqs2 database %s not present. Use download_eggnog_database.py to fetch it' % mmseqs_db, 'red'))
            raise EmapperException()

        if not args.input:
            parser.error('An input fasta file is required (-i)')

        # Output file required
        if not args.output:
            parser.error('An output project name is required (-o)')

        if args.resume == True:
            print(colorify("MMseqs2 jobs cannot be resumed. --resume will be ignored.", 'blue'))
            args.resume = False

        if args.annotate_hits_table is not None:
            print(colorify(f"--annotate_hits_table will be ignored, due to -m {SEARCH_MODE_MMSEQS2}", 'blue'))
            args.annotate_hits_table = None
            
    elif args.mode == SEARCH_MODE_HMMER:

        # if args.usemem == True:
        #     total_workers = args.num_workers * args.num_servers
        #     if args.cpu < total_workers:
        #         parser.error(f"Less cpus ({args.cpu}) than total workers ({total_workers}) were specified.")
        #     if args.cpu % total_workers != 0:
        #         parser.error(f"Number of cpus ({args.cpu}) must be a multiple of total workers ({total_workers}).")        

        #     args.cpus_per_worker = int(args.cpu / total_workers)
        #     sys.stderr.write(f"CPUs per worker: {args.cpus_per_worker}\n")
        # else:
        #     args.cpus_per_worker = args.cpu
        
        if not args.input:
            parser.error('An input file is required (-i)')

        # Output file required
        if not args.output:
            parser.error('An output project name is required (-o)')

        # Hmmer database
        # NOTE: hmmer database format, name and checking if exists is done within hmmer module
        if not args.db:
            parser.error('HMMER mode requires a target database (-d, --database).')

        if args.itype == ITYPE_CDS:
            args.translate = True

        if (args.itype == ITYPE_GENOME or args.itype == ITYPE_META) and args.genepred == GENEPRED_MODE_SEARCH:
            parser.error('HMMER mode is not compatible with "--genepred search" option.')            

        if args.annotate_hits_table is not None:
            print(colorify(f"--annotate_hits_table will be ignored, due to -m {SEARCH_MODE_HMMER}", 'blue'))
            args.annotate_hits_table = None

        if args.clean_overlaps is not None:
            if args.clean_overlaps == "none":
                args.clean_overlaps = None

    elif args.mode == SEARCH_MODE_CACHE:
        if args.cache_file is None:
                parser.error('A file with annotations and md5 of queries is required (-c FILE)')
                
        if args.no_annot == True:
            parser.error(f'Cache mode (-m {SEARCH_MODE_CACHE}) should be used to annotate.')
            
    elif args.mode == SEARCH_MODE_NO_SEARCH:
        if args.no_annot == False and not args.annotate_hits_table:
            parser.error(f'No search mode (-m {SEARCH_MODE_NO_SEARCH}) requires a hits table to annotate (--annotate_hits_table FILE.seed_orthologs)')
        if args.md5 == True and args.input is None:
            parser.error(f'--md5 requires an input FASTA file (-i FASTA).')            
        # if args.no_annot == True and args.report_orthologs == False:
        #     parser.error(f'Nothing to do if running in no search mode (-m {SEARCH_MODE_NO_SEARCH}), with --no_annot and without --report_orthologs.')
            
    else:
        parser.error(f'unrecognized search mode (-m {args.mode})')


    # Search thresholds
    args.dmnd_evalue = args.mmseqs_evalue = args.hmm_evalue = args.evalue
    args.dmnd_score = args.mmseqs_score = args_hmm_score = args.score
    args.qcov = args.query_cover
    
    # Annotation options
    if args.no_annot == False or args.report_orthologs == True:
        if not pexists(get_eggnogdb_file()):
            print(colorify('Annotation database data/eggnog.db not present. Use download_eggnog_database.py to fetch it', 'red'))
            raise EmapperException()

        args.tax_scope_mode, args.tax_scope_id = __parse_tax_scope(args.tax_scope)
        if args.target_taxa is not None:
            args.target_taxa = args.target_taxa.split(",")
        if args.excluded_taxa is not None:
            args.excluded_taxa = args.excluded_taxa.split(",")
        
    # Sets GO evidence bases
    if args.go_evidence == 'experimental':
        args.go_evidence = set(["EXP","IDA","IPI","IMP","IGI","IEP"])
        args.go_excluded = set(["ND", "IEA"])

    elif args.go_evidence == 'non-electronic':
        args.go_evidence = None
        args.go_excluded = set(["ND", "IEA"])

    elif args.go_evidence == 'all':
        args.go_evidence = None
        args.go_excluded = None
        
    else:
        raise ValueError('Invalid --go_evidence value')

    # PFAM annotation options
    if args.pfam_transfer in [PFAM_TRANSFER_BEST_OG, PFAM_TRANSFER_NARROWEST_OG, PFAM_TRANSFER_SEED_ORTHOLOG]:
        pass
    else:
        raise ValueError(f'Invalid --pfam_transfer option {args.pfam_transfer}')
    
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

    try:

        parser = create_arg_parser()
        args = parse_args(parser)

        _total_time = time.time()
        
        print('# ', get_version())
        print('# emapper.py ', ' '.join(sys.argv[1:]))
            
        emapper = Emapper(args.itype, args.genepred, args.mode, (not args.no_annot), args.report_orthologs, args.output, args.output_dir, args.scratch_dir, args.resume, args.override)
        emapper.run(args, args.input, args.annotate_hits_table, args.cache_file)

        print(get_citation([args.mode, args.genepred]))
        print('Total time: %g secs' % (time.time()-_total_time))
        
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
