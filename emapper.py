#!/usr/bin/env python3

import os, sys, time, traceback
import argparse, multiprocessing

if sys.version_info < (3,7):
    sys.exit('Sorry, Python < 3.7 is not supported')
    
# get the path of this script and add it to the "pythonpath"
SCRIPT_PATH = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.insert(0, SCRIPT_PATH)

from eggnogmapper.emapperException import EmapperException
from eggnogmapper.utils import colorify

from eggnogmapper.emapper import Emapper

from eggnogmapper.genepred.genepred_modes import GENEPRED_MODE_SEARCH, GENEPRED_MODE_PRODIGAL

from eggnogmapper.search.search_modes import \
    SEARCH_MODE_NO_SEARCH, SEARCH_MODE_DIAMOND, \
    SEARCH_MODE_HMMER, SEARCH_MODE_MMSEQS2, SEARCH_MODE_CACHE, \
    SEARCH_MODE_NOVEL_FAMS, get_eggnog_dmnd_db

from eggnogmapper.search.diamond.diamond import SENSMODES, SENSMODE_SENSITIVE, \
    ALLOW_OVERLAPS_NONE, ALLOW_OVERLAPS_ALL, ALLOW_OVERLAPS_DIFF_FRAME, ALLOW_OVERLAPS_OPPOSITE_STRAND, \
    DMND_ITERATE_YES, DMND_ITERATE_NO, DMND_ITERATE_DEFAULT, \
    DMND_ALGO_AUTO, DMND_ALGO_0, DMND_ALGO_1, DMND_ALGO_CTG, DMND_ALGO_DEFAULT

from eggnogmapper.search.hmmer.hmmer_search import \
    QUERY_TYPE_SEQ, QUERY_TYPE_HMM, \
    DB_TYPE_SEQ, DB_TYPE_HMM

from eggnogmapper.search.hmmer.hmmer_setup import DEFAULT_PORT, DEFAULT_END_PORT

from eggnogmapper.annotation.pfam.pfam_modes import PFAM_REALIGN_NONE, PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO

from eggnogmapper.deco.decoration import \
    DECORATE_GFF_NONE, DECORATE_GFF_GENEPRED, DECORATE_GFF_FIELD_DEFAULT

from eggnogmapper.annotation.tax_scopes.tax_scopes import \
    parse_tax_scope, print_taxa, \
    TAX_SCOPE_MODE_BROADEST, TAX_SCOPE_MODE_INNER_BROADEST, \
    TAX_SCOPE_MODE_INNER_NARROWEST, TAX_SCOPE_MODE_NARROWEST

from eggnogmapper.common import existing_file, existing_dir, get_data_path, set_data_path, pexists, \
    get_eggnogdb_file, get_eggnog_mmseqs_db, \
    get_version, get_full_version_info, get_citation, get_call_info, \
    ITYPE_CDS, ITYPE_PROTS, ITYPE_GENOME, ITYPE_META, \
    MP_START_METHOD_DEFAULT, MP_START_METHOD_FORK, MP_START_METHOD_SPAWN, MP_START_METHOD_FORKSERVER


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
                        help="List taxa available for --tax_scope/--tax_scope_mode, and exit")

    ##
    pg_exec = parser.add_argument_group('Execution Options')
    
    pg_exec.add_argument('--cpu', type=int, default=1, metavar='NUM_CPU',
                        help="Number of CPUs to be used. --cpu 0 to run with all available CPUs.")

    pg_exec.add_argument('--mp_start_method', type=str, default=MP_START_METHOD_DEFAULT,
                         choices = [MP_START_METHOD_FORK, MP_START_METHOD_SPAWN, MP_START_METHOD_FORKSERVER], 
                         help="Sets the python multiprocessing start method. Check https://docs.python.org/3/library/multiprocessing.html. Only use if the default method is not working properly in your OS.")

    pg_exec.add_argument('--resume', action="store_true",
                          help=("Resumes a previous emapper run, skipping results in existing output files."))

    pg_exec.add_argument('--override', action="store_true",
                        help=(
                            "Overwrites output files if they exist. "
                            "By default, execution is aborted if conflicting files are detected."))
    
    ##
    pg_input = parser.add_argument_group('Input Data Options')

    pg_input.add_argument('-i', dest="input", metavar='FASTA_FILE', type=existing_file,
                          help=f'Input FASTA file containing query sequences (proteins by default; see --itype and --translate). '
                          f'Required unless -m {SEARCH_MODE_NO_SEARCH}.')

    pg_input.add_argument('--itype', dest="itype", choices = [ITYPE_CDS, ITYPE_PROTS, ITYPE_GENOME, ITYPE_META],
                          default=ITYPE_PROTS,
                          help=f'Type of data in the input (-i) file.')
    
    pg_input.add_argument('--translate', action="store_true",
                          help=('When --itype CDS, translate CDS to proteins before search. '
                                'When --itype genome/metagenome and --genepred search, '
                                'translate predicted CDS from blastx hits to proteins.'))

    pg_input.add_argument('--annotate_hits_table', type=str, metavar='SEED_ORTHOLOGS_FILE',
                          help=f'Annotate TSV formatted table with 4 fields:'
                          f' query, hit, evalue, score. '
                          f' Usually, a .seed_orthologs file from a previous emapper.py run. '
                          f' Requires -m {SEARCH_MODE_NO_SEARCH}.')

    pg_input.add_argument('-c', '--cache', dest="cache_file", metavar='FILE', type=existing_file,
                          help=f'File containing annotations and md5 hashes of queries, to be used as cache. '
                          f'Required if -m {SEARCH_MODE_CACHE}')
        
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
                           choices = [SEARCH_MODE_DIAMOND, SEARCH_MODE_MMSEQS2, SEARCH_MODE_HMMER, SEARCH_MODE_NO_SEARCH, SEARCH_MODE_CACHE, SEARCH_MODE_NOVEL_FAMS],
                           default=SEARCH_MODE_DIAMOND,
                           help=(
                               f'{SEARCH_MODE_DIAMOND}: search seed orthologs using diamond (-i is required). '
                               f'{SEARCH_MODE_MMSEQS2}: search seed orthologs using MMseqs2 (-i is required). '
                               f'{SEARCH_MODE_HMMER}: search seed orthologs using HMMER. (-i is required). '
                               f'{SEARCH_MODE_NO_SEARCH}: skip seed orthologs search (--annotate_hits_table is required, unless --no_annot). '
                               f'{SEARCH_MODE_CACHE}: skip seed orthologs search and annotate based on cached results (-i and -c are required).'
                               f'{SEARCH_MODE_NOVEL_FAMS}: search against the novel families database (-i is required).'
                           ))

    ##
    pg_diamond_mmseqs = parser.add_argument_group('Search filtering common options')

    pg_diamond_mmseqs.add_argument('--pident', dest='pident', type=float, default=None,
                                   help=(
                                       f'Report only alignments above or equal to the given percentage of identity (0-100).'
                                       f'No effect if -m {SEARCH_MODE_HMMER}.'))
    
    pg_diamond_mmseqs.add_argument('--query_cover', dest='query_cover', type=float, default=None,
                                   help=(
                                       f'Report only alignments above or equal the given percentage of query cover (0-100).'
                                       f'No effect if -m {SEARCH_MODE_HMMER}.'))
    
    pg_diamond_mmseqs.add_argument('--subject_cover', dest='subject_cover', type=float, default=None,
                                   help=(
                                       f'Report only alignments above or equal the given percentage of subject cover (0-100).'
                                       f'No effect if -m {SEARCH_MODE_HMMER}.'))
    
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

    pg_diamond.add_argument('--sensmode', dest='sensmode', 
                            choices = SENSMODES, 
                            default=SENSMODE_SENSITIVE,
                            help=(
                                "Diamond's sensitivity mode. "
                                "Note that emapper's default is "+SENSMODE_SENSITIVE+", "
                                "which is different from diamond's default, which can "
                                "be activated with --sensmode default."
                            ))

    pg_diamond.add_argument('--dmnd_iterate', dest='dmnd_iterate', choices = [DMND_ITERATE_YES, DMND_ITERATE_NO],
                            default = DMND_ITERATE_DEFAULT,
                            help=(
                                f"--dmnd_iterate {DMND_ITERATE_YES} --> activates the --iterate option of diamond for iterative searches, "
                                f"from faster, less sensitive modes, up to the sensitivity specified with --sensmode. "
                                f"Available since diamond 2.0.11. --dmnd_iterate {DMND_ITERATE_NO} --> disables the --iterate mode. "
                            ))
        
    pg_diamond.add_argument('--matrix', dest='matrix', 
                            choices = ['BLOSUM62', 'BLOSUM90','BLOSUM80','BLOSUM50','BLOSUM45','PAM250','PAM70','PAM30'], 
                            default=None, help='Scoring matrix')

    pg_diamond.add_argument('--dmnd_frameshift', dest='dmnd_frameshift', type=int, default=None, 
                            help='Diamond --frameshift/-F option. Not used by default. Recommended by diamond: 15.')
    
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
    pg_hmmer = parser.add_argument_group('HMMER Search Options')

    pg_hmmer.add_argument('-d', '--database', dest='db', metavar='HMMER_DB_PREFIX',
                          help=('specify the target database for sequence searches. '
                                'Choose among: euk,bact,arch, or a database loaded in a server, db.hmm:host:port (see hmm_server.py)'))

    pg_hmmer.add_argument('--servers_list', dest="servers_list", metavar="FILE",
                          help="A FILE with a list of remote hmmpgmd servers. "
                                "Each row in the file represents a server, in the format 'host:port'. "
                                "If --servers_list is specified, host and port from -d option will be ignored.")
    
    pg_hmmer.add_argument('--qtype',  choices=[QUERY_TYPE_HMM, QUERY_TYPE_SEQ], default=QUERY_TYPE_SEQ,
                          help="Type of input data (-i).")

    pg_hmmer.add_argument('--dbtype', dest="dbtype",
                          choices=[DB_TYPE_HMM, DB_TYPE_SEQ], default=DB_TYPE_HMM,
                          help="Type of data in DB (-db).")

    pg_hmmer.add_argument('--usemem', action="store_true",
                          help='''Use this option to allocate the whole database (-d) in memory using hmmpgmd.
                          If --dbtype hmm, the database must be a hmmpress-ed database.
                          If --dbtype seqdb, the database must be a HMMER-format database created with esl-reformat.
                          Database will be unloaded after execution.
                          Note that this only works for HMMER based searches.
                          To load the eggnog-mapper annotation DB into memory use --dbmem.''')

    pg_hmmer.add_argument('-p', '--port', dest='port', type=int, default=DEFAULT_PORT, metavar='PORT',
                          help=('Port used to setup HMM server, when --usemem. Also used for --pfam_realign modes.'))
    
    pg_hmmer.add_argument('--end_port', dest='end_port', type=int, default=DEFAULT_END_PORT, metavar='PORT',
                          help=('Last port to be used to setup HMM server, when --usemem. Also used for --pfam_realign modes.'))

    pg_hmmer.add_argument('--num_servers', dest='num_servers', type=int, default=1, metavar="NUM_SERVERS",
                          help=("When using --usemem, specify the number of servers to fire up."
                                "Note that cpus specified with --cpu will be distributed among servers and workers."
                                " Also used for --pfam_realign modes."))
    
    pg_hmmer.add_argument('--num_workers', dest='num_workers', type=int, default=1, metavar="NUM_WORKERS",
                          help=("When using --usemem, specify the number of "
                                "workers per server (--num_servers) to fire up. "
                                "By default, cpus specified with --cpu will be "
                                "distributed among servers and workers. "
                                "Also used for --pfam_realign modes."))

    pg_hmmer.add_argument('--hmm_maxhits', dest='maxhits', type=int, default=1, metavar='MAXHITS',
                          help="Max number of hits to report (0 to report all).")

    pg_hmmer.add_argument('--report_no_hits', action="store_true",
                          help="Whether queries without hits should be included in the output table.")

    pg_hmmer.add_argument('--hmm_maxseqlen', dest='maxseqlen', type=int, default=5000, metavar='MAXSEQLEN',
                          help="Ignore query sequences larger than `maxseqlen`.")

    pg_hmmer.add_argument('--Z', dest='Z', type=float, default=40000000, metavar='DB_SIZE',
                          help='Fixed database size used in phmmer/hmmscan'
                          ' (allows comparing e-values among databases).')

    pg_hmmer.add_argument('--cut_ga', action="store_true",
                          help=("Adds the --cut_ga to hmmer commands (useful for "
                                "Pfam mappings, for example). See hmmer documentation."))

    pg_hmmer.add_argument('--clean_overlaps', dest="clean_overlaps", type=str, default=None, metavar="none|all|clans|hmmsearch_all|hmmsearch_clans",
                          help=('Removes those hits which overlap, keeping only the one with best evalue. '
                                'Use the "all" and "clans" options when performing a '
                                'hmmscan type search (i.e. domains are in the database). '
                                'Use the "hmmsearch_all" and "hmmsearch_clans" options '
                                'when using a hmmsearch type search (i.e. domains are the queries from -i file). '
                                'The "clans" and "hmmsearch_clans" and options will only '
                                'have effect for hits to/from Pfam.'))
    
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
                           ' seed orthologs will not be annotated.')

    pg_annot.add_argument('--seed_ortholog_score', default=None, type=float, metavar='MIN_SCORE',
                           help='Min bit score expected when searching for seed eggNOG ortholog.'
                           ' Queries not having a significant'
                           ' seed orthologs will not be annotated.')
    
    pg_annot.add_argument("--tax_scope", type=str, default='auto', 
                          help=("Fix the taxonomic scope used for annotation, so only speciation events from a "
                                "particular clade are used for functional transfer. "
                                "More specifically, the --tax_scope list is intersected with the "
                                "seed orthologs clades, "
                                "and the resulting clades are used for annotation based on --tax_scope_mode. "
                                "Note that those seed orthologs without clades intersecting with --tax_scope "
                                "will be filtered out, and won't annotated. "
                                "Possible arguments for --tax_scope are: "
                                "1) A path to a file defined by the user, which contains "
                                "a list of tax IDs and/or tax names. "
                                "2) The name of a pre-configured tax scope, whose source is "
                                "a file stored within the 'eggnogmapper/annotation/tax_scopes/' directory "
                                "By default, available ones are: 'auto' ('all'), 'auto_broad' ('all_broad'), "
                                "'all_narrow', 'archaea', "
                                "'bacteria', 'bacteria_broad', 'eukaryota', 'eukaryota_broad' "
                                "and 'prokaryota_broad'."
                                "3) A comma-separated list of taxonomic names and/or taxonomic IDs, "
                                "sorted by preference. "
                                "An example of list of tax IDs would be 2759,2157,2,1 for Eukaryota, "
                                "Archaea, Bacteria and root, in that order of preference. "
                                "4) 'none': do not filter out annotations based on taxonomic scope."))

    pg_annot.add_argument("--tax_scope_mode", type=str, default=TAX_SCOPE_MODE_INNER_NARROWEST,
                          help=("For a seed ortholog which passed the filter imposed by --tax_scope, "
                                "--tax_scope_mode controls which specific clade, to which the "
                                "seed ortholog belongs, will be used for annotation. "
                                f"Options: "
                                f"1) {TAX_SCOPE_MODE_BROADEST}: use the broadest clade. "
                                f"2) {TAX_SCOPE_MODE_INNER_BROADEST}: use the broadest clade "
                                "from the intersection with --tax_scope. "
                                f"3) {TAX_SCOPE_MODE_INNER_NARROWEST}: use the narrowest clade "
                                "from the intersection with --tax_scope. "
                                f"4) {TAX_SCOPE_MODE_NARROWEST}: use the narrowest clade. "
                                f"5) A taxonomic scope as in --tax_scope: use this second list "
                                "to intersect with seed ortholog clades and "
                                f"use the narrowest (as in inner_narrowest) from the intersection to annotate."))

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
    
    pg_annot.add_argument('--go_evidence', type=str, choices=('experimental', 'non-electronic', 'all'),
                          default='non-electronic',
                          help='Defines what type of GO terms should be used for annotation. '
                          'experimental = Use only terms inferred from experimental evidence. '
                          'non-electronic = Use only non-electronically curated terms')

    pg_annot.add_argument('--pfam_realign', type=str,
                          choices=(PFAM_REALIGN_NONE, PFAM_REALIGN_REALIGN, PFAM_REALIGN_DENOVO),
                          default=PFAM_REALIGN_NONE,
                          help=('Realign the queries to the PFAM domains. '
                                f'{PFAM_REALIGN_NONE} = no realignment is performed. PFAM annotation will be '
                                'that transferred as specify in the --pfam_transfer option. '
                                f'{PFAM_REALIGN_REALIGN} = queries will be realigned to the PFAM domains '
                                'found according to the --pfam_transfer option. '
                                f'{PFAM_REALIGN_DENOVO} = queries will be realigned to the whole PFAM database, '
                                'ignoring the --pfam_transfer option. '
                                f'Check hmmer options (--num_servers, --num_workers, --port, --end_port) '
                                'to change how the hmmpgmd server is run. '))
    
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
        print_taxa()
        sys.exit(0)

    if args.cpu == 0:
        args.cpu = multiprocessing.cpu_count()
    multiprocessing.set_start_method(args.mp_start_method)

    if args.resume == True and args.override == True:
        parser.error('Only one of --resume or --override is allowed.')        

    # Gene prediction
    if args.training_genome is not None and args.training_file is None:
        parser.error('"--training_genome requires --training_file"')

    if args.training_genome is None and args.training_file is not None:
        if not os.path.isfile(args.training_file):
            parser.error('"--training_file must point to an existing file, if no --training_genome is provided."')
    
    # Search modes
    if args.mode == SEARCH_MODE_DIAMOND or args.mode == SEARCH_MODE_NOVEL_FAMS:
        # dmnd_db = args.dmnd_db if args.dmnd_db else get_eggnog_dmnd_db(args.dmnd_db, args.mode, get_data_path())
        dmnd_db = get_eggnog_dmnd_db(args.dmnd_db, args.mode, get_data_path())
        if not pexists(dmnd_db):
            print(colorify('DIAMOND database %s not present. Use download_eggnog_database.py to fetch it' % dmnd_db, 'red'))
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
            print(colorify('MMseqs2 database %s not present. Use download_eggnog_database.py to fetch it' % mmseqs_db, 'red'))
            raise EmapperException()

        if not args.input:
            parser.error('An input fasta file is required (-i)')

        # Output file required
        if not args.output:
            parser.error('An output project name is required (-o)')

        if args.annotate_hits_table is not None:
            print(colorify(f"--annotate_hits_table will be ignored, due to -m {SEARCH_MODE_MMSEQS2}", 'blue'))
            args.annotate_hits_table = None
            
    elif args.mode == SEARCH_MODE_HMMER:
        
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
        if args.decorate_gff != DECORATE_GFF_NONE:
            print(colorify("WARNING: no GFF will be created for cache-based annotations. It is not implemented yet, sorry.", 'red'))
                
        if args.no_annot == True:
            parser.error(f'Cache mode (-m {SEARCH_MODE_CACHE}) should be used to annotate.')
            
    elif args.mode == SEARCH_MODE_NO_SEARCH:
        if args.no_annot == False and not args.annotate_hits_table:
            parser.error(f'No search mode (-m {SEARCH_MODE_NO_SEARCH}) requires a hits table to annotate (--annotate_hits_table FILE.seed_orthologs)')
        if args.md5 == True and args.input is None:
            parser.error(f'--md5 requires an input FASTA file (-i FASTA).')            
            
    else:
        parser.error(f'unrecognized search mode (-m {args.mode})')


    # Search thresholds
    args.dmnd_evalue = args.mmseqs_evalue = args.hmm_evalue = args.evalue
    args.dmnd_score = args.mmseqs_score = args_hmm_score = args.score
    args.qcov = args.query_cover
    
    # Annotation options
    if args.no_annot == False or args.report_orthologs == True:
        if not pexists(get_eggnogdb_file()) and args.mode != SEARCH_MODE_NOVEL_FAMS:
            print(colorify('Annotation database data/eggnog.db not present. Use download_eggnog_database.py to fetch it', 'red'))
            raise EmapperException()

        args.tax_scope_ids = parse_tax_scope(args.tax_scope)
        
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
        
        n, elapsed_time = emapper.run(args, args.input, args.annotate_hits_table, args.cache_file)

        elapsed_time = time.time() - start_time

        addons = [args.mode, args.genepred]
        # when using novel_fams, diamond is also used
        if args.mode == SEARCH_MODE_NOVEL_FAMS:
            addons.append(SEARCH_MODE_DIAMOND)
            
        print(get_citation(addons))
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
