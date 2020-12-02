#!/usr/bin/env python3

import os, sys, time, traceback
import argparse, multiprocessing

# get the path of this script and add it to the "pythonpath"
SCRIPT_PATH = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.insert(0, SCRIPT_PATH)

from eggnogmapper.common import existing_file, existing_dir, set_data_path, pexists
from eggnogmapper.utils import colorify

from eggnogmapper.emapperException import EmapperException
from eggnogmapper.hmm_mapper import HmmMapper

from eggnogmapper.search.hmmer.hmmer_search import QUERY_TYPE_SEQ, QUERY_TYPE_HMM, DB_TYPE_SEQ, DB_TYPE_HMM
from eggnogmapper.search.hmmer.hmmer_setup import DEFAULT_PORT, DEFAULT_END_PORT


__description__ = ('A program wrapping HMM in-memory searches')
__author__ = 'Jaime Huerta Cepas'
__license__ = "GPL v2"

def create_arg_parser():
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--version', action='store_true',
                        help="show version and exit.")

    ##
    pg_exec = parser.add_argument_group('Execution Options')
    
    pg_exec.add_argument('--cpu', type=int, default=1, metavar='NUM_CPU',
                        help="Number of CPUs to be used. --cpu 0 to run with all available CPUs. Default: 2")
    
    ##
    pg_input = parser.add_argument_group('Input Data Options')

    pg_input.add_argument('-i', dest="input", metavar='FASTA_FILE', type=existing_file,
                          help=f'Input with queries. Either a FASTA file with sequences (proteins by default; see --translate)'
                          ' or a HMM file with profiles (--qtype hmm)')

    pg_input.add_argument('--translate', action="store_true",
                          help='Assume input sequences are CDS instead of proteins (it has effect only if --qtype seq, and also when -d is a plain FASTA file)')

    ##
    pg_hmmer = parser.add_argument_group('HMMER Search Options')

    pg_hmmer.add_argument('-d', '--database', dest='db', metavar='DB_PATH',
                       help=('specify the target database for sequence searches. '
                            'Choose among: db:host:port, or a local database.'))

    pg_hmmer.add_argument('--servers_list', dest="servers_list", metavar="FILE",
                          help="A FILE with a list of remote hmmpgmd servers. "
                                "Each row in the file represents a server, in the format 'host:port'. "
                                "If --servers_list is specified, host and port from -d option will be ignored.")
        
    pg_hmmer.add_argument('--qtype',  choices=[QUERY_TYPE_HMM, QUERY_TYPE_SEQ], default=QUERY_TYPE_SEQ,
                       help="Type of input data (-i). "
                          f"Default: {QUERY_TYPE_SEQ}")

    pg_hmmer.add_argument('--dbtype', dest="dbtype",
                       choices=[DB_TYPE_HMM, DB_TYPE_SEQ], default=DB_TYPE_HMM,
                       help="Type of data in DB (-d). "
                          f"Default: {DB_TYPE_HMM}")

    pg_hmmer.add_argument('--usemem', action="store_true",
                    help='''Use this option to allocate the whole database (-d) in memory using hmmpgmd.
                    If --dbtype hmm, the database must be a hmmpress-ed database.
                    If --dbtype seqdb, the database must be a HMMER-format database created with esl-reformat.
                    Database will be unloaded after execution.''')

    pg_hmmer.add_argument('-p', '--port', dest='port', type=int, default=DEFAULT_PORT, metavar='PORT',
                          help=('Port used to setup HMM server, when --usemem'))
    
    pg_hmmer.add_argument('--end_port', dest='end_port', type=int, default=DEFAULT_END_PORT, metavar='PORT',
                          help=('Last port to be used to setup HMM server, when --usemem'))

    pg_hmmer.add_argument('--num_servers', dest='num_servers', type=int, default=1, metavar="NUM_SERVERS",
                          help="When using --usemem, specify the number of servers to fire up."
                          " By default, cpus specified with --cpu will be distributed among servers and workers.")
    
    pg_hmmer.add_argument('--num_workers', dest='num_workers', type=int, default=1, metavar="NUM_WORKERS",
                          help="When using --usemem, specify the number of workers per server (--num_servers) to fire up."
                          " By default, cpus specified with --cpu will be distributed among servers and workers.")

    pg_hmmer.add_argument('--hmm_maxhits', dest='maxhits', type=int, default=1, metavar='MAXHITS',
                        help="Max number of hits to report (0 to report all). Default=1.")

    pg_hmmer.add_argument('--report_no_hits', action="store_true",
                        help="Whether queries without hits should be included in the output table.")

    pg_hmmer.add_argument('--hmm_maxseqlen', dest='maxseqlen', type=int, default=5000, metavar='MAXSEQLEN',
                        help="Ignore query sequences larger than `maxseqlen`. Default=5000")
        
    pg_hmmer.add_argument('--hmm_evalue', dest='evalue', default=None, type=float, metavar='MIN_E-VALUE',
                          help="E-value threshold. For example, -hmm_evalue 0.001. Default=10")

    pg_hmmer.add_argument('--hmm_score', dest='score', default=None, type=float, metavar='MIN_SCORE',
                          help="Bit score threshold. For example, --hmm_score 20. Default=None")

    pg_hmmer.add_argument('--hmm_qcov', dest='qcov', type=float, metavar='MIN_QCOV',
                        help="min query coverage (from 0 to 1). Default=(disabled)")

    pg_hmmer.add_argument('--Z', dest='Z', type=float, default=40000000, metavar='DB_SIZE',
                        help='Fixed database size used in phmmer/hmmscan'
                        ' (allows comparing e-values among databases). Default=40,000,000')

    pg_hmmer.add_argument('--cut_ga', action="store_true",
                          help="Adds the --cut_ga to hmmer commands (useful for Pfam mappings, for example). See hmmer documentation.")

    pg_hmmer.add_argument('--clean_overlaps', dest="clean_overlaps", type=str, default=None, metavar="none|all|clans|hmmsearch_all|hmmsearch_clans",
                          help='Removes those hits which overlap, keeping only the one with best evalue. '
                          'Use the "all" and "clans" options when performing a hmmscan type search (i.e. domains are in the database). '
                          'Use the "hmmsearch_all" and "hmmsearch_clans" options when using a hmmsearch type search (i.e. domains are the queries from -i file). '
                          'The "clans" and "hmmsearch_clans" and options will only have effecto for hits to/from Pfam.')

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
                        help="Resumes a previous execution skipping reported hits in the output file.")
        
    pg_out.add_argument('--override', action="store_true",
                    help="Overwrites output files if they exist.")

    pg_out.add_argument("--temp_dir", default=os.getcwd(), type=existing_dir, metavar='DIR',
                    help="Where temporary files are created. Better if this is a local disk.")

    pg_out.add_argument('--no_file_comments', action="store_true",
                        help="No header lines nor stats are included in the output files")
        
    return parser

def parse_args(parser):
    
    args = parser.parse_args()

    if args.version:
        print(get_version())
        sys.exit(0)

    if args.cpu == 0:
        args.cpu = multiprocessing.cpu_count()

    if args.usemem == True:
        total_workers = args.num_workers * args.num_servers
        if args.cpu < total_workers:
            parser.error(f"Less cpus ({args.cpu}) than total workers ({total_workers}) were specified.")
        if args.cpu % total_workers != 0:
            parser.error(f"Number of cpus ({args.cpu}) must be a multiple of total workers ({total_workers}).")        

        args.cpus_per_worker = int(args.cpu / total_workers)
        sys.stderr.write(f"CPUs per worker: {args.cpus_per_worker}\n")
    else:
        args.cpus_per_worker = args.cpu    

    # Required files
    if not args.input:
        parser.error('An input file is required (-i)')
        
    if not args.output:
        parser.error('An output project name is required (-o)')
        
    if not args.db:
        parser.error('hmm_mapper requires a target database (-d, --database).')

    if args.clean_overlaps is not None:
        if args.clean_overlaps == "none":
            args.clean_overlaps = None
            
    return args

def get_version():
    return "1.0"

def get_call_info():
    text = []
    text.append('# ' + time.ctime())
    text.append('# ' + get_version())
    text.append('# ' + ' '.join(sys.argv))
    text.append('#')
    return '\n'.join(text)

def get_citation():
    return __author__+" "+__license__+" : "+__description__

if __name__ == "__main__":

    parser = create_arg_parser()
    args = parse_args(parser)

    _total_time = time.time()
    try:
        
        print('# ', get_version())
        print('# hmm_mapper.py ', ' '.join(sys.argv[1:]))

        args.call_info = get_call_info()
        hmm_mapper = HmmMapper(args.output, args.output_dir, args.scratch_dir, args.resume, args.override)
        hmm_mapper.run(args, args.input)
        
    except EmapperException as ee:
        print(ee)
        sys.exit(1)
    except Exception:
        traceback.print_exc()
        sys.exit(1)
    else:
        print("FINISHED")
        sys.exit(0)
    finally:
        print(get_citation())
        print('Total time: %g secs' % (time.time()-_total_time))
        
## END
