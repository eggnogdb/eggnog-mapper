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

from eggnogmapper.search.hmmer_search import QUERY_TYPE_SEQ, QUERY_TYPE_HMM, DB_TYPE_SEQ, DB_TYPE_HMM



__description__ = ('A program wrapping HMM in-memory searches')
__author__ = 'Jaime Huerta Cepas'
__license__ = "GPL v2"

def create_arg_parser():
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--version', action='store_true',
                        help="show version and exit.")

    ##
    pg_exec = parser.add_argument_group('Execution Options')
    
    pg_exec.add_argument('--cpu', type=int, default=2, metavar='NUM_CPU',
                        help="Number of CPUs to be used. --cpu 0 to run with all available CPUs. Default: 2")
    
    ##
    pg_input = parser.add_argument_group('Input Data Options')

    pg_input.add_argument('-i', dest="input", metavar='FASTA_FILE', type=existing_file,
                          help=f'Input FASTA file containing query sequences (proteins by default; see --translate).')

    pg_input.add_argument('--translate', action="store_true",
                          help='Assume input sequences are CDS instead of proteins')

    ##
    pg_hmmer = parser.add_argument_group('HMMER Search Options')

    pg_hmmer.add_argument('-d', '--database', dest='db', metavar='HMMER_DB_PREFIX',
                       help=('specify the target database for sequence searches. '
                            'Choose among: euk,bact,arch, host:port, or a local "hmmpress-ed" database'))
        
    pg_hmmer.add_argument('--qtype',  choices=[QUERY_TYPE_HMM, QUERY_TYPE_SEQ], default=QUERY_TYPE_SEQ,
                       help="Type of input data (-i). "
                          f"Default: {QUERY_TYPE_SEQ}")

    pg_hmmer.add_argument('--dbtype', dest="dbtype",
                       choices=[DB_TYPE_HMM, DB_TYPE_SEQ], default=DB_TYPE_HMM,
                       help="Type of data in DB (-db). "
                          f"Default: {DB_TYPE_HMM}")

    pg_hmmer.add_argument('--usemem', action="store_true",
                    help='''If a local "hmmpress-ed" database is provided as target using --database,
                    --usemem will allocate the whole database in memory using hmmpgmd.
                    Database will be unloaded after execution.''')

    pg_hmmer.add_argument('--hmm_maxhits', dest='maxhits', type=int, default=1, metavar='MAXHITS',
                        help="Max number of hits to report (0 to report all). Default=1.")

    pg_hmmer.add_argument('--hmm_maxseqlen', dest='maxseqlen', type=int, default=5000, metavar='MAXSEQLEN',
                        help="Ignore query sequences larger than `maxseqlen`. Default=5000")
        
    pg_hmmer.add_argument('--hmm_evalue', dest='evalue', default=0.001, type=float, metavar='MIN_E-VALUE',
                        help="E-value threshold. Default=0.001")

    pg_hmmer.add_argument('--hmm_score', dest='score', default=20, type=float, metavar='MIN_SCORE',
                        help="Bit score threshold. Default=20")

    pg_hmmer.add_argument('--hmm_qcov', dest='qcov', type=float, metavar='MIN_QCOV',
                        help="min query coverage (from 0 to 1). Default=(disabled)")

    pg_hmmer.add_argument('--Z', dest='Z', type=float, default=40000000, metavar='DB_SIZE',
                        help='Fixed database size used in phmmer/hmmscan'
                        ' (allows comparing e-values among databases). Default=40,000,000')

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

    if not args.input:
        parser.error('An input file is required (-i)')

    # Output file required
    if not args.output:
        parser.error('An output project name is required (-o)')

    # Hmmer database
    # NOTE: hmmer database format, name and checking if exists is done within hmmer module
    if not args.db:
        parser.error('HMMER mode requires a target database (-d, --database).')
    
    return args

def get_version():
    return "1.0"

def get_citation():
    return __author__+" "+__license__+" : "+__description__

if __name__ == "__main__":

    parser = create_arg_parser()
    args = parse_args(parser)

    _total_time = time.time()
    try:
        
        print('# ', get_version())
        print('# hmm_mapper.py ', ' '.join(sys.argv[1:]))

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
