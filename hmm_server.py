#!/usr/bin/env python3

import os, sys, time, traceback
import argparse, multiprocessing

# get the path of this script and add it to the "pythonpath"
SCRIPT_PATH = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.insert(0, SCRIPT_PATH)

from eggnogmapper.emapperException import EmapperException
from eggnogmapper.common import existing_dir
from eggnogmapper.utils import colorify
from eggnogmapper.search.hmmer_search import QUERY_TYPE_SEQ, QUERY_TYPE_HMM, DB_TYPE_SEQ, DB_TYPE_HMM, SCANTYPE_MEM
from eggnogmapper.search.hmmer_setup import setup_custom_db, start_server


__description__ = ('A program serving HMM in-memory searches')
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
    pg_hmmer = parser.add_argument_group('HMMER Search Options')

    pg_hmmer.add_argument('-d', '--database', dest='db', metavar='HMMER_DB_PREFIX',
                       help=('specify the target database for sequence searches. '
                            'Choose among: euk,bact,arch, host:port, or a local "hmmpress-ed" database'))

    pg_hmmer.add_argument('--dbtype', dest="dbtype",
                       choices=[DB_TYPE_HMM, DB_TYPE_SEQ], default=DB_TYPE_HMM,
                       help="Type of data in DB (-db). "
                          f"Default: {DB_TYPE_HMM}")
        
    return parser

def parse_args(parser):
    
    args = parser.parse_args()

    if args.version:
        print(get_version())
        sys.exit(0)

    if args.cpu == 0:
        args.cpu = multiprocessing.cpu_count()

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

        dbname, dbpath, host, port, end_port, idmap_file = setup_custom_db(args.db, scantype = SCANTYPE_MEM)

        print(f"DB setup: {dbpath} --> {host}:{port}-{end_port}")
        
        dbpath, host, port = start_server(dbpath, host, port, end_port, args.cpu, args.dbtype)

        print(colorify("Server ready listening at %s:%s and using %d CPU cores" % (host, port, args.cpu), 'green'))
        print(colorify("Use `emapper.py -d %s:%s:%s (...)` to search against this server" % (args.db, host, port), 'lblue'))
        while True:
            time.sleep(10)
        raise emapperException("Server {db}:{host}:{port} stopped.")
    
        print(get_citation())
        print('Total time: %g secs' % (time.time()-_total_time))
        
    except EmapperException as ee:
        print(ee)
        sys.exit(1)
    except Exception:
        traceback.print_exc()
        sys.exit(1)
    else:
        print("FINISHED")
        sys.exit(0)

## END
