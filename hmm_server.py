#!/usr/bin/env python3

import os, sys, time, traceback
import argparse, multiprocessing

# get the path of this script and add it to the "pythonpath"
SCRIPT_PATH = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.insert(0, SCRIPT_PATH)

from eggnogmapper.emapperException import EmapperException
from eggnogmapper.common import TIMEOUT_LOAD_SERVER
from eggnogmapper.utils import colorify

from eggnogmapper.search.hmmer.hmmer_search import DB_TYPE_SEQ, DB_TYPE_HMM, SCANTYPE_MEM
from eggnogmapper.search.hmmer.hmmer_setup import setup_custom_db, DEFAULT_PORT, DEFAULT_END_PORT
from eggnogmapper.search.hmmer.hmmer_server import load_server, server_functional, create_servers

__description__ = ('A server for HMMER3 in-memory searches')
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
    pg_server = parser.add_argument_group('HMM Server Options')

    pg_server.add_argument('-d', '--database', dest='db', metavar='HMMER_DB_PREFIX',
                       help=('specify the target database for sequence searches. '
                            'Choose among: euk,bact,arch, host:port, or a local "hmmpress-ed" database'))

    pg_server.add_argument('--dbtype', dest="dbtype",
                       choices=[DB_TYPE_HMM, DB_TYPE_SEQ], default=DB_TYPE_HMM,
                       help="Type of data in DB (-db). "
                          f"Default: {DB_TYPE_HMM}")

    pg_server.add_argument('-p', '--port', dest='port', type=int, default=DEFAULT_PORT, metavar='PORT',
                          help=('Port used by clients to connect to this HMM master server'))

    pg_server.add_argument('--end_port', dest='end_port', type=int, default=DEFAULT_END_PORT, metavar='PORT',
                          help=('Last port to be used to setup HMM server, when --usemem'))

    pg_server.add_argument('-w', '--wport', dest='wport', type=int, default=53001, metavar='PORT',
                          help=('Port used by workers to connect to this HMM master server'))

    pg_server.add_argument('--num_servers', dest='num_servers', type=int, default=1, metavar="NUM_SERVERS",
                          help="When using --usemem, specify the number of servers to fire up."
                          " By default, cpus specified with --cpu will be distributed among servers and workers.")
    
    pg_server.add_argument('--num_workers', dest='num_workers', type=int, default=1, metavar="NUM_WORKERS",
                          help="When using --usemem, specify the number of workers per server (--num_servers) to fire up."
                          " By default, cpus specified with --cpu will be distributed among servers and workers.")

    pg_server.add_argument('-o', '--output_servers_list', dest="output_servers_list", type=str, default=None, metavar="FILE",
                           help='Output the list of running servers to FILE.')
        
    return parser


def parse_args(parser):
    
    args = parser.parse_args()

    if args.version:
        print(get_version())
        sys.exit(0)

    if args.cpu == 0:
        args.cpu = multiprocessing.cpu_count()

    # cpus per worker
    total_workers = args.num_workers * args.num_servers
    if args.cpu < total_workers:
        parser.error(f"Less cpus ({args.cpu}) than total workers ({total_workers}) were specified.")
    if args.cpu % total_workers != 0:
        parser.error(f"Number of cpus ({args.cpu}) must be a multiple of total workers ({total_workers}).")        

    args.cpus_per_worker = int(args.cpu / total_workers)
    sys.stderr.write(f"CPUs per worker: {args.cpus_per_worker}\n")
        
    # Hmmer database
    # NOTE: hmmer database format, name and checking if exists is done within hmmer module
    if not args.db:
        parser.error('The HMMER server requires a target database (-d, --database).')

    if os.path.exists(args.output_servers_list):
        parser.error(f"File {args.output_servers_list} already exists, and won't be overwritten."
                     "Please, remove it and run again to create it.")
    
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
        print('# hmm_server.py ', ' '.join(sys.argv[1:]))

        dbpath, host, idmap_file = setup_custom_db(args.db, scantype = SCANTYPE_MEM, dbtype = args.dbtype)

        host = 'localhost'
        port = args.port
        end_port = args.end_port
        wport = args.wport

        dbpath, host, port, servers = create_servers(args.dbtype, dbpath, host, port, end_port,
                                                     args.num_servers, args.num_workers, args.cpus_per_worker)

        print(colorify("All servers ready and listening", 'green'))
        if args.output_servers_list is not None:
            print(f"Creating servers list file: {args.output_servers_list}")
            with open(args.output_servers_list, 'w') as outfn:
                for server in servers:
                    print(f"{server[0]}:{server[1]}", file=outfn)
            print(f"File {args.output_servers_list} created successfully.")
            
            print(colorify(f"Use `emapper.py (-d db:host:port or --servers_list {args.output_servers_list}) to search against these servers", 'lblue'))
        else:
            print(colorify("Use `emapper.py (-d db:host:port or --servers_list FILE) to search against these servers", 'lblue'))                    
        
        while True:
            time.sleep(10)
        raise emapperException("Servers stopped.")
        
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
        print()
        print(get_citation())
        print('Total time: %g secs' % (time.time()-_total_time))        

## END
