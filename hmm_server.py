#!/usr/bin/env python3

import os, sys, time, traceback
import argparse, multiprocessing

# get the path of this script and add it to the "pythonpath"
SCRIPT_PATH = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.insert(0, SCRIPT_PATH)

from eggnogmapper.emapperException import EmapperException
from eggnogmapper.common import TIMEOUT_LOAD_SERVER
from eggnogmapper.utils import colorify

from eggnogmapper.search.hmmer_search import DB_TYPE_SEQ, DB_TYPE_HMM, SCANTYPE_MEM
from eggnogmapper.search.hmmer_setup import setup_custom_db
from eggnogmapper.search.hmmer_server import load_server, server_functional

__description__ = ('A server for HMMER3 in-memory searches')
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
    pg_server = parser.add_argument_group('HMM Server Options')

    pg_server.add_argument('-d', '--database', dest='db', metavar='HMMER_DB_PREFIX',
                       help=('specify the target database for sequence searches. '
                            'Choose among: euk,bact,arch, host:port, or a local "hmmpress-ed" database'))

    pg_server.add_argument('--dbtype', dest="dbtype",
                       choices=[DB_TYPE_HMM, DB_TYPE_SEQ], default=DB_TYPE_HMM,
                       help="Type of data in DB (-db). "
                          f"Default: {DB_TYPE_HMM}")

    pg_server.add_argument('-p', '--port', dest='port', type=int, default=53000, metavar='PORT',
                          help=('Port used by clients to connect to this HMM master server'))

    pg_server.add_argument('-w', '--wport', dest='wport', type=int, default=53001, metavar='PORT',
                          help=('Port used by workers to connect to this HMM master server'))

    pg_server.add_argument('--is_worker', action="store_true",
                           help='In addition to the master, create a worker also in this host.')
        
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
        parser.error('The HMMER server requires a target database (-d, --database).')
    
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

        dbname, dbpath, host, port, end_port, idmap_file = setup_custom_db(args.db, scantype = SCANTYPE_MEM, dbtype = args.dbtype)

        host = 'localhost'
        port = args.port
        wport = args.wport
            
        print(colorify(f"Loading server at {host}:{port}; workers port:{wport}", 'green'))
        
        dbpath, master_db, worker_db = load_server(dbpath, port, wport, args.cpu, dbtype = args.dbtype, is_worker = args.is_worker)
        
        ready = False
        for _ in range(TIMEOUT_LOAD_SERVER):
            print(f"Waiting for server to become ready at {host}:{port} ...")
            time.sleep(1)
            if master_db.is_alive() and (not args.is_worker or worker_db.is_alive()):
                if not args.is_worker:
                    print(colorify("master is UP", 'green'))
                    break
                else: # worker_db.is_alive
                    if server_functional(host, port, args.dbtype):
                        print(colorify("master and worker are UP", 'green'))
                        break
                
            elif not master_db.is_alive():
                master_db.terminate()
                master_db.join()
                print(colorify("master not alive"), 'red')
                break
            
            elif args.is_worker and not worker_db.is_alive():
                worker_db.terminate()
                worker_db.join()
                print(colorify("worker not alive"), 'red')
                break

        print(colorify("Server ready listening at %s:%s and using %d CPU cores" % (host, port, args.cpu), 'green'))
        print(colorify("Use `emapper.py -d %s:%s:%s (...)` to search against this server" % (args.db, host, port), 'lblue'))
        
        while True:
            time.sleep(10)
        raise emapperException("Server {db}:{host}:{port} stopped.")
        
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
