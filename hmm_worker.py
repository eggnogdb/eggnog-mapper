#!/usr/bin/env python3

import os, sys, time, traceback
import argparse, multiprocessing

# get the path of this script and add it to the "pythonpath"
SCRIPT_PATH = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.insert(0, SCRIPT_PATH)

from eggnogmapper.emapperException import EmapperException
from eggnogmapper.common import TIMEOUT_LOAD_SERVER
from eggnogmapper.utils import colorify

from eggnogmapper.search.hmmer.hmmer_server import load_worker
from eggnogmapper.search.hmmer.hmmer_setup import DEFAULT_PORT

__description__ = ('A worker for HMMER3 in-memory searches')
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
    pg_master = parser.add_argument_group('HMM Master Server Options')
    
    pg_master.add_argument('-@', '--host', dest='host', metavar='HOST',
                       help=('IP address or hostname of HMM master server'))

    pg_master.add_argument('-p', '--port', type=int, dest='port', default=DEFAULT_PORT, metavar='PORT',
                          help=('Port used to connect to the HMM master server'))
        
    return parser

def parse_args(parser):
    
    args = parser.parse_args()

    if args.version:
        print(get_version())
        sys.exit(0)

    if args.cpu == 0:
        args.cpu = multiprocessing.cpu_count()
    
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
        print('# hmm_worker.py ', ' '.join(sys.argv[1:]))

        worker_db = None
    
        print(colorify(f"Loading worker at localhost, port {args.port}, connecting to {args.host}", 'green'))
        worker_db = load_worker(args.host, args.port, args.cpu)
    
        ready = False
        for _ in range(TIMEOUT_LOAD_SERVER):
            print(f"Waiting for worker to become ready at localhost:{args.port} ...")
            time.sleep(1)
            if worker_db.is_alive():
                break
            else:
                worker_db.terminate()
                worker_db.join()
                print(colorify("worker not alive"), 'red')
                break
        
        print(colorify("Worker of master %s ready listening at localhost:%s and using %d CPU cores" % (args.host, args.port, args.cpu), 'lblue'))
        while True:
            time.sleep(10)
        raise emapperException("Worker localhost:{port} connected to {host} stopped.")
        
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
