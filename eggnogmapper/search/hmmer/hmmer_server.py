##
## JHCepas

import sys
import os
import time
import subprocess
from multiprocessing import Process
import signal
import traceback

from ...common import HMMPGMD, TIMEOUT_LOAD_SERVER
from ...utils import colorify

from .hmmer_search_hmmpgmd import get_hits
from .hmmer_search import DB_TYPE_SEQ, DB_TYPE_HMM, QUERY_TYPE_SEQ, QUERY_TYPE_HMM

CHILD_PROC = None
MASTER = None
WORKERS = None


##
def check_servers(dbtype, qtype, dbpath, host, port, servers_list):
    # Get list of servers
    servers = []
    functional = 0
    if servers_list is not None:
        with open(servers_list, 'r') as servers_fn:
            for line in servers_fn:
                host, port = map(str.strip, line.split(":"))
                port = int(port)
                servers.append([host, port, -1, -1]) # set -1 to master and worker PIDs, since they are not needed here
                if server_functional(host, port, dbtype, qtype):
                    functional += 1
                else:
                    print(colorify("warning: eggnog-mapper server not found at %s:%s" % (host, port), 'orange'))
                    
    else:
        servers = [[host, port, -1, -1]] # set -1 to master and worker PIDs, since they are not needed here
        if server_functional(host, port, dbtype, qtype):
            functional += 1
        else:
            print(colorify("eggnog-mapper server not found at %s:%s" % (host, port), 'red'))
            exit(1)
            
    if functional == 0:
        print(colorify("No functional server was found", 'red'))
        exit(1)

    return dbpath, host, port, servers


##
def create_servers(dbtype, dbpath, host, port, end_port, num_servers, num_workers, cpus_per_worker, silent = False):
    if silent == False:
        print(f"create_servers: {dbtype}:{dbpath}:{host}:{port}-{end_port}")
    servers = []
    sdbpath = dbpath
    shost = host
    sport = port
    MAX_CREATE_SERVER_FAILS = 3 # max number of cummulative fails creating servers before aborting creating new ones
    fails = 0
    
    for num_server in range(num_servers):
        if sport >= end_port:
            printf(colorify(f"start port ({sport}) equal or greater than end port ({end_port})", 'red'))
            break

        if silent == False:
            print(f"Creating server number {num_server+1}/{num_servers}")
        try:
            sdbpath, shost, sport, master_pid, workers_pids = start_server(dbpath, host, sport, end_port, cpus_per_worker, num_workers, dbtype, silent = silent)
            servers.append((sdbpath, sport, master_pid, workers_pids))
            fails = 0
        except Exception as e:
            fails += 1
            print(colorify(f"Could not create server number {num_server+1}/{num_servers}. Fails: {fails}", 'red'))
            if fails >= MAX_CREATE_SERVER_FAILS:
                break
            
        sport = sport + 2
        
    dbpath = sdbpath
    host = shost
    port = sport

    if silent == False:
        print(f"Created {len(servers)} out of {num_servers}")
    if len(servers) == 0:
        raise Exception("Could not create hmmpgmd servers")

    return dbpath, host, port, servers


##
def start_server(dbpath, host, port, end_port, cpus_per_worker, num_workers, dbtype, qtype = QUERY_TYPE_SEQ, silent = False):
    master_db = worker_db = workers = None
    MAX_PORTS_TO_TRY = 3
    ports_tried = 0
    ready = False
    for try_port in range(port, end_port, 2):
        if silent == False:
            print(colorify("Loading server at localhost, port %s-%s" %
                           (try_port, try_port + 1), 'lblue'))

        dbpath, master_db, workers = load_server(dbpath, try_port, try_port + 1,
                                                 cpus_per_worker, num_workers = num_workers, dbtype = dbtype,
                                                 silent = silent)
        port = try_port
        # ready = False
        if silent == False:
            print(f"Waiting for server to become ready at {host}:{port} ...")
        for attempt in range(TIMEOUT_LOAD_SERVER):
            time.sleep(attempt+1)
            if not master_db.is_alive() or not any([worker_db.is_alive() for worker_db in workers]):
                master_db.terminate()
                master_db.join()
                for worker_db in workers:
                    worker_db.terminate()
                    worker_db.join()                        
                break
            elif server_functional(host, port, dbtype, qtype):
                if silent == False:
                    print(f"Server ready at {host}:{port}")
                ready = True
                break
            else:
                if silent == False:
                    sys.stdout.write(".")
                    sys.stdout.flush()

        ports_tried += 1
        if ready:
            dbpath = host
            break
        else:
            if ports_tried >= MAX_PORTS_TO_TRY:
                raise Exception("Could not start server after trying {ports_tried} ports.")
            
    if ready == False:
        raise Exception("Could not start server.")        

    return dbpath, host, port, master_db, workers


##
def server_up(host, port):
    import socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex((host, port))
    sock.close()
    if result == 0:
        return True
    else:
        return False

def server_functional(host, port, dbtype = DB_TYPE_HMM, qtype = QUERY_TYPE_SEQ):
    if server_up(host, port):
        try:
            if qtype == QUERY_TYPE_SEQ:
                get_hits("test", "TESTSEQ", host, port, dbtype, qtype=qtype)
            elif qtype == QUERY_TYPE_HMM:

                testhmm = ""
                with open("tests/fixtures/hmmer_custom_dbs/bact.short.hmm", 'r') as hmmfile:
                    for line in hmmfile:
                        testhmm += line
                        
                get_hits("test", testhmm, host, port, dbtype, qtype=qtype)
        except Exception as e:
            # traceback.print_exc()
            # print(e)
            return False
        else:
            return True
    # else:
    #     print(colorify("Server is still down", 'red'))
    return False

def safe_exit(a, b):
    if CHILD_PROC:
        CHILD_PROC.kill()
    sys.exit(0)

def load_worker(master_host, worker_port, cpu, output=None, silent = False):
    global CHILD_PID, WORKERS
    if not output:
        OUT = open(os.devnull, 'w')
    else:
        OUT = output
        
    signal.signal(signal.SIGINT, safe_exit)
    signal.signal(signal.SIGTERM, safe_exit)
              
    def start_worker():
        cmd = HMMPGMD + f' --worker {master_host} --wport {worker_port} --cpu {cpu}'
        if silent == False:
            print(colorify(f"Loading worker: {cmd}", 'orange'))
        CHILD_PROC = subprocess.Popen(cmd.split(), shell=False, stderr=OUT, stdout=OUT)
        while 1:
            time.sleep(60)
            
    worker = Process(target=start_worker)
    worker.start()
    WORKERS = [worker]
    
    return worker

def load_server(dbpath, client_port, worker_port, cpus_per_worker, num_workers=1, output=None, dbtype=DB_TYPE_HMM, is_worker = True, silent = False):
    global CHILD_PID, MASTER, WORKERS
    if not output:
        OUT = open(os.devnull, 'w')
    else:
        OUT = output
        
    signal.signal(signal.SIGINT, safe_exit)
    signal.signal(signal.SIGTERM, safe_exit)

    def start_master():
        cmd = HMMPGMD + f' --master --cport {client_port} --wport {worker_port} --{dbtype} {dbpath}'
        if silent == False:
            print(colorify(f"Loading master: {cmd}", 'orange'))
        CHILD_PROC = subprocess.Popen(cmd.split(), shell=False, stderr=OUT, stdout=OUT)
        while 1:
            time.sleep(60)

    def start_worker():
        cmd = HMMPGMD + f' --worker localhost --wport {worker_port} --cpu {cpus_per_worker}'
        if silent == False:
            print(colorify(f"Loading worker: {cmd}", 'orange'))
        CHILD_PROC = subprocess.Popen(cmd.split(), shell=False, stderr=OUT, stdout=OUT)
        while 1:
            time.sleep(60)

    if silent == False:
        print(f"Creating hmmpgmd server at port {client_port} ...")
    MASTER = Process(target=start_master)
    MASTER.start()

    if is_worker == True and num_workers > 0:
        if silent == False:
            print(f"Creating hmmpgmd workers ({num_workers}) at port {worker_port} ...")
        WORKERS = []
        for i in range(num_workers):
            worker = Process(target=start_worker)
            worker.start()
            WORKERS.append(worker)
            
    return dbpath, MASTER, WORKERS


def shutdown_server_by_pid(MASTER, WORKERS):

    import psutil
    
    try:
        # This is killing THIS python script also, and is UNIX dependent
        # os.killpg(os.getpgid(WORKER.pid), signal.SIGTERM)

        for worker in WORKERS:
            parent = psutil.Process(worker.pid)
            for child in parent.children(recursive=True):  # or parent.children() for recursive=False
                child.kill()
            parent.kill()

    except (OSError, AttributeError):
        print("warning: could not kill hmmpgmd worker")
        pass
    
    try:

        parent = psutil.Process(MASTER.pid)
        for child in parent.children(recursive=True):  # or parent.children() for recursive=False
            child.kill()
        parent.kill()
                
    except Exception as e:
        print(e)
    except (OSError, AttributeError):
        print("warning: could not kill hmmpgmd master")
        pass
    
    return

def shutdown_server():
    global MASTER, WORKERS
    shutdown_server_by_pid(MASTER, WORKERS)
    return
 
    
    
def alive(p):
    """ Check For the existence of a unix pid. """
    return p.is_alive()
    
## END
