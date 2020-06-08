##
## JHCepas

import sys
import os
import time
import subprocess
from multiprocessing import Process
import signal
import traceback

from ..common import HMMPGMD, TIMEOUT_LOAD_SERVER
from ..utils import colorify
from .hmmer_search import get_hits, DB_TYPE_SEQ, DB_TYPE_HMM, QUERY_TYPE_SEQ, QUERY_TYPE_HMM

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
def create_servers(dbtype, dbpath, host, port, end_port, num_servers, num_workers, cpus_per_worker):
    servers = []
    sdbpath = dbpath
    shost = host
    sport = port
    for num_server, server in enumerate(range(num_servers)):
        sdbpath, shost, sport, master_pid, workers_pids = start_server(dbpath, host, sport, end_port, cpus_per_worker, num_workers, dbtype)
        servers.append((sdbpath, sport, master_pid, workers_pids))
        sport = sport + 2
    dbpath = sdbpath
    host = shost
    port = sport

    return dbpath, host, port, servers


##
def start_server(dbpath, host, port, end_port, cpus_per_worker, num_workers, dbtype, qtype = QUERY_TYPE_SEQ):
    master_db = worker_db = workers = None
    for try_port in range(port, end_port, 2):
        print(colorify("Loading server at localhost, port %s-%s" %
                       (try_port, try_port + 1), 'lblue'))

        dbpath, master_db, workers = load_server(dbpath, try_port, try_port + 1,
                                                 cpus_per_worker, num_workers = num_workers, dbtype = dbtype)
        port = try_port
        ready = False
        for _ in range(TIMEOUT_LOAD_SERVER):
            print(f"Waiting for server to become ready at {host}:{port} ...")
            time.sleep(1)
            if not master_db.is_alive() or not any([worker_db.is_alive() for worker_db in workers]):
                master_db.terminate()
                master_db.join()
                for worker_db in workers:
                    worker_db.terminate()
                    worker_db.join()                        
                break
            elif server_functional(host, port, dbtype, qtype):
                print(f"Server ready at {host}:{port}")
                ready = True
                break
            else:
                print(f"Waiting for server to become ready at {host}:{port} ...")

        if ready:
            dbpath = host
            break

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
    else:
        print(colorify("Server is still down", 'red'))
    return False

def safe_exit(a, b):
    if CHILD_PROC:
        CHILD_PROC.kill()
    sys.exit(0)

def load_worker(master_host, worker_port, cpu, output=None):
    global CHILD_PID, WORKERS
    if not output:
        OUT = open(os.devnull, 'w')
    else:
        OUT = output
        
    signal.signal(signal.SIGINT, safe_exit)
    signal.signal(signal.SIGTERM, safe_exit)
              
    def start_worker():
        cmd = HMMPGMD + f' --worker {master_host} --wport {worker_port} --cpu {cpu}'
        print(colorify(f"Loading worker: {cmd}", 'orange'))
        CHILD_PROC = subprocess.Popen(cmd.split(), shell=False, stderr=OUT, stdout=OUT)
        while 1:
            time.sleep(60)
            
    worker = Process(target=start_worker)
    worker.start()
    WORKERS = [worker]
    
    return worker

def load_server(dbpath, client_port, worker_port, cpus_per_worker, num_workers=1, output=None, dbtype=DB_TYPE_HMM, is_worker = True):
    global CHILD_PID, MASTER, WORKERS
    if not output:
        OUT = open(os.devnull, 'w')
    else:
        OUT = output
        
    signal.signal(signal.SIGINT, safe_exit)
    signal.signal(signal.SIGTERM, safe_exit)

    def start_master():
        cmd = HMMPGMD + f' --master --cport {client_port} --wport {worker_port} --{dbtype} {dbpath}'
        print(colorify(f"Loading master: {cmd}", 'orange'))
        CHILD_PROC = subprocess.Popen(cmd.split(), shell=False, stderr=OUT, stdout=OUT)
        while 1:
            time.sleep(60)

    def start_worker():
        cmd = HMMPGMD + f' --worker localhost --wport {worker_port} --cpu {cpus_per_worker}'
        print(colorify(f"Loading worker: {cmd}", 'orange'))
        CHILD_PROC = subprocess.Popen(cmd.split(), shell=False, stderr=OUT, stdout=OUT)
        while 1:
            time.sleep(60)

    print(f"Creating hmmpgmd server at {host}:{port} ...")
    MASTER = Process(target=start_master)
    MASTER.start()

    if is_worker == True and num_workers > 0:
        print(f"Creating hmmpgmd workers ({num_workers}) at {host}:{port} ...")
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
