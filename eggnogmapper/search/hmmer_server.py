##
## JHCepas

import sys
import os
import time
import subprocess
from multiprocessing import Process
import signal
import traceback

from ..common import *
from ..utils import colorify
from .hmmer_search import get_hits, DB_TYPE_SEQ, DB_TYPE_HMM, QUERY_TYPE_SEQ, QUERY_TYPE_HMM

CHILD_PROC = None
MASTER = None
WORKERS = None

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

def load_server(dbpath, client_port, worker_port, cpus_per_worker, workers=1, output=None, dbtype=DB_TYPE_HMM, is_worker = True):
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
    
    MASTER = Process(target=start_master)
    MASTER.start()
        
    if is_worker == True and workers > 0:
        WORKERS = []
        for i in range(workers):
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
    
