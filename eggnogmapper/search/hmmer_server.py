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
from .hmmer_search import get_hits, DB_TYPE_SEQ, DB_TYPE_HMM

CHILD_PROC = None
MASTER = None
WORKER = None

def server_up(host, port):
    import socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex((host, port))
    sock.close()
    if result == 0:
        return True
    else:
        return False

def server_functional(host, port, dbtype = DB_TYPE_HMM):
    if server_up(host, port):
        try:
            get_hits("test", "TESTSEQ", host, port, dbtype)
        except Exception as e:
            traceback.print_exc()
            print(e)
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
    global CHILD_PID, WORKER
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

    WORKER = Process(target=start_worker)
    WORKER.start()
    
    return WORKER

def load_server(dbpath, client_port, worker_port, cpu, output=None, dbtype=DB_TYPE_HMM, is_worker = True):
    global CHILD_PID, MASTER, WORKER
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
        cmd = HMMPGMD + f' --worker localhost --wport {worker_port} --cpu {cpu}'
        print(colorify(f"Loading worker: {cmd}", 'orange'))
        CHILD_PROC = subprocess.Popen(cmd.split(), shell=False, stderr=OUT, stdout=OUT)
        while 1:
            time.sleep(60)
    
    MASTER = Process(target=start_master)
    MASTER.start()

    if is_worker == True:
        WORKER = Process(target=start_worker)
        WORKER.start()
    
    return dbpath, MASTER, WORKER

def shutdown_server():
    global MASTER, WORKER
    try:
        # This is killing THIS python script also, and is UNIX dependent
        # os.killpg(os.getpgid(WORKER.pid), signal.SIGTERM)
        
        import psutil
        parent = psutil.Process(WORKER.pid)
        for child in parent.children(recursive=True):  # or parent.children() for recursive=False
            child.kill()
            parent.kill()
            
        # WORKER.terminate()
        # if WORKER.is_alive():
        #     WORKER.kill()
        #     if WORKER.is_alive:
        #         os.kill(WORKER.pid, signal.SIGKILL)

    except (OSError, AttributeError):
        print("warning: could not kill hmmpgmd worker")
        pass
    
    try:
        # This is killing THIS python script also, and is UNIX dependent
        # os.killpg(os.getpgid(MASTER.pid), signal.SIGTERM)

        import psutil
        parent = psutil.Process(MASTER.pid)
        for child in parent.children(recursive=True):  # or parent.children() for recursive=False
            child.kill()
            parent.kill()
            
        # MASTER.terminate()
        # if MASTER.is_alive():
        #     MASTER.kill()
        #     if MASTER.is_alive():
        #         os.kill(MASTER.pid, signal.SIGKILL)
                
    except Exception as e:
        print(e)
    except (OSError, AttributeError):
        print("warning: could not kill hmmpgmd master")
        pass
    
    return
 
    
    
def alive(p):
    """ Check For the existence of a unix pid. """
    return p.is_alive()
    
