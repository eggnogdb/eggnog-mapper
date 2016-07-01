from __future__ import absolute_import


import sys
import os
import time
import subprocess
from multiprocessing import Process
import signal

from .common import *
from .utils import colorify
from . import search

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

def server_functional(host, port, dbtype):
    if server_up(host, port):
        try:
            search.get_hits("test", "TESTSEQ", host, port, dbtype)
        except Exception, e:
            #print 'Server not ready', e
            return False
        else:
            return True
    return False

def safe_exit(a, b):
    if CHILD_PROC:
        CHILD_PROC.kill()
    sys.exit(0)
        
def load_server(dbpath, client_port, worker_port, cpu, output=None):
    global CHILD_PID, MASTER, WORKER
    if not output:
        OUT = open(os.devnull, 'w')
    else:
        OUT = output
        
    signal.signal(signal.SIGINT, safe_exit)
    signal.signal(signal.SIGTERM, safe_exit)
    
    def start_master():
        cmd = HMMPGMD +' --master --cport %d --wport %s --hmmdb %s' %(client_port, worker_port, dbpath)
        CHILD_PROC = subprocess.Popen(cmd.split(), shell=False, stderr=OUT, stdout=OUT)
        while 1:
            time.sleep(60)
              
    def start_worker():
        cmd = HMMPGMD +' --worker localhost --wport %s --cpu %d' %(worker_port, cpu)
        CHILD_PROC = subprocess.Popen(cmd.split(), shell=False, stderr=OUT, stdout=OUT)
        while 1:
            time.sleep(60)
        
    MASTER = Process(target=start_master)
    MASTER.start()

    WORKER = Process(target=start_worker)
    WORKER.start()
    
    return dbpath, MASTER, WORKER

def shutdown_server():
    global MASTER, WORKER
    try:
        os.killpg(os.getpgid(MASTER.pid), signal.SIGTERM)
    except (OSError, AttributeError):
        pass
    try:
        os.killpg(os.getpgid(WORKER.pid), signal.SIGTERM)
    except (OSError, AttributeError):
        pass
 
    
    
def alive(p):
    """ Check For the existence of a unix pid. """
    return p.is_alive()


def generate_idmap(dbpath):
    if dbpath.endswith(".h3f"):
        dbpath = dbpath.replace(".h3f", "")
    cmd = """%s %s |grep -v '#'|awk '{print $1" "$2}' > %s""" %(HMMSTAT, dbpath, dbpath+'.idmap')
    print colorify(cmd, "cyan")
    print('Generating idmap in '+dbpath+'.idmap')
    return os.system(cmd) == 0
    
