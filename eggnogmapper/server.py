from __future__ import absolute_import


import sys
import os
import time
import subprocess
from multiprocessing import Process
import signal

from common import *
import search

CHILD_PROC = None

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
    try:
        search.get_hits("test", "TESTSEQ", host, port, dbtype)
    except Exception, e:
        print 'Server not ready', e
        return False
    return True

def safe_exit(a, b):
    if CHILD_PROC:
        CHILD_PROC.kill()
    sys.exit(0)
        
def load_server(dbpath, client_port, worker_port, cpu, output=None):
    global CHILD_PID
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
        
    master = Process(target=start_master)
    master.start()

    worker = Process(target=start_worker)
    worker.start()
    
    return dbpath, master, worker
    
def alive(p):
    """ Check For the existence of a unix pid. """
    return p.is_alive()


def generate_idmap(dbpath):
    cmd = """%s %s |grep -v '#'|awk '{print $1" "$2}' > %s""" %(HMMSTAT, dbpath, dbpath+'.idmap')
    print('Generating idmap in '+dbpath+'.idmap')
    return os.system(cmd) == 0



    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', dest='db', nargs='+', choices=['euk', 'bact', 'arch'])
    parser.add_argument('--cpu', dest='cpu', type=int, default=4)
    
    args = parser.parse_args()

    checker = []
    
    for dbname in args.db:
        db = DBDATA[dbname]
        checker.append(load_server(db['db_path'], db['client_port'], db['worker_port'], args.cpu, sys.stdout))

    # Keep the server runninf
    while 1:
        print 'Eggnog-mapper serving databases...', time.ctime()
        for dbpath, master, worker in checker:
            print "  % 30s - Master (% 5d):%s - Worker (% 5d):%s" (dbpath, master.pid, alive(master), worker.pid, alive(worker))        
        time.sleep(60)
