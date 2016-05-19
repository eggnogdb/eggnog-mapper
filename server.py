#!/usr/bin/env python

import sys
import argparse
import os
import time
import subprocess
from multiprocessing import Process

HMMPGMD_BIN = 'hmmpgmd'
BASEPATH = os.path.split(os.path.abspath(__file__))[0]

# Predefined eggnog databases
DBDATA = {
    'euk': { 'name': 'euk_500', 'db_path':os.path.join(BASEPATH, 'hmmdb/euk_500/euk_500.hmm'), 'client_port':51400, 'worker_port':51401, 'idmap':os.path.join(BASEPATH, 'hmmdb/euk_500/euk_500.pkl'), 'cpu':20},
    'bact':{ 'name': 'bact_50', 'db_path':os.path.join(BASEPATH, 'hmmdb/bact_50/bact_50.hmm'), 'client_port':51500, 'worker_port':51501, 'idmap':os.path.join(BASEPATH, 'hmmdb/bact_50/bact_50.pkl'), 'cpu':20},
    'arch':{ 'name': 'arch_1', 'db_path':os.path.join(BASEPATH, 'hmmdb/arch_1/arch_1.hmm'), 'client_port':51600, 'worker_port':51601, 'idmap':os.path.join(BASEPATH, 'hmmdb/arch_1/arch_1.pkl'), 'cpu':20},
    }

def load_server(dbpath, client_port, worker_port, cpu):        
    FNULL = open(os.devnull, 'w')
    def start_master():
        cmd = HMMPGMD_BIN +' --master --cport %d --wport %s --hmmdb %s' %(client_port, worker_port, dbpath)
        subprocess.call(cmd, shell=True, stderr=FNULL, stdout=FNULL)
        
    def start_worker():
        cmd = HMMPGMD_BIN +' --worker localhost --wport %s --cpu %d' %(worker_port, cpu)
        subprocess.call(cmd, shell=True, stderr=FNULL, stdout=FNULL)
        
    
    master = Process(target=start_master)
    master.start()

    worker = Process(target=start_worker)
    worker.start()
    return dbpath, master, worker

    
def check_pid(pid):
    """ Check For the existence of a unix pid. """
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    else:
        return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', dest='db', nargs='+', choices=['euk', 'bact', 'arch'])
    parser.add_argument('--cpu', dest='cpu', type=int, default=20)
        
    args = parser.parse_args()

    checker = []
    
    for dbname in args.db:
        db = DBDATA[dbname]
        checker.append(load_server(db['db_path'], db['client_port'], db['worker_port'], db['cpu']))

    # Keep the server runninf
    while 1:
        print 'Eggnog-mapper serving databases...', time.ctime()
        for dbpath, master, worker in checker:
            print "  % 30s - Master (% 5d):%s - Worker (% 5d):%s" (dbpath, master.pid, check_pid(master.pid), worker.pid, check_pid(worker.pid))        
        time.sleep(60)
