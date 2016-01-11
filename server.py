#!/usr/bin/env python

import sys
import argparse
import os
import time

HMMPGMD = 'hmmpgmd'
BASEPATH = os.path.split(os.path.abspath(__file__))[0]

DBDATA = {
    'euk': { 'name': 'euk_500', 'db_path':os.path.join(BASEPATH, 'hmmdb/euk_500/euk_500.hmm'), 'client_port':51400, 'worker_port':51401, 'idmap':os.path.join(BASEPATH, 'hmmdb/euk_500/euk_500.pkl'), 'cpu':20},
    'bact':{ 'name': 'bact_50', 'db_path':os.path.join(BASEPATH, 'hmmdb/bact_50/bact_50.hmm'), 'client_port':51500, 'worker_port':51501, 'idmap':os.path.join(BASEPATH, 'hmmdb/bact_50/bact_50.pkl'), 'cpu':20},
    'arch':{ 'name': 'arch_1', 'db_path':os.path.join(BASEPATH, 'hmmdb/arch_1/arch_1.hmm'), 'client_port':51600, 'worker_port':51601, 'idmap':os.path.join(BASEPATH, 'hmmdb/arch_1/arch_1.pkl'), 'cpu':20},
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', dest='db', nargs='+', choices=['euk', 'bact', 'arch'])
    #parser.add_argument('--master', dest='master', nargs='+', choices=['euk', 'bact', 'arch'])
    #parser.add_argument('--worker', dest='worker', nargs='+', choices=['euk', 'bact', 'arch'])
    parser.add_argument('--cpu', dest='cpu', type=int, default=20)
    args = parser.parse_args()
    
    for dbname in args.db:
        pid = os.fork()
        if pid == 0:
            cmd = HMMPGMD + ' --master --cport %(client_port)d --wport %(worker_port)s --hmmdb %(db_path)s --pid %(name)s.m.pid' %(DBDATA[dbname])
            print cmd
            os.system(cmd)
            sys.exit(0)
            

    for dbname in args.db:
        DBDATA[dbname]["cpu"] = args.cpu
        pid = os.fork()
        if pid == 0:
            cmd = HMMPGMD + ' --worker localhost --wport %(worker_port)s --pid %(name)s.w.pid --cpu %(cpu)s' %(DBDATA[dbname])
            print cmd
            os.system(cmd)
            sys.exit(0)
            
    while 1:
        print 'Serving ', time.ctime()
        time.sleep(60)

        
        
