##
## CPCantalapiedra 2020

import os, sys
import time

from os.path import exists as pexists
from os.path import join as pjoin

from ..common import EGGNOG_DATABASES, get_db_info, TIMEOUT_LOAD_SERVER, get_data_path
from ..utils import colorify

from .hmmer_server import server_functional, load_server
from .hmmer_search import SCANTYPE_MEM
from .hmmer_idmap import generate_idmap

SETUP_TYPE_REMOTE = "remote"
SETUP_TYPE_EGGNOG = "eggnog"
SETUP_TYPE_CUSTOM = "custom"

##
def setup_hmm_search(db, scantype, dbtype, no_refine, cpu, servermode):

    setup_type = None
    
    if ":" in db:
        dbname, dbpath, host, port, idmap_file = setup_remote_db(db, dbtype)
        setup_type = SETUP_TYPE_REMOTE

    else: # setup_local_db --> dbpath, host, port, idmap_file
        if db in EGGNOG_DATABASES:
            dbname, dbpath, host, port, end_port, idmap_file = setup_eggnog_db(db, no_refine, scantype)
            setup_type = SETUP_TYPE_EGGNOG
        else:
            dbname, dbpath, host, port, end_port, idmap_file = setup_custom_db(db, scantype)
            setup_type = SETUP_TYPE_CUSTOM

        if scantype == SCANTYPE_MEM:
            dbpath, host, port = start_server(dbpath, host, port, end_port, cpu, dbtype)

    # If server mode, just listen for connections and exit when interrupted
    if servermode:
        while True:
            print(colorify("Server ready listening at %s:%s and using %d CPU cores" % (host, port, cpu), 'green'))
            print(colorify("Use `emapper.py -d %s:%s:%s (...)` to search against this server" % (db, host, port), 'lblue'))
            time.sleep(10)
        raise emapperException("Server {db}:{host}:{port} stopped.")
    else:
        return dbname, dbpath, host, port, idmap_file, setup_type
    
##
def setup_eggnog_db(db, no_refine, scantype):

    dbpath, port = get_db_info(db)

    db_present = [pexists(dbpath + "." + ext)
                  for ext in 'h3f h3i h3m h3p idmap'.split()]

    if False in db_present:
        print(db_present)
        print(colorify('Database %s not present. Use download_eggnog_database.py to fetch it' % (db), 'red'))
        raise ValueError('Database not found')

    if not no_refine:
        if not pexists(pjoin(get_data_path(), 'OG_fasta')):
            print(colorify('Database data/OG_fasta/ not present. Use download_eggnog_database.py to fetch it', 'red'))
            raise ValueError('Database fasta sequences not found')

    if scantype == SCANTYPE_MEM:
        host = 'localhost'
        end_port = 53200
        idmap_file = dbpath + '.idmap'
    else:
        host = None
        port = None
        end_port = None
        idmap_file = None

    return db, dbpath, host, port, end_port, idmap_file

##
def setup_custom_db(db, scantype):

    if not os.path.isfile(db + '.h3f'):
        print(colorify('Database %s not present. Use hmmpress to create the h3* files' % (db), 'red'))
        raise ValueError('Database not found')

    print(colorify(f"Preparing to query custom database {db}", 'green'))
    dbpath = db

    if scantype == SCANTYPE_MEM:
        host = 'localhost'
        port = 53000
        end_port = 53200

        idmap_file = dbpath + ".idmap"
        if not pexists(idmap_file):
            if generate_idmap(db):
                idmap_file = db + ".idmap"
                print("idmap succesfully created!", file=sys.stderr)
            else:
                raise ValueError("idmap could not be created!")
    else:
        host = None
        port = None
        end_port = None
        idmap_file = None

    return db, dbpath, host, port, end_port, idmap_file

##
def setup_remote_db(db, dbtype):
    dbname, host, port = map(str.strip, db.split(":"))
    port = int(port)
    if dbname in EGGNOG_DATABASES:
        dbfile, port = get_db_info(dbname)
        db = dbname
    else:
        dbfile = dbname

    idmap_file = dbfile + '.idmap'
    if not pexists(idmap_file):
        raise ValueError("idmap file not found: %s" % idmap_file)

    dbpath = host
    if not server_functional(host, port, dbtype):
        print(colorify("eggnog-mapper server not found at %s:%s" % (host, port), 'red'))
        exit(1)

    return dbname, dbpath, host, port, idmap_file

##
def start_server(dbpath, host, port, end_port, cpu, dbtype):
    master_db, worker_db = None, None
    for try_port in range(port, end_port, 2):
        print(colorify("Loading server at localhost, port %s-%s" %
                       (try_port, try_port + 1), 'lblue'))
        dbpath, master_db, worker_db = load_server(
            dbpath, try_port, try_port + 1, cpu, dbtype = dbtype)
        port = try_port
        ready = False
        for _ in range(TIMEOUT_LOAD_SERVER):
            print(f"Waiting for server to become ready at {host}:{port} ...")
            time.sleep(1)
            if not master_db.is_alive() or not worker_db.is_alive():
                master_db.terminate()
                master_db.join()
                worker_db.terminate()
                worker_db.join()
                break
            elif server_functional(host, port, dbtype):
                ready = True
                break
            
        if ready:
            dbpath = host
            break

    return dbpath, host, port

## END
