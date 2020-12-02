##
## CPCantalapiedra 2020

import os, sys
import time

from os.path import exists as pexists
from os.path import join as pjoin

from ...common import get_db_info, get_data_path, get_hmmer_databases, get_hmmer_dbpath
from ...utils import colorify
from ...emapperException import EmapperException

# from .hmmer_server import server_functional
from .hmmer_search import SCANTYPE_MEM, DB_TYPE_SEQ, DB_TYPE_HMM, QUERY_TYPE_SEQ
from .hmmer_idmap import generate_idmap

SETUP_TYPE_REMOTE = "remote"
SETUP_TYPE_EGGNOG = "eggnog"
SETUP_TYPE_CUSTOM = "custom"

DEFAULT_PORT = 51700
DEFAULT_END_PORT = 53200

##
def setup_hmm_search(db, scantype, dbtype, qtype = QUERY_TYPE_SEQ, port = DEFAULT_PORT, end_port = DEFAULT_END_PORT, servers_list = None, silent = False):

    setup_type = None
    
    if ":" in db or servers_list is not None:
        
        dbname, dbpath, host, port, idmap_file = setup_remote_db(db, dbtype, qtype)
        end_port = port
        setup_type = SETUP_TYPE_REMOTE

    else: # setup_local_db --> dbpath, host, port, idmap_file
        
        if db in get_hmmer_databases():
            dbpath, host, idmap_file = setup_eggnog_db(db, scantype)
            dbname = db
            setup_type = SETUP_TYPE_EGGNOG
        else:
            dbpath, host, idmap_file = setup_custom_db(db, scantype, dbtype, silent)
            dbname = db
            setup_type = SETUP_TYPE_CUSTOM

    return dbname, dbpath, host, port, end_port, idmap_file, setup_type


##
def setup_eggnog_db(db, scantype):
    dbpath = get_hmmer_dbpath(db)
    print(dbpath)
    db_present = [pexists(dbpath + "." + ext) for ext in 'h3f h3i h3m h3p idmap'.split()]

    if False in db_present:
        raise EmapperException(f'Error: files (.h3* or .idmap) for database {db} not present. Use "download_eggnog_database.py" to fetch it')

    if scantype == SCANTYPE_MEM:
        host = 'localhost'
        idmap_file = dbpath + '.idmap'
    else:
        host = None
        idmap_file = None

    return dbpath, host, idmap_file


##
def setup_custom_db(db, scantype, dbtype = DB_TYPE_HMM, silent = False):
    
    if dbtype == DB_TYPE_HMM:
        dbpath, host, idmap_file = setup_custom_hmmdb(db, scantype, silent)
    elif dbtype == DB_TYPE_SEQ:
        dbpath, host, idmap_file = setup_custom_seqdb(db, scantype, silent)
    else:
        raise EmapperException(f"Unrecognized dbtype {dbtype}.")

    return dbpath, host, idmap_file
    
        
##
def setup_custom_hmmdb(db, scantype, silent = False):
    if not os.path.isfile(db + '.h3f'):
        print(colorify(f'Database {db} not present. Use "hmmpress" to create the .h3* files', 'red'))
        raise EmapperException(f'Database {db} not found')

    if silent == False:
        print(colorify(f"Preparing to query custom database {db}", 'green'))

    if scantype == SCANTYPE_MEM:
        host = 'localhost'
        idmap_file = db + ".idmap"
        if not pexists(idmap_file):
            if generate_idmap(db):
                if silent == False:
                    print("idmap succesfully created!", file=sys.stderr)
            else:
                raise EmapperException(f"idmap file {idmap_file} could not be created!")
    else:
        host = None
        idmap_file = None

    dbpath = db
    
    return dbpath, host, idmap_file


##
def setup_custom_seqdb(db, scantype, silent = False):

    if silent == False:
        print(colorify(f"Preparing to query custom database {db}", 'green'))

    if scantype == SCANTYPE_MEM:
        if not os.path.isfile(db + '.seqdb'):
            print(colorify(f'esl-reformat database (with name {db}.seqdb) not present. Use esl-reformat to create the .seqdb file', 'red'))
            raise EmapperException(f'esl-reformat database {db}.seqdb not found')
        dbpath = db + ".seqdb"
        host = 'localhost'

        idmap_file = db + ".map"
        if not pexists(idmap_file):
            print(colorify(f'ID map file (with name {db}.map) not present. Use esl-reformat to create the .map file', 'red'))
            raise ValueError(f'esl-reformat {db}.map file not found')
    else:
        dbpath = db
        host = None
        idmap_file = None
    
    return dbpath, host, idmap_file


##
def setup_remote_db(db, dbtype, qtype):
    
    if dbtype == DB_TYPE_HMM:
        dbname, dbpath, host, port, idmap_file = setup_remote_hmmdb(db, dbtype, qtype)
    elif dbtype == DB_TYPE_SEQ:
        dbname, dbpath, host, port, idmap_file = setup_remote_seqdb(db, dbtype, qtype)
    else:
        raise EmapperException(f"Unrecognized dbtype {dbtype}.")

    return dbname, dbpath, host, port, idmap_file

##
def setup_remote_hmmdb(db, dbtype, qtype):
    if ":" in db:
        dbname, host, port = map(str.strip, db.split(":"))
        dbpath = host
        port = int(port)
        
    else:
        dbname = db
        dbpath = None
        host = None
        port = None
        
    if dbname in get_hmmer_databases():
        dbfile, port = get_db_info(dbname)
        db = dbname
    else:
        dbfile = dbname

    idmap_file = dbfile + '.idmap'
    if not pexists(idmap_file):
        raise EmapperException(f"idmap file {idmap_file} not found")

    return dbname, dbpath, host, port, idmap_file


##
def setup_remote_seqdb(db, dbtype, qtype):
    if ":" in db:
        dbname, host, port = map(str.strip, db.split(":"))
        dbpath = host
        port = int(port)
        
    else:
        dbname = db
        dbpath = None
        host = None
        port = None

    dbfile = dbname + ".seqdb"
    if not pexists(dbfile):
        raise EmapperException(f"db file {dbfile} not found")
    
    idmap_file = dbname + '.map'
    if not pexists(idmap_file):
        raise EmapperException(f"db idmap file {idmap_file} file not found")

    return dbname, dbpath, host, port, idmap_file

## END
