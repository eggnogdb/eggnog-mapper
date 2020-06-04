##
## CPCantalapiedra 2020

import os, sys
import time

from os.path import exists as pexists
from os.path import join as pjoin

from ..common import EGGNOG_DATABASES, get_db_info, get_data_path
from ..utils import colorify

# from .hmmer_server import server_functional
from .hmmer_search import SCANTYPE_MEM, DB_TYPE_SEQ, DB_TYPE_HMM, QUERY_TYPE_SEQ
from .hmmer_idmap import generate_idmap

SETUP_TYPE_REMOTE = "remote"
SETUP_TYPE_EGGNOG = "eggnog"
SETUP_TYPE_CUSTOM = "custom"

##
def setup_hmm_search(db, scantype, dbtype, qtype = QUERY_TYPE_SEQ, servers_list = None):

    setup_type = None
    
    if ":" in db or servers_list is not None:
        
        dbname, dbpath, host, port, idmap_file = setup_remote_db(db, dbtype, qtype)
        end_port = port
        setup_type = SETUP_TYPE_REMOTE

    else: # setup_local_db --> dbpath, host, port, idmap_file
        if db in EGGNOG_DATABASES:
            dbname, dbpath, host, port, end_port, idmap_file = setup_eggnog_db(db, scantype)
            setup_type = SETUP_TYPE_EGGNOG
        else:
            dbname, dbpath, host, port, end_port, idmap_file = setup_custom_db(db, scantype, dbtype)
            setup_type = SETUP_TYPE_CUSTOM

    return dbname, dbpath, host, port, end_port, idmap_file, setup_type


##
def setup_eggnog_db(db, scantype):

    dbpath, port = get_db_info(db)

    db_present = [pexists(dbpath + "." + ext)
                  for ext in 'h3f h3i h3m h3p idmap'.split()]

    if False in db_present:
        print(db_present)
        print(colorify('Database %s not present. Use download_eggnog_database.py to fetch it' % (db), 'red'))
        raise ValueError('Database not found')

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
def setup_custom_db(db, scantype, dbtype = DB_TYPE_HMM):
    
    if dbtype == DB_TYPE_HMM:
        db, dbpath, host, port, end_port, idmap_file = setup_custom_hmmdb(db, scantype)
    elif dbtype == DB_TYPE_SEQ:
        db, dbpath, host, port, end_port, idmap_file = setup_custom_seqdb(db, scantype)
    else:
        raise Exception("Unrecognized dbtype {dbtype}.")

    return db, dbpath, host, port, end_port, idmap_file
    
        
##
def setup_custom_hmmdb(db, scantype):
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
                # idmap_file = db + ".idmap"
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
def setup_custom_seqdb(db, scantype):

    # if not os.path.isfile(db + '.seqdb'):
    #     print(colorify('esl-reformat database (with name %s.seqdb) not present. Use esl-reformat to create the .seqdb file' % (db), 'red'))
    #     raise ValueError('esl-reformat database not found')

    print(colorify(f"Preparing to query custom database {db}", 'green'))

    if scantype == SCANTYPE_MEM:
        if not os.path.isfile(db + '.seqdb'):
            print(colorify('esl-reformat database (with name %s.seqdb) not present. Use esl-reformat to create the .seqdb file' % (db), 'red'))
            raise ValueError('esl-reformat database not found')
        dbpath = db + ".seqdb"
        host = 'localhost'
        port = 53000
        end_port = 53200

        idmap_file = db + ".map"
        if not pexists(idmap_file):
            print(colorify('ID map file (with name %s.map) not present. Use esl-reformat to create the .map file' % (db), 'red'))
            raise ValueError('esl-reformat .map file not found')
    else:
        dbpath = db
        host = None
        port = None
        end_port = None
        idmap_file = None
    
    return db, dbpath, host, port, end_port, idmap_file


##
def setup_remote_db(db, dbtype, qtype):
    
    if dbtype == DB_TYPE_HMM:
        dbname, dbpath, host, port, idmap_file = setup_remote_hmmdb(db, dbtype, qtype)
    elif dbtype == DB_TYPE_SEQ:
        dbname, dbpath, host, port, idmap_file = setup_remote_seqdb(db, dbtype, qtype)
    else:
        raise Exception("Unrecognized dbtype {dbtype}.")

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

    if dbname in EGGNOG_DATABASES:
        dbfile, port = get_db_info(dbname)
        db = dbname
    else:
        dbfile = dbname

    idmap_file = dbfile + '.idmap'
    if not pexists(idmap_file):
        raise ValueError("idmap file not found: %s" % idmap_file)

    # if not server_functional(host, port, dbtype, qtype):
    #     print(colorify("eggnog-mapper server not found at %s:%s" % (host, port), 'red'))
    #     exit(1)

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
        raise ValueError("%s file not found" % dbfile)
    
    idmap_file = dbname + '.map'
    if not pexists(idmap_file):
        raise ValueError("%s file not found" % idmap_file)

    # if not server_functional(host, port, dbtype, qtype):
    #     print(colorify("eggnog-mapper server not found at %s:%s" % (host, port), 'red'))
    #     exit(1)

    return dbname, dbpath, host, port, idmap_file

## END
