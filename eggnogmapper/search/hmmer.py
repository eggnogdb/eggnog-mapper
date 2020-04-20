##
## CPCantalapiedra 2020


class HmmerSearcher:

    ##
    def __init__(self, args):
        return

    ##
    def search(self, in_file, seed_orthologs_file, hmm_hits_file):
        host, port, dbpath, scantype, idmap = setup_hmm_search(args)

        # Start HMM SCANNING sequences (if requested)                                                                                                                              
        if not pexists(hmm_hits_file) or args.override:
            dump_hmm_matches(in_file, hmm_hits_file, dbpath, port, scantype, idmap, args)

        if not args.no_refine and (not pexists(seed_orthologs_file) or args.override):
            if args.db == 'viruses':
                print 'Skipping seed ortholog detection in "viruses" database'
            elif args.db in EGGNOG_DATABASES:
                refine_matches(in_file, seed_orthologs_file, hmm_hits_file, args)
            else:
                print 'refined hits not available for custom hmm databases.'

        return

    ##
    def setup_hmm_search(args):
        host = 'localhost'
        idmap = None
        if args.usemem:
            scantype = 'mem'
        else:
           scantype = 'disk'

        connecting_to_server = False
        # If searching against a predefined database name
        if args.db in EGGNOG_DATABASES:
            dbpath, port = get_db_info(args.db)
            print dbpath
            db_present = [pexists(dbpath + "." + ext)
                          for ext in 'h3f h3i h3m h3p idmap'.split()]

            if False in db_present:
                print db_present
                print colorify('Database %s not present. Use download_eggnog_database.py to fetch it' % (args.db), 'red')
                raise ValueError('Database not found')

            if not args.no_refine:
                if not pexists(pjoin(get_data_path(), 'OG_fasta')):
                    print colorify('Database data/OG_fasta/ not present. Use download_eggnog_database.py to fetch it', 'red')
                    raise ValueError('Database not found')

            if scantype == 'mem':
                idmap_file = dbpath + '.idmap'
                end_port = 53200

        # If searching against a custom hmm database
        elif os.path.isfile(args.db + '.h3f'):
            dbpath = args.db
            if scantype == 'mem':
                idmap_file = args.db + ".idmap"
                if not pexists(idmap_file):
                    if generate_idmap(args.db):
                        idmap_file = args.db + ".idmap"
                        print >>sys.stderr, "idmap succesfully created!"
                    else:
                        raise ValueError("idmap could not be created!")
                port = 53000
                end_port = 53200
            else:
                idmap_file = None
                port = None

        # If searching against a emapper hmm server
        elif ":" in args.db:
            dbname, host, port = map(str.strip, args.db.split(":"))
            scantype = 'mem'
            port = int(port)
            if dbname in EGGNOG_DATABASES:
                dbfile, port = get_db_info(dbname)
                args.db = dbname
            else:
                dbfile = dbname

            idmap_file = dbfile + '.idmap'
            if not pexists(idmap_file):
                raise ValueError("idmap file not found: %s" % idmap_file)

            dbpath = host
            if not server_functional(host, port, args.dbtype):
                print colorify("eggnog-mapper server not found at %s:%s" % (host, port), 'red')
                exit(1)
            connecting_to_server = True
        else:
            raise ValueError('Invalid database name/server')


        # If memory based searches requested, start server
        if scantype == "mem" and not connecting_to_server:
            master_db, worker_db = None, None
            for try_port in range(port, end_port, 2):
                print colorify("Loading server at localhost, port %s-%s" %
                               (try_port, try_port + 1), 'lblue')
                dbpath, master_db, worker_db = load_server(
                    dbpath, try_port, try_port + 1, args.cpu)
                port = try_port
                ready = False
                for _ in xrange(TIMEOUT_LOAD_SERVER):
                    print "Waiting for server to become ready...", host, try_port
                    time.sleep(1)
                    if not master_db.is_alive() or not worker_db.is_alive():
                        master_db.terminate()
                        master_db.join()
                        worker_db.terminate()
                        worker_db.join()
                        break
                    elif server_functional(host, port, args.dbtype):
                        ready = True
                        break
                if ready:
                    dbpath = host
                    break
        elif scantype == "mem":
            print colorify("DB Server already running or not needed!", 'yellow')
            dbpath = host

        # Preload seqid map to translate hits from hmmpgmd
        if scantype == "mem":
            print colorify("Reading idmap %s" % idmap_file, color='lblue')
            idmap = {}
            for _lnum, _line in enumerate(open(idmap_file)):
                if not _line.strip():
                    continue
                try:
                    _seqid, _seqname = map(str, _line.strip().split(' '))
                except ValueError:
                    if _lnum == 0:
                        # idmap generated by esl_reformat has info line at beginning
                        continue  
                    else:
                        raise
                _seqid = int(_seqid)
                idmap[_seqid] = [_seqname]
            print len(idmap), "names loaded"

        # If server mode, just listen for connections and exit when interrupted
        if args.servermode:
            while True:
                print colorify("Server ready listening at %s:%s and using %d CPU cores" % (host, port, args.cpu), 'green')
                print colorify("Use `emapper.py -d %s:%s:%s (...)` to search against this server" % (args.db, host, port), 'lblue')
                time.sleep(10)
            raise emapperException()
        else:
            return host, port, dbpath, scantype, idmap
    

## END
