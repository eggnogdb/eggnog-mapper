#!/usr/bin/env python
import os
from argparse import ArgumentParser
from eggnogmapper.common import EGGNOG_DATABASES, DATA_PATH, HMMDB_PATH, pexists, pjoin, get_level_base_path
from eggnogmapper.utils import ask, colorify

def run(cmd):
    print colorify(cmd, 'cyan')
    if not args.simulate:
        os.system(cmd)

def download_hmm_database(level):
    level_base_path = get_level_base_path(level)
    url = 'http://eggnogdb.embl.de/download/eggnog_4.5/eggnog-mapper-data/hmmdb_levels/%s/' %level_base_path
    if not args.force:
        flag = '-N'
    else:
        flag = ''
    cmd = 'mkdir -p %s; cd %s; wget %s -nH --user-agent=Mozilla/5.0 --relative -r --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off %s' %(HMMDB_PATH, HMMDB_PATH, flag, url)
    run(cmd)

def download_annotations():
    url = 'http://eggnogdb.embl.de/download/eggnog_4.5/eggnog-mapper-data/eggnog.db.gz'
    cmd = 'cd %s && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.db.gz %s && echo Decompressing... && gunzip eggnog.db.gz' %(DATA_PATH, url)
    run(cmd)

def download_groups():
    url = 'http://eggnogdb.embl.de/download/eggnog_4.5/eggnog-mapper-data/OG_fasta.tar.gz'
    cmd = 'cd %s && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O OG_fasta.tar.gz  %s && echo Decompressing... && tar -zxf OG_fasta.tar.gz && rm OG_fasta.tar.gz' %(DATA_PATH,  url)
    run(cmd)

def download_diamond_db():
    url = 'http://eggnogdb.embl.de/download/eggnog_4.5/eggnog-mapper-data/eggnog_proteins.dmnd.gz'
    cmd = 'cd %s && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog_proteins.dmnd.gz  %s && echo Decompressing... && gunzip eggnog_proteins.dmnd.gz' %(DATA_PATH,  url)
    run(cmd)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('dbs', metavar='dbs', nargs='+', choices=sorted(EGGNOG_DATABASES.keys()+['all', 'none']),
                        help='list of eggNOG HMM databases to download. Choose "none" if only diamond will be used')

    parser.add_argument('-y', action="store_true", dest='allyes',
                        help='assume "yes" to all questions')

    parser.add_argument('-f', action="store_true", dest='force',
                        help='forces download even if the files exist')

    parser.add_argument('-s', action="store_true", dest='simulate',
                        help='simulate and print commands. Nothing is downloaded')


    args = parser.parse_args()
    if 'all' in args.dbs:
        args.dbs = EGGNOG_DATABASES

    if args.force or not pexists(pjoin(DATA_PATH, 'eggnog.db')):
        if args.allyes or ask("Download main annotation database?") == 'y':
            print colorify('Downloading "eggnog.db" at %s...' %DATA_PATH, 'green')
            download_annotations()
        else:
            print 'Skipping'

    else:
        print colorify('Skipping eggnog.db database (already present). Use -f to force download', 'lblue')

    if args.force or not pexists(pjoin(DATA_PATH, 'OG_fasta')):
        if args.allyes or ask("Download OG fasta files for annotation refinement (~20GB after decompression)?") == 'y':
            print colorify('Downloading fasta files " at %s/OG_fasta...' %DATA_PATH, 'green')
            download_groups()
        else:
            print 'Skipping'

    else:
        print colorify('Skipping OG_fasta/ database (already present). Use -f to force download', 'lblue')

    if args.force or not pexists(pjoin(DATA_PATH, 'eggnog_proteins.dmnd')):
        if args.allyes or ask("Download diamond database (~4GB after decompression)?") == 'y':
            print colorify('Downloading fasta files " at %s/eggnog_proteins.dmnd...' %DATA_PATH, 'green')
            download_diamond_db()
        else:
            print 'Skipping'
    else:
        print colorify('Skipping diamond database (or already present). Use -f to force download', 'lblue')

    if set(args.dbs) != set(['none']):
        if args.allyes or ask("Download %d HMM database(s): %s?"%(len(args.dbs), ','.join(args.dbs))) == 'y':
            for db in args.dbs:
                print colorify('Downloading %s HMM database " at %s/%s\_hmm ...' %(db, HMMDB_PATH, db), 'green')
                download_hmm_database(db)
        else:
            print 'Skipping'
