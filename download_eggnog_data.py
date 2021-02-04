#!/usr/bin/env python3

import os, sys
from argparse import ArgumentParser

from eggnogmapper.common import get_eggnogdb_file, get_ncbitaxadb_file, get_eggnog_dmnd_db, get_eggnog_mmseqs_dbpath, get_pfam_dbpath, get_hmmer_base_dbpath
from eggnogmapper.common import pexists, set_data_path, get_data_path, existing_dir, HMMPRESS
from eggnogmapper.utils import ask, ask_name, colorify
from eggnogmapper.version import __DB_VERSION__

if sys.version_info < (3,7):
    sys.exit('Sorry, Python < 3.7 is not supported')
    
BASE_URL = f'http://eggnogdb.embl.de/download/emapperdb-{__DB_VERSION__}'
EGGNOG_URL = f'http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level'
EGGNOG_DOWNLOADS_URL = 'http://eggnog5.embl.de/#/app/downloads'

def run(cmd):
    print(colorify(cmd, 'cyan'))
    if not args.simulate:
        os.system(cmd)

def gunzip_flag():
    if args.force:
        return '-f'
    return ''

##
# Annotation DBs
def download_annotations(data_path):
    url = BASE_URL + '/eggnog.db.gz'
    cmd = (
        f'cd {data_path} && '
        f'wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.db.gz {url} && '
        f'echo Decompressing... && '
        f'gunzip eggnog.db.gz {gunzip_flag()}'
    )
    run(cmd)

##
# Taxa DBs
def download_taxa(data_path):
    url = BASE_URL + '/eggnog.taxa.tar.gz'
    cmd = (
        f'cd {data_path} && '
        f'wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.taxa.tar.gz {url} && '
        f'echo Decompressing... && '
        f'tar -zxf eggnog.taxa.tar.gz && '
        f'rm eggnog.taxa.tar.gz'
    )
    run(cmd)
    
##
# Diamond DBs
def download_diamond_db(data_path):
    url = BASE_URL + '/eggnog_proteins.dmnd.gz'
    cmd = (
        f'cd {data_path} && '
        f'wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog_proteins.dmnd.gz {url} && '
        f'echo Decompressing... && '
        f'gunzip eggnog_proteins.dmnd.gz {gunzip_flag()}'
    )
    run(cmd)

##
# MMseqs2 DB
def download_mmseqs_db(data_path):
    url = BASE_URL + '/mmseqs.tar.gz'
    cmd = (
        f'cd {data_path} && '
        f'wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O mmseqs.tar.gz {url} && '
        f'echo Decompressing... && '
        f'tar -zxf mmseqs.tar.gz && '
        f'rm mmseqs.tar.gz'
    )
    run(cmd)

##
# PFAM DB
def download_pfam_db(data_path):
    url = BASE_URL + '/pfam.tar.gz'
    cmd = (
        f'cd {data_path} && '
        f'wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O pfam.tar.gz {url} && '
        f'echo Decompressing... && '
        f'tar -zxf pfam.tar.gz && '
        f'rm pfam.tar.gz'
    )
    run(cmd)
    
##
# HMMER mode DBs
def download_hmm_database(level, dbname, dbpath):
    if not os.path.exists(dbpath):
        os.makedirs(dbpath)

    baseurl = f'{EGGNOG_URL}/{level}/'
    hmmsurl = f'{baseurl}/{level}_hmms.tar.gz'
    seqsurl = f'{baseurl}/{level}_raw_algs.tar'
    
    if not args.force:
        flag = '-N'
    else:
        flag = ''

    # # Download HMM and FASTA files
    # cmd = (
    #     f'cd {dbpath}; '


    # )
    
    # run(cmd)

    # Create HMMER database
    cmd = (
        f'cd {dbpath}; '
        f'echo Downloading HMMs... && '
        f'wget {flag} -nH --user-agent=Mozilla/5.0 --relative -r --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off {hmmsurl} && '
        f'echo Decompressing HMMs... && '
        f'tar zxf {level}_hmms.tar.gz && '
        f'echo {level}/* | xargs mv -t ./ && rm -r {level} && '
        f'rm {level}_hmms.tar.gz; '
        'numf=$(find ./ | grep -c ".hmm$"); '
        'curr=0; '
        f'cat /dev/null > {dbname}.hmm_tmp; '
        'for file in $(find ./ | grep ".hmm$"); do '
        'curr=$((curr+1)); '
        'echo "merging HMMs... ${file} (${curr}/${numf})"; '
        f'cat "${{file}}" | sed -e "s/.faa.final_tree.fa//" -e "s/.faa.final_tree//" >> {dbname}.hmm_tmp; '
        'rm "${file}"; '
        'done; '
        f'mv {dbname}.hmm_tmp {dbname}.hmm; '
        f'(if [ -f {dbname}.hmm.h3i ]; then rm {dbname}.hmm.h3*; fi) && '
        'echo "hmmpress-ing HMMs... " && '
        f'{HMMPRESS} {dbname}.hmm && '
        'echo "generating idmap file... " && '
        f'cat {dbname}.hmm | grep "^NAME" | sed -e "s/^NAME *//" | awk \'{{print NR"\t"$0}}\' > {dbname}.hmm.idmap && '
        'echo "removing single OG hmm files... " && '
        f'echo ./*hmm | xargs rm; '
    )
    
    run(cmd)

    # Transform alignment files to fasta files
    cmd = (
        f'cd {dbpath}; '
        f'echo Downloading FASTAs... && '
        f'wget {flag} -nH --user-agent=Mozilla/5.0 --relative -r --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off {seqsurl} && '
        f'echo Decompressing FASTAs... && '
        f'tar xf {level}_raw_algs.tar && '
        f'echo {level}/* | xargs mv -t ./ && rm -r {level} && '
        f'rm {level}_raw_algs.tar; '
        'numf=$(find ./ | grep -c ".faa.gz$"); '
        'curr=0; '
        'for file in $(find ./ | grep ".faa.gz$"); do '
        'curr=$((curr+1)); '
        'echo "processing FASTAs...  ${file} (${curr}/${numf})"; '
        'outf=$(echo "$file" | sed "s/\.raw_alg\.faa\.gz/\.fa/"); '
        'zcat "$file" | awk \'/^>/{print; next}{gsub("-", ""); print}\' > "$outf" && '
        'rm "$file"; '
        'done'
    )
    
    run(cmd)
    
    return


##
# MAIN
if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument('-D', action="store_true", dest='skip_diamond',
                        help='Do not install the diamond database')

    parser.add_argument('-P', action="store_true", dest='pfam',
                        help='Install the Pfam database, required for de novo annotation or realignment')

    parser.add_argument('-M', action="store_true", dest='mmseqs',
                        help='Install the MMseqs2 database, required for "emapper.py -m mmseqs"')

    parser.add_argument('-H', action="store_true", dest='hmmer',
                        help='Install the HMMER database specified with "-d TAXID". Required for "emapper.py -m hmmer -d TAXID"')

    parser.add_argument('-d', type=str, dest="hmmer_dbs",
                        help=(
                            f'Tax ID of eggNOG HMM database to download. e.g. "-H -d 2" for Bacteria. Required if "-H". '
                            'Available tax IDs can be found at {EGGNOG_DOWNLOADS_URL}.'
                        ))
    
    parser.add_argument('-y', action="store_true", dest='allyes',
                        help='assume "yes" to all questions')

    parser.add_argument('-f', action="store_true", dest='force',
                        help='forces download even if the files exist')

    parser.add_argument('-s', action="store_true", dest='simulate',
                        help='simulate and print commands. Nothing is downloaded')

    parser.add_argument('-q', action="store_true", dest='quiet',
                        help='quiet_mode')

    parser.add_argument("--data_dir", metavar='', type=existing_dir,
                        help='Directory to use for DATA_PATH.')

    args = parser.parse_args()

    if "EGGNOG_DATA_DIR" in os.environ:
        set_data_path(os.environ["EGGNOG_DATA_DIR"])

    if args.data_dir:
        set_data_path(args.data_dir)

    data_path = get_data_path()

    ##
    # Annotation DB
    
    if args.force or not pexists(get_eggnogdb_file()):
        if args.allyes or ask("Download main annotation database?") == 'y':
            print(colorify(f'Downloading "eggnog.db" at {data_path}...', 'green'))
            download_annotations(data_path)
        else:
            print('Skipping')
    else:
        if not args.quiet:
            print(colorify('Skipping eggnog.db database (already present). Use -f to force download', 'lblue'))

    ##
    # NCBI taxa
    
    if args.force or not pexists(get_ncbitaxadb_file()):
        if args.allyes or ask("Download taxa database?") == 'y':
            print(colorify(f'Downloading "eggnog.taxa.db" at {data_path}...', 'green'))
            download_taxa(data_path)
        else:
            print('Skipping')
    else:
        if not args.quiet:
            print(colorify('Skipping eggnog.taxa.db database (already present). Use -f to force download', 'lblue'))
            

    ##
    # Diamond DB
    
    if not args.skip_diamond and (args.force or not pexists(get_eggnog_dmnd_db())):
        if args.allyes or ask("Download diamond database (~4GB after decompression)?") == 'y':
            print(colorify(f'Downloading fasta files " at {data_path}...', 'green'))
            download_diamond_db(data_path)
        else:
            print('Skipping')
    else:
        if not args.quiet:
            print(colorify('Skipping diamond database (or already present). Use -f to force download', 'lblue'))


    ## PFAM
    if args.pfam and (args.force or not pexists(get_pfam_dbpath())):
        if args.allyes or ask("Download pfam database (~3GB after decompression)?") == 'y':
            print(colorify(f'Downloading Pfam files " at {data_path}...', 'green'))
            download_pfam_db(data_path)
        else:
            print('Skipping')
    else:
        if not args.quiet:
            print(colorify('Skipping Pfam database (or already present). Use -P and -f to force download', 'lblue'))


    ## MMseqs
    if args.mmseqs and (args.force or not pexists(get_eggnog_mmseqs_dbpath())):
        if args.allyes or ask("Download MMseqs2 database (~10GB after decompression)?") == 'y':
            print(colorify(f'Downloading MMseqs2 files " at {data_path}...', 'green'))
            download_mmseqs_db(data_path)
        else:
            print('Skipping')
    else:
        if not args.quiet:
            print(colorify('Skipping MMseqs2 database (or already present). Use -M and -f to force download', 'lblue'))

    ## HMMER
    if args.hmmer == True:
        if args.allyes or ask(f"Download HMMER database of tax ID {args.hmmer_dbs}?") == 'y':
            
            dbname = args.hmmer_dbs if args.allyes == True else ask_name('Please, specify a non-empty name for the database (e.g. Bacteria)', args.hmmer_dbs)

            dbspath = get_hmmer_base_dbpath(dbname)
            if args.force or not pexists(dbspath):                    
                print(colorify(f'Downloading HMMER database of tax ID {args.hmmer_dbs} as "{dbname}" to {dbspath}', 'green'))
                print(colorify(f'Note that this can take a long time for large taxonomic levels', 'red'))
                download_hmm_database(args.hmmer_dbs, dbname, dbspath)
            else:
                if not args.quiet:
                    print(colorify(f'HMMER database {dbname} already present at {dbspath}. Use "-f" to force download', 'lblue'))                    
        else:
            print(colorify(f'Skipping HMMER database', 'lblue'))
    else:
        print(colorify('No HMMER database requested. Use "-H -d taxid" to download the hmmer database for taxid', 'lblue'))

    print(colorify("Finished.", "green"))
    
## END
