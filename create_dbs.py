#!/usr/bin/env python3
## CPCantalapiedra 2021

import os, sys, shutil
from argparse import ArgumentParser

from eggnogmapper.common import set_data_path, get_data_path, pexists, pjoin, existing_dir
from eggnogmapper.utils import ask, ask_name, colorify
from eggnogmapper.search.diamond.diamond import create_diamond_db
from eggnogmapper.search.mmseqs.mmseqs import create_mmseqs_db, create_mmseqs_index

if sys.version_info < (3,7):
    sys.exit('Sorry, Python < 3.7 is not supported')


def get_eggnog_proteins_file(): return pjoin(get_data_path(), "e5.proteomes.faa")
def get_eggnog_taxid_info_file(): return pjoin(get_data_path(), "e5.taxid_info.tsv")

BASE_URL = f'http://eggnog5.embl.de/download/eggnog_5.0'

def run(cmd):
    print(colorify(cmd, 'cyan'))
    if not args.simulate:
        os.system(cmd)
        
def download_proteins(data_path):
    url = BASE_URL + '/e5.proteomes.faa'
    cmd = (
        f'cd {data_path} && '
        f'wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O e5.proteomes.faa {url} && '
        f'echo Downloaded.'
    )
    run(cmd)
    return

def download_taxid_info(data_path):
    url = BASE_URL + '/e5.taxid_info.tsv'
    cmd = (
        f'cd {data_path} && '
        f'wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O e5.taxid_info.tsv {url} && '
        f'echo Downloaded.'
    )
    run(cmd)
    return

def parse_taxa_table(taxid_info_file, list_taxids, list_taxa):

    taxa_set = set()

    print("Parsing taxa table ...", file=sys.stderr)
    
    if not args.simulate:
        is_list_of_ids = False
        
        if list_taxids is not None and list_taxids != "":
            taxlist = set(map(int, list_taxids.split(",")))
            is_list_of_ids = True
        else:
            taxlist = set(map(str, list_taxa.split(",")))
            is_list_of_ids = False

        with open(taxid_info_file, 'r') as taxtable:
            for i, line in enumerate(taxtable):
                if i==0: continue

                line_data = line.strip().split("\t")
                if is_list_of_ids == True:
                    line_taxlist = set(map(int, line_data[4].strip().split(",")))
                else:
                    line_taxlist = set(map(str, line_data[3].strip().split(",")))
                inters = taxlist & line_taxlist
                if len(inters) > 0:
                    taxa_set.add(line_data[0])

    print(f"Total taxa selected: {len(taxa_set)}", file=sys.stderr)
                     
    return taxa_set

def parse_proteins(out_file, proteins_file, taxa_set):

    print("Parsing fasta file ...", file=sys.stderr)
    total_selected = 0
    total_discarded = 0
    
    if not args.simulate:
        with open(out_file, 'w') as out_fasta:
            accept = False
            with open(proteins_file, 'r') as infasta:
                for line in infasta:
                    if line.startswith(">"):
                        taxid = line[1:].split(".")[0]
                        if taxid in taxa_set:
                            accept = True
                            print(line.strip(), file=out_fasta)
                        else:
                            accept = False
                    else:
                        if accept == True:
                            print(line.strip(), file=out_fasta)
                            total_selected += 1
                        else:
                            total_discarded += 1

    print(f"Total sequences selected: {total_selected}", file=sys.stderr)
    print(f"Total sequences discarded: {total_discarded}", file=sys.stderr)

    return


##
# MAIN
if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument('-m', dest='mode', 
                        choices = ['diamond', 'mmseqs'],
                        default='diamond',
                        help=(
                            'diamond: search seed orthologs using diamond (-i is required). '
                            'mmseqs: search seed orthologs using MMseqs2 (-i is required). '
                            'Default:diamond'
                        ))

    parser.add_argument('-x', action="store_true", dest='skip_mmseqs_index',
                        help='Skip the MMseqs2 DB index.')
    
    parser.add_argument('--dbname', type=str, dest="dbname",
                        help=(
                            'A name prefix for the DB to be created. e.g. "--dbname eggnog_proteins.bact_arch". '
                        ))
    
    parser.add_argument('--taxids', type=str, dest="taxids",
                        help=(
                            'Comma-separated list of tax IDs to which the proteins belong. e.g. "--taxids 2,2157" is equivalent to "--taxa Bacteria,Archaea. '
                        ))

    parser.add_argument('--taxa', type=str, dest="taxa",
                        help=(
                            'Comma-separated list of tax names to which the proteins belong. e.g. "--taxa Bacteria,Archaea". '
                        ))

    parser.add_argument("--data_dir", metavar='', type=existing_dir,
                        help='Directory to use for DATA_PATH.')

    parser.add_argument('-y', action="store_true", dest='allyes',
                        help='assume "yes" to all questions')

    parser.add_argument('-s', action="store_true", dest='simulate',
                        help='simulate and print commands. Nothing is downloaded')

    ##
    
    args = parser.parse_args()

    if args.dbname is None or args.dbname == "":
        print(colorify(f'A prefix name for the DB to be created is required. Use the --dbname option.', 'red'))
        sys.exit(1)
        
    if (args.taxids is None or args.taxids == "") and (args.taxa is None or args.taxa == ""):
        print(colorify(f'Either --taxids or --taxa parameter is required', 'red'))
        sys.exit(1)
        
    if (args.taxids is not None and args.taxids != "") and (args.taxa is not None and args.taxa != ""):
        print(colorify(f'Use either --taxids or --taxa, not both', 'red'))
        sys.exit(1)

    ##

    if "EGGNOG_DATA_DIR" in os.environ:
        set_data_path(os.environ["EGGNOG_DATA_DIR"])

    if args.data_dir:
        set_data_path(args.data_dir)

    data_path = get_data_path()

    # http://eggnog5.embl.de/download/eggnog_5.0/e5.proteomes.faa
    if not pexists(get_eggnog_proteins_file()):
        if args.allyes or ask(f"Download eggnog5 proteins to {data_path}? ~9GB (It is required to create new databases)") == 'y':
            print(colorify(f'Downloading eggnog5 proteins file to {data_path}...', 'green'))
            download_proteins(data_path)
        else:
            print(colorify(f'eggnog5 proteins file was not found. Use --data_dir to specify another data path, or allow the download', 'red'))
            sys.exit(1)
    else:
        print(colorify(f'Using existing eggnog5 proteins file found at {get_eggnog_proteins_file()}', 'green'))

    # http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv
    if not pexists(get_eggnog_taxid_info_file()):
        if args.allyes or ask(f"Download eggnog5 taxid info table to {data_path}? (It is required to create new databases)") == 'y':
            print(colorify(f'Downloading eggnog5 taxid info table to {data_path}...', 'green'))
            download_taxid_info(data_path)
        else:
            print(colorify(f'eggnog5 taxid info table was not found. Use --data_dir to specify another data path, or allow the download', 'red'))
            sys.exit(1)
    else:
        print(colorify(f'Using existing eggnog5 taxid info table found at {get_eggnog_taxid_info_file()}', 'green'))

    taxa_set = parse_taxa_table(get_eggnog_taxid_info_file(), args.taxids, args.taxa)
    proteins_out_file = pjoin(get_data_path(), f"{args.dbname}.faa")
    parse_proteins(proteins_out_file, get_eggnog_proteins_file(), taxa_set)

    # generate DB
    dbprefix = None
    if args.mode == "diamond":
        dbprefix = pjoin(get_data_path(), f'{args.dbname}.dmnd')
        create_diamond_db(dbprefix, proteins_out_file)
        
    elif args.mode == "mmseqs":
        dbdir = pjoin(get_data_path(), f'{args.dbname}.mmseqs')
        dbprefix = pjoin(get_data_path(), f'{args.dbname}.mmseqs', f'{args.dbname}.mmseqs')
        if not os.path.isdir(os.path.realpath(dbdir)):
            os.mkdir(dbdir)
            
        create_mmseqs_db(dbprefix, proteins_out_file)
        
        if not args.skip_mmseqs_index:
            tmp_dir = pjoin(get_data_path(), f'{args.dbname}.mmseqs.tmp')
            if not os.path.isdir(os.path.realpath(tmp_dir)):
                os.mkdir(tmp_dir)        
            create_mmseqs_index(dbprefix, tmp_dir)
            # remove tmp dir
            if os.path.isdir(os.path.realpath(tmp_dir)):
                shutil.rmtree(tmp_dir)
            
    else:
        print(f"Unrecognized mode '-m {args.mode}'", file=sys.stderr)
        sys.exit(1)
        
    # remove files
    if os.path.isfile(proteins_out_file):
        os.remove(proteins_out_file)

    print(f"Created {dbprefix} DB.", file=sys.stderr)
    print("Finished", file=sys.stderr)
    
    sys.exit(0)


## END
