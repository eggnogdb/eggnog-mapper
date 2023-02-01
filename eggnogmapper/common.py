##

import sys
import os
import time
from distutils.spawn import find_executable
from os.path import join as pjoin
from os.path import exists as pexists
import gzip
import shutil
import errno
from subprocess import Popen, PIPE, run, CalledProcessError

from .utils import colorify

try:
    from .version import __VERSION__
    from .version import __DB_VERSION__
except ImportError:
    __VERSION__ = 'unknown'
    __DB_VERSION__ = 'unknown'

# Multiprocessing start method
# check
# https://docs.python.org/3/library/multiprocessing.html
# section "Contexts and start methods"
MP_START_METHOD_SPAWN = "spawn" # Available on Unix and Windows. The default on Windows and macOS.
MP_START_METHOD_FORK = "fork" # Available on Unix only. The default on Unix.
MP_START_METHOD_FORKSERVER = "forkserver" # Available on Unix platforms which support passing file descriptors over Unix pipes.
MP_START_METHOD_DEFAULT = MP_START_METHOD_SPAWN
    
# Input types
ITYPE_CDS = "CDS"
ITYPE_PROTS = "proteins"
ITYPE_GENOME = "genome"
ITYPE_META = "metagenome"

TIMEOUT_LOAD_SERVER = 10

BASE_PATH = os.path.abspath(os.path.split(os.path.abspath(__file__))[0] + '/..')

HMMSEARCH = find_executable('hmmsearch') or pjoin(BASE_PATH, 'eggnogmapper', 'bin', 'hmmsearch')
HMMSCAN = find_executable('hmmscan') or pjoin(BASE_PATH, 'eggnogmapper', 'bin', 'hmmscan')
HMMSTAT = find_executable('hmmstat') or pjoin(BASE_PATH, 'eggnogmapper', 'bin', 'hmmstat')
HMMPGMD = find_executable('hmmpgmd') or pjoin(BASE_PATH, 'eggnogmapper', 'bin', 'hmmpgmd')
PHMMER = find_executable('phmmer') or pjoin(BASE_PATH, 'eggnogmapper', 'bin', 'phmmer')
HMMFETCH = find_executable('hmmfetch') or pjoin(BASE_PATH, 'eggnogmapper', 'bin', 'hmmfetch')
HMMPRESS = find_executable('hmmpress') or pjoin(BASE_PATH, 'eggnogmapper', 'bin', 'hmmpress')
ESL_REFORMAT = find_executable('esl-reformat') or pjoin(BASE_PATH, 'eggnogmapper', 'bin', 'esl-reformat')
LOCAL_DIAMOND = pjoin(BASE_PATH, 'eggnogmapper', 'bin', 'diamond')
DIAMOND = find_executable('diamond') or LOCAL_DIAMOND
LOCAL_MMSEQS2 = pjoin(BASE_PATH, 'eggnogmapper', 'bin', 'mmseqs')
MMSEQS2 = find_executable('mmseqs') or LOCAL_MMSEQS2
PRODIGAL = find_executable('prodigal') or pjoin(BASE_PATH, 'eggnogmapper', 'bin', 'prodigal.linux')

DATA_PATH = pjoin(BASE_PATH, "data")
def get_data_path(): return DATA_PATH
def get_eggnogdb_file(): return pjoin(DATA_PATH, "eggnog.db")
def get_ncbitaxadb_file(): return pjoin(DATA_PATH, "eggnog.taxa.db")
def get_eggnog_mmseqs_dbpath(): return pjoin(DATA_PATH, "mmseqs")
def get_eggnog_mmseqs_db(): return pjoin(DATA_PATH, "mmseqs", "mmseqs.db")
def get_pfam_dbpath(): return pjoin(DATA_PATH, "pfam")
def get_pfam_db(): return pjoin(DATA_PATH, "pfam", "Pfam-A.hmm")
def get_pfam_clans_file(): return pjoin(DATA_PATH, "pfam", "Pfam-A.clans.tsv.gz")
def get_hmmer_dbpath(dbname): return pjoin(DATA_PATH, 'hmmer', dbname, dbname+".hmm")
def get_hmmer_base_dbpath(dbname): return pjoin(DATA_PATH, 'hmmer', dbname)
def get_hmmdb_path(): return pjoin(DATA_PATH, "hmmer")
def get_OG_fasta_path(dbname, og): return pjoin(DATA_PATH, 'hmmer', dbname, f"{og}.fa")
def get_hmmer_databases(): return os.listdir(get_hmmdb_path()) if os.path.isdir(os.path.realpath(get_hmmdb_path())) else []

def get_tax_scopes_path(): return pjoin(BASE_PATH, "eggnogmapper", "annotation", "tax_scopes")

def get_oglevels_file(): return pjoin(DATA_PATH, "og2level.tsv.gz")

def set_data_path(data_path):
    global DATA_PATH
    DATA_PATH = existing_dir(data_path)

##
def cleanup_og_name(name):
    import re
    # names in the hmm databases are sometiemes not clean eggnog OG names
    # m = re.search('\w+\.((ENOG41|COG|KOG|arCOG)\w+)\.', name)
    m = re.search('.*((ENOG41|COG|KOG|arCOG)\w+)\.', name)
    if m:
        name = m.groups()[0]
    name = re.sub("^ENOG41", "", name)
    return name

def get_call_info():
    text = []
    text.append('## ' + time.ctime())
    text.append('## ' + get_version())
    text.append('## ' + ' '.join(sys.argv))
    text.append('##')
    return '\n'.join(text)


def get_full_version_info():

    version = get_version()
    
    exp_db_version = __DB_VERSION__
    if exp_db_version is not None:
        version = f"{version} / Expected eggNOG DB version: {exp_db_version}"

    db_version = None
    try:
        db_version = get_db_version()
    except Exception as e:
        print(colorify(f"There was an error retrieving eggnog-mapper DB data: {e}", 'red'))
        print(colorify("Maybe you need to run download_eggnog_data.py", 'white'))
        db_version = "unknown"
        
    if db_version is not None:
        version = f"{version} / Installed eggNOG DB version: {db_version}"

    if exp_db_version is not None and db_version is not None:
        if exp_db_version != db_version and db_version != "unknown":
            print(colorify(f"Warning: expected DB version ({exp_db_version}) is different than the one found ({db_version}).", 'red'))

    dmnd_version = get_diamond_version()
    if dmnd_version is not None:
        version = f"{version} / {dmnd_version}"

    mmseqs_version = get_mmseqs_version()
    if mmseqs_version is not None:
        version = f"{version} / {mmseqs_version}"
        
    return version


def get_version():
    _version = ''
    try:
        # If on a git repository and tags are available
        # Use a tag based code (e.g. 3.1.1b2-8-gb2d12f4)
        p = Popen(["git", "describe", "--tags"], stdout=PIPE, stderr=PIPE, cwd=BASE_PATH)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            # Git not installed or package not under git
            pass
        else:
            raise
    else:
        if p.returncode == 0:
            _version += f"-{bytes.decode(out).rstrip()}"
        else:
            # If tags were not available
            # Use a short hash for the current commit (e.g. b2d12f4)
            p = Popen(["git", "rev-parse", "--short", "HEAD"], stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()

            if p.returncode == 0:
                # prefix also with __VERSION__
                # https://github.com/eggnogdb/eggnog-mapper/issues/302
                _version += f"-{__VERSION__}-{bytes.decode(out).rstrip()}"
    if _version != '':
        version = 'emapper' + _version
    else:
        version = 'emapper-'+ __VERSION__
        
    return version

def get_db_version():
    from .annotation.db_sqlite import get_fresh_eggnog_db
    eggnog_db = get_fresh_eggnog_db(usemem = False)
    return eggnog_db.get_db_version()

def get_diamond_version():
    dmnd_version = None
    cmd = f"{DIAMOND} --version"
    try:
        completed_process = run(cmd, capture_output=True, check=True, shell=True)

        if completed_process is not None:
            dmnd_version = f"Diamond version found: {completed_process.stdout.decode('utf-8').strip()}"
            
    except CalledProcessError as cpe:
        print("Couldn't find diamond: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1], file = sys.stderr)
        dmnd_version = "Diamond was not found."
    
    return dmnd_version

def get_mmseqs_version():
    mmseqs_version = None
    cmd = f"{MMSEQS2} version"
    try:
        completed_process = run(cmd, capture_output=True, check=True, shell=True)

        if completed_process is not None:
            mmseqs_version = f"MMseqs2 version found: {completed_process.stdout.decode('utf-8').strip()}"
            
    except CalledProcessError as cpe:
        print("Couldn't find MMseqs2: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1], file = sys.stderr)
        mmseqs_version = "MMseqs2 was not found."

    return mmseqs_version


EGGNOG_DATABASES = {k:51700+(i*2) for i, k in enumerate('NOG,aciNOG,acidNOG,acoNOG,actNOG,agaNOG,agarNOG,apiNOG,aproNOG,aquNOG,arNOG,arcNOG,artNOG,arthNOG,ascNOG,aveNOG,bacNOG,bactNOG,bacteNOG,basNOG,bctoNOG,biNOG,bproNOG,braNOG,carNOG,chaNOG,chlNOG,chlaNOG,chloNOG,chlorNOG,chloroNOG,chorNOG,chrNOG,cloNOG,cocNOG,creNOG,cryNOG,cyaNOG,cytNOG,debNOG,defNOG,dehNOG,deiNOG,delNOG,dipNOG,dotNOG,dproNOG,droNOG,eproNOG,eryNOG,euNOG,eurNOG,euroNOG,eurotNOG,fiNOG,firmNOG,flaNOG,fuNOG,fusoNOG,gproNOG,haeNOG,halNOG,homNOG,hymNOG,hypNOG,inNOG,kinNOG,lepNOG,lilNOG,maNOG,magNOG,meNOG,metNOG,methNOG,methaNOG,necNOG,negNOG,nemNOG,onyNOG,opiNOG,perNOG,plaNOG,pleNOG,poaNOG,prNOG,proNOG,rhaNOG,roNOG,sacNOG,saccNOG,sorNOG,sordNOG,sphNOG,spiNOG,spriNOG,strNOG,synNOG,tenNOG,thaNOG,theNOG,therNOG,thermNOG,treNOG,veNOG,verNOG,verrNOG,virNOG,viruses'.split(','))}
EGGNOG_DATABASES.update({'euk':51400, 'bact':51500, 'arch':51600})

def get_db_info(level):
    return (get_hmmer_dbpath(level), EGGNOG_DATABASES[level])

def get_db_present(level):
    dbpath, port = get_db_info(level)
    db_present = all([pexists(dbpath + "." + ext) for ext in 'h3f h3i h3m h3p idmap'.split()])
    return db_present

def get_citation(addons=['hmmer']):
    EXAMPLE = """
e.g. Functional annotation was performed using eggNOG-mapper (version %s) [1]
"""%get_version()
    
    CITATION = """
================================================================================
CITATION:
If you use this software, please cite:

[1] eggNOG-mapper v2: functional annotation, orthology assignments, and domain 
      prediction at the metagenomic scale. Carlos P. Cantalapiedra, 
      Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021.
      Molecular Biology and Evolution, msab293, https://doi.org/10.1093/molbev/msab293
"""

    if 'novel_fams' in addons:
        CITATION += """
[2] Functional and evolutionary significance of unknown genes from uncultivated taxa. 
        Álvaro Rodríguez del Río, Joaquín Giner-Lamia, Carlos P. Cantalapiedra, 
        Jorge Botas, Ziqi Deng, Ana Hernández-Plaza, Lucas Paoli, Thomas S.B. Schmidt, 
        Shinichi Sunagawa, Peer Bork, Luis Pedro Coelho, Jaime Huerta-Cepas. 
        2022. bioRxiv 2022.01.26.477801. https://doi.org/10.1101/2022.01.26.477801
"""
        EXAMPLE += " based on novel families from [2]."
    else:
        CITATION += """
[2] eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated
      orthology resource based on 5090 organisms and 2502 viruses. Jaime
      Huerta-Cepas, Damian Szklarczyk, Davide Heller, Ana Hernandez-Plaza,
      Sofia K Forslund, Helen Cook, Daniel R Mende, Ivica Letunic, Thomas
      Rattei, Lars J Jensen, Christian von Mering and Peer Bork. Nucleic Acids
      Research, Volume 47, Issue D1, 8 January 2019, Pages D309-D314,
      https://doi.org/10.1093/nar/gky1085 
"""
        EXAMPLE += " based on eggNOG orthology data [2]."

    if 'hmmer' in addons:
        CITATION += """
[3] Accelerated Profile HMM Searches. 
       Eddy SR. 2011. PLoS Comput. Biol. 7:e1002195. 
"""
    elif 'diamond' in addons:
        CITATION += """
[3] Sensitive protein alignments at tree-of-life scale using DIAMOND.
       Buchfink B, Reuter K, Drost HG. 2021.
       Nature Methods 18, 366–368 (2021). https://doi.org/10.1038/s41592-021-01101-x
"""
    elif 'mmseqs' in addons:
        CITATION += """
[3] MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.
       Steinegger M & Söding J. 2017. Nat. Biotech. 35, 1026–1028. https://doi.org/10.1038/nbt.3988
"""

    EXAMPLE += " Sequence searches were performed using [3]."

    if 'prodigal' in addons:
        CITATION += """
[4] Prodigal: prokaryotic gene recognition and translation initiation site identification.
       Hyatt et al. 2010. BMC Bioinformatics 11, 119. https://doi.org/10.1186/1471-2105-11-119.
"""
        EXAMPLE += " Gene prediction was performed using [4]."

    CITATION += EXAMPLE
    CITATION += """

================================================================================
"""

    return CITATION


LICENSE = """
LICENSE:
[1] eggNOG-mapper is free software distributed under the GPL v2 terms.
Built-in databases (e.g. eggNOG data) might be subjected to different licensing.

[2] eggNOG v5.0 data are distributed under the terms of the Creative Commons Non-Commercial Attribution
License (http://creativecommons.org/licenses/by-nc/4.0/), which permits
unrestricted reuse, distribution, and reproduction in any medium, provided the
original work is properly cited.
"""

def gopen(fname):
    if fname.endswith('.gz'):
        return gzip.open(fname, 'rt')
    else:
        return open(fname)

def silent_rm(f):
    if pexists(f):
        os.remove(f)

def silent_cp(f, dst):
    if pexists(f):
        shutil.copy(f, dst)


def existing_file(fname):
    fname = os.path.realpath(fname)
    if os.path.isfile(fname):
        return fname
    else:
        raise TypeError('not a valid file "%s"' %fname)

def existing_dir(dname):
    dname = os.path.realpath(dname)
    if os.path.isdir(dname):
        return dname
    else:
        raise TypeError('not a valid directory "%s"' %dname)

## END
