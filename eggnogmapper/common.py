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
from .emapperException import EmapperException

try:
    from .version import __VERSION__
    from .version import __DB_VERSION__
except ImportError:
    __VERSION__ = 'unknown'
    __DB_VERSION__ = 'unknown'

# Input types
ITYPE_CDS = "CDS"
ITYPE_PROTS = "proteins"
ITYPE_GENOME = "genome"
ITYPE_META = "metagenome"

# ANNOTATIONS_HEADER = list(map(str.strip, 'Preferred_name, GOs, EC, KEGG_ko, KEGG_Pathway, KEGG_Module, KEGG_Reaction, KEGG_rclass, BRITE, KEGG_TC, CAZy, BiGG_Reaction'.split(',')))

TIMEOUT_LOAD_SERVER = 10

LEVEL_CONTENT = {
"NOG":["arNOG","bactNOG","euNOG","thaNOG","eurNOG","creNOG","synNOG","spiNOG","firmNOG","fusoNOG","aquNOG","aciNOG","therNOG","tenNOG","proNOG","defNOG","plaNOG","actNOG","chloNOG","cyaNOG","deiNOG","bctoNOG","chlNOG","chlaNOG","verNOG","kinNOG","perNOG","virNOG","apiNOG","opiNOG","arcNOG","metNOG","methNOG","thermNOG","methaNOG","halNOG","theNOG","negNOG","cloNOG","eryNOG","bacNOG","acidNOG","delNOG","gproNOG","aproNOG","bproNOG","chlorNOG","dehNOG","cytNOG","bacteNOG","sphNOG","flaNOG","verrNOG","strNOG","chloroNOG","acoNOG","cocNOG","meNOG","fuNOG","dproNOG","eproNOG","braNOG","lilNOG","haeNOG","cryNOG","biNOG","basNOG","ascNOG","poaNOG","nemNOG","artNOG","chorNOG","agarNOG","treNOG","saccNOG","euroNOG","sordNOG","dotNOG","chrNOG","inNOG","veNOG","agaNOG","sacNOG","debNOG","eurotNOG","onyNOG","hypNOG","magNOG","sorNOG","pleNOG","rhaNOG","lepNOG","dipNOG","hymNOG","fiNOG","aveNOG","maNOG","arthNOG","necNOG","chaNOG","droNOG","spriNOG","carNOG","prNOG","roNOG","homNOG"],
"arNOG":["thaNOG","eurNOG","creNOG","arcNOG","metNOG","methNOG","thermNOG","methaNOG","halNOG","theNOG"],
"bactNOG":["synNOG","spiNOG","firmNOG","fusoNOG","aquNOG","aciNOG","therNOG","tenNOG","proNOG","defNOG","plaNOG","actNOG","chloNOG","cyaNOG","deiNOG","bctoNOG","chlNOG","chlaNOG","verNOG","negNOG","cloNOG","eryNOG","bacNOG","acidNOG","delNOG","gproNOG","aproNOG","bproNOG","chlorNOG","dehNOG","cytNOG","bacteNOG","sphNOG","flaNOG","verrNOG","dproNOG","eproNOG"],
"euNOG":["kinNOG","perNOG","virNOG","apiNOG","opiNOG","strNOG","chloroNOG","acoNOG","cocNOG","meNOG","fuNOG","braNOG","lilNOG","haeNOG","cryNOG","biNOG","basNOG","ascNOG","poaNOG","nemNOG","artNOG","chorNOG","agarNOG","treNOG","saccNOG","euroNOG","sordNOG","dotNOG","chrNOG","inNOG","veNOG","agaNOG","sacNOG","debNOG","eurotNOG","onyNOG","hypNOG","magNOG","sorNOG","pleNOG","rhaNOG","lepNOG","dipNOG","hymNOG","fiNOG","aveNOG","maNOG","arthNOG","necNOG","chaNOG","droNOG","spriNOG","carNOG","prNOG","roNOG","homNOG"],
"eurNOG":["arcNOG","metNOG","methNOG","thermNOG","methaNOG","halNOG","theNOG"],
"firmNOG":["negNOG","cloNOG","eryNOG","bacNOG"],
"aciNOG":["acidNOG"],
"proNOG":["delNOG","gproNOG","aproNOG","bproNOG","dproNOG","eproNOG"],
"chloNOG":["chlorNOG","dehNOG"],
"bctoNOG":["cytNOG","bacteNOG","sphNOG","flaNOG"],
"verNOG":["verrNOG"],
"virNOG":["strNOG","chloroNOG","braNOG","lilNOG","poaNOG"],
"apiNOG":["acoNOG","cocNOG","haeNOG","cryNOG"],
"opiNOG":["meNOG","fuNOG","biNOG","basNOG","ascNOG","nemNOG","artNOG","chorNOG","agarNOG","treNOG","saccNOG","euroNOG","sordNOG","dotNOG","chrNOG","inNOG","veNOG","agaNOG","sacNOG","debNOG","eurotNOG","onyNOG","hypNOG","magNOG","sorNOG","pleNOG","rhaNOG","lepNOG","dipNOG","hymNOG","fiNOG","aveNOG","maNOG","arthNOG","necNOG","chaNOG","droNOG","spriNOG","carNOG","prNOG","roNOG","homNOG"],
"delNOG":["dproNOG","eproNOG"],
"strNOG":["braNOG","lilNOG","poaNOG"],
"acoNOG":["haeNOG"],
"cocNOG":["cryNOG"],
"meNOG":["biNOG","nemNOG","artNOG","chorNOG","chrNOG","inNOG","veNOG","rhaNOG","lepNOG","dipNOG","hymNOG","fiNOG","aveNOG","maNOG","droNOG","spriNOG","carNOG","prNOG","roNOG","homNOG"],
"fuNOG":["basNOG","ascNOG","agarNOG","treNOG","saccNOG","euroNOG","sordNOG","dotNOG","agaNOG","sacNOG","debNOG","eurotNOG","onyNOG","hypNOG","magNOG","sorNOG","pleNOG","arthNOG","necNOG","chaNOG"],
"lilNOG":["poaNOG"],
"biNOG":["nemNOG","artNOG","chorNOG","chrNOG","inNOG","veNOG","rhaNOG","lepNOG","dipNOG","hymNOG","fiNOG","aveNOG","maNOG","droNOG","spriNOG","carNOG","prNOG","roNOG","homNOG"],
"basNOG":["agarNOG","treNOG","agaNOG"],
"ascNOG":["saccNOG","euroNOG","sordNOG","dotNOG","sacNOG","debNOG","eurotNOG","onyNOG","hypNOG","magNOG","sorNOG","pleNOG","arthNOG","necNOG","chaNOG"],
"nemNOG":["chrNOG","rhaNOG"],
"artNOG":["inNOG","lepNOG","dipNOG","hymNOG","droNOG"],
"chorNOG":["veNOG","fiNOG","aveNOG","maNOG","spriNOG","carNOG","prNOG","roNOG","homNOG"],
"agarNOG":["agaNOG"],
"saccNOG":["sacNOG","debNOG"],
"euroNOG":["eurotNOG","onyNOG","arthNOG"],
"sordNOG":["hypNOG","magNOG","sorNOG","necNOG","chaNOG"],
"dotNOG":["pleNOG"],
"chrNOG":["rhaNOG"],
"inNOG":["lepNOG","dipNOG","hymNOG","droNOG"],
"veNOG":["fiNOG","aveNOG","maNOG","spriNOG","carNOG","prNOG","roNOG","homNOG"],
"onyNOG":["arthNOG"],
"hypNOG":["necNOG"],
"sorNOG":["chaNOG"],
"dipNOG":["droNOG"],
"maNOG":["spriNOG","carNOG","prNOG","roNOG","homNOG"],
"spriNOG":["prNOG","roNOG","homNOG"],
"prNOG":["homNOG"]
}


LEVEL_HIERARCHY = '((thaNOG,(arcNOG,metNOG,methNOG,thermNOG,methaNOG,halNOG,theNOG)eurNOG,creNOG)arNOG,(synNOG,spiNOG,(negNOG,cloNOG,eryNOG,bacNOG)firmNOG,fusoNOG,aquNOG,(acidNOG)aciNOG,therNOG,tenNOG,((cytNOG,bacteNOG,sphNOG,flaNOG)bctoNOG,chlNOG)NoName,((dproNOG,eproNOG)delNOG,gproNOG,aproNOG,bproNOG)proNOG,defNOG,plaNOG,actNOG,(chlorNOG,dehNOG)chloNOG,cyaNOG,(chlaNOG,(verrNOG)verNOG)NoName,deiNOG)bactNOG,(kinNOG,perNOG,(((braNOG,(poaNOG)lilNOG)NoName)strNOG,chloroNOG)virNOG,((haeNOG)acoNOG,(cryNOG)cocNOG)apiNOG,(((((((lepNOG,(droNOG)dipNOG,hymNOG)NoName)inNOG)artNOG,(((fiNOG,(aveNOG,((((homNOG)prNOG,roNOG)spriNOG,carNOG)NoName)maNOG)NoName)NoName)veNOG)chorNOG)NoName,((rhaNOG)chrNOG)nemNOG)biNOG)meNOG,(((((agaNOG)agarNOG,treNOG)NoName)basNOG,(((((eurotNOG,(arthNOG)onyNOG)NoName)euroNOG,((magNOG,(chaNOG)sorNOG)NoName,(necNOG)hypNOG)sordNOG,(pleNOG)dotNOG)NoName,((sacNOG,debNOG)NoName)saccNOG)NoName)ascNOG)NoName)fuNOG)opiNOG)euNOG)NOG;'
EGGNOG_DATABASES = {k:51700+(i*2) for i, k in enumerate('NOG,aciNOG,acidNOG,acoNOG,actNOG,agaNOG,agarNOG,apiNOG,aproNOG,aquNOG,arNOG,arcNOG,artNOG,arthNOG,ascNOG,aveNOG,bacNOG,bactNOG,bacteNOG,basNOG,bctoNOG,biNOG,bproNOG,braNOG,carNOG,chaNOG,chlNOG,chlaNOG,chloNOG,chlorNOG,chloroNOG,chorNOG,chrNOG,cloNOG,cocNOG,creNOG,cryNOG,cyaNOG,cytNOG,debNOG,defNOG,dehNOG,deiNOG,delNOG,dipNOG,dotNOG,dproNOG,droNOG,eproNOG,eryNOG,euNOG,eurNOG,euroNOG,eurotNOG,fiNOG,firmNOG,flaNOG,fuNOG,fusoNOG,gproNOG,haeNOG,halNOG,homNOG,hymNOG,hypNOG,inNOG,kinNOG,lepNOG,lilNOG,maNOG,magNOG,meNOG,metNOG,methNOG,methaNOG,necNOG,negNOG,nemNOG,onyNOG,opiNOG,perNOG,plaNOG,pleNOG,poaNOG,prNOG,proNOG,rhaNOG,roNOG,sacNOG,saccNOG,sorNOG,sordNOG,sphNOG,spiNOG,spriNOG,strNOG,synNOG,tenNOG,thaNOG,theNOG,therNOG,thermNOG,treNOG,veNOG,verNOG,verrNOG,virNOG,viruses'.split(','))}
EGGNOG_DATABASES.update({'euk':51400, 'bact':51500, 'arch':51600})


TAX_LVL_6 = ["84992", "225057", "422676", "135624", "91061", "976", "213481", "204458", "1090", "32061", "135613", "186801", "84998", "1117", "301297", "28221", "29547", "526524", "119069", "118969", "135618", "909932", "206351", "206350", "135619", "414999", "135625", "4776", "121069", "206389", "204441", "766", "84995", "204457", "189775", "72273", "203494", "135623", "135614"]
TAX_LVL_5 = ["204432", "201174", "28211", "5794", "183980", "2836", "28216", "204428", "200795", "3041", "5878", "28889", "1297", "1239", "4751", "1236", "142182", "183963", "5653", "11989", "10477", "33208", "183925", "183939", "224756", "11157", "10662", "76804", "464095", "203682", "10744", "157897", "10699", "35493", "544448", "651137", "183968", "183967", "675063", "74201"]
TAX_LVL_4 = ["57723", "554915", "200783", "423358", "28883", "200930", "28890", "10474", "32066", "10404", "548681", "10860", "1511857", "10841", "40117", "33154", "1224", "11632", "203691", "508458", "10656", "200940", "200918", "33090"]
TAX_LVL_321 = ["2157", "2", "2759", "10239", "1"]
TAX_SCOPE_AUTO = TAX_LVL_6 + TAX_LVL_5 + TAX_LVL_4 + TAX_LVL_321

TAX_SCOPE_AUTO_BROAD = ["10239", # viruses
                        "5794", #apicomplexa
                        "33090", # plants
                        "6231", "6656", # nematods, arthopods
                        "40674", "78", "8782", # mammals, fishes, avian
                        "33208", # metazoa
                        "4751", # fungi
                        "33154", # opithokonta
                        '2759', '2157', '2', # euk, arch, bact
                        '1']

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
def get_eggnog_dmnd_db(): return pjoin(DATA_PATH, "eggnog_proteins.dmnd")
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

def get_oglevels_file(): return pjoin(DATA_PATH, "og2level.tsv.gz")

def set_data_path(data_path):
    global DATA_PATH
    DATA_PATH = existing_dir(data_path)
    # show_binaries()

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

def show_binaries():
    for e in (HMMSEARCH, HMMSCAN, HMMSTAT, HMMPGMD, PHMMER, DIAMOND, DATA_PATH,
              get_fasta_path(), get_hmmdb_path(), get_eggnogdb_file(), get_oglevels_file(), get_eggnog_dmnd_db()):
        print("% 65s" %e, pexists(e))

def get_call_info():
    text = []
    text.append('# ' + time.ctime())
    text.append('# ' + get_version())
    text.append('# ' + ' '.join(sys.argv))
    text.append('#')
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
            _version += "-{}".format(bytes.decode(out).rstrip())
        else:
            # If tags were not available
            # Use a short hash for the current commit (e.g. b2d12f4)
            p = Popen(["git", "rev-parse", "--short", "HEAD"], stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()

            if p.returncode == 0:
                _version += "-{}".format(bytes.decode(out).rstrip())
    if _version != '':
        version = 'emapper' + _version
    else:
        version = 'emapper-'+ __VERSION__
        
    return version

def get_db_version():
    from .annotation import db_sqlite
    db_sqlite.connect()
    return db_sqlite.get_db_version()

def get_diamond_version():
    dmnd_version = None
    cmd = f"{LOCAL_DIAMOND} --version"
    try:
        completed_process = run(cmd, capture_output=True, check=True, shell=True)
    except CalledProcessError as cpe:
        raise EmapperException("Error running local diamond: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])

    if completed_process is not None:
        dmnd_version = f"Local diamond version: {completed_process.stdout.decode('utf-8').strip()}"
    
    return dmnd_version

def get_mmseqs_version():
    mmseqs_version = None
    cmd = f"{LOCAL_MMSEQS2} version"
    try:
        completed_process = run(cmd, capture_output=True, check=True, shell=True)
    except CalledProcessError as cpe:
        raise EmapperException("Error running local mmseqs: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])

    if completed_process is not None:
        mmseqs_version = f"Local MMseqs2 version: {completed_process.stdout.decode('utf-8').strip()}"

    return mmseqs_version


def get_db_info(level):
    return (get_hmmer_dbpath(level), EGGNOG_DATABASES[level])

def get_db_present(level):
    dbpath, port = get_db_info(level)
    db_present = all([pexists(dbpath + "." + ext) for ext in 'h3f h3i h3m h3p idmap'.split()])
    return db_present

def get_citation(addons=['hmmer']):
    EXAMPLE = """
e.g. Functional annotation was performed using %s [1]
 based on eggNOG orthology data [2]. Sequence searches were performed using [3].
"""%get_version()
    
    CITATION = """
================================================================================
CITATION:
If you use this software, please cite:

[1] Fast genome-wide functional annotation through orthology assignment by
      eggNOG-mapper. Jaime Huerta-Cepas, Kristoffer Forslund, Luis Pedro Coelho,
      Damian Szklarczyk, Lars Juhl Jensen, Christian von Mering and Peer Bork.
      Mol Biol Evol (2017). doi: https://doi.org/10.1093/molbev/msx148

[2] eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated
      orthology resource based on 5090 organisms and 2502 viruses. Jaime
      Huerta-Cepas, Damian Szklarczyk, Davide Heller, Ana Hernandez-Plaza,
      Sofia K Forslund, Helen Cook, Daniel R Mende, Ivica Letunic, Thomas
      Rattei, Lars J Jensen, Christian von Mering and Peer Bork. Nucleic Acids
      Research, Volume 47, Issue D1, 8 January 2019, Pages D309-D314,
      https://doi.org/10.1093/nar/gky1085 """

    if 'hmmer' in addons:
        CITATION += """

[3] Accelerated Profile HMM Searches. 
       Eddy SR. 2011. PLoS Comput. Biol. 7:e1002195. 
"""
    elif 'diamond' in addons:
        CITATION += """

[3] Fast and Sensitive Protein Alignment using DIAMOND. Buchfink B, Xie C,
       Huson DH. 2015. Nat. Methods 12, 59–60. https://doi.org/10.1038/nmeth.3176
"""
    elif 'mmseqs' in addons:
        CITATION += """

[3] MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.
       Steinegger M & Söding J. 2017. Nat. Biotech. 35, 1026–1028. https://doi.org/10.1038/nbt.3988
"""

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

# def load_nog_lineages():
#     if os.path.exists('NOG_hierarchy.pkl'):
#         nog2lineage = cPickle.load(open('NOG_hierarchy.pkl'))
#     else:
#         nog2lineage = {}
#         for line in open('../build_db/NOG_hierarchy.tsv'):
#             fields = line.strip().split('\t')
#             nog2lineage[fields[0].split('@')[0]] = map(lambda x: tuple(x.split('@')), fields[2].split(','))
#         cPickle.dump(nog2lineage, open('NOG_hierarchy.pkl', 'wb'), protocol=2)
#     return nog2lineage

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
