from __future__ import absolute_import

import sys
import os
import time
from distutils.spawn import find_executable
from os.path import join as pjoin
from os.path import exists as pexists
import gzip
import shutil
import errno
from subprocess import Popen, PIPE

try:
    from .version import __VERSION__
except ImportError:
    __VERSION__ = 'unknown'

TIMEOUT_LOAD_SERVER = 1800

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
TAXID2LEVEL = {6656: 'artNOG', 1: 'NOG', 2: 'bactNOG', 5125: 'hypNOG', 5139: 'sorNOG', 5653: 'kinNOG', 110618: 'necNOG', 7711: 'chorNOG', 508458: 'synNOG', 639021: 'magNOG', 7214: 'droNOG', 28211: 'aproNOG', 28216: 'bproNOG', 28221: 'dproNOG', 7742: 'veNOG', 1090: 'chlNOG', 8782: 'aveNOG', 200783: 'aquNOG', 34384: 'arthNOG', 5204: 'basNOG', 147541: 'dotNOG', 6231: 'nemNOG', 147545: 'euroNOG', 200795: 'chloNOG', 6236: 'rhaNOG', 1117: 'cyaNOG', 147550: 'sordNOG', 92860: 'pleNOG', 909932: 'negNOG', 2157: 'arNOG', 5234: 'treNOG', 3699: 'braNOG', 183925: 'metNOG', 183939: 'methNOG', 204428: 'chlaNOG', 4751: 'fuNOG', 204432: 'acidNOG', 183963: 'halNOG', 183967: 'theNOG', 183968: 'thermNOG', 5794: 'apiNOG', 5796: 'cocNOG', 35493: 'strNOG', 4776: 'perNOG', 183980: 'arcNOG', 4893: 'sacNOG', 5819: 'haeNOG', 526524: 'eryNOG', 544448: 'tenNOG', 2759: 'euNOG', 1224: 'proNOG', 1236: 'gproNOG', 200918: 'therNOG', 1239: 'firmNOG', 28889: 'creNOG', 7898: 'fiNOG', 40674: 'maNOG', 9443: 'prNOG', 203494: 'verrNOG', 7399: 'hymNOG', 301297: 'dehNOG', 9989: 'roNOG', 35082: 'cryNOG', 1297: 'deiNOG', 33554: 'carNOG', 422676: 'acoNOG', 4890: 'ascNOG', 4891: 'saccNOG', 28890: 'eurNOG', 5338: 'agaNOG', 314146: 'spriNOG', 766764: 'debNOG', 119089: 'chrNOG', 32061: 'chlorNOG', 32066: 'fusoNOG', 200930: 'defNOG', 4447: 'lilNOG', 29547: 'eproNOG', 57723: 'aciNOG', 50557: 'inNOG', 651137: 'thaNOG', 33154: 'opiNOG', 9604: 'homNOG', 35718: 'chaNOG', 33090: 'virNOG', 33183: 'onyNOG', 203682: 'plaNOG', 38820: 'poaNOG', 203691: 'spiNOG', 68525: 'delNOG', 7088: 'lepNOG', 186801: 'cloNOG', 5042: 'eurotNOG', 91061: 'bacNOG', 33208: 'meNOG', 33213: 'biNOG', 200643: 'bacteNOG', 976: 'bctoNOG', 201174: 'actNOG', 74201: 'verNOG', 3041: 'chloroNOG', 155619: 'agarNOG', 7147: 'dipNOG', 117743: 'flaNOG', 117747: 'sphNOG', 224756: 'methaNOG', 768503: 'cytNOG'}
TAXONOMIC_RESOLUTION = ["apiNOG",
                        "virNOG",
                        "nemNOG", "artNOG",
                        "maNOG","fiNOG", "aveNOG",
                        "meNOG",
                        "fuNOG",
                        "opiNOG",
                        'euNOG', 'arNOG', 'bactNOG',
                        'NOG']

BASE_PATH = os.path.abspath(os.path.split(os.path.abspath(__file__))[0] + '/..')

HMMSEARCH = find_executable('hmmsearch') or pjoin(BASE_PATH, 'bin', 'hmmsearch')
HMMSCAN = find_executable('hmmscan') or pjoin(BASE_PATH, 'bin', 'hmmscan')
HMMSTAT = find_executable('hmmstat') or pjoin(BASE_PATH, 'bin', 'hmmstat')
HMMPGMD = find_executable('hmmpgmd') or pjoin(BASE_PATH, 'bin', 'hmmpgmd')
PHMMER = find_executable('phmmer') or pjoin(BASE_PATH, 'bin', 'phmmer')
DIAMOND = find_executable('diamond') or pjoin(BASE_PATH, 'bin', 'diamond')

DATA_PATH = pjoin(BASE_PATH, "data")
FASTA_PATH = pjoin(DATA_PATH, "OG_fasta")
HMMDB_PATH = pjoin(DATA_PATH, "hmmdb_levels")
EGGNOGDB_FILE = pjoin(DATA_PATH, "eggnog.db")
OGLEVELS_FILE = pjoin(DATA_PATH, "og2level.tsv.gz")
EGGNOG_DMND_DB = pjoin(DATA_PATH, "eggnog_proteins.dmnd")

def show_binaries():
    for e in (HMMSEARCH, HMMSCAN, HMMSTAT, HMMPGMD, PHMMER, DIAMOND, DATA_PATH,
              FASTA_PATH, HMMDB_PATH, EGGNOGDB_FILE, OGLEVELS_FILE, EGGNOG_DMND_DB):
        print "% 65s" %e, pexists(e)

def get_call_info():
    text = []
    text.append('# ' + time.ctime())
    text.append('# ' + get_version())
    text.append('# ' + ' '.join(sys.argv))
    text.append('#')
    return '\n'.join(text)

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
    if _version:
        version = 'emapper' + _version
    else:
        version = 'emapper-'+ __VERSION__

    return version

def get_level_base_path(level):
    if level == 'euk':
        level = 'euk_500'
    elif level == 'bact':
        level = 'bact_50'
    elif level == 'arch':
        level = 'arch_1'
    else:
        level = level+"_hmm"
    return level

def get_db_info(level):
    if level == 'euk':
        return (pjoin(HMMDB_PATH,"euk_500/euk_500.hmm"), EGGNOG_DATABASES[level])
    elif level == 'bact':
        return (pjoin(HMMDB_PATH,"bact_50/bact_50.hmm"), EGGNOG_DATABASES[level])
    elif level == 'arch':
        return (pjoin(HMMDB_PATH,"arch_1/arch_1.hmm"), EGGNOG_DATABASES[level])
    else:
        return (pjoin(HMMDB_PATH, level+"_hmm", level + "_hmm.all_hmm"), EGGNOG_DATABASES[level])


def get_citation(addons=['hmmer']):
    CITATION = """
================================================================================
CITATION:
If you use this software, please cite:

[1] Fast genome-wide functional annotation through orthology assignment by
      eggNOG-mapper. Jaime Huerta-Cepas, Kristoffer Forslund, Luis Pedro Coelho,
      Damian Szklarczyk, Lars Juhl Jensen, Christian von Mering and Peer Bork.
      Mol Biol Evol (2017). doi: https://doi.org/10.1093/molbev/msx148

[2] eggNOG 4.5: a hierarchical orthology framework with improved functional
      annotations for eukaryotic, prokaryotic and viral sequences. Jaime
      Huerta-Cepas, Damian Szklarczyk, Kristoffer Forslund, Helen Cook, Davide
      Heller, Mathias C. Walter, Thomas Rattei, Daniel R. Mende, Shinichi
      Sunagawa, Michael Kuhn, Lars Juhl Jensen, Christian von Mering, and Peer
      Bork. Nucl. Acids Res. (04 January 2016) 44 (D1): D286-D293. doi:
      https://doi.org/10.1093/nar/gkv1248
"""

    if 'hmmer' in addons:
        CITATION += """
[3] Accelerated Profile HMM Searches. PLoS Comput. Biol. 7:e1002195. Eddy SR.
       2011.
"""
    elif 'diamond' in addons:
        CITATION += """
[3] Fast and Sensitive Protein Alignment using DIAMOND. Buchfink B, Xie C,
       Huson DH. 2015. Nat. Methods [Internet] 12.
"""

    CITATION += """

(e.g. Functional annotation was performed using %s [1]
 based on eggNOG orthology data [2]. Sequence searches were performed
 using [3].)

================================================================================
""" %get_version()

    return CITATION


LICENSE = """
LICENSE:
[1] eggNOG-mapper is free software distributed under the GPL v2 terms.

[2] eggNOG data are distributed under the terms of the Creative Commons Attribution
License (http://creativecommons.org/licenses/by/4.0/), which permits
unrestricted reuse, distribution, and reproduction in any medium, provided the
original work is properly cited.
"""

def gopen(fname):
    if fname.endswith('.gz'):
        return gzip.open(fname, 'r:gz')
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
