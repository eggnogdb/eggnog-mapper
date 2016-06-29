from __future__ import absolute_import

import os
from distutils.spawn import find_executable
from os.path import join as pjoin
from os.path import exists as pexists

EGGNOG_DATABASES = {k:51700+(i*2) for i, k in enumerate('NOG,aciNOG,acidNOG,acoNOG,actNOG,agaNOG,agarNOG,apiNOG,aproNOG,aquNOG,arNOG,arcNOG,artNOG,arthNOG,ascNOG,aveNOG,bacNOG,bactNOG,bacteNOG,basNOG,bctoNOG,biNOG,bproNOG,braNOG,carNOG,chaNOG,chlNOG,chlaNOG,chloNOG,chlorNOG,chloroNOG,chorNOG,chrNOG,cloNOG,cocNOG,creNOG,cryNOG,cyaNOG,cytNOG,debNOG,defNOG,dehNOG,deiNOG,delNOG,dipNOG,dotNOG,dproNOG,droNOG,eproNOG,eryNOG,euNOG,eurNOG,euroNOG,eurotNOG,fiNOG,firmNOG,flaNOG,fuNOG,fusoNOG,gproNOG,haeNOG,halNOG,homNOG,hymNOG,hypNOG,inNOG,kinNOG,lepNOG,lilNOG,maNOG,magNOG,meNOG,metNOG,methNOG,methaNOG,necNOG,negNOG,nemNOG,onyNOG,opiNOG,perNOG,plaNOG,pleNOG,poaNOG,prNOG,proNOG,rhaNOG,roNOG,sacNOG,saccNOG,sorNOG,sordNOG,sphNOG,spiNOG,spriNOG,strNOG,synNOG,tenNOG,thaNOG,theNOG,therNOG,thermNOG,treNOG,veNOG,verNOG,verrNOG,virNOG'.split(','))}
EGGNOG_DATABASES.update({'euk':51400, 'bact':51500, 'arch':51600})

BASE_PATH = os.path.abspath(os.path.split(os.path.abspath(__file__))[0]+'/..')

HMMSEARCH = find_executable('hmmsearch')
HMMSCAN = find_executable('hmmscan')
HMMPGMD = find_executable('hmmpgmd')
PHMMER = find_executable('phmmer')

DATA_PATH = pjoin(BASE_PATH, "data")
FASTA_PATH = pjoin(DATA_PATH, "OG_fasta")
HMMDB_PATH = pjoin(DATA_PATH, "hmmdb_levels")
EGGNOGDB_PATH = pjoin(BASE_PATH, "db", "eggnog.db")

def get_db_info(level):
    if level == 'euk':
        return (pjoin(HMMDB_PATH,"euk_500/euk_500.hmm"), EGGNOG_DATABASES[level]) 
    elif level == 'bact':
        return (pjoin(HMMDB_PATH,"bact_50/bact_50.hmm"), EGGNOG_DATABASES[level])
    elif level == 'arch':
        return (pjoin(HMMDB_PATH,"arch_1/arch_1.hmm"), EGGNOG_DATABASES[level])
    else:
        return (pjoin(HMMDB_PATH, level+"_hmm", level + "_hmm.all_hmm"), EGGNOG_DATABASES[level])


CITATION = """
CITATION:
If you use this software, please cite:

[1] eggNOG-mapper: eggNOG-mapper: fast proteome-scale functional annotation through orthology assignments.
aime Huerta-Cepas, Kristoffer Forslund, Damian Szklarczyk, Lars Juhl Jensen, Christian von Mering and Peer Bork. In preparation.

[2] eggNOG 4.5: a hierarchical orthology framework with improved functional annotations for eukaryotic, prokaryotic and viral sequences.
Jaime Huerta-Cepas, Damian Szklarczyk, Kristoffer Forslund, Helen Cook, Davide Heller, Mathias C. Walter, Thomas Rattei, Daniel R. Mende, Shinichi Sunagawa, Michael Kuhn, Lars Juhl Jensen, Christian von Mering, and Peer Bork.
Nucl. Acids Res. (04 January 2016) 44 (D1): D286-D293. doi: 10.1093/nar/gkv1248
"""

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

def load_nog_lineages():
    if os.path.exists('NOG_hierarchy.pkl'):
        nog2lineage = cPickle.load(open('NOG_hierarchy.pkl'))
    else:
        nog2lineage = {}
        for line in open('../build_db/NOG_hierarchy.tsv'):
            fields = line.strip().split('\t')
            nog2lineage[fields[0].split('@')[0]] = map(lambda x: tuple(x.split('@')), fields[2].split(','))
        cPickle.dump(nog2lineage, open('NOG_hierarchy.pkl', 'wb'), protocol=2)
    return nog2lineage
