import sys
from eggnogmapper.common import EGGNOG_DATABASES, DATA_PATH

def download_hmm_data9base(level):
    if level == 'euk':
        level = 'euk_500'
    elif level == 'bact':
        level = 'euk_50'
    elif level == 'arch':
        level = 'arch_1'
        
    url = 'http://beta-eggnogdb.embl.de/download/eggnog_4.5/eggnog-mapper-data/hmmdb_levels/%s_hmm/' %level
    cmd = 'mkdir -p %s; cd %s; wget -nH --user-agent=Mozilla/5.0 --relative -r --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off %s' %(HMMDB_PATH, HMMDB_PATH, url)
    return cmd

def download_annotations():
    url = 'http://beta-eggnogdb.embl.de/download/eggnog_4.5/eggnog-mapper-data/eggnog.db.gz' 
    cmd = 'cd %s && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.db.gz %s && gunzip eggnog.db.gz' %(DATA_PATH, url)
    return cmd
    
def download_groups():
    url = 'http://beta-eggnogdb.embl.de/download/eggnog_4.5/eggnog-mapper-data/OG_fasta.tar.gz' 
    cmd = 'cd %s && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O OG_fasta.tar.gz  %s && tar -zxf OG_fasta.tar.gz && rm OG_fasta.tar.gz' %(DATA_PATH,  url)
    return cmd

    
print download_groups()
print download_annotations()
