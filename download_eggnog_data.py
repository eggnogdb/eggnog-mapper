import sys
from eggnogmapper.common import EGGNOG_DATABASES, DATA_PATH

def download_hmm_database(level):
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
    url = 'http://beta-eggnogdb.embl.de/download/eggnog_4.5/eggnog-mapper-data/eggnog.db' 
    cmd = 'cd %s; wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.db -f %s' %(DATA_PATH, url)
    return cmd
    
def download_groups():
    url = 'http://beta-eggnogdb.embl.de/download/eggnog_4.5/eggnog-mapper-data/OG_data.tar.gz' 
    cmd = 'cd %s; wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O OG_data.tar.gz -f %s && tar -zxf OG_data.tar.gz' %(DATA_PATH,  url)
    return cmd

    
print download_groups()
print download_annotations()
