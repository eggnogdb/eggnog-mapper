##
## CPCantalapiedra 2020

from os.path import join as pjoin

from ..emapperException import EmapperException

from .diamond.diamond import DiamondSearcher
from .hmmer.hmmer import HmmerSearcher
from .mmseqs.mmseqs import MMseqs2Searcher

SEARCH_MODE_NO_SEARCH = "no_search"
SEARCH_MODE_CACHE = "cache"
SEARCH_MODE_DIAMOND = "diamond"
SEARCH_MODE_HMMER = "hmmer"
SEARCH_MODE_MMSEQS2 = "mmseqs"
SEARCH_MODE_NOVEL_FAMS = "novel_fams"

##
def get_searcher(args, mode, data_path):
    searcher = None
    
    if mode in [SEARCH_MODE_NO_SEARCH, SEARCH_MODE_CACHE]:
        searcher = None

    elif mode == SEARCH_MODE_DIAMOND:
        searcher = DiamondSearcher(args, get_eggnog_dmnd_db(args.dmnd_db, mode, data_path))

    elif mode == SEARCH_MODE_HMMER:
        searcher = HmmerSearcher(args)

    elif mode == SEARCH_MODE_MMSEQS2:
        searcher = MMseqs2Searcher(args)

    elif mode == SEARCH_MODE_NOVEL_FAMS:
        searcher = DiamondSearcher(args, get_eggnog_dmnd_db(args.dmnd_db, mode, data_path))

    else:
        raise EmapperException("Unknown search mode %s" % mode)
    
    return searcher


#
def get_eggnog_dmnd_db(dmnd_db, mode, data_path):
    ret_dmnd_db = None

    if dmnd_db is not None:
        ret_dmnd_db = dmnd_db
        
    else:
        if mode is None or mode == SEARCH_MODE_DIAMOND:
            ret_dmnd_db = pjoin(data_path, "eggnog_proteins.dmnd")
        
        elif mode == SEARCH_MODE_NOVEL_FAMS:
            ret_dmnd_db = pjoin(data_path, "novel_fams.dmnd")
        
        else:
            raise EmapperException(f"Unrecognized mode (-m) {mode} for get_eggnog_dmnd_db")
    
    return ret_dmnd_db

## END
