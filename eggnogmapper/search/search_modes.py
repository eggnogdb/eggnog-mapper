##
## CPCantalapiedra 2020

from ..emapperException import EmapperException

from .diamond.diamond import DiamondSearcher
from .hmmer.hmmer import HmmerSearcher
from .mmseqs.mmseqs import MMseqs2Searcher

SEARCH_MODE_NO_SEARCH = "no_search"
SEARCH_MODE_CACHE = "cache"
SEARCH_MODE_DIAMOND = "diamond"
SEARCH_MODE_HMMER = "hmmer"
SEARCH_MODE_MMSEQS2 = "mmseqs"

##
def get_searcher(args, mode):
    searcher = None
    
    if mode in [SEARCH_MODE_NO_SEARCH, SEARCH_MODE_CACHE]:
        searcher = None

    elif mode == SEARCH_MODE_DIAMOND:
        searcher = DiamondSearcher(args)

    elif mode == SEARCH_MODE_HMMER:
        searcher = HmmerSearcher(args)

    elif mode == SEARCH_MODE_MMSEQS2:
        searcher = MMseqs2Searcher(args)

    else:
        raise EmapperException("Unknown search mode %s" % mode)
    
    return searcher



## END
