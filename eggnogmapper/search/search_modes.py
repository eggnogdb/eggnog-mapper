##
## CPCantalapiedra 2020

from .diamond import DiamondSearcher
from .hmmer import HmmerSearcher

SEARCH_MODE_NO_SEARCH = "no_search"
SEARCH_MODE_DIAMOND = "diamond"
SEARCH_MODE_HMMER = "hmmer"

##
def get_searcher(args, mode):
    searcher = None
    
    if mode == SEARCH_MODE_NO_SEARCH:
        searcher = None

    elif mode == SEARCH_MODE_DIAMOND:
        searcher = DiamondSearcher(args)

    elif mode == SEARCH_MODE_HMMER:
        searcher = HmmerSearcher(args)

    else:
        raise EmapperException("Unknown search mode %s" % mode)
    
    return searcher



## END
