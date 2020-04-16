##
## CPCantalapiedra 2020

from .diamond import DiamondSearcher

SEARCH_MODE_NO_SEARCH = "no_search"
SEARCH_MODE_DIAMOND = "diamond"

##
def get_searcher(args, mode):
    searcher = None
    
    if mode == SEARCH_MODE_NO_SEARCH:
        searcher = None

    elif mode == SEARCH_MODE_DIAMOND:
        searcher = DiamondSearcher(args)

    else:
        raise EmapperException("Unknown search mode %s" % mode)
    
    return searcher



## END
