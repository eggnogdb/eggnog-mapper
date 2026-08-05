##
## CPCantalapiedra 2020

from os.path import join as pjoin

from ..emapperException import EmapperException

from .diamond.diamond import DiamondSearcher
from .mmseqs.mmseqs import MMseqs2Searcher

SEARCH_MODE_NO_SEARCH = "no_search"
SEARCH_MODE_DIAMOND = "diamond"
SEARCH_MODE_MMSEQS2 = "mmseqs"


def get_searcher(args, mode, data_path):
    if mode == SEARCH_MODE_NO_SEARCH:
        return None
    if mode == SEARCH_MODE_DIAMOND:
        return DiamondSearcher(args, get_eggnog_dmnd_db(args.dmnd_db, mode, data_path))
    if mode == SEARCH_MODE_MMSEQS2:
        return MMseqs2Searcher(args)
    raise EmapperException("Unknown search mode %s" % mode)


def get_eggnog_dmnd_db(dmnd_db, mode, data_path):
    if dmnd_db is not None:
        return dmnd_db
    if mode is None or mode == SEARCH_MODE_DIAMOND:
        return pjoin(data_path, "eggnog_proteins.dmnd")
    raise EmapperException(f"Unrecognized mode (-m) {mode} for get_eggnog_dmnd_db")

## END
