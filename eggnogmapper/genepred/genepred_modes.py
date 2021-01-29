##
## CPCantalapiedra 2020

from ..emapperException import EmapperException

from .prodigal import ProdigalPredictor

GENEPRED_MODE_SEARCH = "search"
GENEPRED_MODE_PRODIGAL = "prodigal"
# GENEPRED_MODE_ORFS = "prodigal_all_orfs"

##
def get_predictor(args, mode):
    predictor = None
    
    if mode == GENEPRED_MODE_SEARCH:
        predictor = None

    elif mode == GENEPRED_MODE_PRODIGAL:
        predictor = ProdigalPredictor(args)

    else:
        raise EmapperException("Unknown gene prediction mode %s" % mode)
    
    return predictor

## END
