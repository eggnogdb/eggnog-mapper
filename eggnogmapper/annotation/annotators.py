##
# CPCantalapiedra 2020

from .cache_annotator import CacheAnnotator
from .annotator import Annotator

##
def get_cache_annotator(args, annot, report_orthologs):
    annotator = None

    annotator = CacheAnnotator(args, annot, report_orthologs)
    
    return annotator

##
def get_annotator(args, annot, report_orthologs):
    annotator = None

    annotator = Annotator(args, annot, report_orthologs)
    
    return annotator

## END
