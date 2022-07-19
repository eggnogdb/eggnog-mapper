##
# CPCantalapiedra 2020

from .cache_annotator import CacheAnnotator
from .annotator import Annotator
from .annotator_novel_fams import AnnotatorNovelFams

##
def get_cache_annotator(args):
    annotator = None

    annotator = CacheAnnotator(args)
    
    return annotator

##
def get_annotator(args, annot, excel, report_orthologs):
    annotator = None

    annotator = Annotator(args, annot, excel, report_orthologs)
    
    return annotator

##
def get_annotator_novel_fams(args, annot, excel, report_orthologs):
    annotator = None

    annotator = AnnotatorNovelFams(args, annot, excel, report_orthologs)
    
    return annotator

## END
