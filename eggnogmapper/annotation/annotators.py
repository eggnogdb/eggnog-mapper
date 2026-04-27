"""Annotator factory.

In v3 the only supported annotator is the v7-batch `Annotator`. The legacy
`AnnotatorNovelFams` and `CacheAnnotator` were removed in Phase 2.
"""

from .annotator import Annotator


def get_annotator(args, annot, excel, report_orthologs):
    return Annotator(args, annot, excel, report_orthologs)
