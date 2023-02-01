##
## CPCantalapiedra 2021

from ..emapperException import EmapperException

# from . import orthologs as ortho
# from . import annota

##
# Retrieve annotations and orthologs of a hit
##
def annotate_hit_line(arguments):
    
    hit, annot, seed_ortholog_score, seed_ortholog_evalue, \
        tax_scope_mode, tax_scope_ids, \
        target_taxa, target_orthologs, excluded_taxa, \
        go_evidence, go_excluded, data_dir, annotation = arguments

    # For hits already annotated (this is important when --resume),
    # just return the annotation back
    if annotation is not None:
        return ((hit, annotation), True) # and mark as already present in output file
    
    try:
        query_name = hit[0]
        best_hit_name = "_".join(hit[1].split("_")[1:])
        best_hit_evalue = float(hit[2])
        best_hit_score = float(hit[3])
        
        ##
        # Filter by empty hit, error, evalue and/or score
        if filter_out(best_hit_name, best_hit_evalue, best_hit_score, seed_ortholog_evalue, seed_ortholog_score):
            pass
        else:
            novel_fam = hit[1].split("_")[0]
            annotation = (query_name, best_hit_name, best_hit_evalue, best_hit_score, novel_fam)
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise EmapperException(f'Error: annotation went wrong for hit {hit}. '+str(e))

    # False: new annotation, not present in output file (this is important when --resume)
    return ((hit, annotation), False)
        
##
def filter_out(hit_name, hit_evalue, hit_score, threshold_evalue, threshold_score):
    """
    Filter hit if ERROR, by score or by evalue
    """
    if hit_name == '-' or hit_name == 'ERROR':
        return True

    if threshold_evalue is not None and hit_evalue is not None and hit_evalue > threshold_evalue:
        return True

    if threshold_score is not None and hit_score is not None and hit_score < threshold_score:
        return True
    
    return False

## END
