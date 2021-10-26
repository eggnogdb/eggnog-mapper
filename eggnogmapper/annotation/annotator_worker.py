##
## CPCantalapiedra 2021

from collections import Counter

from ..emapperException import EmapperException
from ..common import set_data_path

from .tax_scopes.tax_scopes import parse_nogs

from .db_sqlite import get_eggnog_db

from . import orthologs as ortho
from . import annota
from . import output

ANNOTATIONS_HEADER = output.ANNOTATIONS_HEADER

##
def annotate_hit_line_mem(arguments):
    eggnog_db = get_eggnog_db(usemem = True)

    return annotate_hit_line(arguments, eggnog_db)

##
def annotate_hit_line_ondisk(arguments):
    data_dir = arguments[-2]
    set_data_path(data_dir)
    eggnog_db = get_eggnog_db(usemem = False)
    
    return annotate_hit_line(arguments, eggnog_db)

##
# Retrieve annotations and orthologs of a hit
##
def annotate_hit_line(arguments, eggnog_db):
    
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
        best_hit_name = hit[1]
        best_hit_evalue = float(hit[2])
        best_hit_score = float(hit[3])
        
        ##
        # Filter by empty hit, error, evalue and/or score
        if filter_out(best_hit_name, best_hit_evalue, best_hit_score, seed_ortholog_evalue, seed_ortholog_score):
            pass
        else:
            ##
            # Retrieve OGs (orthologs groups) the hit belongs to
            match_nogs = get_member_ogs(best_hit_name, eggnog_db)
            if not match_nogs:
                pass
            else:
                ##
                # Obtain OGs sorted by depth (aka tax level) and with name added
                # and also the narrowest OG, and the best OG according to tax scope ids list and mode

                match_nogs, match_nogs_names, narr_ogs, best_ogs = parse_nogs(match_nogs,
                                                                              tax_scope_mode,
                                                                              tax_scope_ids)
                
                if best_ogs is None:
                    pass
                else:
                    ##
                    # Retrieve co-orthologs of seed ortholog
                    try:
                        all_orthologies, best_OG = ortho.get_member_orthologs(best_hit_name,
                                                                              best_ogs,
                                                                              match_nogs,
                                                                              eggnog_db)
                        if best_OG is not None:
                            best_ogs = [best_OG]

                        # filter co-orthologs to keep only target_orthologs: "all", "one2one", ...
                        annot_orthologs = _filter_orthologs(all_orthologies,
                                                            target_orthologs,
                                                            target_taxa,
                                                            excluded_taxa)

                    except Exception as e:
                        # import traceback
                        # traceback.print_exc()
                        raise EmapperException(f'Error: orthology retrieval went wrong for hit {hit}. '+str(e))

                    ##
                    # Retrieve annotations of co-orthologs
                    if annot == True and annot_orthologs is not None and len(annot_orthologs) > 0:
                        annotations = annota.summarize_annotations(annot_orthologs,
                                                                   annotations_fields = ANNOTATIONS_HEADER,
                                                                   target_go_ev = go_evidence,
                                                                   excluded_go_ev = go_excluded,
                                                                   eggnog_db = eggnog_db)

                    else:
                        annotations = {}
                    
                    match_nogs_descriptions = get_ogs_descriptions(match_nogs, eggnog_db)

                    # best_ogs[0] because all best_ogs MUST HAVE the same tax lvl
                    max_annot_lvl = best_ogs[0][2].split("@")[1]

                    annotation = (query_name, best_hit_name, best_hit_evalue, best_hit_score,
                                  annotations,
                                  match_nogs_descriptions,
                                  max_annot_lvl,
                                  match_nogs_names,
                                  all_orthologies, annot_orthologs)
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise EmapperException(f'Error: annotation went wrong for hit {hit}. '+str(e))

    # False: new annotation, not present in output file (this is important when --resume)
    return ((hit, annotation), False)


def _filter_orthologs(all_orthologies, target_orthologs, target_taxa, excluded_taxa):
    
    orthologs = sorted(all_orthologies[target_orthologs])
    
    if excluded_taxa is not None:
        orthologs = [o for o in orthologs if int(o.split(".")[0]) not in excluded_taxa]
    if target_taxa is not None:
        orthologs = [o for o in orthologs if int(o.split(".")[0]) in target_taxa]
    return orthologs
        
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


def get_member_ogs(name, eggnog_db):
    ogs = None
    match = eggnog_db.get_member_ogs(name)
    if match is not None and match[0] is not None:
        ogs = [str(x).strip() for x in match[0].split(',')]
    return ogs


def get_ogs_descriptions(nogs, eggnog_db):
    og_name = None
    cat = None
    desc = None
    
    for nog in reversed(nogs):
        nog_id, nog_tax_id, nog_name, nog_depth = nog
        nog_cat, nog_desc = get_og_description(nog_id, nog_tax_id, eggnog_db)

        if nog_desc != "-":
            og_name = nog_name
            cat = nog_cat.replace('\n', '')
            desc = nog_desc.replace('\n', '')
            break
            
    return (og_name, cat, desc)

def get_og_description(og, tax_id, eggnog_db):
    best = ['-', '-', '-']
    
    for og, nm, desc, cat in eggnog_db.get_ogs_description(og, tax_id):
        desc = desc.strip()
        if desc and desc != 'N/A' and desc != 'NA':
            best = [nm, cat, desc]
            break
    
    return best[1], best[2]

## END
