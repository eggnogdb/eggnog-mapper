##
## CPCantalapiedra 2021

import os

from ...emapperException import EmapperException
from ...common import get_tax_scopes_path

TAX_SCOPE_MODE_BROADEST = "broadest"
TAX_SCOPE_MODE_INNER_BROADEST = "inner_broadest"
TAX_SCOPE_MODE_INNER_NARROWEST = "inner_narrowest"
TAX_SCOPE_MODE_NARROWEST = "narrowest"

from .vars import LEVEL_DEPTH, LEVEL_DICT, LEVEL_NAMES, LEVEL_PARENTS

def print_taxa():
    print("tax_name\ttax_id\tdepth\tparents\tparents_names")
    for tax_name, tax_id in LEVEL_DICT.items():
        depth = LEVEL_DEPTH.get(tax_id, "-")
        parents = LEVEL_PARENTS.get(tax_id, "-")
        parents_names = [LEVEL_NAMES.get(x, "-").replace(",", "") for x in parents]
        print(f"{tax_name}\t{tax_id}\t{depth}\t{','.join(parents)}\t{','.join(parents_names)}")
    return

##
# Parse a tax scope file
def parse_tax_scope_file(tax_scope_file):
    tax_scope_ids = []
    with open(tax_scope_file, 'r') as f:
        for line in f:
            tax_scope_ids.append(line.strip())
    return tax_scope_ids

##
# Parses tax_scope command line argument
# to define the list of tax_scope_ids (one or more tax IDs)
def parse_tax_scope(tax_scope):
    tax_scope_ids = None

    # 1st option, check if there is an actual file with that path
    # e.g. tests/fixtures/tax_scope.bacteria_broad
    if os.path.exists(tax_scope) and os.path.isfile(tax_scope):
        # parse tax scope file
        tax_scope_ids = parse_tax_scope_file(tax_scope)        
    else:
        # 2nd option, a file within eggnogmapper/annotation/tax_scope
        tax_scope_file = os.path.join(get_tax_scopes_path(), tax_scope)
        if os.path.exists(tax_scope_file) and os.path.isfile(tax_scope_file):
            # parse tax scope file
            tax_scope_ids = parse_tax_scope_file(tax_scope_file)
            
        # 3rd option, tax scope is defined in the argument itself
        elif tax_scope is not None and tax_scope != "none":
            tax_scope_ids = tax_scope.strip().split(",")

        # 4th option, no tax scope
        elif tax_scope is None or tax_scope == "none":
            tax_scope_ids = None
            
        else:
            raise EmapperException(f"Unrecognized tax scope {tax_scope}")

    # Create a list which contains only tax IDs (tax names are translated)
    # and check that those tax IDs are recognized by eggNOG-mapper
    if tax_scope_ids is not None and len(tax_scope_ids) > 0:
        tax_scope_ids_int = []

        for tax_id in tax_scope_ids:
            if tax_id in LEVEL_NAMES:
                tax_scope_ids_int.append(tax_id)
            elif tax_id in LEVEL_DICT:
                tax_scope_ids_int.append(LEVEL_DICT[tax_id])
            else:
                raise EmapperException(f"Unrecognized tax ID, tax name or tax_scope mode: '{tax_id}'.")

        tax_scope_ids = tax_scope_ids_int
    
    return tax_scope_ids


##
def parse_nogs(match_nogs, tax_scope_mode, tax_scope_ids):
    match_nogs_full = []
    match_nogs_names = None
    narr_ogs = None
    best_ogs = None

    for nog in match_nogs:
        nog_id, nog_tax_id = nog.split("@")
        nog_name = f"{nog}|{LEVEL_NAMES.get(nog_tax_id, nog_tax_id).replace(',', ' ')}"
        
        match_nogs_full.append((nog_id, nog_tax_id, nog_name, LEVEL_DEPTH[nog_tax_id]))

    # sort by depth and then by full name
    match_nogs_full = sorted(match_nogs_full, key=lambda x: (x[3], x[2]))

    match_nogs_depths = [nog[3] for nog in match_nogs_full]
    
    ##
    # Obtain narrowest OGs
    
    max_depth = max(match_nogs_depths)
    max_depth_ogs = [nog for nog in match_nogs_full if nog[3] == max_depth]
    if len(max_depth_ogs) > 0:
        narr_ogs = max_depth_ogs

    ##
    # Obtain best OG when there is no tax scope

    # if no tax scope, base OG selection on tax scope mode alone
    if tax_scope_ids is None:
        # choose broadest OGs
        if tax_scope_mode in {TAX_SCOPE_MODE_BROADEST, TAX_SCOPE_MODE_INNER_BROADEST}:
            min_depth = min(match_nogs_depths)
            min_depth_ogs = [nog for nog in match_nogs_full if nog[3] == min_depth]
            if len(min_depth_ogs) > 0:
                best_ogs = min_depth_ogs

        # choose narrowest OGs
        elif tax_scope_mode in {TAX_SCOPE_MODE_INNER_NARROWEST, TAX_SCOPE_MODE_NARROWEST}:
            best_ogs = narr_ogs
            
        else:
            # If tax_scope_mode is another tax_scope
            tax_scope_mode_ids = parse_tax_scope(tax_scope_mode)
            if tax_scope_mode_ids is not None:
                inters = set(tax_scope_mode_ids) & set([nog[1] for nog in match_nogs_full])
                if len(inters) > 0:
                    max_depth = max([LEVEL_DEPTH[tax_id] for tax_id in inters])
                    candidate_best_ogs = [nog for nog in match_nogs_full if nog[3] == max_depth]
                    if len(candidate_best_ogs) > 0:
                        best_ogs = candidate_best_ogs
            else:
                raise EmapperException(f"Could not recognize tax scope mode {tax_scope_mode}.")

    ##
    # Obtain best OG based on tax scope
    
    else:
        # Intersection of tax scope and OGs
        candidate_best_ogs = None
        inters = set(tax_scope_ids) & set([nog[1] for nog in match_nogs_full])
        if len(inters) > 0:
            # pick best OG tax id based on tax scope mode
            if tax_scope_mode == TAX_SCOPE_MODE_BROADEST:
                min_depth = min(match_nogs_depths)
                candidate_best_ogs = [nog for nog in match_nogs_full if nog[3] == min_depth]
                
            elif tax_scope_mode == TAX_SCOPE_MODE_INNER_BROADEST:
                min_depth = min([LEVEL_DEPTH[tax_id] for tax_id in inters])
                candidate_best_ogs = [nog for nog in match_nogs_full if nog[3] == min_depth]
                
            elif tax_scope_mode == TAX_SCOPE_MODE_INNER_NARROWEST:
                max_depth = max([LEVEL_DEPTH[tax_id] for tax_id in inters])
                candidate_best_ogs = [nog for nog in match_nogs_full if nog[3] == max_depth]
                
            elif tax_scope_mode == TAX_SCOPE_MODE_NARROWEST:
                candidate_best_ogs = max_depth_ogs
            else:
                # If tax_scope_mode is another tax_scope
                tax_scope_mode_ids = parse_tax_scope(tax_scope_mode)
                if tax_scope_mode_ids is not None:
                    inters = set(tax_scope_mode_ids) & set([nog[1] for nog in match_nogs_full])
                    if len(inters) > 0:
                        max_depth = max([LEVEL_DEPTH[tax_id] for tax_id in inters])
                        candidate_best_ogs = [nog for nog in match_nogs_full if nog[3] == max_depth]
                else:                
                    raise EmapperException(f"Unrecognized tax scope mode {tax_scope_mode}.")

            if candidate_best_ogs is not None and len(candidate_best_ogs) > 0:
                best_ogs = candidate_best_ogs

    match_nogs_names = [f"{nog[2]}" for nog in match_nogs_full]

    return match_nogs_full, match_nogs_names, narr_ogs, best_ogs

## END
