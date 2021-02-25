##
## CPCantalapiedra 2021

import os

from ...emapperException import EmapperException
from ...common import get_tax_scopes_path

from .vars import LEVEL_DEPTH, LEVEL_DICT, LEVEL_NAMES, LEVEL_PARENTS

def print_taxa():
    print("tax_name\ttax_id\tdepth\tparents\tparents_names")
    for tax_name, tax_id in LEVEL_DICT.items():
        depth = LEVEL_DEPTH.get(tax_id, "-")
        parents = LEVEL_PARENTS.get(tax_id, "-")
        parents_names = [LEVEL_NAMES.get(x, "-") for x in parents]
        print(f"{tax_name}\t{tax_id}\t{depth}\t{','.join(parents)}\t{','.join(parents_names)}")
    return

##
# Parse a tax scope file
def parse_tax_scope_file(tax_scope_file):
    tax_scope_fields = []
    with open(tax_scope_file, 'r') as f:
        for line in f:
            tax_scope_fields.append(line.strip())
    return tax_scope_fields

##
# Parses tax_scope command line argument
# to define tax_scope_mode and tax_scope_id (one or more tax IDs)
def parse_tax_scope(tax_scope):
    tax_scope_mode = None
    tax_scope_id = None

    # 1st option, check if there is an actual file with that path
    # e.g. tests/fixtures/tax_scope.bacteria_broad
    if os.path.exists(tax_scope) and os.path.isfile(tax_scope):
        # parse tax scope file
        tax_scope_mode = "none"
        tax_scope_fields = parse_tax_scope_file(tax_scope)        
    else:
        # 2nd option, a file within eggnogmapper/annotation/tax_scope
        tax_scope_file = os.path.join(get_tax_scopes_path(), tax_scope)
        if os.path.exists(tax_scope_file) and os.path.isfile(tax_scope_file):
            # parse tax scope file
            tax_scope_mode = "none"
            tax_scope_fields = parse_tax_scope_file(tax_scope_file)
            
        # 3rd option, tax scope is defined in the argument itself
        else:
            tax_scope_fields = tax_scope.strip().split(",")
            tax_scope_mode = tax_scope_fields[0]
        
    # Narrowest
    if tax_scope_mode == "narrowest":
        tax_scope_id = None

    # Comma-separated list of tax IDs or list of tax IDs from file
    else:
        # Only the specified tax ID
        if len(tax_scope_fields) == 1:
            tax_scope_mode = "none"
            tax_scope_id = [tax_scope_fields[0]]

        # Tax ID list, with or without mode for those not found in the list
        elif len(tax_scope_fields) > 1:
            last_pos = tax_scope_fields[-1]
            if last_pos in ["narrowest", "auto", "none"]:
                tax_scope_mode = last_pos
                tax_scope_id = tax_scope_fields[:-1]
            else:
                tax_scope_mode = "none"
                tax_scope_id = tax_scope_fields
        else:
            raise EmapperException(f"Error: unrecognized tax scope format {tax_scope}.")

    # Create a list which contains only tax IDs (tax names are translated)
    # and check that those tax IDs are recognized by eggNOG-mapper
    if tax_scope_id is not None and len(tax_scope_id) > 0:
        tax_scope_id_int = []

        for tax_id in tax_scope_id:
            if tax_id in LEVEL_NAMES:
                tax_scope_id_int.append(tax_id)
            elif tax_id in LEVEL_DICT:
                tax_scope_id_int.append(LEVEL_DICT[tax_id])
            else:
                raise EmapperException(f"Unrecognized tax ID, tax name or tax_scope mode: '{tax_id}'.")

        tax_scope_id = tax_scope_id_int
    
    return tax_scope_mode, tax_scope_id


##
def parse_nogs(match_nogs, tax_scope_mode, tax_scope_id):        
    match_nogs_names = []
    best_og = None
    narr_og = None
    
    # best_og_id = None
    # best_og_tax_id = None
    # best_og_depth = None
    # narr_og_id = None
    # narr_og_tax_id = None
    # narr_og_depth = None
    
    tax_scope_id_2 = tax_scope_id

    match_nogs_sorted = sorted(match_nogs, key=lambda x: LEVEL_DEPTH[x.split("@")[1]])
    print(match_nogs_sorted)

    # Obtain narrowest OG
    narr_og_id, narr_og_tax_id = match_nogs_sorted[-1].split("@")
    narr_og_name = f"{narr_og_id}@{narr_og_tax_id}|{LEVEL_NAMES.get(narr_og_tax_id, narr_og_tax_id)}"
    narr_og = (narr_og_id, narr_og_tax_id, narr_og_name) if narr_og_id is not None else None
    
    # If use narrowest, best_og will be just the narrowest
    if tax_scope_id is None and tax_scope_mode == "narrowest":
        best_og = narr_og
        match_nog_names = [f"{nog}|{LEVEL_NAMES.get(nog.split('@')[1], nog.split('@')[1])}" for nog in match_nogs_sorted]
            
    elif tax_scope_id is None:
            raise EmapperException(f"Unrecognized tax scope mode {tax_scope_mode}")
    # Obtain best OG based on tax scope
    else:
        match_nogs_tax_ids = set([nog.split("@")[1] for nog in match_nogs_sorted])
        print(f"NOGs tax ids: {match_nogs_tax_ids}")
        print(f"Tax scope: {tax_scope_id}")
        
        for nog in match_nogs_sorted:
            nog_id = nog.split("@")[0]
            nog_tax_id = nog.split("@")[1]

            nog_name = f"{nog}|{LEVEL_NAMES.get(nog_tax_id, nog_tax_id)}"
            match_nogs_names.append(nog_name)

            nog_depth = LEVEL_DEPTH[nog_tax_id]

                    
            for i, filter_tax_id in enumerate(tax_scope_id_2):
                if filter_tax_id == nog_tax_id:
                    best_og_id = nog_id
                    best_og_level = nog_tax_id
                    best_og_depth = nog_depth
                    # only leave IDs of equal or more priority (left-most in the list)
                    tax_scope_id_2 = tax_scope_id_2[:i+1]
                    break
                    

        if best_og_id is None:
            best_og = None
        else:
            best_og_name = f"{best_og_id}@{best_og_level}|{LEVEL_NAMES.get(best_og_level, best_og_level)}"
            best_og = (best_og_id, best_og_level, best_og_name)
            
            # # Obtain narrowest OG
            # if narr_og_depth is None or nog_depth >= narr_og_depth:
            #     narr_og_id = nog_id
            #     narr_og_level = nog_tax_id
            #     narr_og_depth = nog_depth
            #     if tax_scope_id is None and tax_scope_mode == "narrowest":
            #         best_og_id = narr_og_id
            #         best_og_level = narr_og_level
            #         best_og_depth = narr_og_depth

            # # Obtain best OG based on tax scope
            # if tax_scope_id is None:
            #     if tax_scope_mode == "narrowest":
            #         pass # Already processed above

            #     else:
            #         raise EmapperException(f"Unrecognized tax scope mode {tax_scope_mode}")    

            # else: # tax_scope_id is not None:
                # for i, filter_tax_id in enumerate(tax_scope_id_2):
                #     if filter_tax_id == nog_tax_id:
                #         best_og_id = nog_id
                #         best_og_level = nog_tax_id
                #         best_og_depth = nog_depth
                #         # only leave IDs of equal or more priority (left-most in the list)
                #         tax_scope_id_2 = tax_scope_id_2[:i+1]
                #         break

    # print(match_nogs_names)
    # print(f"Best OG: {best_og_id}-{best_og_level}")

    return match_nogs_names, narr_og, best_og

## END
