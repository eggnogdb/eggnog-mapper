##

def get_member_orthologs(member, best_ogs, all_nogs, eggnog_db):
    
    # Try to setup orthology using best OG
    
    orthology = __setup_orthology(member, best_ogs, eggnog_db)
    if orthology is not None and len(orthology) > 0:
        all_orthologs = __load_orthology(orthology)
        best_OG = None
    else:

        ##
        # If no orthology was obtained from best OG, try using other NOGs in its vicinity

        # 1) Split the list of OGs from the best OGs, in 2 parts:
        # - bot_top_nogs: bottom-top NOGs from the best OGs
        # - top_bot_nogs: top-bottom NOGs from the best OGs
        
        best_ogs_indexes = [all_nogs.index(best_og) for best_og in best_ogs]
        
        if min(best_ogs_indexes) <= len(all_nogs):
            bot_top_nogs = all_nogs[:min(best_ogs_indexes)]
            bot_top_nogs.reverse()
        else:
            bot_top_nogs = []
        
        if max(best_ogs_indexes) + 1 <= len(all_nogs) - 1:
            top_bot_nogs = all_nogs[max(best_ogs_indexes) + 1:]
        else:
            top_bot_nogs = []

        # 2) Iterate alternatively from one list to the other trying to find
        # the closesest NOG to the best OGs which yields a valid orthology
        if len(top_bot_nogs) > 0:
            curr_list = top_bot_nogs
            prev_list = bot_top_nogs
        else:
            curr_list = bot_top_nogs
            prev_list = []
            
        while len(curr_list) > 0:
            nog = curr_list[0]
            orthology = __setup_orthology(member, [nog], eggnog_db)

            # If a valid orthology is found, the new best OG will be this NOG
            if orthology is not None and len(orthology) > 0:
                all_orthologs = __load_orthology(orthology)
                best_OG = nog
                break

            # If no valid orthology was found, keep looping alternatively the 2 lists
            if len(prev_list) > 0:
                tmp_list = prev_list
                prev_list = curr_list[1:]
                curr_list = tmp_list
            # If one list was already depleted, iterate over the remaining list
            else:
                curr_list = curr_list[1:]

        ##
        # If no orthology was found after iterating over all the NOGs, use the seed ortholog for annotation
        
        if orthology is None or len(orthology) == 0:
            all_orthologs = {
                "one2one": {member},
                "one2many": set(),
                "many2many": set(),
                "many2one": set(),
                "all": {member},
            }
            best_OG = ("seed_ortholog", "-", f"seed_ortholog@{member}|-", None)

    return all_orthologs, best_OG


def __load_orthology(orthology):
    all_orthologs = {
        "one2one": set(),
        "one2many": set(),
        "many2many": set(),
        "many2one": set(),
        "all": set(),
    } # each set contains a list of taxa.sequence items
    
    # k: (species, list_of_co-orthologs_species)
    # v: set of species with orthologs:  set((species1, list_orths), (species2, list_orths), ...)
    for k, v in orthology.items():

        all_orthologs['all'].update(k[1])

        if len(k[1]) == 1:
            otype_prefix = "one2"
        else:
            otype_prefix = "many2"

        for t2, co2 in v:

            all_orthologs['all'].update(co2)

            if len(co2) == 1:
                otype = otype_prefix + "one"
            else:
                otype = otype_prefix + "many"

            all_orthologs[otype].update(k[1])
            all_orthologs[otype].update(co2)
            
    return all_orthologs

def __setup_orthology(member, best_ogs, eggnog_db):
    orthology = {}
    
    member_as_set = set([member])

    best_ogs_tax_ids = set([best_og[1] for best_og in best_ogs])
    
    for level, _side1, _side2 in eggnog_db.get_member_events(member.strip(), best_ogs_tax_ids):
        side1 = [m.split('.', 1) for m in _side1.split(',')]
        side2 = [m.split('.', 1) for m in _side2.split(',')]

        # filter by taxa (by species)
        by_sp1 = __by_species(side1)#, query_taxa)
        by_sp2 = __by_species(side2)#, query_taxa)
        
        # merge by coorthologs
        __set_coorthologs(by_sp1, by_sp2, member_as_set, orthology)
        __set_coorthologs(by_sp2, by_sp1, member_as_set, orthology)
    
    return orthology


def __set_coorthologs(by_sp1, by_sp2, target_members, orthology):
    # spX: taxa (species); coX: set of sequences (co-orthologs)
    for sp1, co1 in by_sp1.items():
        if target_members & co1:
            key1 = (sp1, tuple(sorted((co1))))

            for sp2, co2 in by_sp2.items():
                key2 = (sp2, tuple(sorted(co2)))
                orthology.setdefault(key1, set()).add(key2)            
    return


def __by_species(side):
    by_sp = {}
    # t: taxa; s: sequence
    for t, s in side:            
        mid = "%s.%s" % (t, s)
        by_sp.setdefault(t, set()).add(mid)
    return by_sp

## END
