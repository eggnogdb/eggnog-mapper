##

def get_member_orthologs(member, best_ogs, all_nogs, eggnog_db):

    ##
    # Try to setup orthology using best OG

    # First, obtain all OGs from bestOG to the narrowestOG
    best_ogs_pos = min([all_nogs.index(best_og) for best_og in best_ogs])
    annot_ogs = all_nogs[best_ogs_pos:]
    
    orthology = __setup_orthology(member, annot_ogs, eggnog_db)
    if orthology is not None and len(orthology) > 0:
        all_orthologs = __load_orthology(member, orthology)
        best_OG = None
    else:

        ##
        # If no orthology was obtained from best OG and its children OG,
        # try looking for a valid ortholog from best OG to root (bottom-top)

        if best_ogs_pos > 0:
            bot_top_nogs = all_nogs[:best_ogs_pos]
            bot_top_nogs.reverse() # bottom-top order

            while len(bot_top_nogs) > 0:
                nog = bot_top_nogs[0]
                orthology = __setup_orthology(member, [nog], eggnog_db)

                # If a valid orthology is found, the new best OG will be this NOG
                if orthology is not None and len(orthology) > 0:
                    all_orthologs = __load_orthology(member, orthology)
                    best_OG = nog
                    break

                # If no valid orthology was found, keep looping bottom-top
                bot_top_nogs = bot_top_nogs[1:]

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


def __load_orthology(member, orthology):
    all_orthologs = {
        "one2one": set(),
        "one2many": set(),
        "many2many": set(),
        "many2one": set(),
        "all": set(),
    } # each set contains a list of taxa.sequence items

    # member
    # e.g. 1041607.K0KSN3

    # orthology
    # k: (species, (tuple_of_co-orthologs))
    # v: set of species with orthologs:  set((species1, list_orths), (species2, list_orths), ...)

    # e.g. of many2one relationship
    # {
    #     ('1041607', ('1041607.K0KPV8', '1041607.K0KSN3')): {
    #         ('27289', ('27289.XP_003672691.1',)),
    #         ('4956', ('4956.XP_002498479.1',)),
    #         ('381046', ('381046.XP_002553862.1',)),
    #         ('28985', ('28985.XP_453001.1',)),
    #         ('113608', ('113608.XP_003687705.1',)),
    #         ('588726', ('588726.J7S748',)),
    #         ('27288', ('27288.XP_003673468.1',)),
    #         ('36033', ('36033.XP_001645185.1',)),
    #         ('4932', ('4932.YGR015C',)),
    #         ('5478', ('5478.XP_446354.1',)),
    #         ('432096', ('432096.XP_003958079.1',)),
    #         ('4950', ('4950.XP_003681095.1',))
    #     }
    # }
    
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

def __setup_orthology(member, ogs, eggnog_db):
    orthology = {}
    
    member_as_set = set([member])

    ogs_tax_ids = set([og[1] for og in ogs])
    
    for level, _side1, _side2 in eggnog_db.get_member_events(member.strip(), ogs_tax_ids):
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
