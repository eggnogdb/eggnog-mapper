##

# from ..utils import colorify

from . import db_sqlite as db_sqlite

def get_member_orthologs(member, target_levels, all_nogs):
    
    # Try to setup orthology using best OG
    
    orthology = __setup_orthology(member, target_levels)
    if orthology is not None and len(orthology) > 0:
        all_orthologs = __load_orthology(orthology)
        best_OG = None
    else:
        
        # If no orthology from best OG, try using other NOGs from narrowest to widest

        for nog in reversed(all_nogs):
            nog_level = nog.split("|")[0].split("@")[1]
            orthology = __setup_orthology(member, [nog_level])
            
            if orthology is not None and len(orthology) > 0:
                # print(colorify(f"Warning: we found no orthologs for auto best OG for {member}. Using OG {nog} for annotation.", 'orange'))
                all_orthologs = __load_orthology(orthology)
                best_OG = nog
                break

        # If no orthology found, use seed ortholog for annotation
        
        if orthology is None or len(orthology) == 0:
            # print(colorify(f"Warning: we found no orthologs for {member}. Using seed ortholog for annotation.", 'orange'))
            all_orthologs = {
                "one2one": {member},
                "one2many": set(),
                "many2many": set(),
                "many2one": set(),
                "all": {member},
            }
            best_OG = f"seed_ortholog@{member}|-"

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

def __setup_orthology(member, target_levels):
    orthology = {}
    
    member_as_set = set([member])
    
    for level, _side1, _side2 in db_sqlite.get_member_events(member.strip(), target_levels):

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
