##

from . import db_sqlite as db_sqlite

def get_member_orthologs(member, target_taxa=None, target_levels=None):

    all_orthologs = {
        "one2one": set(),
        "one2many": set(),
        "many2many": set(),
        "many2one": set(),
        "all": set(),
    }
    # each set contains a list of taxa.sequence items
    
    orthology = __setup_orthology(member, target_taxa, target_levels)

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


def __setup_orthology(member, target_taxa, target_levels):
    orthology = {}
    
    # query_taxa = member.split('.', 1)[0]
    if target_taxa:
        target_taxa = list(map(str, target_taxa))
    member_as_set = set([member])
    
    for level, _side1, _side2 in db_sqlite.get_member_events(member.strip(), target_levels):

        side1 = [m.split('.', 1) for m in _side1.split(',')]
        side2 = [m.split('.', 1) for m in _side2.split(',')]

        # filter by taxa (by species)
        by_sp1 = __by_species(side1, target_taxa)#, query_taxa)
        by_sp2 = __by_species(side2, target_taxa)#, query_taxa)
        
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


def __by_species(side, target_taxa):#, query_taxa):
    by_sp = {}
    # t: taxa; s: sequence
    for t, s in side:            
        if not target_taxa or t in target_taxa:# or t == query_taxa:
            mid = "%s.%s" % (t, s)
            by_sp.setdefault(t, set()).add(mid)
    return by_sp

## END
