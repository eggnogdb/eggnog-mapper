##

def get_member_orthologs(member, best_ogs, all_nogs, eggnog_db):

    ##
    # Try to setup orthology using best OG

    # First, obtain all OGs from bestOG to the narrowestOG
    best_ogs_pos = min([all_nogs.index(best_og) for best_og in best_ogs])
    annot_ogs = all_nogs[best_ogs_pos:]

    # In int_mode, member is an integer ID string from DIAMOND (e.g. "42")
    if eggnog_db._int_mode:
        member_key = int(member)
    else:
        member_key = member

    orthology = __setup_orthology(member_key, annot_ogs, eggnog_db)
    if orthology is not None and len(orthology) > 0:
        all_orthologs = __load_orthology(member_key, orthology)
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
                orthology = __setup_orthology(member_key, [nog], eggnog_db)

                # If a valid orthology is found, the new best OG will be this NOG
                if orthology is not None and len(orthology) > 0:
                    all_orthologs = __load_orthology(member_key, orthology)
                    best_OG = nog
                    break

                # If no valid orthology was found, keep looping bottom-top
                bot_top_nogs = bot_top_nogs[1:]

        ##
        # If no orthology was found after iterating over all the NOGs, use the seed ortholog for annotation

        if orthology is None or len(orthology) == 0:
            all_orthologs = {
                "one2one": {member_key},
                "one2many": set(),
                "many2many": set(),
                "many2one": set(),
                "all": {member_key},
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
    }

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

    if eggnog_db._int_mode:
        # Integer mode: side1/side2 are Python int lists (decoded from BLOB)
        _get_taxid = eggnog_db.get_taxid
        for level, side1_ids, side2_ids in eggnog_db.get_member_events(member, ogs_tax_ids):
            by_sp1 = __by_species_int(side1_ids, _get_taxid)
            by_sp2 = __by_species_int(side2_ids, _get_taxid)

            __set_coorthologs(by_sp1, by_sp2, member_as_set, orthology)
            __set_coorthologs(by_sp2, by_sp1, member_as_set, orthology)
    else:
        # Legacy mode: side1/side2 contain "TAXID.PROTNAME" CSV strings
        for level, _side1, _side2 in eggnog_db.get_member_events(member.strip(), ogs_tax_ids):
            side1 = [m.split('.', 1) for m in _side1.split(',')]
            side2 = [m.split('.', 1) for m in _side2.split(',')]

            by_sp1 = __by_species(side1)
            by_sp2 = __by_species(side2)

            __set_coorthologs(by_sp1, by_sp2, member_as_set, orthology)
            __set_coorthologs(by_sp2, by_sp1, member_as_set, orthology)

    return orthology


def __set_coorthologs(by_sp1, by_sp2, target_members, orthology):
    # spX: taxa (species); coX: set of sequences/IDs (co-orthologs)
    for sp1, co1 in by_sp1.items():
        if target_members & co1:
            key1 = (sp1, tuple(sorted((co1))))

            for sp2, co2 in by_sp2.items():
                key2 = (sp2, tuple(sorted(co2)))
                orthology.setdefault(key1, set()).add(key2)
    return


def __by_species(side):
    """Legacy: group [taxid, protname] pairs by species. Values are "TAXID.PROTNAME" strings."""
    by_sp = {}
    for t, s in side:
        mid = "%s.%s" % (t, s)
        by_sp.setdefault(t, set()).add(mid)
    return by_sp


def __by_species_int_arr(ids, taxid_arr):
    """Integer mode: group protein integer IDs by species using numpy array.

    Optimized version: direct array indexing, no function call overhead.
    """
    from collections import defaultdict
    if not ids:
        return {}

    by_sp = defaultdict(set)
    max_idx = len(taxid_arr) - 1
    for pid in ids:
        idx = pid if pid <= max_idx else 0
        taxid = taxid_arr[idx]
        by_sp[taxid].add(pid)
    return {str(k): v for k, v in by_sp.items()}


def __by_species_int(ids, get_taxid):
    """Integer mode: group protein integer IDs by species (taxid from in-memory array).

    Legacy wrapper - calls optimized version if taxid array is available.
    """
    from collections import defaultdict
    by_sp = defaultdict(set)
    for pid in ids:
        taxid = get_taxid(pid)
        by_sp[taxid].add(pid)
    return {str(k): v for k, v in by_sp.items()}

def resolve_orthologs_cached(member_id, best_ogs, all_nogs,
                              cached_events, event_ids, get_taxid):
    """Same logic as get_member_orthologs but using pre-fetched events.

    cached_events: {event_i: (level, side1_ids, side2_ids)}
    event_ids: list of event IDs for this member (from event_index)
    get_taxid: callable(protein_id) → taxid (int)
    """
    best_ogs_pos = min([all_nogs.index(best_og) for best_og in best_ogs])
    annot_ogs = all_nogs[best_ogs_pos:]

    orthology = _setup_orthology_cached(member_id, annot_ogs,
                                         cached_events, event_ids, get_taxid)
    if orthology is not None and len(orthology) > 0:
        all_orthologs = __load_orthology(member_id, orthology)
        best_OG = None
    else:
        if best_ogs_pos > 0:
            bot_top_nogs = all_nogs[:best_ogs_pos]
            bot_top_nogs.reverse()

            while len(bot_top_nogs) > 0:
                nog = bot_top_nogs[0]
                orthology = _setup_orthology_cached(member_id, [nog],
                                                     cached_events, event_ids, get_taxid)
                if orthology is not None and len(orthology) > 0:
                    all_orthologs = __load_orthology(member_id, orthology)
                    best_OG = nog
                    break
                bot_top_nogs = bot_top_nogs[1:]

        if orthology is None or len(orthology) == 0:
            all_orthologs = {
                "one2one": {member_id},
                "one2many": set(),
                "many2many": set(),
                "many2one": set(),
                "all": {member_id},
            }
            best_OG = ("seed_ortholog", "-", f"seed_ortholog@{member_id}|-", None)

    return all_orthologs, best_OG


_orthology_stats = {'calls': 0, 'events': 0}

def _setup_orthology_cached(member, ogs, cached_events, event_ids, get_taxid):
    """Ortholog resolution using pre-fetched and pregrouped events.

    cached_events format: {event_id: (level, by_sp1, by_sp2)}
    where by_spX is {taxid_str: set(protein_ids)} - already grouped by species.
    """
    orthology = {}
    member_as_set = set([member])
    ogs_tax_ids = set([og[1] for og in ogs])
    _orthology_stats['calls'] += 1

    for eid in event_ids:
        ev = cached_events.get(eid)
        if ev is None:
            continue
        level, by_sp1, by_sp2 = ev
        if level not in ogs_tax_ids:
            continue

        _orthology_stats['events'] += 1
        __set_coorthologs(by_sp1, by_sp2, member_as_set, orthology)
        __set_coorthologs(by_sp2, by_sp1, member_as_set, orthology)

    return orthology


def get_orthology_stats():
    """Return and reset orthology processing stats."""
    stats = dict(_orthology_stats)
    _orthology_stats['calls'] = _orthology_stats['events'] = _orthology_stats['members'] = 0
    return stats


def collect_orthologs_simple(member_id, best_ogs, all_nogs, cached_events, event_ids):
    """Simplified ortholog collection: flat set of IDs, no species grouping.

    This is much faster than resolve_orthologs_cached() because it skips:
    - Species grouping (the main bottleneck - millions of dict operations)
    - Orthology type classification (one2one, one2many, etc.)

    Returns:
        (orthologs, best_OG) where orthologs is a flat set of protein IDs
        and best_OG is None or a fallback tuple if no events matched.

    The tradeoff: no orthology type info, just "all" orthologs.
    """
    best_ogs_pos = min([all_nogs.index(best_og) for best_og in best_ogs])
    annot_ogs = all_nogs[best_ogs_pos:]

    orthologs = _collect_from_events(member_id, annot_ogs, cached_events, event_ids)

    if orthologs:
        return orthologs, None

    # Try bottom-top if no orthologs found at best level
    if best_ogs_pos > 0:
        for nog in reversed(all_nogs[:best_ogs_pos]):
            orthologs = _collect_from_events(member_id, [nog], cached_events, event_ids)
            if orthologs:
                return orthologs, nog

    # Fallback: just the seed ortholog
    return {member_id}, ("seed_ortholog", "-", f"seed_ortholog@{member_id}|-", None)


def _collect_from_events(member_id, ogs, cached_events, event_ids):
    """Collect all ortholog IDs from events matching the target OG levels.

    No species grouping - just union all protein IDs from both sides
    of matching events where the member appears.
    """
    ogs_tax_ids = set(og[1] for og in ogs)
    orthologs = set()

    for eid in event_ids:
        ev = cached_events.get(eid)
        if ev is None:
            continue
        level, side1_ids, side2_ids = ev
        if level not in ogs_tax_ids:
            continue

        # Check if member is on either side and collect the "other" side
        if member_id in side1_ids:
            orthologs.update(side1_ids)
            orthologs.update(side2_ids)
        elif member_id in side2_ids:
            orthologs.update(side1_ids)
            orthologs.update(side2_ids)

    return orthologs

## END
