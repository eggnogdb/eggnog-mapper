##
## CPCantalapiedra 2020

import gzip, os, sys

from ...common import get_pfam_clans_file

CLEAN_OVERLAPS_ALL = "all"
CLEAN_OVERLAPS_CLANS = "clans"
CLEAN_OVERLAPS_NONE = "none"
CLEAN_OVERLAPS_HMMSEARCH_ALL = "hmmsearch_all"
CLEAN_OVERLAPS_HMMSEARCH_CLANS = "hmmsearch_clans"

##
def process_overlaps(hits, clean_overlaps, idmap_idx = None):
    if clean_overlaps == CLEAN_OVERLAPS_ALL:
        hits = process_overlaps_all(hits)
    elif clean_overlaps == CLEAN_OVERLAPS_CLANS:
        hits = process_overlaps_clans(hits, idmap_idx)
    elif clean_overlaps == CLEAN_OVERLAPS_HMMSEARCH_ALL:
        hits = process_overlaps_all_queries(hits)
    elif clean_overlaps == CLEAN_OVERLAPS_HMMSEARCH_CLANS:
        hits = process_overlaps_clans_queries(hits)
    else:
        sys.stderr.write(f"Warning: unrecognized clean_overlaps option ({clean_overlaps}). No overlaps will be processed.\n")

    return hits


##
def process_overlaps_clans(hits, idmap_idx = None):
    CLANS_FILE = get_pfam_clans_file()

    if not os.path.exists(CLANS_FILE) or not os.path.isfile(CLANS_FILE):
        raise Exception(f"Couldn't find PFAM clans file at path {CLANS_FILE}, or it is not a file.")

    # sys.stderr.write("Loading clans data...\n")
    
    clans_dict = {}
    with gzip.open(CLANS_FILE, 'rt') as clans_f:
        for line in clans_f:
            data = line.strip().split("\t")
            pfname = data[3]
            clan = data[1]
            if clan is not None and clan != "":
                clans_dict[pfname] = clan

    # sys.stderr.write("Checking clans overlaps...\n")
    clean_doms = []
    total_range = set()

    for hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore in hits:
        hmmfrom, hmmto, sqfrom, sqto = map(int, [hmmfrom, hmmto, sqfrom, sqto])
        new_span = set(range(sqfrom, sqto+1))

        total_overlap = new_span & total_range
        if len(total_overlap) > 0:
            best = True
            tmp_clean_doms = []
            tmp_overlapping = []

            hitname = hid
            if idmap_idx:
                hitname = idmap_idx[hid][0]
                
            hitclan = None

            for phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore in clean_doms:
                prev_span = set(range(psqfrom, psqto+1))
                overlap = new_span & prev_span

                phitname = phid
                if idmap_idx:
                    phitname = idmap_idx[phid][0]
                    
                if hitclan is None:
                    hitclan = clans_dict.get(hitname)
                phitclan = clans_dict.get(phitname)

                # print(f"hmmer_overlaps.py: {hitname}-{hitclan} / {phitname}-{phitclan}")
                
                if len(overlap) > 0 and best == True and hitclan is not None and hitclan == phitclan:
                    if heval > pheval:
                        best = False
                    tmp_overlapping.append([phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore])
                else:
                    tmp_clean_doms.append([phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore])

            if best == True:
                tmp_clean_doms.append([hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore])
            else:
                tmp_clean_doms.extend(tmp_overlapping)

            # update clean_doms and total_range
            clean_doms = tmp_clean_doms
            for phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore in clean_doms:
                clean_span = set(range(psqfrom, psqto+1))
                total_range.update(clean_span)
        else:
            clean_doms.append([hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore])
            total_range.update(new_span)

    return clean_doms


##
def process_overlaps_all(hits):
    CLANS_FILE = get_pfam_clans_file()
    
    clean_doms = []
    total_range = set()

    for hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore in hits:
        hmmfrom, hmmto, sqfrom, sqto = map(int, [hmmfrom, hmmto, sqfrom, sqto])
        new_span = set(range(sqfrom, sqto+1))

        total_overlap = new_span & total_range
        if len(total_overlap) > 0:
            best = True
            tmp_clean_doms = []
            tmp_overlapping = []

            for phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore in clean_doms:
                prev_span = set(range(psqfrom, psqto+1))
                overlap = new_span & prev_span
                if len(overlap) > 0 and best == True:
                    if heval > pheval:
                        best = False
                    tmp_overlapping.append([phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore])
                else:
                    tmp_clean_doms.append([phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore])

            if best == True:
                tmp_clean_doms.append([hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore])
            else:
                tmp_clean_doms.extend(tmp_overlapping)

            # update clean_doms and total_range
            clean_doms = tmp_clean_doms
            for phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore in clean_doms:
                clean_span = set(range(psqfrom, psqto+1))
                total_range.update(clean_span)
        else:
            clean_doms.append([hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore])
            total_range.update(new_span)

    return clean_doms


##
def process_overlaps_all_queries(namedhits):
    CLANS_FILE = get_pfam_clans_file()
    
    targets_hits = {}

    for name, querylen, hits in namedhits:
    
        for hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore in hits:
            hmmfrom, hmmto, sqfrom, sqto = map(int, [hmmfrom, hmmto, sqfrom, sqto])
            new_span = set(range(sqfrom, sqto+1))
            
            if hid in targets_hits:
                total_range, clean_doms = targets_hits[hid]
                
                total_overlap = new_span & total_range
                if len(total_overlap) > 0:
                    best = True
                    tmp_clean_doms = []
                    tmp_overlapping = []

                    for pname, pquerylen, phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore in clean_doms:
                        prev_span = set(range(psqfrom, psqto+1))
                        overlap = new_span & prev_span
                        if len(overlap) > 0 and best == True:
                            if heval > pheval:
                                best = False
                            tmp_overlapping.append([pname, pquerylen, phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore])
                        else:
                            tmp_clean_doms.append([pname, pquerylen, phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore])

                    if best == True:
                        tmp_clean_doms.append([name, querylen, hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore])
                    else:
                        tmp_clean_doms.extend(tmp_overlapping)

                    # update clean_doms and total_range
                    for pname, pquerylen, phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore in clean_doms:
                        clean_span = set(range(psqfrom, psqto+1))
                        total_range.update(clean_span)
                    targets_hits[hid] = (total_range, tmp_clean_doms)
                else:
                    clean_doms.append([name, querylen, hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore])
                    total_range.update(new_span)
            else:
                clean_doms = [[name, querylen, hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore]]
                total_range = set()
                total_range.update(new_span)
                targets_hits[hid] = (total_range, clean_doms)

    clean_doms = []
    for hid, (total_range, hclean_doms) in targets_hits.items():
        for name, querylen, hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore in sorted(hclean_doms, key=lambda x: x[7]):
            clean_doms.append((name, querylen, [[hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore]]))

    return clean_doms


##
def process_overlaps_clans_queries(namedhits):
    CLANS_FILE = get_pfam_clans_file()
    
    if not os.path.exists(CLANS_FILE) or not os.path.isfile(CLANS_FILE):
        raise Exception(f"Couldn't find PFAM clans file at path {CLANS_FILE}, or it is not a file.")

    # sys.stderr.write("Loading clans data...\n")
    
    clans_dict = {}
    with gzip.open(CLANS_FILE, 'rt') as clans_f:
        for line in clans_f:
            data = line.strip().split("\t")
            pfname = data[3]
            clan = data[1]
            if clan is not None and clan != "":
                clans_dict[pfname] = clan
                
    targets_hits = {}

    for name, querylen, hits in namedhits:

        hitclan = None
        
        for hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore in hits:
            hmmfrom, hmmto, sqfrom, sqto = map(int, [hmmfrom, hmmto, sqfrom, sqto])
            new_span = set(range(sqfrom, sqto+1))
            
            if hid in targets_hits:
                total_range, clean_doms = targets_hits[hid]
                
                total_overlap = new_span & total_range
                if len(total_overlap) > 0:
                    best = True
                    tmp_clean_doms = []
                    tmp_overlapping = []

                    for pname, pquerylen, phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore in clean_doms:
                        prev_span = set(range(psqfrom, psqto+1))
                        overlap = new_span & prev_span

                        if hitclan is None:
                            hitclan = clans_dict.get(name)
                        phitclan = clans_dict.get(pname)

                        # print(name+"\t"+str(hitclan))
                        # print(pname+"\t"+str(phitclan))
                        # print("//")
                    
                        if len(overlap) > 0 and best == True and hitclan is not None and hitclan == phitclan:
                            if heval > pheval:
                                best = False
                            tmp_overlapping.append([pname, pquerylen, phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore])
                        else:
                            tmp_clean_doms.append([pname, pquerylen, phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore])

                    if best == True:
                        tmp_clean_doms.append([name, querylen, hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore])
                    else:
                        tmp_clean_doms.extend(tmp_overlapping)

                    # update clean_doms and total_range
                    for pname, pquerylen, phid, pheval, phscore, phmmfrom, phmmto, psqfrom, psqto, pdomscore in clean_doms:
                        clean_span = set(range(psqfrom, psqto+1))
                        total_range.update(clean_span)
                    targets_hits[hid] = (total_range, tmp_clean_doms)
                else:
                    clean_doms.append([name, querylen, hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore])
                    total_range.update(new_span)
            else:
                clean_doms = [[name, querylen, hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore]]
                total_range = set()
                total_range.update(new_span)
                targets_hits[hid] = (total_range, clean_doms)

    clean_doms = []
    for hid, (total_range, hclean_doms) in targets_hits.items():
        for name, querylen, hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore in sorted(hclean_doms, key=lambda x: x[7]):
            clean_doms.append((name, querylen, [[hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore]]))

    return clean_doms

## END
