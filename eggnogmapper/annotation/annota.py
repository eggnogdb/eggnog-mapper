#
## JHCepas

from collections import Counter, defaultdict

# from . import db_sqlite


def summarize_annotations(seq_names, annotations_fields, target_go_ev, excluded_go_ev, eggnog_db):

    MIN_COUNT_PFAMS = 1
    MIN_RATE_PFAMS = 0.05
    annotations = defaultdict(Counter)
    
    # for fields in db_sqlite.get_annotations(','.join(['"%s"' % n for n in seq_names])):
    for fields in eggnog_db.get_annotations(','.join(['"%s"' % n for n in seq_names])):
        for i, h in enumerate(annotations_fields):
            if not fields[i]:
                continue
            if h == 'GOs':
                gos = fields[i]
                annotations[h].update(parse_gos(gos, target_go_ev, excluded_go_ev))
            elif h == 'Preferred_name':
                annotations[h].update([fields[i].strip()])
            else:
                values = [str(x).strip() for x in fields[i].split(',')]
                annotations[h].update(values)
                
    for h in annotations:
        del annotations[h]['']

    if annotations:
        try:
            pname = annotations['Preferred_name'].most_common(1)
            if pname: 
                name_candidate, freq = annotations['Preferred_name'].most_common(1)[0]
            else:
                freq =  0
        except:
            print(annotations)
            raise 
        if freq >= 2:
            annotations['Preferred_name'] = [name_candidate]
        else:
            annotations['Preferred_name'] = ['']

    if "PFAMs" in annotations:
        total_c = len(seq_names)
        annotations["PFAMs"] = Counter({k: c for k, c in annotations["PFAMs"].items()
                                        if c > MIN_COUNT_PFAMS and c / total_c > MIN_RATE_PFAMS})

    return annotations


def parse_gos(gos, target_go_ev, excluded_go_ev):
    selected_gos = set()
    for g in gos.strip().split(','):
        if not g:
            continue
        gocat, gid, gevidence = list(map(str, g.strip().split('|')))
        if not target_go_ev or gevidence in target_go_ev:
            if not excluded_go_ev or gevidence not in excluded_go_ev:
                selected_gos.add(gid)
    return selected_gos
    
## END
