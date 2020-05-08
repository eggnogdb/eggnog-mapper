##
## JHCepas
## CPCantalapiedra 2020

from . import orthology
import multiprocessing

def iter_hit_lines(filename, args):
    for line in open(filename):
        if line.startswith('#') or not line.strip():
            continue
        yield (line, args)


def dump_orthologs(seed_orthologs_file, orthologs_file, args):
#Copy from predict_orthologs.py
    OUT = open(orthologs_file, "w")

    if args.predict_output_format == "per_query":
        ortholog_header = ("#Query", "Orthologs")
    elif args.predict_output_format == "per_species":
        ortholog_header = ("#Query", "Species", "Orthologs")

    print("\t".join(ortholog_header), file=OUT)

    if args.target_taxa != 'all':
        args._expanded_target_taxa = orthology.normalize_target_taxa(args.target_taxa)
    else:
        # report orthologs from any species by default
        args._expanded_target_taxa = None

    pool = multiprocessing.Pool(args.cpu)
    for result in pool.imap(find_orthologs_per_hit, iter_hit_lines(seed_orthologs_file, args)):
        if result:
            write_orthologs_in_file(result, OUT, args)

    pool.terminate()

def find_orthologs_per_hit(arguments):
    #Copy from predict_orthologs.py

    orthology.connect()
    line, args = arguments

    if not line.strip() or line.startswith('#'):
        return None
    r = list(map(str.strip, line.split('\t')))

    query_name = r[0]
    best_hit_name = r[1]
    if best_hit_name == '-' or best_hit_name == 'ERROR':
        return None

    best_hit_evalue = float(r[2])
    best_hit_score = float(r[3])

    if best_hit_score < args.seed_ortholog_score or best_hit_evalue > args.seed_ortholog_evalue:
        return None

    target_taxa = args._expanded_target_taxa
    
    orthologs_pred = orthology.predict_orthologs_by_seed(best_hit_name, target_taxa=target_taxa, target_levels = None)
    return (query_name, best_hit_name, orthologs_pred)


def write_orthologs_in_file(result_line, ORTHOLOGS, args):
#Copy from predict_orthologs.py

    """
    Writes orthologs in file for all output formats except json
    """
    query_name, best_hit_name, orthologs_pred = result_line
    
    if args._expanded_target_taxa:  
        target_taxa = list(args._expanded_target_taxa)
    else:
        target_taxa = None
    
    if args.predict_output_format == "per_query":
        orthologs = []
        for key in orthologs_pred:
            if target_taxa is not None:
                    if key in target_taxa:
                        members = (','.join(orthologs_pred[key]))
                        orthologs.append(members)
                        print('\t'.join(map(str, (query_name, ','.join(orthologs)))), file=ORTHOLOGS)

            else:
                members = (','.join(orthologs_pred[key]))
                orthologs.append(members)
                print('\t'.join(map(str, (query_name, ','.join(orthologs)))), file=ORTHOLOGS)

    elif args.predict_output_format == "per_species":
        for key in orthologs_pred:
            sp_taxid = int(key)
            if target_taxa is not None: 
                if sp_taxid in target_taxa:
                    print('\t'.join(map(str, (query_name, key,
                                                    ','.join(orthologs_pred[key])))), file=ORTHOLOGS) 

            else:
                print('\t'.join(map(str, (query_name, key,
                                                    ','.join(orthologs_pred[key])))), file=ORTHOLOGS)
    ORTHOLOGS.flush()
