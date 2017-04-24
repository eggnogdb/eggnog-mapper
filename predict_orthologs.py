import argparse
import multiprocessing
import os
from collections import defaultdict
import json

from eggnogmapper.common import *
from eggnogmapper import search
from eggnogmapper import annota
from emapper import setup_hmm_search, dump_hmm_matches, get_seq_hmm_matches, refine_matches
from ete3 import NCBITaxa

# Example run: python predict_orthologs.py -i test/testCOG0515.fa -d bact

# working options:  --output_format --orthology_type
# not working: --target_taxa

# TO DO:
# - debug error for target_taxa option:   
# File "/g/bork1/huerta/miniconda/lib/python2.7/multiprocessing/pool.py", line 668, in next raise value
# NameError: global name 'query_taxa' is not defined
# - test if output_formats are ok (especially the separation of the target_taxa)

def main(args):
    fasta_file = args.input
    fields = fasta_file.split("/")
    file_id = fields[len(fields)-1].split(".fa")[0]
    orthologs_file = "%s.orthologs" %file_id
    hmm_hits_file = "tmp.emapper.hmm_hits"
    seed_orthologs_file = "tmp.emapper.seed_orthologs"
    if args.target_taxa:
        target_taxa = args.target_taxa

    # Sequence search with hmmer
    host, port, dbpath, scantype, idmap = setup_hmm_search(args)
    # Start HMM SCANNING sequences (if requested)
    if not pexists(hmm_hits_file) or args.override:
        dump_hmm_matches(args.input, hmm_hits_file, dbpath, port, scantype, idmap, args)

        if not args.no_refine and (not pexists(seed_orthologs_file) or args.override):
            if args.db == 'viruses':
                print 'Skipping seed ortholog detection in "viruses" database'
            elif args.db in EGGNOG_DATABASES:
                refine_matches(args.input, seed_orthologs_file, hmm_hits_file, args)
            else:
                print 'refined hits not available for custom hmm databases.'

    # Orthologs search
    annota.connect()
    find_orthologs(seed_orthologs_file, orthologs_file, hmm_hits_file, args)

    #os.system("rm %s %s" % (hmm_hits_file, seed_orthologs_file))
    
    print "done"

def find_orthologs(seed_orthologs_file, orthologs_file, hmm_hits_file, args):

    OUT = open(orthologs_file, "w")

    if args.output_format in ("per_query", "per_orthology_type"):
        if args.output_format == "per_query":
            ortholog_header = ("#QUERY", "TARGET_TAXON", "TAXID", "ORTHOLOGS") 
        else:
            ortholog_header = ("#QUERY", "ORTHO_TYPE","IN-PARALOGS","TARGET_TAXON", "TAXID", "ORTHOLOGS")

        print >> OUT, "\t".join(ortholog_header)
        
    if pexists(hmm_hits_file):
        seq2bestOG = get_seq_hmm_matches(hmm_hits_file)
         
    if args.target_taxa != "all":
        args.target_taxa = normalize_target_taxa(args.target_taxa)
    
    pool = multiprocessing.Pool(args.cpu)
    json_dict = {}
    for query_result in pool.imap(find_orthologs_per_hit, iter_hit_lines(seed_orthologs_file, args)):
        if query_result:
            for result_line in query_result:
                if args.output_format == "json":
                    json_dict = build_json_format(result_line, json_dict)
                else:
                    write_orthologs_in_file(result_line, OUT, args)

    pool.terminate()

    if args.output_format == "json":
        json.dump(json_dict, OUT)


def find_orthologs_per_hit(arguments):
    """
    Returns a list of results, one result per query_taxon or 
    one unique result if query_taxa is "all"
    Each result is in the form:
    (query_name, all_orthologs, target_taxon_name, target_taxon_taxid, best_hit_name)
    """
    annota.connect()
    line, args = arguments

    if not line.strip() or line.startswith('#'):
        return None
    r = map(str.strip, line.split('\t'))

    query_name = r[0]
    best_hit_name = r[1]
    if best_hit_name == '-' or best_hit_name == 'ERROR':
        return None

    best_hit_evalue = float(r[2])
    best_hit_score = float(r[3])
    
    if best_hit_score < args.seed_ortholog_score or best_hit_evalue > args.seed_ortholog_evalue:
        return None

    target_taxa = args.target_taxa
    if target_taxa == "all":
        all_orthologs = annota.get_member_orthologs(best_hit_name)
        return [(query_name, all_orthologs, "All", "All", best_hit_name)]
    else:
        results = []
        for taxid, name in target_taxa:
            species = target_taxa[(taxid,name)]
            all_orthologs = annota.get_member_orthologs(best_hit_name, target_taxa=species)
            results.append(query_name, all_orthologs, name, taxid, best_hit_name)
        return results
    
def write_orthologs_in_file(result_line, OUT, args):
    """
    Writes orthologs in file for the output formats "per_query" and "per_orthology_type"
    """
    query_name, all_orthologs, target_taxon_name, target_taxid, best_hit_name = result_line

    if args.output_format == "per_query":
        orthologs = sorted(all_orthologs[args.orthology_type])
        print >> OUT, '\t'.join(map(str, (query_name, target_taxon_name, target_taxid, ','.join(orthologs))))

    elif args.output_format == "per_orthology_type" :
        for ortho_type in all_orthologs:
            if len(all_orthologs[ortho_type])==0 or ortho_type=="all":
                continue
            in_paralogs = []
            orthologs = []
            if ortho_type == "one2one" or ortho_type == "one2many":
                in_paralogs.append("N/A")
                orthologs = all_orthologs[ortho_type]
            else:
                for protein in all_orthologs[ortho_type]:
                    if protein.split(".")[0] == best_hit_name.split(".")[0] and protein != best_hit_name:
                        in_paralogs.append(protein)
                    elif protein != best_hit_name:
                        orthologs.append(protein)
            print >> OUT, '\t'.join(map(str, (query_name, ortho_type, ",".join(in_paralogs), target_taxon_name, target_taxid, ','.join(all_orthologs[ortho_type]))))
                                
    OUT.flush()

def build_json_format(result_line, json_dict):
   
    query_name, all_orthologs, target_taxon_name, target_taxid, best_hit_name = result_line
    json_dict[query_name] = {}
    json_dict[query_name][target_taxid] = {}
    for ortho_type in all_orthologs:
        if ortho_type == "all":
            continue
        if ortho_type == "one2one" or ortho_type == "one2many":
            json_dict[query_name][target_taxid][ortho_type] = list(all_orthologs[ortho_type])
        else:
            in_paralogs = []
            orthologs = []
            for protein in all_orthologs[ortho_type]:
                if protein.split(".")[0] == best_hit_name.split(".")[0] and protein != best_hit_name:
                    in_paralogs.append(protein)
                elif protein != best_hit_name:
                    orthologs.append(protein)

            json_dict[query_name][target_taxid][ortho_type] = [in_paralogs, orthologs]
             
    return json_dict

    
def iter_hit_lines(filename, args):
    for line in open(filename):
        yield (line, args)    
                             
# def normalize_target_taxa(target_taxa):
#     """
#     Receives a list of taxa IDs and/or taxa names and returns this structure:
#     {
#        (taxon_taxid1, name): [species_taxid1, species_taxid2, ...]],
#        (taxon_taxid2, name): [species_taxid1, species_taxid2, ...]],
#        ...
#     }
#     """
#     final_taxa = defaultdict(list)
#     TAXIDS_F = open("eggnog4.core_species_info.txt", "r")
#     for line in TAXIDS_F:
#         fields = line.split("\t")
#         taxid = fields[0]
#         name = fields[1]
#         name_lineage = fields[3].split(",")
#         taxid_lineage = fields[4].split(",")
#         if taxid in target_taxa or name in target_taxa:
#             final_taxa[(taxid, name)].append(taxid)
#             continue
#         for target_taxon in target_taxa:
#             if target_taxon in taxid_lineage or target_taxon in name_lineage:
#                 idx = 0
#                 if target_taxon in taxid_lineage:
#                     idx = taxid_lineage.index(target_taxon)
#                 else:
#                     idx = name_lineage.index(target_taxon)
#                 final_taxa[(taxid_lineage[idx], name_lineage[idx])].append(taxid_lineage.pop()) 
#     print final_taxa
#     return final_taxa
    

def normalize_target_taxa(target_taxa):

    """
    Receives a list of taxa IDs and/or taxa names and returns this structure:
    {
       (taxon_taxid1, name): [species_taxid1, species_taxid2, ...]],
       (taxon_taxid2, name): [species_taxid1, species_taxid2, ...]],
       ...
    }

    """
    final_taxa = defaultdict(list)
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()
    for taxon in target_taxa:
        taxid = ""
        try:
            taxid = int(taxon)
            taxon = ncbi.get_taxid_translator([taxid])[taxid]
        except ValueError:
            taxid = ncbi.get_name_translator([taxon])[taxon]        
        
        species = ncbi.get_descendant_taxa(taxid, collapse_subspecies=True)
        for sp in species: 
            final_taxa[(taxid, taxon)].append(sp)

    return final_taxa

    
def parse_args(parser):
    args = parser.parse_args()

    if not args.input:
        parser.error('An input fasta file is required (-i)')

    return args

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # my options
    parser.add_argument('--orthology_type', choices=["one2one", "many2one",
                                                         "one2many","many2many", "all"],
                          default="all",
                          help='defines what type of orthologs should be found')

    parser.add_argument('--target_taxa', type=str,
                          default="all", nargs="+",
                          help='taxa that will be searched for orthologs')
    parser.add_argument('--output_format', choices=["per_query","per_orthology_type", "json"],
                         default="per_query", help="Choose the output format among: per_query, per_orthology_type")
    # server
    pg_db = parser.add_argument_group('Target HMM Database Options')

    pg_db.add_argument('--guessdb', type=int, metavar='',
                       help='guess eggnog db based on the provided taxid')

    pg_db.add_argument('--database', '-d', dest='db', metavar='',
                       help=('specify the target database for sequence searches'
                             '. Choose among: euk,bact,arch, host:port, or a local hmmpressed database'))

    pg_db.add_argument('--dbtype', dest="dbtype",
                    choices=["hmmdb", "seqdb"], default="hmmdb")

    pg_db.add_argument('--qtype',  choices=["hmm", "seq"], default="seq")



    pg_hmm = parser.add_argument_group('HMM search_options')

    pg_hmm.add_argument('--hmm_maxhits', dest='maxhits', type=int, default=1, metavar='',
                    help="Max number of hits to report. Default=1")

    pg_hmm.add_argument('--hmm_evalue', dest='evalue', default=0.001, type=float, metavar='',
                    help="E-value threshold. Default=0.001")

    pg_hmm.add_argument('--hmm_score', dest='score', default=20, type=float, metavar='',
                    help="Bit score threshold. Default=20")

    pg_hmm.add_argument('--hmm_maxseqlen', dest='maxseqlen', type=int, default=5000, metavar='',
                    help="Ignore query sequences larger than `maxseqlen`. Default=5000")

    pg_hmm.add_argument('--hmm_qcov', dest='qcov', type=float, metavar='',
                    help="min query coverage (from 0 to 1). Default=(disabled)")

    pg_hmm.add_argument('--Z', dest='Z', type=float, default=40000000, metavar='',
                    help='Fixed database size used in phmmer/hmmscan'
                        ' (allows comparing e-values among databases). Default=40,000,000')

    pg_out = parser.add_argument_group('Output options')

#    pg_out.add_argument('--output', '-o', type=str, metavar='',
#                    help="base name for output files")

    pg_out.add_argument('--resume', action="store_true",
                    help="Resumes a previous execution skipping reported hits in the output file.")

    pg_out.add_argument('--override', action="store_true",
                    help="Overwrites output files if they exist.")

    pg_out.add_argument("--no_refine", action="store_true",
                    help="Skip hit refinement, reporting only HMM hits.")

    pg_out.add_argument("--no_annot", action="store_true",
                    help="Skip functional annotation, reporting only hits")

    pg_out.add_argument("--no_search", action="store_true",
                    help="Skip HMM search mapping. Use existing hits file")

    pg_out.add_argument("--report_orthologs", action="store_true",
                    help="The list of orthologs used for functional transferred are dumped into a separate file")

    pg_out.add_argument("--scratch_dir", metavar='', type=existing_dir,
                    help='Write output files in a temporary scratch dir, move them to final the final'
                        ' output dir when finished. Speed up large computations using network file'
                        ' systems.')

    pg_out.add_argument("--output_dir", default=os.getcwd(), type=existing_dir, metavar='',
                    help="Where output files should be written")

    pg_out.add_argument("--temp_dir", default=os.getcwd(), type=existing_dir, metavar='',
                    help="Where temporary files are created. Better if this is a local disk.")

    pg_out.add_argument('--no_file_comments', action="store_true",
                        help="No header lines nor stats are included in the output files")

    pg_out.add_argument('--keep_mapping_files', action='store_true',
                        help='Do not delete temporary mapping files used for annotation (i.e. HMMER and'
                        ' DIAMOND search outputs)')
    
    g4 = parser.add_argument_group('Execution options')
    g4.add_argument('-m', dest='mode', choices = ['hmmer', 'diamond'], default='hmmer',
                    help='Default:hmmer')

    g4.add_argument('-i', dest="input", metavar='', type=existing_file,
                    help='Computes annotations for the provided FASTA file')

    g4.add_argument('--translate', action="store_true",
                    help='Assume sequences are genes instead of proteins')

    g4.add_argument("--servermode", action="store_true",
                    help='Loads target database in memory and keeps running in server mode,'
                    ' so another instance of eggnog-mapper can connect to this sever.'
                    ' Auto turns on the --usemem flag')

    g4.add_argument('--usemem', action="store_true",
                    help="""If a local hmmpressed database is provided as target using --db,
                    this flag will allocate the whole database in memory using hmmpgmd.
                    Database will be unloaded after execution.""")

    g4.add_argument('--cpu', type=int, default=2, metavar='')

    g4.add_argument('--annotate_hits_table', type=str, metavar='',
                    help='Annotatate TSV formatted table of query->hits. 4 fields required:'
                    ' query, hit, evalue, score. Implies --no_search and --no_refine.')

    pg_annot = parser.add_argument_group('Annotation Options')

    pg_annot.add_argument("--tax_scope", type=str, choices=TAXID2LEVEL.values()+["auto"],
                    default='auto', metavar='',
                    help=("Fix the taxonomic scope used for annotation, so only orthologs from a "
                          "particular clade are used for functional transfer. "
                          "By default, this is automatically adjusted for every query sequence."))

    pg_annot.add_argument('--excluded_taxa', type=int, metavar='',
                          help='(for debugging and benchmark purposes)')

    pg_annot.add_argument('--go_evidence', type=str, choices=('experimental', 'non-electronic'),
                          default='non-electronic',
                          help='Defines what type of GO terms should be used for annotation:'
                          'experimental = Use only terms inferred from experimental evidence'
                          'non-electronic = Use only non-electronically curated terms')

    pg_seed = parser.add_argument_group('Seed ortholog search option')

    pg_seed.add_argument('--seed_ortholog_evalue', default=0.001, type=float, metavar='',
                    help='Min E-value expected when searching for seed eggNOG ortholog.'
                         ' Applies to phmmer/diamond searches. Queries not having a significant'
                         ' seed orthologs will not be annotated. Default=0.001')

    pg_seed.add_argument('--seed_ortholog_score', default=60, type=float, metavar='',
                    help='Min bit score expected when searching for seed eggNOG ortholog.'
                         ' Applies to phmmer/diamond searches. Queries not having a significant'
                         ' seed orthologs will not be annotated. Default=60')

    
    args = parse_args(parser)

    try:
        main(args)
    except:
        raise
    else:
        sys.exit(0)
