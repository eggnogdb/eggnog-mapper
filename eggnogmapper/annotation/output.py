##
## CPCantalapiedra 2020

from ..common import get_call_info

#############
# Orthologs

##
def output_orthologs_header(out, no_file_comments):
    if not no_file_comments:        
        # Call info
        print(get_call_info(), file=out)

    # Header
    header = ["#query", "orth_type", "species", "orthologs"]
    print('\t'.join(header), file=out)

    return

##
def output_orthologs_closure(out, ncbi):
    def output_orthologs_rows(rows):
        # Rows
        for query_name in rows:
            query_orthologs = rows[query_name]

            if "annot_orthologs" in query_orthologs:
                annot_orthologs = query_orthologs["annot_orthologs"]
            else:
                annot_orthologs = set()

            for target in query_orthologs:
                if target == "all": continue
                if target == "annot_orthologs": continue
                query_target_orths = query_orthologs[target]
                if query_target_orths is None or len(query_target_orths) == 0:
                    continue

                orthologs_taxids = set([int(x.split(".")[0]) for x in query_target_orths])
                orthologs_taxnames = sorted(ncbi.get_taxid_translator(orthologs_taxids).items(), key=lambda x: x[1])

                for taxid, taxname in orthologs_taxnames:
                    orth_names = []
                    for orth in [x for x in query_target_orths if int(x.split(".")[0]) == taxid]:
                        orth_name = orth.split(".")[1]
                        if orth in annot_orthologs:
                            orth_name = f"*{orth_name}"
                        orth_names.append(orth_name)

                    row = [query_name, target, f"{taxname}({taxid})", ",".join(sorted(orth_names))]
                    print('\t'.join(row), file=out)

        return
    return output_orthologs_rows

##
def output_orthologs_footer(out, no_file_comments, qn, elapsed_time):
    # Timings
    if not no_file_comments:
        print('# %d queries scanned' % (qn), file=out)
        print('# Total time (seconds):', elapsed_time, file=out)
        print('# Rate:', "%0.2f q/s" % ((float(qn) / elapsed_time)), file=out)

    return

##############
# Annotations

HIT_HEADER = ["#query_name",
              "seed_eggNOG_ortholog",
              "seed_ortholog_evalue",
              "seed_ortholog_score",
              "eggNOG OGs",
              "narr_og_name",
              "narr_og_cat",
              "narr_og_desc"]

BEST_OG_HEADER = ["best_og_name",
                  "best_og_cat",
                  "best_og_desc"]

ANNOTATIONS_HEADER = ['Preferred_name',
                      'GOs',
                      'EC',
                      'KEGG_ko',
                      'KEGG_Pathway',
                      'KEGG_Module',
                      'KEGG_Reaction',
                      'KEGG_rclass',
                      'BRITE',
                      'KEGG_TC',
                      'CAZy',
                      'BiGG_Reaction',
                      'PFAMs']

##
def output_annotations_header(out, no_file_comments, md5_field):

    if not no_file_comments:
        print(get_call_info(), file=out)

    print('\t'.join(HIT_HEADER), end="\t", file=out)

    print('\t'.join(BEST_OG_HEADER), end="\t", file=out)

    annot_header = ANNOTATIONS_HEADER
    if md5_field == True:
        annot_header.append("md5")

    print('\t'.join(annot_header), file=out)

    return

##
def output_annotations_closure(out, md5_field, md5_queries):
    def output_annotations_rows(rows):
        for annot_columns in rows:
            if md5_field == True:
                query_name = annot_columns[0]
                if query_name in md5_queries:
                    annot_columns.append(md5_queries[query_name])
                else:
                    annot_columns.append("-")
            print('\t'.join([x if x is not None and x.strip() != "" else "-" for x in annot_columns]), file=out)

        return

    return output_annotations_rows

##
def output_annotations_footer(out, no_file_comments, qn, elapsed_time):
    if not no_file_comments:
        print('# %d queries scanned' % (qn), file=out)
        print('# Total time (seconds):', elapsed_time, file=out)
        print('# Rate:', "%0.2f q/s" % ((float(qn) / elapsed_time)), file=out)

    return

## END
