##
## CPCantalapiedra 2020

import time

from ..common import get_call_info

from .ncbitaxa.ncbiquery import get_ncbi

#############
# Orthologs

##
def output_orthologs(annots, orthologs_file, resume, no_file_comments):
    start_time = time.time()

    ncbi = get_ncbi(usemem = True)

    if resume == True:
        file_mode = 'a'
    else:
        file_mode = 'w'
        
    with open(orthologs_file, file_mode) as ORTHOLOGS_OUT:
        output_orthologs_header(ORTHOLOGS_OUT, no_file_comments, not resume)

        qn = 0
        for ((hit, annotation), exists) in annots:

            # exists == False (--resume)
            
            if exists == False and annotation is not None:
                output_orthologs_row(ORTHOLOGS_OUT, annotation, ncbi)
                
            yield (hit, annotation), exists
            qn += 1

        elapsed_time = time.time() - start_time
        output_orthologs_footer(ORTHOLOGS_OUT, no_file_comments, qn, elapsed_time)

    if ncbi is not None: ncbi.close()
    
    return

##
def output_orthologs_row(out, annotation, ncbi):
    (query_name, best_hit_name, best_hit_evalue, best_hit_score,
     annotations,
     (og_name, og_cat, og_desc),
     max_annot_lvl,
     match_nog_names,
     all_orthologies, annot_orthologs) = annotation

    best_hit_name_id = best_hit_name.split(".")[1]

    all_orthologies["annot_orthologs"] = annot_orthologs

    seed_shown = False # show seed ortholog only once for each query
    
    for target in all_orthologies:
        if target == "all": continue
        if target == "annot_orthologs": continue
        
        query_target_orths = all_orthologies[target]
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

                # if it is the seed, show it separately, and only once
                if orth_name in {best_hit_name_id, f"*{best_hit_name_id}"}:
                    
                    if seed_shown == False:
                        row = [query_name, "seed", f"{taxname}({taxid})", orth_name]
                        print('\t'.join(row), file=out)
                        seed_shown = True
                        
                    # else: DON'T SHOW AGAIN THE SEED
                    #     pass                    
                else:
                    orth_names.append(orth_name)

            if len(orth_names) > 0:
                row = [query_name, target, f"{taxname}({taxid})", ",".join(sorted(orth_names))]
                print('\t'.join(row), file=out)
                
    return

##
def output_orthologs_header(out, no_file_comments, print_header):
    if not no_file_comments:        
        # Call info
        print(get_call_info(), file=out)

    # Header
    if print_header == True:
        header = ["#query", "orth_type", "species", "orthologs"]
        print('\t'.join(header), file=out)

    return

##
def output_orthologs_footer(out, no_file_comments, qn, elapsed_time):
    # Timings
    if not no_file_comments:
        print('## %d queries scanned' % (qn), file=out)
        print('## Total time (seconds):', elapsed_time, file=out)
        print('## Rate:', "%0.2f q/s" % ((float(qn) / elapsed_time)), file=out)

    return

##############
# Annotations

##

HIT_HEADER = ["query",
              "seed_ortholog",
              "evalue",
              "score",
              "eggNOG_OGs",
              "max_annot_lvl",
              "COG_category",
              "Description"]

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

ANNOTATIONS_WHOLE_HEADER = HIT_HEADER + ANNOTATIONS_HEADER

##
def output_annotations(annots, annot_file, resume, no_file_comments, md5_field, md5_queries):

    if resume == True:
        file_mode = 'a'
    else:
        file_mode = 'w'
        
    start_time = time.time()
    
    with open(annot_file, file_mode) as ANNOTATIONS_OUT:
        output_annotations_header(ANNOTATIONS_OUT, no_file_comments, md5_field, not resume)

        qn = 0
        for (hit, annotation), exists in annots:

            # exists == False (--resume)
            
            if exists == False and annotation is not None:
                output_annotations_row(ANNOTATIONS_OUT, annotation, md5_field, md5_queries)
                
            yield (hit, annotation), exists
            qn += 1
        
        elapsed_time = time.time() - start_time
        output_annotations_footer(ANNOTATIONS_OUT, no_file_comments, qn, elapsed_time)
    return

##
def output_annotations_row(out, annotation, md5_field, md5_queries):

    (query_name, best_hit_name, best_hit_evalue, best_hit_score,
     annotations,
     (og_name, og_cat, og_desc),
     max_annot_lvl,
     match_nog_names,
     all_orthologies, annot_orthologs) = annotation

    annot_columns = [query_name, best_hit_name, str(best_hit_evalue), str(best_hit_score),
                     ",".join(match_nog_names), str(max_annot_lvl),
                     og_cat, og_desc]
    
    for h in ANNOTATIONS_HEADER:
        if h in annotations and annotations[h] is not None:
            annot_columns.append(",".join(sorted(list(annotations[h]))))
        else:
            annot_columns.append('-')
                    
    if md5_field == True:
        query_name = annot_columns[0]
        if query_name in md5_queries:
            annot_columns.append(md5_queries[query_name])
        else:
            annot_columns.append("-")
            
    print('\t'.join([x if x is not None and x.strip() != "" else "-" for x in annot_columns]), file=out)
    
    return

##
def output_annotations_header(out, no_file_comments, md5_field, print_header):

    if not no_file_comments:
        print(get_call_info(), file=out)

    if print_header == True:
        print("#", end="", file=out)
        
        annot_header = ANNOTATIONS_WHOLE_HEADER
        if md5_field == True:
            annot_header.append("md5")

        print('\t'.join(annot_header), file=out)

    return



##
def output_annotations_footer(out, no_file_comments, qn, elapsed_time):
    if not no_file_comments:
        print('## %d queries scanned' % (qn), file=out)
        print('## Total time (seconds):', elapsed_time, file=out)
        print('## Rate:', "%0.2f q/s" % ((float(qn) / elapsed_time)), file=out)

    return


######################
# Annotations as excel

def output_excel(annots, excel_filen, resume, no_file_comments, md5_field, md5_queries):
    import xlsxwriter                                           

    workbook = xlsxwriter.Workbook(excel_filen)                   
    worksheet = workbook.add_worksheet()

    start_time = time.time()

    # header
    row = output_excel_header(worksheet, no_file_comments, md5_field)

    # rows
    qn = 0
    for (hit, annotation), exists in annots:

        if annotation is not None:
            output_excel_row(worksheet, row, annotation, md5_field, md5_queries)
            row+=1

        yield (hit, annotation), exists
        qn += 1
        
        elapsed_time = time.time() - start_time
        output_excel_footer(worksheet, row, no_file_comments, qn, elapsed_time)

    workbook.close()
        
    return

##
def output_excel_header(worksheet, no_file_comments, md5_field):

    # row = 0, col = 0
    row = 0
    if not no_file_comments:
        worksheet.write(row, 0, get_call_info())
        row+=1
        
        # row = 1, col = 0
        worksheet.write(1, 0, "#")
        row+=1

    # Header fields
    annot_header = ANNOTATIONS_WHOLE_HEADER
    if md5_field == True:
        annot_header.append("md5")

    # row = 2
    col = 0
    for field in annot_header:
        worksheet.write(row, col, field)
        col+=1

    row+=1

    return row

##
def output_excel_row(worksheet, row, annotation, md5_field, md5_queries):

    (query_name, best_hit_name, best_hit_evalue, best_hit_score,
     annotations,
     (og_name, og_cat, og_desc),
     max_annot_lvl,
     match_nog_names,
     all_orthologies, annot_orthologs) = annotation

    annot_columns = [query_name, best_hit_name, str(best_hit_evalue), str(best_hit_score),
                     ",".join(match_nog_names), str(max_annot_lvl),
                     og_cat, og_desc]
    
    for h in ANNOTATIONS_HEADER:
        if h in annotations and annotations[h] is not None:
            annot_columns.append(",".join(sorted(list(annotations[h]))))
        else:
            annot_columns.append('-')
                    
    if md5_field == True:
        query_name = annot_columns[0]
        if query_name in md5_queries:
            annot_columns.append(md5_queries[query_name])
        else:
            annot_columns.append("-")

    col = 0
    for x in annot_columns:
        field = x if x is not None and x.strip() != "" else "-"
        worksheet.write(row, col, field)
        col+=1
    
    return

##
def output_excel_footer(worksheet, row, no_file_comments, qn, elapsed_time):
    if not no_file_comments:
        worksheet.write(row, 0, f'## {qn} queries scanned')
        worksheet.write(row+1, 0, f'## Total time (seconds): {elapsed_time}')
        worksheet.write(row+2, 0, f'## Rate: {float(qn) / elapsed_time:.2f} q/s')
    return

## END
