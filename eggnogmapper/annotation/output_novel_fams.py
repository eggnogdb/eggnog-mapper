##
## CPCantalapiedra 2020

import time

from ..common import get_call_info

#############
# Orthologs
## not available for novel fams

##############
# Annotations

##

HIT_HEADER = ["query",
              "target",
              "evalue",
              "score",
              "novel_fam"]

ANNOTATIONS_HEADER = []

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

    (query_name, best_hit_name, best_hit_evalue, best_hit_score, novel_fam) = annotation

    annot_columns = [query_name, best_hit_name, str(best_hit_evalue), str(best_hit_score), novel_fam]
                    
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
# no excel output for novel fams yet

## END
