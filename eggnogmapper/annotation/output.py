##
## CPCantalapiedra 2020

import time

from ..common import get_call_info

from .ncbitaxa.ncbiquery import get_ncbi
from .db_sqlite import get_eggnog_db

#############
# Orthologs

##
def output_orthologs(annots, orthologs_file, resume, no_file_comments, eggnog_db=None):
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
                output_orthologs_row(ORTHOLOGS_OUT, annotation, ncbi, eggnog_db)

            yield (hit, annotation), exists
            qn += 1

        elapsed_time = time.time() - start_time
        output_orthologs_footer(ORTHOLOGS_OUT, no_file_comments, qn, elapsed_time)

    if ncbi is not None: ncbi.close()

    return

##
def output_orthologs_row(out, annotation, ncbi, eggnog_db=None):
    # Tuple shapes accepted: 10 (legacy), 11 (v3 cascade), 12 (v3.1
    # adds tax_scope_used). The orthologs row only needs the first
    # 10 fields; later additions are unpacked but unused here.
    n = len(annotation)
    if n == 10:
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs) = annotation
    elif n == 11:
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs, _confidence) = annotation
    else:
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs,
         _confidence, _tax_scope_used) = annotation

    int_mode = eggnog_db is not None and eggnog_db._int_mode

    if int_mode:
        # Decode all protein IDs in all_orthologies + annot_orthologs
        all_ids = set()
        all_ids.add(int(best_hit_name))
        for target, orths in all_orthologies.items():
            if orths:
                all_ids.update(orths)
        if annot_orthologs:
            all_ids.update(annot_orthologs)
        id_to_name = eggnog_db.decode_protein_ids(all_ids)

        best_hit_real_name = id_to_name.get(int(best_hit_name), str(best_hit_name))
        best_hit_name_id = best_hit_real_name.split(".", 1)[1] if "." in best_hit_real_name else best_hit_real_name
    else:
        id_to_name = None
        best_hit_name_id = best_hit_name.split(".")[1]

    all_orthologies["annot_orthologs"] = annot_orthologs

    seed_shown = False # show seed ortholog only once for each query

    for target in all_orthologies:
        if target == "all": continue
        if target == "annot_orthologs": continue

        query_target_orths = all_orthologies[target]
        if query_target_orths is None or len(query_target_orths) == 0:
            continue

        if int_mode:
            # Decode integer IDs to get taxids and names
            orthologs_taxids = set()
            orth_by_taxid = {}
            for o in query_target_orths:
                name = id_to_name.get(o, "0.unknown")
                taxid = int(name.split(".")[0])
                orthologs_taxids.add(taxid)
                orth_by_taxid.setdefault(taxid, []).append((o, name))
        else:
            orthologs_taxids = set([int(x.split(".")[0]) for x in query_target_orths])

        orthologs_taxnames = sorted(ncbi.get_taxid_translator(orthologs_taxids).items(), key=lambda x: x[1])

        for taxid, taxname in orthologs_taxnames:
            orth_names = []

            if int_mode:
                orths_for_taxid = orth_by_taxid.get(taxid, [])
                for orth_id, orth_full_name in orths_for_taxid:
                    orth_name = orth_full_name.split(".", 1)[1] if "." in orth_full_name else orth_full_name
                    if orth_id in (annot_orthologs if annot_orthologs else []):
                        orth_name = f"*{orth_name}"
                    if orth_name in {best_hit_name_id, f"*{best_hit_name_id}"}:
                        if seed_shown == False:
                            row = [query_name, "seed", f"{taxname}({taxid})", orth_name]
                            print('\t'.join(row), file=out)
                            seed_shown = True
                    else:
                        orth_names.append(orth_name)
            else:
                for orth in [x for x in query_target_orths if int(x.split(".")[0]) == taxid]:
                    orth_name = orth.split(".")[1]
                    if orth in annot_orthologs:
                        orth_name = f"*{orth_name}"
                    if orth_name in {best_hit_name_id, f"*{best_hit_name_id}"}:
                        if seed_shown == False:
                            row = [query_name, "seed", f"{taxname}({taxid})", orth_name]
                            print('\t'.join(row), file=out)
                            seed_shown = True
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

# Phase 3c: per-source confidence is appended as a single column. Format
# is `field=tier;field=tier;...` listing only the fields that were
# emitted (and therefore got a confidence label). Reading code can split
# on `;` and `=` without ever quoting tabs.
#
# Phase 7.1b: `tax_scope_used` records the resolved per-seed taxonomic
# scope decision (e.g. "Metazoa", "Bacteria,Archaea", "explicit:33090",
# "none"). Closes the v3 docs gap where the auto-scope per-seed outcome
# was visible only on stderr.
ANNOTATIONS_WHOLE_HEADER = HIT_HEADER + ANNOTATIONS_HEADER + [
    'annotation_confidence', 'tax_scope_used',
]

##
def output_annotations(annots, annot_file, resume, no_file_comments, md5_field, md5_queries,
                       applied_filters=None):

    if resume == True:
        file_mode = 'a'
    else:
        file_mode = 'w'

    try:
        eggnog_db = get_eggnog_db()
    except Exception:
        eggnog_db = None

    start_time = time.time()

    with open(annot_file, file_mode) as ANNOTATIONS_OUT:
        output_annotations_header(
            ANNOTATIONS_OUT, no_file_comments, md5_field, not resume,
            applied_filters=applied_filters,
        )

        qn = 0
        for (hit, annotation), exists in annots:

            # exists == False (--resume)

            if exists == False and annotation is not None:
                output_annotations_row(ANNOTATIONS_OUT, annotation, md5_field, md5_queries,
                                       eggnog_db=eggnog_db)

            yield (hit, annotation), exists
            qn += 1

        elapsed_time = time.time() - start_time
        output_annotations_footer(ANNOTATIONS_OUT, no_file_comments, qn, elapsed_time)
    return

##
def output_annotations_row(out, annotation, md5_field, md5_queries, eggnog_db=None):

    # Annotation tuple shapes:
    #   10 elem  — pre-v3 legacy (no confidence, no tax_scope)
    #   11 elem  — v3 cascade   (annotations_confidence appended)
    #   12 elem  — v3.1 (Phase 7.1b: tax_scope_used appended)
    # The two earlier shapes are accepted as a graceful fallback for older
    # test fixtures and any in-process callers from before Phase 7.1.
    annotations_confidence = {}
    tax_scope_used = "none"
    if len(annotation) == 10:
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs) = annotation
    elif len(annotation) == 11:
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs,
         annotations_confidence) = annotation
    else:
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs,
         annotations_confidence, tax_scope_used) = annotation

    # In int_mode, translate integer seed_ortholog to real protein name
    if eggnog_db is not None and eggnog_db._int_mode:
        seed_display = eggnog_db.get_protein_name(best_hit_name)
    else:
        seed_display = best_hit_name

    annot_columns = [query_name, seed_display, str(best_hit_evalue), str(best_hit_score),
                     ",".join(match_nog_names), str(max_annot_lvl),
                     og_cat, og_desc]

    for h in ANNOTATIONS_HEADER:
        if h in annotations and annotations[h] is not None:
            annot_columns.append(",".join(sorted(list(annotations[h]))))
        else:
            annot_columns.append('-')

    # Per-source confidence column. Only the fields that actually appear
    # in `annotations_confidence` are written; "-" if there's nothing.
    if annotations_confidence:
        conf_str = ";".join(
            f"{field}={tier}"
            for field, tier in sorted(annotations_confidence.items())
        )
    else:
        conf_str = "-"
    annot_columns.append(conf_str)

    # Phase 7.1b: per-seed resolved tax_scope decision.
    annot_columns.append(tax_scope_used or "-")

    if md5_field == True:
        query_name = annot_columns[0]
        if query_name in md5_queries:
            annot_columns.append(md5_queries[query_name])
        else:
            annot_columns.append("-")
            
    print('\t'.join([x if x is not None and x.strip() != "" else "-" for x in annot_columns]), file=out)
    
    return

##
def output_annotations_header(out, no_file_comments, md5_field, print_header,
                              applied_filters=None):

    if not no_file_comments:
        print(get_call_info(), file=out)
        # Phase 7.1a — record the resolved values of every annotation-stage
        # filter and threshold the run actually used. Closes the v3 docs gap
        # where reproducibility relied on parsing the raw command line on
        # line 3 (defaults were never expanded).
        if applied_filters:
            print(format_applied_filters(applied_filters), file=out)

    if print_header == True:
        print("#", end="", file=out)

        annot_header = ANNOTATIONS_WHOLE_HEADER
        if md5_field == True:
            annot_header.append("md5")

        print('\t'.join(annot_header), file=out)

    return


def format_applied_filters(filters: dict) -> str:
    """Render the resolved filter / threshold values as a comment block.

    Output is a multi-line string starting with ``## applied filters:``
    followed by one ``##   key=value`` per filter, sorted by key. ``None``
    becomes the literal string ``null`` so the block is unambiguous when a
    threshold is intentionally unset.
    """
    lines = ["## applied filters:"]
    for key in sorted(filters):
        value = filters[key]
        if value is None:
            rendered = "null"
        elif isinstance(value, (list, tuple, set, frozenset)):
            rendered = ",".join(map(str, sorted(value))) if value else "null"
        else:
            rendered = str(value)
        lines.append(f"##   {key}={rendered}")
    return "\n".join(lines)



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
