##
## CPCantalapiedra 2020

import os
import sys
import time

from ..common import get_call_info

from .ncbitaxa.ncbiquery import get_ncbi
from .db_sqlite import get_eggnog_db

# Translation table: replace tab, CR, LF, and other ASCII control chars
# with a single space so no field value can break TSV structure.
_CTRL_CHARS = "".join(chr(i) for i in range(32) if i not in (9, 10, 13)) + chr(127)
_CTRL_TABLE = str.maketrans(
    '\t\n\r' + _CTRL_CHARS,
    ' ' * (3 + len(_CTRL_CHARS))
)


def _tsv_field(x) -> str:
    """Sanitize a value for safe TSV output.

    Replaces embedded tabs, newlines, carriage returns, and other control
    characters with a space. Returns "-" for None or blank values.
    """
    if x is None:
        return "-"
    s = x if isinstance(x, str) else str(x)
    s = s.translate(_CTRL_TABLE).strip()
    return s if s else "-"


def _fix_truncated_tail(path):
    """Truncate an incomplete last line before appending on --resume.

    A process killed mid-row leaves the file without a trailing newline.
    Truncate to the last complete line so the resumed run appends cleanly
    and the reader does not choke on a malformed partial row.
    """
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return
    with open(path, "r+b") as fh:
        fh.seek(-1, 2)
        if fh.read(1) == b"\n":
            return  # file ends cleanly — nothing to do
        fh.seek(0, 2)
        size = fh.tell()
        chunk_size = min(65536, size)
        fh.seek(size - chunk_size)
        data = fh.read()
        idx = data.rfind(b"\n")
        if idx == -1 and chunk_size < size:
            # Last newline not in the sampled chunk; scan the whole file.
            fh.seek(0)
            data = fh.read()
            idx = data.rfind(b"\n")
            truncate_at = (idx + 1) if idx != -1 else 0
        elif idx == -1:
            truncate_at = 0  # no newline at all — truncate to empty
        else:
            truncate_at = (size - chunk_size) + idx + 1
        removed = size - truncate_at
        fh.truncate(truncate_at)
    print(
        f"WARNING --resume: removed {removed} byte(s) of incomplete last"
        f" line from {path}",
        file=sys.stderr,
    )


#############
# Orthologs

##
def output_orthologs(annots, orthologs_file, resume, no_file_comments, eggnog_db=None,
                     applied_filters=None):
    start_time = time.time()

    ncbi = get_ncbi(usemem = True)

    if resume == True:
        file_mode = 'a'
        _fix_truncated_tail(orthologs_file)
    else:
        file_mode = 'w'

    with open(orthologs_file, file_mode, encoding="utf-8", newline="\n") as ORTHOLOGS_OUT:
        output_orthologs_header(ORTHOLOGS_OUT, no_file_comments, not resume,
                                applied_filters=applied_filters)

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
    """Write one ortholog row to ``out``.

    The orthologs row only needs the first 10 fields; trailing elements
    (confidence, ceiling, farthest_donor) are unpacked but unused here.

    Handles multiple annotation tuple shapes (11, 12, or 14 elements)
    for resume compatibility.

    Args:
        out: File-like output object.
        annotation: Annotation tuple from ``annotate_batch``.
        ncbi: NCBI taxa accessor for taxid → taxname lookups.
        eggnog_db: Open annotation DB for protein-name decoding.
    """
    n = len(annotation)
    if n == 11:
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs, _confidence) = annotation
    elif n == 12:
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs,
         _confidence, _tax_ceiling) = annotation
    else:
        # 14-element current shape.
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs,
         _confidence, _tax_ceiling,
         _farthest_donor_taxid, _farthest_donor_lineage) = annotation

    # Decode all integer protein IDs (best hit + all orthologs) up front.
    all_ids = {int(best_hit_name)}
    for orths in all_orthologies.values():
        if orths:
            all_ids.update(orths)
    if annot_orthologs:
        all_ids.update(annot_orthologs)
    id_to_name = eggnog_db.decode_protein_ids(all_ids)

    best_hit_real_name = id_to_name.get(int(best_hit_name), str(best_hit_name))
    best_hit_name_id = best_hit_real_name.split(".", 1)[1] if "." in best_hit_real_name else best_hit_real_name

    all_orthologies["annot_orthologs"] = annot_orthologs

    seed_shown = False  # show seed ortholog only once for each query

    for target in all_orthologies:
        if target in ("all", "annot_orthologs"):
            continue

        query_target_orths = all_orthologies[target]
        if not query_target_orths:
            continue

        # Decode integer IDs to get taxids and names
        orthologs_taxids = set()
        orth_by_taxid = {}
        for o in query_target_orths:
            name = id_to_name.get(o, "0.unknown")
            taxid = int(name.split(".")[0])
            orthologs_taxids.add(taxid)
            orth_by_taxid.setdefault(taxid, []).append((o, name))

        orthologs_taxnames = sorted(
            ncbi.get_taxid_translator(orthologs_taxids).items(),
            key=lambda x: x[1],
        )

        for taxid, taxname in orthologs_taxnames:
            orth_names = []
            for orth_id, orth_full_name in orth_by_taxid.get(taxid, []):
                orth_name = orth_full_name.split(".", 1)[1] if "." in orth_full_name else orth_full_name
                if orth_id in (annot_orthologs or ()):
                    orth_name = f"*{orth_name}"
                if orth_name in {best_hit_name_id, f"*{best_hit_name_id}"}:
                    if not seed_shown:
                        row = [query_name, "seed", f"{taxname}({taxid})", orth_name]
                        print("\t".join(_tsv_field(c) for c in row), file=out)
                        seed_shown = True
                else:
                    orth_names.append(orth_name)

            if orth_names:
                row = [query_name, target, f"{taxname}({taxid})", ",".join(sorted(orth_names))]
                print("\t".join(_tsv_field(c) for c in row), file=out)

    return

##
def output_orthologs_header(out, no_file_comments, print_header, applied_filters=None):
    if print_header:
        if not no_file_comments:
            print(get_call_info(), file=out)
            if applied_filters:
                print(format_applied_filters(applied_filters), file=out)
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
# `tax_ceiling`: resolved per-seed ev_lca ceiling clade name (replaces
# the old `tax_scope_used`).  `farthest_donor_taxid` / `farthest_donor_lineage`:
# the taxid and full semicolon-separated lineage of the evolutionarily
# most-distant donor ortholog used.
ANNOTATIONS_WHOLE_HEADER = HIT_HEADER + ANNOTATIONS_HEADER + [
    'annotation_confidence',
    'tax_ceiling',
    'farthest_donor_taxid',
    'farthest_donor_lineage',
]

##
def output_annotations(annots, annot_file, resume, no_file_comments, md5_field, md5_queries,
                       applied_filters=None):

    if resume == True:
        file_mode = 'a'
        _fix_truncated_tail(annot_file)
    else:
        file_mode = 'w'

    try:
        eggnog_db = get_eggnog_db()
    except Exception:
        eggnog_db = None

    start_time = time.time()

    with open(annot_file, file_mode, encoding="utf-8", newline="\n") as ANNOTATIONS_OUT:
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
    """Write one annotation row to ``out``.

    Handles multiple annotation tuple shapes for resume compatibility:

    - 11 elements — v3 cascade (annotations_confidence appended, no ceiling).
    - 12 elements — v3.1 (tax_scope_used appended; mapped to tax_ceiling).
    - 14 elements — current (tax_ceiling + farthest_donor_taxid +
      farthest_donor_lineage).

    Args:
        out: File-like output object.
        annotation: Annotation tuple from ``annotate_batch``.
        md5_field: Whether to append an md5 column.
        md5_queries: Mapping from query name to md5 hex string.
        eggnog_db: Open annotation DB for protein-name decoding.
    """
    # Default values for fields that may be absent in older tuple shapes.
    tax_ceiling = "-"
    farthest_donor_taxid = "-"
    farthest_donor_lineage = "-"

    n = len(annotation)
    if n == 11:
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs,
         annotations_confidence) = annotation
    elif n == 12:
        # Legacy shape: tax_scope_used → reuse as tax_ceiling.
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs,
         annotations_confidence, tax_ceiling) = annotation
    else:
        # Current 14-element shape.
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs,
         annotations_confidence,
         tax_ceiling,
         farthest_donor_taxid,
         farthest_donor_lineage) = annotation

    # Translate integer seed_ortholog to display name.
    seed_display = eggnog_db.get_protein_name(best_hit_name)

    if annotations_confidence:
        conf_str = ";".join(
            f"{field}={tier}"
            for field, tier in sorted(annotations_confidence.items())
        )
    else:
        conf_str = "-"

    annot_columns = [
        query_name,
        seed_display,
        str(best_hit_evalue),
        str(best_hit_score),
        ",".join(match_nog_names),
        str(max_annot_lvl),
        og_cat,
        og_desc,
    ]

    for h in ANNOTATIONS_HEADER:
        if h in annotations and annotations[h] is not None:
            annot_columns.append(",".join(sorted(list(annotations[h]))))
        else:
            annot_columns.append("-")

    annot_columns.append(conf_str)
    annot_columns.append(tax_ceiling or "-")
    annot_columns.append(farthest_donor_taxid or "-")
    annot_columns.append(farthest_donor_lineage or "-")

    if md5_field is True:
        annot_columns.append(md5_queries.get(query_name, "-"))

    print("\t".join(_tsv_field(x) for x in annot_columns), file=out)

    return

##
def output_annotations_header(out, no_file_comments, md5_field, print_header,
                              applied_filters=None):
    if print_header:
        if not no_file_comments:
            print(get_call_info(), file=out)
            # Phase 7.1a — record the resolved values of every annotation-stage
            # filter and threshold the run actually used. Closes the v3 docs gap
            # where reproducibility relied on parsing the raw command line on
            # line 3 (defaults were never expanded).
            if applied_filters:
                print(format_applied_filters(applied_filters), file=out)
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


def read_recorded_sort_entries(path):
    """Read the ``sort_entries`` mode recorded in a file's applied-filters
    header block.

    Scans only the leading comment header (stops at the first data line).
    Returns ``True``/``False`` when ``##   sort_entries=...`` is present, or
    ``None`` when it is absent (e.g. the file was written with
    ``--no_file_comments``, or predates this marker) — in which case the
    caller cannot verify the mode and should not block resume.
    """
    try:
        with open(path, "r") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break  # header ended, data begins
                stripped = line.strip()
                if stripped.startswith("##   sort_entries="):
                    return stripped.split("=", 1)[1].strip() == "True"
    except OSError:
        return None
    return None



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
    """Write one annotation row to ``worksheet``.

    Handles multiple annotation tuple shapes for resume compatibility:

    - 10 elements — pre-refactor (no confidence / ceiling fields).
    - 11 elements — v3 cascade (annotations_confidence appended).
    - 12 elements — v3.1 (tax_scope_used / tax_ceiling appended).
    - 14 elements — current (tax_ceiling + farthest_donor_taxid +
      farthest_donor_lineage).
    """
    # Default values for fields absent in older tuple shapes.
    annotations_confidence = None
    tax_ceiling = "-"
    farthest_donor_taxid = "-"
    farthest_donor_lineage = "-"

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
         all_orthologies, annot_orthologs,
         annotations_confidence) = annotation
    elif n == 12:
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs,
         annotations_confidence, tax_ceiling) = annotation
    else:
        # Current 14-element shape.
        (query_name, best_hit_name, best_hit_evalue, best_hit_score,
         annotations,
         (og_name, og_cat, og_desc),
         max_annot_lvl,
         match_nog_names,
         all_orthologies, annot_orthologs,
         annotations_confidence,
         tax_ceiling,
         farthest_donor_taxid,
         farthest_donor_lineage) = annotation

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
