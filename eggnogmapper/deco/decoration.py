##
## CPCantalapiedra 2021

import shutil
from sys import stderr as serr

from ..common import get_version, colorify
from ..annotation.output import ANNOTATIONS_WHOLE_HEADER

DECORATE_GFF_NONE = "no"
DECORATE_GFF_GENEPRED = "yes"

DECORATE_GFF_FIELD_DEFAULT = "ID"
        
##
def run_gff_decoration(mode, resume, gff_ID_field, is_prodigal, is_blastx,
                       gff_outfile, predictor, searcher_name, annotated_hits):

    annot_generator = None
    
    ##
    # Generate GFF data to decorate
    if mode == DECORATE_GFF_NONE:
        
        # if prodigal, just rename its GFF file, since it won't decorated
        if is_prodigal:
            if resume == False and predictor is not None:
                shutil.move(predictor.outfile, gff_outfile)
            annot_generator = annotated_hits
            

        # if blastx, create a gff with the hits
        elif is_blastx:
            # rm_suffix is to remove the "_int" added for
            # gene prediction hits (to recover the contig name)
            rm_suffix = True
            annot_generator = create_gff(searcher_name, get_version(),
                                         annotated_hits, gff_outfile, resume,
                                         rm_suffix, gff_ID_field)

        else:
            annot_generator = annotated_hits

    elif mode == DECORATE_GFF_GENEPRED:

        if annotated_hits is None and is_prodigal:
            if resume == False and predictor is not None:
                shutil.move(predictor.outfile, gff_outfile)
            annot_generator = annotated_hits
            
        elif annotated_hits is None:
            print("Hits are required to create a GFF.")
            annot_generator = annotated_hits
            
        else:
            # rm_suffix is to remove the "_int" added for gene prediction hits (to recover the contig name)
            if is_prodigal or is_blastx:
                rm_suffix = True
            else:
                rm_suffix = False

            annot_generator = create_gff(searcher_name, get_version(),
                                         annotated_hits, gff_outfile, resume,
                                         rm_suffix, gff_ID_field)

    else: # decorate user specified file

        if annotated_hits is None:
            print("No GFF will be created, since there are no annotated hits.")
            
        annot_generator = decorate_gff(mode, gff_ID_field,
                                       gff_outfile, annotated_hits,
                                       get_version(), searcher_name)

    return annot_generator

##
# Parse a GFF and create a new one adding hits and/or annotations
def decorate_gff(gff_file, gff_ID_field, outfile, annotated_hits, version, searcher_name):

    print(colorify(f"Decorating gff file {gff_file}...", 'lgreen'), file=serr)

    # 1) Parse GFF
    gff_comments = []
    gff_dict = {}
    with open(gff_file, 'r') as gff_f:
        for line in gff_f:
            if line.startswith("##gff-version"): continue
            if line.startswith("#"):
                gff_comments.append(line.strip())
                continue
            
            (g_seqid, g_source, g_type, g_start, g_end,
             g_score, g_strand, g_phase, g_attrs) = list(map(str.strip, line.split("\t")))

            g_start = int(g_start)
            g_end = int(g_end)
            g_score = float(g_score)

            attrs_list = [attr for attr in g_attrs.split(";") if attr is not None and attr != ""]
            record_key = [attr.split("=")[1] for attr in attrs_list
                          if attr.split("=")[0] == gff_ID_field][0]

            gff_dict[record_key] = [g_seqid, g_source, g_type, g_start, g_end,
                                    g_score, g_strand, g_phase, attrs_list]

    # 2) Parse annotated hits and yield them again
    for hit, annotation in parse_annotations(annotated_hits):
        query = hit[0]
        if query in gff_dict:
            attrs_list = gff_dict[query][8]
            # include hit
            if hit is not None:
                (query, target, evalue, score,
                 qstart, qend, sstart, send,
                 pident, qcov, scov,
                 strand, phase, attrs) = hit_to_gff(hit, gff_ID_field)

                attrs_list.extend(attrs[1:]) # excluding ID which already exists in the GFF attrs
                
            # include annotations
            if annotation is not None:
                attrs = annotation_to_gff(annotation)
                
                attrs_list.extend(attrs)

        yield hit, annotation

    # 3) Output GFF
    with open(outfile, 'w') as OUT:
        print("##gff-version 3", file=OUT)
        print(f"## decorated with {version}", file=OUT)

        for comment in gff_comments:
            print(comment, file=OUT)

        for v in sorted(gff_dict.values(), key=lambda x: (x[0], x[3], x[4], x[5])):
            fields = v[:-1]
            attrs_list = ";".join(v[-1])
            fields = fields + [attrs_list]
            fields = "\t".join(list(map(str, fields)))
            print(fields, file=OUT)
            
    return

##
# Create GFF file by parsing the hits
##
def create_gff(searcher_name, version, annotated_hits, outfile, resume, rm_suffix, gff_ID_field):

    print(colorify(f"Decorating gff file {outfile}...", 'lgreen'), file=serr)
    
    last_resumed_query = None
    if resume == True:
        # find last query in existing file
        with open(outfile, 'r') as gff_f:
            for line in gff_f:
                if line.startswith("#"): continue
                attrs = line.split("\t")[8].strip()
                attrs = {f.split("=")[0]:f.split("=")[1] for f in attrs.split(";")}
                if gff_ID_field in attrs:
                    last_resumed_query = attrs[gff_ID_field]
                    
        file_mode = 'a'
    else:
        file_mode = 'w'

    # semaphore to start processing new hits
    last_resumed_query_found = False if last_resumed_query is not None else True
    
    with open(outfile, file_mode) as OUT:

        print("##gff-version 3", file=OUT)
        print(f"## created with {version}", file=OUT)

        # The sorted function breaks the generators flow
        # but here it is necessary to sort the gff records by position
        for hit, annotation in sorted(parse_annotations(annotated_hits),
                                      key=lambda hit: sort_annotated_hits(hit, rm_suffix)):

            (query, target, evalue, score,
             qstart, qend, sstart, send,
             pident, qcov, scov,
             strand, phase, attrs) = hit_to_gff(hit, gff_ID_field)

            # --resume from last query found
            if last_resumed_query is not None:
                if query == last_resumed_query:
                    last_resumed_query_found = True
                    yield hit, annotation
                    continue
                else:
                    if last_resumed_query_found == False:
                        yield hit, annotation
                        continue
                    else:
                        last_resumed_query = None # start parsing new queries

            if searcher_name is None:
                attrs.append(f"em_searcher=unk")
            else:
                attrs.append(f"em_searcher={searcher_name}")

            # include annotations
            if annotation is not None:
                attrs.extend(annotation_to_gff(annotation))

            if rm_suffix:
                contig = query[:query.rfind("_")]
            else:
                contig = query
                
            fields = "\t".join((str(x) for x in [contig, "eggNOG-mapper", "CDS", qstart, qend,
                                                 score, strand, phase, ";".join(attrs)]))
            
            print(fields, file=OUT)

            yield hit, annotation
    return

#
def parse_annotations(annotated_hits):
    for hit, annotation in annotated_hits:
        if len(hit) == 4:
            hit = hit + [-1, -1, -1, -1,
                         ".", ".", "."]
        yield hit, annotation
    return

#
def sort_annotated_hits(annotated_hit, rm_suffix):
    hit, annotation = annotated_hit
    query = hit[0]
    if rm_suffix == True:
        contig = query[:query.rfind("_")]
    else:
        contig = query
        
    # reverse strand
    if hit[5] < hit[4]:
        ret = (contig, hit[5], hit[4], hit[3])
    else:
        ret = (contig, hit[4], hit[5], hit[3])
        
    return ret

#
def hit_to_gff(hit, gff_ID_field):
    (query, target, evalue, score,
     qstart, qend, sstart, send,
     pident, qcov, scov) = hit
    if qstart <= qend:
        strand = "+"
    else:
        strand = "-"
        qend = hit[4]
        qstart = hit[5]
        
    phase = "." # we cannot know the phase as we align against proteins

    attrs = [f"{gff_ID_field}={query}", f"em_target={target}",
             f"em_score={score}", f"em_evalue={evalue}",
             f"em_tcov={scov}"]
                
    return (query, target, evalue, score, qstart, qend,
            sstart, send, pident, qcov, scov, strand, phase, attrs)

def annotation_to_gff(annotation):
    attrs = []
    (query_name, best_hit_name, best_hit_evalue, best_hit_score,
     annotations,
     (og_name, og_cat, og_desc),
     max_annot_lvl,
     match_nog_names,
     all_orthologies, annot_orthologs) = annotation
    
    match_nog_names = ",".join(match_nog_names)
    attrs.append(f"em_OGs={match_nog_names}")
    attrs.append(f"em_COG_cat={og_cat}")
    attrs.append(f"em_desc={og_desc}")
    attrs.append(f"em_max_annot_lvl={max_annot_lvl}")

    for k, v in annotations.items():
        tag = f"em_{k.replace(' ', '_')}"
        if v is not None:
            value = ",".join(sorted(list(v)))
            attrs.append(f"{tag}={value}")
        
    return attrs

## END
