##
## CPCantalapiedra 2021

import shutil
from sys import stderr as serr
from collections import defaultdict

from ..common import get_version, colorify
from ..annotation.output import ANNOTATIONS_WHOLE_HEADER

DECORATE_GFF_NONE = "no"
DECORATE_GFF_GENEPRED = "yes"

DECORATE_GFF_FIELD_DEFAULT = "ID"
        
##
def run_gff_decoration(mode, gff_ID_field,
                       is_prodigal, is_blastx,
                       gff_genepred_file, gff_genepred_fasta, gff_outfile,
                       predictor, searcher_name, annotated_hits):

    annot_generator = None
    
    ##
    # Generate GFF data to decorate
    if mode == DECORATE_GFF_NONE:
        # Do nothing, just return the hits and annotations
        annot_generator = annotated_hits

    elif mode == DECORATE_GFF_GENEPRED:

        if is_prodigal:
            annot_generator = decorate_gff(gff_genepred_file, gff_ID_field,
                                           gff_outfile, annotated_hits,
                                           get_version(), searcher_name, gff_genepred_fasta)

        elif is_blastx:
            annot_generator = decorate_blastx_gff(annotated_hits, gff_outfile, searcher_name, gff_ID_field)

        else: # proteins, CDS, seeds
            rm_suffix = False

            annot_generator = create_gff(searcher_name, get_version(),
                                         annotated_hits, gff_outfile, 
                                         rm_suffix, gff_ID_field)

    else: # decorate user specified file

        if annotated_hits is None:
            print("No GFF will be created, since there are no annotated hits.")
            
        annot_generator = decorate_gff(mode, gff_ID_field,
                                       gff_outfile, annotated_hits,
                                       get_version(), searcher_name, None)

    return annot_generator



##
# Parse a GFF and create a new one adding hits and/or annotations
def decorate_gff(gff_file, gff_ID_field, outfile, annotated_hits, version, searcher_name, translation_file):

    print(colorify(f"Decorating gff file {gff_file}...", 'lgreen'), file=serr)

    # 1) Parse translation file, if any
    translation_table = None
    # REFACTOR: this should go outside this function
    # which should receive already a translation dictionary
    if translation_file is not None:
        translation_table = {}
        with open(translation_file, 'r') as t_f:
            for line in t_f:
                if line.startswith(">"):
                    annot_id = line.split(" ")[0][1:]
                    gff_id = line.split("ID=")[1].split(";")[0]
                    translation_table[annot_id] = gff_id                    
                # else: continue

    # 2) Parse GFF
    gff_comments = []
    gff_dict = defaultdict(list)
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
            g_score = "." if g_score == "." else float(g_score)

            attrs_list = [attr for attr in g_attrs.split(";") if attr is not None and attr != ""]
            
            record_key = [attr.split("=")[1] for attr in attrs_list
                          if attr.split("=")[0] == gff_ID_field]

            if len(record_key) > 0:
                record_key = record_key[0]
            else:
                continue

            gff_dict[record_key].append([g_seqid, g_source, g_type, g_start, g_end,
                                         g_score, g_strand, g_phase, attrs_list])

    # 3) Parse annotated hits and yield them again
    for hit, annotation in parse_annotations(annotated_hits):
        query = hit[0]
        if translation_table is None:
            gff_id = query
        else:
            gff_id = translation_table[query]
            
        if gff_id in gff_dict:
            for gff_record in gff_dict[gff_id]:
                attrs_list = gff_record[8]
                # include hit
                if hit is not None:
                    (query, target, evalue, score,
                     qstart, qend, sstart, send,
                     pident, qcov, scov,
                     strand, phase, attrs) = hit_to_gff(hit, gff_ID_field)

                    # if there is a translation table, add the emapper query ID also as
                    # an attribute
                    if translation_table is not None:
                        attrs_list.append(f"em_ID={query}")
                        
                    attrs_list.extend(attrs[1:]) # excluding ID

                # include annotations
                if annotation is not None:
                    attrs = annotation_to_gff(annotation)

                    attrs_list.extend(attrs)

        yield hit, annotation

    # 4) Output GFF
    with open(outfile, 'w') as OUT:
        print("##gff-version 3", file=OUT)
        print(f"## decorated with {version}", file=OUT)

        for comment in gff_comments:
            print(comment, file=OUT)
            
        # create a generator to transform the dict of lists to a list of lists
        # then sort and print to output file
        for v in sorted((v for val in gff_dict.values() for v in val),
                        key=lambda x: (x[0], x[3], x[4], 0 if x[5] == "." else x[5], x[-1])):
            fields = v[:-1]
            attrs_list = ";".join(v[-1])
            fields = fields + [attrs_list]
            fields = "\t".join(list(map(str, fields)))
            print(fields, file=OUT)
            
    return

##
# Create GFF file by parsing the hits
##
def create_gff(searcher_name, version, annotated_hits, outfile, rm_suffix, gff_ID_field):

    print(colorify(f"Decorating gff file {outfile}...", 'lgreen'), file=serr)
    
    with open(outfile, 'w') as OUT:

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


# Create the gff for the hits
def create_blastx_hits_gff(hits_generator, outfile, searcher_name, gff_ID_field):

    print(colorify(f"Creating gff file {outfile}...", 'lgreen'), file=serr)

    version = get_version()

    with open(outfile, 'w') as OUT:
        print("##gff-version 3", file=OUT)
        print(f"## created with {version}", file=OUT)

        # The sorted function breaks the generators flow
        # but here it is necessary to sort the gff records by position
        for hit in sorted(parse_hits(hits_generator),
                          key=lambda hit: sort_hits(hit, True)):
            (query, target, evalue, score,
             qstart, qend, sstart, send,
             pident, qcov, scov,
             strand, phase, attrs) = hit_to_gff(hit, gff_ID_field)

            if searcher_name is None:
                attrs.append(f"em_searcher=unk")
            else:
                attrs.append(f"em_searcher={searcher_name}")
                
            contig = query[:query.rfind("_")]

            fields = "\t".join((str(x) for x in [contig, "eggNOG-mapper", "CDS", qstart, qend,
                                                 score, strand, phase, ";".join(attrs)]))

            print(fields, file=OUT)

            yield hit

    return


def decorate_blastx_gff(annotated_hits, outfile, searcher_name, gff_ID_field):
    
    print(colorify(f"Creating {searcher_name} gff file {outfile}...", 'lgreen'), file=serr)

    version = get_version()
    
    with open(outfile, 'w') as OUT:
        print("##gff-version 3", file=OUT)
        print(f"## created with {version}", file=OUT)

        # No need to sort, hits were already sorted in create_blastx_hits_gff
        for hit, annotation in annotated_hits:
        # The sorted function breaks the generators flow
        # but here it is necessary to sort the gff records by position
        # for hit, annotation in sorted(parse_annotations(annotated_hits),
        #                               key=lambda hit: sort_annotated_hits(hit, True)):

            (query, target, evalue, score,
             qstart, qend, sstart, send,
             pident, qcov, scov,
             strand, phase, attrs) = hit_to_gff(hit, gff_ID_field)

            if searcher_name is None:
                attrs.append(f"em_searcher=unk")
            else:
                attrs.append(f"em_searcher={searcher_name}")

            # include annotations
            if annotation is not None:
                attrs.extend(annotation_to_gff(annotation))

            contig = query[:query.rfind("_")]
                
            fields = "\t".join((str(x) for x in [contig, "eggNOG-mapper", "CDS", qstart, qend,
                                                 score, strand, phase, ";".join(attrs)]))
            
            print(fields, file=OUT)

            yield hit, annotation            
    
    return

#
def parse_hits(hits):
    for hit in hits:
        if len(hit) == 4:
            hit = hit + [-1, -1, -1, -1,
                         ".", ".", "."]
        yield hit
    return

def sort_hits(currhit, rm_suffix):
    hit = currhit
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
