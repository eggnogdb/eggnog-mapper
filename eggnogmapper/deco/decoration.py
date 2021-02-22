##
## CPCantalapiedra 2021

import shutil

from ..common import get_version
from ..annotation.output import ANNOTATIONS_WHOLE_HEADER

DECORATE_GFF_NONE = "no"
DECORATE_GFF_GENEPRED = "yes"
DECORATE_GFF_FILE = "FILE:FIELD"

DECORATE_GFF_FIELD_DEFAULT = "ID"

##
def run_gff_decoration(mode, is_prodigal, is_blastx, gff_outfile, predictor, searcher, annotator):
    
    ##
    # Generate GFF data to decorate
    if mode == DECORATE_GFF_NONE:
        # if prodigal, just rename its GFF file, since it will be no decorated
        if is_prodigal:
            shutil.move(predictor.outfile, gff_outfile)

        # if blastx, create a gff with the hits
        elif is_blastx:
            hits = searcher.get_hits()
            annotations = None
            rm_suffix = True # rm_suffix is to remove the "_int" added for gene prediction hits (to recover the contig name)
            create_gff(searcher.name, get_version(), hits, annotations, gff_outfile, rm_suffix)

        # else: DO NOTHING

    elif mode == DECORATE_GFF_GENEPRED:
            
        if searcher is None: # -m no_search (IT SHOULD NEVER HAPPEN FOR is_blastx)
            hits = None
            searcher_name = None
        else:
            hits = searcher.get_hits()
            searcher_name = searcher.name
            
        if annotator is None: # --no_annot
            annotations = None
        else:
            annotations = annotator.get_annotations_dict()
            hits = annotator.get_hits()

        if hits is None and annotations is None and is_prodigal:
            shutil.move(predictor.outfile, gff_outfile)
        elif hits is None:
            print("Hits are required to create a GFF.")            
        elif len(hits[0]) == 4: # short hits
            print("No GFF will be created from short hits.")
        else:
            # rm_suffix is to remove the "_int" added for gene prediction hits (to recover the contig name)
            if is_prodigal or is_blastx:
                rm_suffix = True
            else:
                rm_suffix = False
            create_gff(searcher_name, get_version(), hits, annotations, gff_outfile, rm_suffix)

    else: # DECORATE_GFF_FILE --> decorate user specified file
        if ":" in mode:
            parts = mode.split(":")
            decorate_gff_file = parts[0]
            decorate_gff_field = parts[1]
        else:
            decorate_gff_file = mode
            decorate_gff_field = DECORATE_GFF_FIELD_DEFAULT

        if searcher is None: # -m no_search (IT SHOULD NEVER HAPPEN FOR is_blastx)
            hits = None
            searcher_name = None
        else:
            hits = searcher.get_hits_dict()
            searcher_name = searcher.name
            
        if annotator is None: # --no_annot
            annotations = None
        else:
            annotations = annotator.get_annotations_dict()
            hits = annotator.get_hits_dict()

        if hits is None and annotations is None:
            print("No GFF will be created, since there are no hits nor annotations.")
            
        decorate_gff(decorate_gff_file, decorate_gff_field, gff_outfile, hits, annotations, get_version(), searcher_name) # TODO

    return

##
# Parse a GFF and create a new one adding hits and/or annotations
def decorate_gff(gff_file, gff_field, outfile, hits_dict, annotations_dict, version, searcher_name):

    if annotations_dict is not None:
        # print([(name,col) for col,name in enumerate(ANNOTATIONS_WHOLE_HEADER)])
        annotations_fields = [(name,col) for col,name in enumerate(ANNOTATIONS_WHOLE_HEADER)][8:]
        
    with open(outfile, 'w') as OUT, open(gff_file, 'r') as gff_f:
        print("##gff-version 3", file=OUT)
        print(f"## decorated with {version}", file=OUT)
        
        for line in gff_f:
            if line.startswith("##gff-version"): continue
            if line.startswith("#"):
                print(line.strip(), file=OUT)
                continue

            (g_seqid, g_source, g_type, g_start, g_end, g_score, g_strand, g_phase, g_attrs) = list(map(str.strip, line.split("\t")))

            attrs_list = [attr for attr in g_attrs.split(";") if attr is not None and attr != ""]
            attrs_dict = {attr.split("=")[0]:attr.split("=")[1] for attr in attrs_list}
            
            if gff_field in attrs_dict:
                g_field_v = attrs_dict[gff_field]
                if g_field_v in hits_dict:
                    hit = hits_dict[g_field_v]
                    if len(hit) == 4:
                        query, target, evalue, score = short_hit_to_gff(hit)
                        attrs_list.extend([f"em_target={target}", f"em_score={score}", f"em_evalue={evalue}"])
                    elif len(hit) == 11:
                        query, target, evalue, score, qstart, qend, sstart, send, pident, qcov, scov, strand, phase = hit_to_gff(hit)
                        attrs_list.extend([f"em_target={target}", f"em_score={score}",
                                           f"em_evalue={evalue}", f"em_tcov={scov}",
                                           f"em_sstart={sstart}", f"em_send={send}"])
                    if searcher_name is None:
                        attrs_list.append(f"em_searcher=unk")
                    else:
                        attrs_list.append(f"em_searcher={searcher_name}")
                        
                if g_field_v in annotations_dict:
                    annotations = annotations_dict[query]
                    for name,col in annotations_fields:
                        attrs_list.append(f"em_{name}={annotations[col]}")
            else:
                pass
            
            fields = "\t".join((str(x) for x in [g_seqid, g_source, g_type, g_start, g_end, g_score, g_strand, g_phase, ";".join(attrs_list)]))
            
            print(fields, file=OUT)
            
    return

##
# Create GFF file by parsing the hits
#
def create_gff(searcher_name, version, hits, annotations_dict, outfile, rm_suffix):
    
    # preload annotations fields to be output in the GFF
    if annotations_dict is not None:
        # print([(name,col) for col,name in enumerate(ANNOTATIONS_WHOLE_HEADER)])
        annotations_fields = [(name,col) for col,name in enumerate(ANNOTATIONS_WHOLE_HEADER)][8:]
        
    with open(outfile, 'w') as OUT:

        print("##gff-version 3", file=OUT)
        print(f"## created with {version}", file=OUT)

        for hit in sorted(hits, key=lambda x: (x[0],x[4],x[5],x[3])):

            query, target, evalue, score, qstart, qend, sstart, send, pident, qcov, scov, strand, phase = hit_to_gff(hit)
            
            attrs = [f"ID={query}", f"em_target={target}", f"em_score={score}", f"em_evalue={evalue}",
                     f"em_tcov={scov}", f"em_sstart={sstart}", f"em_send={send}"]

            if searcher_name is None:
                attrs.append(f"em_searcher=unk")
            else:
                attrs.append(f"em_searcher={searcher_name}")

            # include annotations
            if annotations_dict is not None and query in annotations_dict:
                annotations = annotations_dict[query]
                for name,col in annotations_fields:
                    attrs.append(f"em_{name}={annotations[col]}")

            if rm_suffix:
                contig = query[:query.rfind("_")]
            else:
                contig = query
                
            fields = "\t".join((str(x) for x in [contig, "eggNOG-mapper", "CDS", qstart, qend, score, strand, phase, ";".join(attrs)]))
            
            print(fields, file=OUT)

    return

#
def hit_to_gff(hit):
    (query, target, evalue, score, qstart, qend, sstart, send, pident, qcov, scov) = hit
    if qstart <= qend:
        strand = "+"
    else:
        strand = "-"
        qend = hit[4]
        qstart = hit[5]
        
    phase = "." # we cannot know the phase as we align against proteins
                
    return query, target, evalue, score, qstart, qend, sstart, send, pident, qcov, scov, strand, phase


def short_hit_to_gff(hit):
    (query, target, evalue, score) = hit
                
    return query, target, evalue, score

## END
