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
def decorate_gff(mode, is_prodigal, is_blastx, gff_outfile, predictor, searcher, annotator):
    
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
            create_gff(searcher.name, get_version(), hits, annotations, gff_outfile)

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
        elif len(hits[0]) == 4: # short hits
            print("No GFF will be created from short hits.")
        else:
            create_gff(searcher_name, get_version(), hits, annotations, gff_outfile)

    else: # DECORATE_GFF_FILE --> decorate user specified file
        if ":" in mode:
            parts = mode.split(":")
            decorate_gff_file = parts[0]
            decorate_gff_field = parts[1]
        else:
            decorate_gff_file = mode
            decorate_gff_field = DECORATE_GFF_FIELD_DEFAULT
            
        decorate_gff(decorate_gff_file, decorate_gff_field, gff_outfile, queries_results) # TODO

    return


##
# Create GFF file by parsing the hits
#
def create_gff(searcher_name, version, hits, annotations_dict, outfile):
    hits_dict = {}

    # preload annotations fields to be output in the GFF
    if annotations_dict is not None:
        # print([(name,col) for col,name in enumerate(ANNOTATIONS_WHOLE_HEADER)])
        annotations_fields = [(name,col) for col,name in enumerate(ANNOTATIONS_WHOLE_HEADER)][8:]
        
    with open(outfile, 'w') as OUT:

        print("##gff-version 3", file=OUT)
        print(f"## {version}", file=OUT)

        for hit in sorted(hits, key=lambda x: (x[0],x[4],x[5],x[3])):
            query = hit[0]
            target = hit[1]
            evalue = hit[2]
            score = hit[3]
            qstart = hit[4]
            qend = hit[5]
            sstart = hit[6]
            send = hit[7]
            scov = hit[9]
            if qstart <= qend:
                strand = "+"
            else:
                strand = "-"
                qend = hit[4]
                qstart = hit[5]

            phase = "." # we cannot know the phase as we align against proteins

            if query in hits_dict:
                hits_dict[query] += 1
            else:
                hits_dict[query] = 0
            suffix = hits_dict[query]
            
            attrs = [f"ID={query}", f"score={score}", f"evalue={evalue}",
                     f"eggNOG_target={target}", f"target_cov={scov}",
                     f"sstart={sstart}", f"send={send}"]

            if searcher_name is not None:
                attrs.append(f"searcher={searcher_name}")

            # include annotations
            if annotations_dict is not None and query in annotations_dict:
                annotations = annotations_dict[query]
                annot_attrs = []
                for name,col in annotations_fields:
                    annot_attrs.append(f"{name}={annotations[col]}")
                attrs = attrs + annot_attrs
            
            fields = "\t".join((str(x) for x in [query, "eggNOG-mapper", "CDS", qstart, qend, score, strand, phase, ";".join(attrs)]))
            
            print(fields, file=OUT)

    return

## END
