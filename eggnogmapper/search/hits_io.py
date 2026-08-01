##
## CPCantalapiedra 2019

import os
import sys
import time

from ..common import get_call_info
from ..emapperException import EmapperException

##
# Generator of hits from filename
def parse_seeds(filename):
    for lineno, raw in enumerate(open(filename, 'r'), 1):
        if raw.startswith('#') or not raw.strip():
            continue

        line = list(map(str.strip, raw.split('\t')))

        try:
            # short hits: query, target, evalue, score
            if len(line) == 4:
                hit = [line[0], line[1], float(line[2]), float(line[3])]
            # full hits: query, target, evalue, score, qstart, qend,
            # sstart, send, pident, qcov, scov
            elif len(line) == 11:
                hit = [line[0], line[1], float(line[2]), float(line[3]),
                       int(line[4]), int(line[5]), int(line[6]), int(line[7]),
                       float(line[8]), float(line[9]), float(line[10])]
            else:
                # Fail loudly instead of silently re-yielding the previous
                # line's hit (which would misattribute its annotation).
                raise EmapperException(
                    f"Malformed hits table: {filename} line {lineno} has "
                    f"{len(line)} tab-separated fields, expected 4 or 11: "
                    f"{raw.rstrip()!r}")
        except (ValueError, IndexError) as exc:
            raise EmapperException(
                f"Could not parse hits table {filename} line {lineno}: {exc} "
                f"(line: {raw.rstrip()!r})") from exc

        yield hit
    return

##
# Write filtered hits to out_file eagerly (sequential, no parallelism).
# Returns a list of orig_hits for downstream use (e.g. blastx GFF/FASTA
# creation needs the original contig-relative coordinates, not the
# ORF-relative ones written to the seeds file).
def output_seeds(cmds, hits, out_file, no_file_comments, change_seeds_coords=False):
    start_time = time.time()
    orig_hits = []

    with open(out_file, 'w') as OUT:

        # comments
        if not no_file_comments:
            print(get_call_info(), file=OUT)
            if cmds is not None:
                for cmd in cmds:
                    print('##'+cmd, file=OUT)

        # header
        print('#'+"\t".join(("qseqid sseqid evalue bitscore qstart qend "
                             "sstart send pident qcov scov").split(" ")), file=OUT)

        for hit in hits:
            if change_seeds_coords:
                orig_hit = hit
                hit = change_seed_coords(orig_hit)
            else:
                orig_hit = hit

            print('\t'.join(map(str, hit)), file=OUT)
            orig_hits.append(orig_hit)

        elapsed_time = time.time() - start_time
        qn = len(orig_hits)
        if not no_file_comments:
            print('## %d queries scanned' % qn, file=OUT)
            print('## Total time (seconds):', elapsed_time, file=OUT)
            print('## Rate:', "%0.2f q/s" % (float(qn) / elapsed_time if elapsed_time > 0 else 0.0),
                  file=OUT)

    print(f"  {qn} seed orthologs written to {out_file} ({elapsed_time:.1f}s)", file=sys.stderr)
    return orig_hits


##
# Completeness sentinel for --resume. output_seeds writes the
# "## N queries scanned" footer line immediately after the data loop
# finishes (see below), so a newline-terminated occurrence of that line
# proves every seed row was written. This lets --resume reuse a complete
# .seed_orthologs file instead of re-parsing the whole .hits file.
#
# Returns False when the file is missing, empty, has no footer (e.g. it was
# written with --no_file_comments) or the footer is only partially written —
# in every such case the caller falls back to full regeneration, so the check
# is conservative by design (a false negative only costs a rewrite; it can
# never wrongly reuse an incomplete file).
def seed_orthologs_complete(out_file):
    if not os.path.isfile(out_file):
        return False
    try:
        with open(out_file, 'rb') as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            read_n = min(size, 4096)
            fh.seek(size - read_n)
            tail = fh.read(read_n)
    except OSError:
        return False
    # Require the sentinel AND a trailing newline: if the run died mid-footer
    # the last line has no newline, so we treat the file as incomplete.
    return b"queries scanned" in tail and tail.endswith(b"\n")


def change_seed_coords(hit):
    [query, target, evalue, score, qstart, qend, sstart, send, pident, qcov, scov] = hit
    if qstart <= qend:
        qend = qend - (qstart - 1)
        qstart = 1
    else:
        qstart = qstart - (qend - 1)
        qend = 1
    return [query, target, evalue, score, qstart, qend, sstart, send, pident, qcov, scov]

## END
