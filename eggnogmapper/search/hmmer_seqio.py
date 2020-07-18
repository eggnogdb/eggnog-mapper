##

import os, sys
import re
import gzip

CLEAN_SEQ = re.compile("[\s\-\.]+")
def iter_fasta_seqs(source, translate=False, silent=False):    
    """Iter seq records in a FASTA file"""

    if silent == False:
        sys.stderr.write(f"Parsing fasta file {source}...\n")
    if translate:
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        
    from pathlib import Path
    
    if os.path.isfile(source) or Path(source).is_file():
        if source.endswith('.gz'):
            _source = gzip.open(source, "rt")
        else:
            _source = open(source, "rU")
    else:
        _source = iter(source.split("\n"))

    seq_chunks = []
    seq_name = None
    for line in _source:
        line = line.strip()
        if line.startswith('#') or not line:
            continue       
        elif line.startswith('>'):
            # yield seq if finished
            if seq_name and not seq_chunks:
                raise ValueError("Error parsing fasta file. %s has no sequence" %seq_name)
            elif seq_name:
                if translate:
                    full_seq = ''.join(seq_chunks)
                    prot = Seq(full_seq, generic_dna).translate(to_stop=True)
                    yield seq_name, str(prot)
                else:
                    yield seq_name, ''.join(seq_chunks)
                
            seq_name = line[1:].split()[0].strip()
            seq_chunks = []
        else:
            if seq_name is None:
                raise Exception("Error reading sequences: Wrong format.")
            seq_chunks.append(re.sub(CLEAN_SEQ, '', line))

    # return last sequence
    if seq_name and not seq_chunks:
        raise ValueError("Error parsing fasta file. %s has no sequence" %seq_name)
    elif seq_name:
        if translate:
            full_seq = ''.join(seq_chunks)
            prot = Seq(full_seq, generic_dna).translate(to_stop=True)
            yield seq_name, str(prot)
        else:
            yield seq_name, ''.join(seq_chunks)

    if silent == False:
        sys.stderr.write(f"Fasta file {source} parsing complete.\n")
    
    return
