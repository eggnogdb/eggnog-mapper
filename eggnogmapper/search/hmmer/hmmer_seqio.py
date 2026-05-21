##

import bz2
import gzip
import os
import re
import sys

CLEAN_SEQ = re.compile("[\s\-\.]+")

def _fasta_open(source):
    """Open a FASTA source for text reading, auto-detecting gzip/bzip2."""
    from pathlib import Path
    if not (os.path.isfile(source) or Path(source).is_file()):
        return iter(source.split("\n"))
    try:
        with open(source, 'rb') as fh:
            magic = fh.read(3)
    except OSError:
        return open(source, 'r')
    if magic[:2] == b'\x1f\x8b':
        return gzip.open(source, 'rt')
    if magic[:3] == b'BZh':
        return bz2.open(source, 'rt')
    return open(source, 'r')

def iter_fasta_seqs(source, translate=False, silent=False, trans_table = 1):
    """Iter seq records in a FASTA file"""

    if silent == False:
        sys.stderr.write(f"Parsing fasta file {source}...\n")
    if translate:
        from Bio.Seq import Seq
        try:
            from Bio.Alphabet import generic_dna
        except ImportError:
            generic_dna = None

    _source = _fasta_open(source)

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
                    trans_table = trans_table if trans_table is not None else 1
                    
                    if generic_dna is not None:
                        prot = Seq(full_seq, generic_dna).translate(to_stop=True, table=trans_table)
                    else:
                        prot = Seq(full_seq).translate(to_stop=True, table=trans_table)
                        
                    if prot is None or prot == "":
                        print(f"No translation found for sequence {seq_name}", file=sys.stderr)
                    else:
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
            trans_table = trans_table if trans_table is not None else 1
            if generic_dna is not None:
                prot = Seq(full_seq, generic_dna).translate(to_stop=True, table=trans_table)
            else:
                prot = Seq(full_seq).translate(to_stop=True, table=trans_table)
            yield seq_name, str(prot)
        else:
            yield seq_name, ''.join(seq_chunks)

    if silent == False:
        sys.stderr.write(f"Fasta file {source} parsing complete.\n")
    
    return
