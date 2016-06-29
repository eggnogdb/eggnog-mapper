# Description

eggnog-mapper is a tool designed for efficient functinal annotation of proteins
using precomputed eggNOG-based orthology assignments. Obvious examples include
the annotation of novel genomes, transcriptomes or even metagenomic sequences.

The tool can also be used for fast mapping of any kind of sequences to
precomputed HMM databases.

This tool is also used for the web services currently hosted at
http://beta-eggnogdb.embl.de/#/app/seqmapper.

# Requirements: 
- Python 2.7+
- BioPython
- HMMER 3 available in your path

Depending on your target databases, 
- ~130GB of disk
- ~130GB of RAM memory (~90GB for the eukaryotes DB, ~32 for the bacteria DB) 

# Install: 
- Clone or download this repository (master branch)
```
git clone https://github.com/jhcepas/eggnog-mapper.git
```

- run `bootstrap.py` to download and parse eggNOG HMM and annotation databases. It may take a while and will require about 150GB of disk space.
```
cd eggnog-mapper
python upgrade.py
```

# Usage: 

## Mapping and annotation using eggNOG databases

### Optimized databases:

Three optimized databases are provided covering Eukaryota, Bacteria and
Archea. These databases include all HMM models from their corresponding
taxonomic levels in eggNOG. Additionally, for those HMM models including more
than 50 eukaryotic or 500 bacterial sequences, fine-grained HMM models covering
other taxonmic eggNOG levels are included.

### Example: annotate a set of bacterial sequences 
```
eggnog_mapper/$ python server.py --db bact

# in a different shell
python eggnog_mapper.py test/testCOG0515.fa --output testCOG0515.hits --db bact --evalue 0.001 --maxhits 20
python annotate.py testCOG0515.hits --output testCOG0515_annotations
```

### Taxonomic restricted databases

eggNOG v4.5 provides a granularity of 107 taxomic levels that can be used to
restrict sequence mapping and functional annotation.  Levels are listed at
http://eggnogdb.embl.de/#/app/downloads. You can refer to any level code as
target database (i.e. maNOG for Mammal). The viral dataset is also available
(code 'viruses').

### Example: annotate a set of bacterial sequences 
```
eggnog_mapper/$ python server.py --db bact

# in a different shell
python eggnog_mapper.py test/testCOG0515.fa --output testCOG0515.hits --db bact --evalue 0.001 --maxhits 20
python annotate.py testCOG0515.hits --output testCOG0515_annotations
```

## Mapping to custom databases

Functional annotation is not available when using a custom HMM dataset. 

### Example: annotate a set of bacterial sequences 
```
eggnog_mapper/$ python server.py --db bact

# in a different shell
python eggnog_mapper.py test/testCOG0515.fa --output testCOG0515.hits --db bact --evalue 0.001 --maxhits 20
python annotate.py testCOG0515.hits --output testCOG0515_annotations
```

