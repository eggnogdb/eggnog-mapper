# Description

eggnog-mapper is a tool designed for fast functinal annotation of proteins using
precomputed eggNOG-based orthology assignments. Obvious examples include the
annotation of complete novel genomes, transcriptomes or even metagenomic gene
catalogs.

The tool can also be used for fast mapping of of sequences against custom HMM
databases. 

This tool is also used for the web services currently hosted at
http://beta-eggnogdb.embl.de/#/app/seqmapper.

# Software Requirements: 

- Python 2.7+
- wget 
- HMMER 3 binaries available in your PATH
- BioPython (not required if input file contain protein sequences)

# Storage Requirements:
 
- ~130GB of disk for the three main eggNOG databases (euk, bact, arch).

- from 1 to 40 GBs for each taxonomic-specific eggNOG HMM databases. You don't
  have to download all, just pick the ones you are interested.

(you can check the size of individual datasets at
http://beta-eggnogdb.embl.de/download/eggnog_4.5/hmmdb_levels/)

# Memory requirements:

eggnog-mapper allows to run very fast searches by allocating the HMM databases
into memory (using HMMER3 hmmpgmd program). This is enabled by means of the
`--usemem` flag and will require a lot of RAM memory (depending on the size of
the target database). As a reference:

- ~90GB to load the optimized eukaryotic databases (euk)
- ~32GB to load the optimized bacterial database (bact)
- ~10GB to load the optimized archeal database (arch)

Note: 

- Searches are still possible in low memory systems. However, Disk I/O will be
  the a bottleneck (specially when using multiple CPUs). Place databases in
  the fastest disk posible.


# Installation 

- Download the lastest version of eggnog-mapper from https://github.com/jhcepas/eggnog-mapper/releases

```
wget https://github.com/jhcepas/eggnog-mapper/releases && tar -zxf latest.tar.gz
```

# eggNOG databases retrieval 

- eggNOG mapper provides 107 taxonomically restricted HMM databases (`xxxNOG`),
  three optimized databases (Eukaryota-`euk`, Bacteria-`bact` and Archea-`arch`)
  and a virus specific database (`viruses`).

- The three optimized databases include all HMM models from their corresponding
  taxonomic levels in eggNOG (euNOG, bactNOG, arNOG) plus additional models
  spliting large alignments into taxonomically restricted (smaller) HMM
  models. In particular, HMM models with more than 500 (euk) or 50 (bact)
  sequences are expanded. The `arch` database includes all models from all
  archeal taxonomic levels in eggNOG.

- taxonomically restricted databases are listed here and can be referred by its
  code (i.e. maNOG for Mammals).


To download a given database, execute: 
```
 execute `download_eggnog_databases.py` to retrieve your preferred eggNOG
  database. This may take a while and may require substantial . For example: 

```
download_eggnog_databases.py euk bact arch viruses
```

# Usage: 

## Mapping and annotation using eggNOG databases

### Using the optimized databases:


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

