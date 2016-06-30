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

## Download

- Download the latest version of eggnog-mapper from
  https://github.com/jhcepas/eggnog-mapper/releases. The program does not
  require compilation nor installation.

```
wget https://github.com/jhcepas/eggnog-mapper/releases && tar -zxf latest.tar.gz
```

## eggNOG database retrieval 

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


To download a given database, execute the download script providing a the list
of databases to fetch:

```
download_eggnog_databases.py euk bact arch viruses
```

# Basic Usage

## Mapping and annotation using eggNOG databases

```
python emapper.py -i protsequences.fa --output testCOG0515.hits --db bact
```

## Mapping to custom databases

You can also provide a custom hmmpressed HMMR3 database. For this, just provide
the path and base name of the database (removing the `.h3f` like extension).

```
python emapper.py -i protsequences.fa --output testCOG0515.hits --db custom/database.hmm
```
