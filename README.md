# Description

eggnog-mapper is a tool designed for bulk and efficient mappings of protein
sequences aginst HMM databases. Several preconfigured databases are provided
based on **eggNOG v4.5 orthologous groups** (OGs), allowing to functionally
annotate novel sequences based on orthology predictions. Furthermore, fast orthology
assignments of novel sequences can be executed against 2031 organisms.

Some of the scripts under this package can be used to map protein
sequences against any random HMMER3 database. 

This tool is also used for the web services currently hosted at
http://beta-eggnogdb.embl.de/#/app/seqmapper.

## Current Tools: 

- `server.py`: preloads eggNOG HMM databases into memory using hmmpgmd from [HMMER3](http://hmmer.org/). Three eggNOG
  taxonomy-restricted databases are currently available: `euk`(Eukaryotes), 
  `bact` (Bacteria), `arch` (Archea).

- `eggnog_mapper.py`: maps protein sequences against the above mentioned 
  eggNOG HMM-databases. `eggnog_mapper.py` can also be used
  as for mapping sequences against custom HMMER3 databases. 
 
- `annotate.py`: reads the output produced by eggnog_mapper and provides functional
  predictions based on eggNOG-based orthology annotations. Functional mappings include GO terms, KEGG pathways, COG
  functional categories and a consensus functional description.

- `refine.py` (experimental): reads the ouput of eggnog_mapper and uses eggNOG precomputed fine
  grained orthology assignments to predict orthologs for the query
  sequences. Additionally, functional annotations and predicted protein names
  (sub-families) can be refined using orthology predictions.


# Requirements: 
- Python 2.7+
- BioPython
- HMMER 3 available in your path
- sqlite3 
- ~130GB of disk
- ~130GB of RAM memory (~90GB for the eukaryotes DB, ~32 for the bacteria DB) 

# Install: 
- Clone or download this repository (master branch)
```
git clone https://github.com/jhcepas/eggnog-mapper.git
```

- run `upgrade.py` to download and parse eggNOG HMM and annotation databases. It may take a while and will require about 150GB of disk space.
```
cd eggnog-mapper
python upgrade.py
```

# Usage: 

## Mapping against the main eggNOG databases: 

1) start the HMMER server in a separate terminal, specifying the number of CPUs available and the target databases to load into memory. Note that the full set of databases will require ~130GB of RAM. 
Loading the databases may take a few minutes. 
```
   python server.py --db euk bact arch --cpu 20
```

2) map sequences to eggnog orthologous groups. While `server.py` is running, use a different shell to execute `eggnog_mapper.py`. Basic usage is as follows:
```
   python eggnog_mapper.py --output output.hits --db euk --evalue 0.001 FASTA_File.fa
```

2) annotate sequences based on the resulting hits using `annotate.py`. GO terms, KEGG pathways, SMART/PFAM domains and general functional descriptions are available. x 
```
   python annotate.py output.hits --ouput my_annotations
```


### Example: annotate a set of bacterial sequences 
```
eggnog_mapper/$ python server.py --db bact

# in a different shell
python eggnog_mapper.py test/testCOG0515.fa --output testCOG0515.hits --db bact --evalue 0.001 --maxhits 20
python annotate.py testCOG0515.hits --output testCOG0515_annotations
```

## Mapping against taxonomically restricted eggNOG databases: 

1) Download any of taxonomicly restricted eggNOG HMM-database from: 
http://beta-eggnogdb.embl.de/download/eggnog_4.5/hmmdb_levels/

2) map sequences:
```
python eggnog_mapper.py test/testCOG0515.fa --output testCOG0515.hits --evalue 0.001 --maxhits 20 --db eggnog_dbs/viruses.all_hmm

``` 
If the database is very large, consider using the `--usemem` parameter, which allows to preloads the database in memory for faster mappings. 
```
python eggnog_mapper.py test/testCOG0515.fa --output testCOG0515.hits --evalue 0.001 --maxhits 20 --db eggnog_dbs/primates.all_hmm --usemem
```

3) annotate functionally
...to be completed...

