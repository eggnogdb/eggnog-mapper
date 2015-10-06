This is a tool for bulk and efficient mapping of novel sequences to eggNOG v4.1
orthologous groups (OGs). It uses HMMER hmmpgmd to load several precomputed
databases in memory and perform parallelized searches. 

Three databases are currently available: 
- euk (Eukaryotes)
- bact (Bacteria)
- arch (Archea)

### Requirements: 
- BioPython
- HMMER 3 available in your path
- sqlite3 

- ~130GB of disk
- ~130GB of RAM memory (~90GB for the eukaryote specific db, ~32 for the bacterial specific db) 

### Install: 
- Clone or download this repository
```
git clone https://github.com/jhcepas/eggnog-mapper.git
```

- run `upgrade.sh` to download and parse the hmm and annotation databases. It may take a while.
```
cd eggnog-mapper
sh upgrade.sh
```

### Usage: 
1) start the HMMER server in a separate terminal, specifying the number of CPUs available and the target databases to load into memory. Note that the full set of databases will require ~130GB of RAM. 
```
   python server.py --master euk bact arch --worker euk bact arch --cpu 20
```

2) map sequences to eggnog orthologous groups. Using a different shell, use eggnog_mapper as follows: 
```
   python annotate.py fastaFile.fa --output output.hits --db euk bact
```


  
