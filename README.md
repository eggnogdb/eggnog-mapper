This is a tool for bulk and efficient mapping of novel sequences to eggNOG v4.1
orthologous groups (OGs). It uses HMMER hmmpgmd to load several precomputed
databases in memory and perform parallelized searches. Note that **loading all
databases requires about 150GB of RAM memory.**

Three databases are currently available: 
- euk (Eukaryotes)
- bact (Bacteria)
- arch (Archea)



### Requirements: 
- BioPython
- HMMER 3 available in your path
- sqlite3 

- ~150GB of disk
- ~150GB of RAM memory (~90GB for the eukaryote specific db, ~32 for the bacterial specific db) 

### Install: 
- Clone or download this repository
```
git clone https://github.com/jhcepas/eggnog-mapper.git
```

- run `upgrade.sh` and follow instructions
```
cd eggnog-mapper
sh upgrade.sh
```

### Usage: 
- start the HMMER server, specifying number of CPUs available and databases to load into memory. Note that the full set of databases will require ~130GB of RAM. 
 
   python server.py --master euk bact arch --worker euk bact arch --cpu 20

- Using a different shell, annotate sequences using the annotate.py script: 

   python annotate.py fastaFile.fa --output output.hits --db euk bact



  
