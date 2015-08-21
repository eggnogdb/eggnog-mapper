This software allows for bulk mapping of novel sequences to eggNOG 4.1 orthologous groups. 
Three databases are available: euk (Eukaryotes), bact (Bacteria) and arch (Archea).  

### Installation: 

- Install BioPython
- Install HMMER 3
- Download and decompress eggNOG HMM databases from http://eggnogdb.embl.de/download/eggnog_4.1/hmmdb.euk_bact_arch.tar.gz (it will require ~108GB). 
- If necessary, mv the uncompressed hmmdb/ directory to the root directory of this tool. 

### Usage: 
- start the HMMER server, specifying number of CPUs available and databases to load into memory. Note that the full set of databases will require ~130GB of RAM. 
 
python server.py --master euk bact arch --worker euk bact arch --cpu 20

- start annotating sequences using the annotate.py script: 

python annotate.py fastaFile.fa --output output.hits --db bact

 
  
