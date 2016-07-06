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
- BioPython (not required if input sequences are proteins)

# Storage Requirements:

- ~20GB for the eggNOG annotation database

- ~20GB for eggNOG fasta files

- ~130GB for the three optimized eggNOG databases (euk, bact, arch), and from
  1GB to 35GB for each taxonomic-specific eggNOG HMM database. You don't have to
  download all, just pick the ones you are interested.

(you can check the size of individual datasets at
http://beta-eggnogdb.embl.de/download/eggnog_4.5/hmmdb_levels/)

# Memory requirements:

eggnog-mapper allows to run very fast searches by allocating the target HMM
databases into memory (using the HMMER3 hmmpgmd program). This is enabled when
using of the `--usemem` or `--servermode` flags, and it will require a lot of
RAM memory (depending on the size of the target database). As a reference:

- ~90GB to load the optimized eukaryotic databases (`euk`)
- ~32GB to load the optimized bacterial database (`bact`)
- ~10GB to load the optimized archeal database (`arch`)

Note: 

- Searches are still possible in low memory systems. However, Disk I/O will be
  the a bottleneck (specially when using multiple CPUs). Place databases in
  the fastest disk posible.


# Installation 

## Download

- Download and decompress the latest version of eggnog-mapper from
  https://github.com/jhcepas/eggnog-mapper/releases. The program does not
  require compilation nor installation.

- or clone this git repository:

```
git clone https://github.com/jhcepas/eggnog-mapper.git
```

## eggNOG database retrieval 

- eggNOG mapper provides 107 taxonomically restricted HMM databases (`xxxNOG`),
  three optimized databases [Eukaryota (`euk`), Bacteria (`bact`) and Archea (`arch`)]
  and a virus specific database (`viruses`).

- The three optimized databases include all HMM models from their corresponding
  taxonomic levels in eggNOG (euNOG, bactNOG, arNOG) plus additional models
  spliting large alignments into taxonomically restricted (smaller) HMM
  models. In particular, HMM models with more than 500 (euk) or 50 (bact)
  sequences are expanded. The `arch` database includes all models from all
  archeal taxonomic levels in eggNOG.

- taxonomically restricted databases are listed
  [here](http://beta-eggnogdb.embl.de/#/app/downloads). They can be referred by
  its code (i.e. maNOG for Mammals).


To download a given database, execute the download script providing a the list
of databases to fetch:

```
download_eggnog_data.py euk bact arch viruses
```

This will fetch and decompress all precomputed eggNOG data into the data/ directory. 

# Basic Usage

## Mapping and annotation using eggNOG databases

- Disk based searches on the optimized bacterial database 

```
python emapper.py -i test/polb.fa --output polb_bact -d bact
```

- Disk based searches on the optimized database of viral models

```
python emapper.py -i test/polb.fa --output polb_viruses  -d viruses
```

- Disk based searches on the mammal specific database

```
python emapper.py -i test/p53.fa --output p53_maNOG -d maNOG
```

- Memory based searches on the mammal specific database

```
python emapper.py -i test/p53.fa --output p53_maNOG -d maNOG --usemem
```

- Mapping sequences and skipping functional annotation
```
python emapper.py -i test/p53.fa --output p53_maNOG -d maNOG --hits_only 
```

## Mapping to custom databases

You can also provide a custom hmmpressed HMMR3 database. For this, just provide
the path and base name of the database (removing the `.h3f` like extension).

```
python emapper.py -i test/polb.fa --output polb_pfam -d pfam/pfam.hmm
```

# Output files

eggnog-mapper will produce two text tab-delimited files: 

- `outputname.hits`: A list of HMM hits for each query sequence in the input
  file.

- `outputname.hits.annot`: A list of HMM hits for each query sequence in the
  input file, and their associated eggNOG functional description and COG
  functional category.

- `outputname.annot`: A list of predicted fine-grained orthologs, gene names, GO
  terms and KEGG pathways associated to each query sequence. (*the viruses
  database does not support fine grained annotation*)

## Example Output:

```bash
$ python emapper.py -d euk -i test/polb.fa -o polb_bact --override
Sequence mapping starts now!
processed queries:0 total_time:86.3751051426 rate:0.00 q/s
Functional annotation starts now!
Done


$ head polb_bact.*

==> polb_bact.hits <==
# Fri Jul  1 18:35:52 2016
# emapper.py -d euk -i test/polb.fa -o polb_bact --override
# #query_name   hit     evalue  sum_score       query_length    hmmfrom hmmto   seqfrom seqto   query_coverage
362663.ECP_0061 euNOG.KOG0970.meta_raw  3.1e-87 297.3   783     543     1152    193     782     0.752234993614
362663.ECP_0061 euNOG.KOG0969.meta_raw  7e-56   193.3   783     374     962     204     768     0.72030651341
362663.ECP_0061 euNOG.KOG0968.meta_raw  3.6e-41 144.1   783     1271    1549    148     428     0.357598978289
362663.ECP_0061 euNOG.KOG0968.meta_raw  3.6e-41 144.1   783     1644    1912    490     759     0.343550446999
362663.ECP_0061 virNOG.ENOG411DSZ9.meta_raw     7.1e-14 54.4    783     988     1033    193     238     0.0574712643678
362663.ECP_0061 virNOG.ENOG411DSZ9.meta_raw     7.1e-14 54.4    783     1234    1259    404     428     0.0306513409962
362663.ECP_0061 virNOG.ENOG411DSZ9.meta_raw     7.1e-14 54.4    783     1357    1498    491     643     0.194125159642

==> polb_bact.hits.annot <==
#query_name     hit     level   evalue  sum_score       query_length    hmmfrom hmmto   seqfrom seqto   query_coverage  members_in_og   og_description  og_COG_categories
362663.ECP_0061 KOG0970 euNOG   3.1e-87 297.3   783     543     1152    193     782     0.752234993614  271     DNA polymerase  L
362663.ECP_0061 KOG0969 euNOG   7e-56   193.3   783     374     962     204     768     0.72030651341   224     DNA polymerase  L
362663.ECP_0061 KOG0968 euNOG   3.6e-41 144.1   783     1271    1549    148     428     0.357598978289  345     DNA polymerase  L
362663.ECP_0061 KOG0968 euNOG   3.6e-41 144.1   783     1644    1912    490     759     0.343550446999  345     DNA polymerase  L
362663.ECP_0061 1DSZ9   virNOG  7.1e-14 54.4    783     988     1033    193     238     0.0574712643678 20      DNA polymerase  L
362663.ECP_0061 1DSZ9   virNOG  7.1e-14 54.4    783     1234    1259    404     428     0.0306513409962 20      DNA polymerase  L
362663.ECP_0061 1DSZ9   virNOG  7.1e-14 54.4    783     1357    1498    491     643     0.194125159642  20      DNA polymerase  L
362663.ECP_0061 1BF16   strNOG  2.3e-13 52.6    783     964     1009    193     238     0.0574712643678 19      DNA polymerase  L
362663.ECP_0061 1BF16   strNOG  2.3e-13 52.6    783     1213    1239    404     429     0.0319284802043 19      DNA polymerase  L

==> polb_bact.annot <==
#query_name     best_hit_eggNOG_ortholog        best_hit_evalue best_hit_score  predicted_name  strict_orthologs        GO      KEGG(pathway)
362663.ECP_0061 39946.BGIOSGA039704-PA  2.2e-248        820.7   POLB    1006551.KOX_10750,1028307.EAE_11080,1045856.EcWSU1_00673,155864.Z0068,198214.SF0055,198628.Dda3937_01373,199310.c0071,218491.ECA3852,218493.SBG_0085,220341.STY0112,272620.KPN_00059,290338.CKO_03322,290339.ESA_03280,316407.85674307,362663.ECP_0061,39946.BGIOSGA039704-PA,399742.Ent638_0607,406817.XNC1_4058,406818.XBJ1_1795,465817.ETA_07310,469595.CSAG_03356,469613.HMPREF0864_02496,471874.PROSTU_01143,481805.EcolC_3597,498217.ETAE_0607,500637.PROVRUST_06038,500639.ENTCAN_05148,500640.CIT292_09418,502347.ESCAB7627_3201,511145.b0060,517433.PanABDRAFT_1681,520999.PROVALCAL_02276,521000.PROVRETT_06774,529507.PMI2327,561229.Dd1591_0570,561230.PC1_3628,561231.Pecwa_3819,579405.Dd703_0609,590409.Dd586_3566,592316.Pat9b_0640,634499.EpC_07180,634500.EbC_06980,634503.NT01EI_0703,637910.ROD_00651,640513.Entas_0660,665029.EAMY_2919,701347.Entcl_3665,706191.PANA_0690,712898.Pvag_0098,716541.ECL_00856,741091.Rahaq_3758,745277.Rahaq2_3852,882884.SARI_02907,99287.STM0097  GO:0003674,GO:0003824,GO:0003887,GO:0004518,GO:0004527,GO:0004529,GO:0004536,GO:0005575,GO:0005622,GO:0005623,GO:0005694,GO:0006139,GO:0006259,GO:0006260,GO:0006261,GO:0006281,GO:0006289,GO:0006297,GO:0006301,GO:0006725,GO:0006807,GO:0006950,GO:0006974,GO:0007154,GO:0008150,GO:0008152,GO:0008296,GO:0008408,GO:0009058,GO:0009059,GO:0009432,GO:0009605,GO:0009987,GO:0009991,GO:0016740,GO:0016772,GO:0016779,GO:0016787,GO:0016788,GO:0016796,GO:0016895,GO:0018130,GO:0019438,GO:0019985,GO:0031668,GO:0033554,GO:0034061,GO:0034641,GO:0034645,GO:0034654,GO:0043170,GO:0043226,GO:0043228,GO:0043229,GO:0043232,GO:0044237,GO:0044238,GO:0044249,GO:0044260,GO:0044271,GO:0044424,GO:0044464,GO:0044699,GO:0044763,GO:0045004,GO:0045005,GO:0046483,GO:0050896,GO:0051716,GO:0071496,GO:0071704,GO:0071897,GO:0090304,GO:0090305,GO:1901360,GO:1901362,GO:1901576
# Total time (seconds): 4.29354810715

```



# Advance usage

## Speeding up annotation using memory based multi-threaded based searches.

If only one input file is going to be annotated, simply use the `--usemem` and
`--cpu XX` options. For instance: 

```
python emapper.py -i test/polb.fa --output polb_pfam -d pfam/pfam.hmm --usemem --cpu 10
``` 

If you are planning to use the same database for annotating multiple files, you
can start eggnog-mapper in server mode (this will load the target database in
memory and keep it there until stopped). Then you can use another eggnog-mapper
instance to connect to the server. For instance, 

In terminal 1, execute:

```
python emapper.py -d arch --cpu 10 --servermode
```

This will load the memory and give you the address to connect to the
database. Then, in a different terminal, execute:

```
python emapper.py -d arch:localhost:51600 -i test/polb.fa -o polb_arch
```

# Citation

eggNOG 4.5: a hierarchical orthology framework with improved functional annotations for eukaryotic, prokaryotic and viral sequences.
Jaime Huerta-Cepas, Damian Szklarczyk, Kristoffer Forslund, Helen Cook, Davide Heller, Mathias C. Walter, Thomas Rattei, Daniel R. Mende, Shinichi Sunagawa, Michael Kuhn, Lars Juhl Jensen, Christian von Mering, and Peer Bork.
Nucl. Acids Res. (04 January 2016) 44 (D1): D286-D293. doi: 10.1093/nar/gkv1248
