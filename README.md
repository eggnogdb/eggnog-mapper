[![Build Status](https://travis-ci.com/eggnogdb/eggnog-mapper.svg?branch=master)](https://travis-ci.com/eggnogdb/eggnog-mapper)

# Overview
**eggNOG-mapper** is a tool for fast functional annotation of novel sequences. It uses precomputed orthologous groups and phylogenies from the eggNOG database to transfer functional information from fine-grained orthologs only.

Common uses of eggNOG-mapper include the annotation of novel genomes, transcriptomes or even metagenomic gene catalogs.

The use of orthology predictions for functional annotation permits a higher precision than traditional homology searches (i.e. BLAST searches), as it avoids transferring annotations from close paralogs (duplicate genes with a higher chance of being involved in functional divergence).

Benchmarks comparing different eggNOG-mapper options against BLAST and InterProScan [can be found here](https://github.com/jhcepas/emapper-benchmark/blob/master/benchmark_analysis.ipynb).

EggNOG-mapper is also available as a public online resource: http://eggnog-mapper.embl.de

# Documentation
https://github.com/jhcepas/eggnog-mapper/wiki

# Citation

If you use this software, please cite:
```

[1] Fast genome-wide functional annotation through orthology assignment by
     eggNOG-mapper. Jaime Huerta-Cepas, Kristoffer Forslund, Luis Pedro Coelho,
     Damian Szklarczyk, Lars Juhl Jensen, Christian von Mering and Peer Bork.
     Mol Biol Evol (2017). [doi:
     10.1093/molbev/msx148](https://doi.org/10.1093/molbev/msx148)

[2] eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated
      orthology resource based on 5090 organisms and 2502 viruses. Jaime
      Huerta-Cepas, Damian Szklarczyk, Davide Heller, Ana Hernández-Plaza, Sofia
      K Forslund, Helen Cook, Daniel R Mende, Ivica Letunic, Thomas Rattei, Lars
      J Jensen, Christian von Mering, Peer Bork Nucleic Acids Res. 2019 Jan 8;
      47(Database issue): D309–D314. doi: 10.1093/nar/gky1085 
```

