[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=eggnog_mapper)

# eggNOG-mapper v3

**eggNOG-mapper** is a tool for fast functional annotation of novel sequences using
precomputed orthologous groups and phylogenies from the
[eggNOG database](https://doi.org/10.1093/nar/gkaf1249).
Functional information is transferred exclusively from fine-grained orthologs,
yielding higher precision than homology-based approaches (e.g. BLAST) by avoiding
annotation transfer from close paralogs.

Common uses include annotation of novel genomes, transcriptomes, and metagenomic
gene catalogs.

eggNOG-mapper is also available as a public web server: http://mapper.eggnogdb.org

## What's new in v3

v3 is a major release targeting the **eggNOG v7** database and a completely
redesigned annotation engine.

- **eggNOG v7 database** — integer-encoded orthology, phylogeny-aware speciation
  events, and ~12M proteins across ~10k taxa. eggNOG v5 databases are no longer
  supported.
- **Per-seed taxonomic ceiling** — replaces the old `--tax_scope` predefined
  scope lists. Each query seed gets its own `ev_lca`-based ceiling automatically
  narrowed to the most informative phylogenetic level (`--tax_scope auto`,
  default). Fixed clades (`Metazoa`, `33208`, etc.) are still accepted.
- **Cascade annotation engine** — for each functional source (GO, KEGG, Pfam,
  EC, …) donors are walked from closest+best-typed first, with the seed's own
  curated annotation as the strongest tier-0 donor.
- **No bundled binaries** — DIAMOND, HMMER, MMseqs2, and Prodigal must be
  installed externally (see Requirements below). The wheel shrinks from ~150 MB
  to ~5 MB and cross-platform installs (macOS, Windows) now work.
- **Compressed input** — gzip and bzip2 FASTA inputs are autodetected by magic
  bytes.
- **Parallel annotation** — `--cpu N` parallelises both search and annotation.
- **Cython-accelerated inner loops** — `_codec` and `_collect_inner` extensions
  give ~2–3× speedup on the annotation phase.
- **`--resume`** — safely resumes an interrupted run, reusing the existing hits
  file.
- **Apptainer/Singularity image** — a self-contained HPC image is provided via
  `apptainer/build.sh`.

## Requirements

- Python ≥ 3.9
- At least one search backend:

| Tool | Install |
|------|---------|
| [DIAMOND](https://github.com/bbuchfink/diamond) | `conda install -c bioconda diamond` |
| [HMMER](http://hmmer.org) | `conda install -c bioconda hmmer` |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | `conda install -c bioconda mmseqs2` |
| [Prodigal](https://github.com/hyattpd/Prodigal) | `conda install -c bioconda prodigal` (gene prediction only) |

## Installation

```bash
pip install eggnog-mapper
```

Or from source:

```bash
git clone https://github.com/eggnogdb/eggnog-mapper.git
cd eggnog-mapper
pip install .
```

### Download the eggNOG v7 database

```bash
download_eggnog_data.py --data_dir /path/to/eggnog-data
```

## Quick start

```bash
# Protein sequences against eggNOG v7 using DIAMOND
emapper.py -m diamond -i proteins.fa --itype proteins \
    --data_dir /path/to/eggnog-data \
    -o my_annotation --output_dir results/ --cpu 20

# Two-step: search first, annotate later
emapper.py -m diamond -i proteins.fa --itype proteins \
    --data_dir /path/to/eggnog-data \
    -o my_annotation --output_dir results/ --no_annot --cpu 20

emapper.py -m no_search --annotate_hits_table results/my_annotation.emapper.seed_orthologs \
    --data_dir /path/to/eggnog-data \
    -o my_annotation --output_dir results/
```

## Documentation

https://github.com/eggnogdb/eggnog-mapper/wiki

## Citation

If you use eggNOG-mapper, please cite:

```
[1] eggNOG-mapper v2: functional annotation, orthology assignments, and domain
    prediction at the metagenomic scale. Carlos P. Cantalapiedra,
    Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021.
    Molecular Biology and Evolution, msab293, https://doi.org/10.1093/molbev/msab293

[2] eggNOG v7: phylogeny-based orthology predictions and functional annotations.
    Ana Hernández-Plaza, Ziqi Deng, Fabian Robledo-Yagüe, Damian Szklarczyk,
    Christian von Mering, Peer Bork, Jaime Huerta-Cepas. Nucleic Acids Research,
    Volume 54, Issue D1, 6 January 2026, Pages D402-D408.
    https://doi.org/10.1093/nar/gkaf1249
```

Please also cite the search tool used:

```
[DIAMOND] Sensitive protein alignments at tree-of-life scale using DIAMOND.
          Buchfink B, Reuter K, Drost HG. 2021.
          Nature Methods 18, 366–368. https://doi.org/10.1038/s41592-021-01101-x

[HMMER]   Accelerated Profile HMM Searches.
          Eddy SR. 2011. PLoS Comput. Biol. 7:e1002195.

[MMSEQS2] MMseqs2 enables sensitive protein sequence searching for the analysis
          of massive data sets. Steinegger M & Söding J. 2017.
          Nat. Biotech. 35, 1026–1028. https://doi.org/10.1038/nbt.3988

[PRODIGAL] Prodigal: prokaryotic gene recognition and translation initiation
           site identification. Hyatt et al. 2010.
           BMC Bioinformatics 11, 119. https://doi.org/10.1186/1471-2105-11-119
```

## Legacy v2 (eggNOG v5)

If you are working with eggNOG v5 databases, use the
[v2 branch](https://github.com/eggnogdb/eggnog-mapper/tree/v2)
or install the last v2 release from PyPI:

```bash
pip install eggnog-mapper==2.1.15
```

v2 and v3 databases are not interchangeable. v3 only works with eggNOG v7.
