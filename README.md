[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=eggnog_mapper)

# eggNOG-mapper v3

**eggNOG-mapper** is a tool for fast functional annotation of novel sequences
(proteins, CDS, genomes, metagenomes) using precomputed orthology data from the
[eggNOG database](https://doi.org/10.1093/nar/gkaf1249). Functional terms — GO,
KEGG orthologs/pathways/modules, EC numbers, PFAM domains, CAZy families, BiGG
reactions and more — are transferred from fine-grained orthologs, which is more
precise than best-hit/BLAST transfer because it avoids annotation from close
paralogs.

> **v3 targets the eggNOG 7 database** (~59 million proteins across the tree of
> life). It is **not** compatible with eggNOG 5 data. See **[USAGE.md](USAGE.md)**
> for the complete manual (basic → advanced → HPC → containers).

## Three ways to run

| | Best for | Install |
|---|---|---|
| **① Free web service** — https://eggnog-mapper.cgmlab.org | small/moderate jobs, zero setup | none |
| **② Apptainer/Singularity image** | reproducible runs, HPC clusters | one `.sif` (all tools bundled) |
| **③ pip / from source** | pipelines, development | Python package + databases |

Prebuilt images and databases: **https://data.cgmlab.org/eggnog-mapper/**

## What's new in v3

- **eggNOG 7 database** with integer-encoded orthology and phylogeny-aware
  speciation events (~59M proteins). eggNOG 5 is no longer supported.
- **Curated-only functional donors** — only manually curated terms are used as
  annotation donors, stopping the propagation of automated misannotations, while
  achieving better coverage than v2.
- **Per-seed taxonomic ceiling** (`--tax_scope auto`, default) replaces the old
  predefined scope lists; each seed narrows to its most informative level. Fixed
  clades (`Metazoa`, `33208`, …) are still accepted.
- **Cascade annotation engine** with a **lazy closest-cascade on by default** and
  seed sorting + de-duplication — large speedups on redundant proteome/UniProt
  scale inputs, byte-identical output.
- **Compact 22-column output** with a positional `annotation_confidence` field
  (per-source confidence, documented in the header). See
  [USAGE › Output files](USAGE.md#output-files).
- **`--selftest`** — one command downloads a tiny reference set and verifies the
  install reproduces the expected annotations.
- **Self-contained Apptainer image** bundling DIAMOND, MMseqs2, HMMER and
  Prodigal (`apptainer/build.sh`).
- **Optional MMseqs2 backend**, gzip/bzip2 input autodetection, Cython-accelerated
  inner loops, and `--resume` for interrupted runs.

## Install

```bash
pip install eggnog-mapper
# or from source:
git clone https://github.com/eggnogdb/eggnog-mapper.git
cd eggnog-mapper && pip install .
```

Requires **Python ≥ 3.9** and the search tools on `PATH`:

| tool | role | in v3 |
|---|---|---|
| [DIAMOND](https://github.com/bbuchfink/diamond) | seed search | **default backend** |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | seed search | optional backend (`-m mmseqs`) |
| [Prodigal](https://github.com/hyattpd/Prodigal) | gene prediction | genome/metagenome input |
| [HMMER](http://hmmer.org) | domain realignment | only for `--pfam_realign` |

Install e.g. `conda install -c bioconda diamond mmseqs2 prodigal hmmer`.
**All of these are bundled in the Apptainer image** — nothing to install there.

### Databases

```bash
python download_eggnog_data.py -y --data_dir /path/to/data
```

Data is versioned by **MAJOR.MINOR** (`emapper-3.0/`) and pinned to the eggNOG DB
version: a **minor** bump (`3.0` → `3.1`) means the DB changed and must be
re-downloaded; **patch** releases reuse the same data. Core download ≈ **45 GB**
(annotation DB + DIAMOND DB + taxonomy + prebuilt caches); add `-M` for the
optional MMseqs2 DB, `-P` for Pfam. Details:
[USAGE › Databases](USAGE.md#databases).

## Quick start

```bash
# Proteins (default input type)
emapper.py -i proteins.fa -o my_run --data_dir /path/to/data --cpu 8

# CDS / genome / metagenome
emapper.py -i cds.fna       --itype CDS        -o my_run --data_dir DATA --cpu 8
emapper.py -i genome.fna    --itype genome     -o my_run --data_dir DATA --cpu 8
emapper.py -i contigs.fna   --itype metagenome -o my_run --data_dir DATA --cpu 8

# Verify your installation end-to-end
emapper.py --selftest
```

Output is written as `my_run.emapper.annotations` (22-column TSV) plus optional
`.xlsx` / `.orthologs` / `.gff`. Full options, output-column reference, advanced
usage and HPC performance tuning are in **[USAGE.md](USAGE.md)**.

## Documentation

**[USAGE.md](USAGE.md)** — the complete manual: basic usage, databases, output
format, Apptainer, advanced options, and maximum-performance tips for HPC
clusters.

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

Please also cite the search tool used (the end-of-run message prints the exact
versions):

```
[DIAMOND] Sensitive protein alignments at tree-of-life scale using DIAMOND.
          Buchfink B, Reuter K, Drost HG. 2021.
          Nature Methods 18, 366–368. https://doi.org/10.1038/s41592-021-01101-x

[MMSEQS2] MMseqs2 enables sensitive protein sequence searching for the analysis
          of massive data sets. Steinegger M & Söding J. 2017.
          Nat. Biotech. 35, 1026–1028. https://doi.org/10.1038/nbt.3988

[PRODIGAL] Prodigal: prokaryotic gene recognition and translation initiation
           site identification. Hyatt et al. 2010.
           BMC Bioinformatics 11, 119. https://doi.org/10.1186/1471-2105-11-119
```

## Legacy v2 (eggNOG 5)

For eggNOG 5 databases, use the
[v2 branch](https://github.com/eggnogdb/eggnog-mapper/tree/v2) or the last v2
release:

```bash
pip install eggnog-mapper==2.1.15
```

v2 and v3 databases are **not** interchangeable — v3 only works with eggNOG 7.
