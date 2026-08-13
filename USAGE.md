# eggNOG-mapper v3 — Usage

Fast functional annotation of novel sequences (proteins, CDS, genomes,
metagenomes) using precomputed **eggNOG 7** orthology assignments. eggNOG-mapper
transfers GO terms, KEGG orthologs/pathways/modules, EC numbers, PFAM domains,
CAZy families, BiGG reactions and more from fine-grained orthologs, which is more
precise than transferring annotations from best-match homologs.

> **Version:** this document tracks eggNOG-mapper **v3** (eggNOG DB **v7**).
> `v3` requires v7 databases; it does **not** work with eggNOG 5 data.

## Contents

- [Three ways to run](#three-ways-to-run)
- [Basic usage](#basic-usage)
- [Databases](#databases)
- [Output files](#output-files)
- [Verify your install (`--selftest`)](#verify-your-install---selftest)
- [Apptainer / Singularity](#apptainer--singularity)
- [Advanced usage](#advanced-usage)
- [Maximum performance on HPC clusters](#maximum-performance-on-hpc-clusters)
- [Data & image downloads](#data--image-downloads)
- [Citation](#citation)

---

## Three ways to run

| | Best for | Install | Compute |
|---|---|---|---|
| **① Free web service** | small jobs, no setup | none | on our servers |
| **② Apptainer image** | reproducible / HPC | one `.sif` file | your machine/cluster |
| **③ pip / source** | integration, dev | Python + deps + DBs | your machine/cluster |

### ① Free cloud service — no install

For small to moderate jobs, run online for free (no download, no install):

**➡️ https://eggnog-mapper.cgmlab.org**

Upload your FASTA, pick options, and download results. Ideal for a handful of
genomes/proteomes when you don't want to host the ~45 GB database.

### ② Apptainer image (recommended for clusters)

A single self-contained image bundles emapper **and** all search tools
(DIAMOND, MMseqs2, HMMER, Prodigal). See
[Apptainer / Singularity](#apptainer--singularity).

```bash
# get the image (see the downloads section for the exact URL)
apptainer run eggnog-mapper-3.0.0.sif emapper.py --version
```

### ③ pip / from source

```bash
# from a source checkout
python -m pip install .
# or editable, for development
python -m pip install -e .
```
Requires Python ≥ 3.9 and the search binaries on `PATH` (`diamond`, and
optionally `mmseqs`, `prodigal`, `hmmer`). Then fetch the databases (below).

---

## Basic usage

Every run needs an **input FASTA** (`-i`), an **output prefix** (`-o`), and the
**databases** (via `--data_dir` or the `EGGNOG_DATA_DIR` environment variable).

```bash
# Proteins (the default input type)
emapper.py -i proteins.faa -o my_run --data_dir /path/to/data --cpu 8
```

The input type is set with `--itype`:

```bash
# Coding sequences (CDS) — auto-translated before the protein search
emapper.py -i cds.fna --itype CDS -o my_run --data_dir DATA --cpu 8

# A genome — genes are predicted, then annotated
emapper.py -i genome.fna --itype genome -o my_run --data_dir DATA --cpu 8

# A metagenome assembly
emapper.py -i contigs.fna --itype metagenome -o my_run --data_dir DATA --cpu 8
```

Common options:

| option | meaning |
|---|---|
| `-i FILE` | input FASTA (proteins by default) |
| `--itype {proteins,CDS,genome,metagenome}` | input type (default `proteins`) |
| `-o PREFIX` | output file prefix (required) |
| `--output_dir DIR` | where outputs go (default: current dir) |
| `--data_dir DIR` | database directory (or set `EGGNOG_DATA_DIR`) |
| `--cpu N` | threads (`--cpu 0` = all cores; default 1) |
| `--resume` | continue a previous run, skipping finished work |
| `--override` | overwrite existing output files |

> **Tip — a genome or contigs, not proteins?** Use `--itype genome` (one
> organism) or `--itype metagenome` (mixed community). Gene prediction runs
> automatically.

---

## Databases

Databases are versioned by **MAJOR.MINOR** on the server and pinned to the
eggNOG DB version:

- **MAJOR** = eggNOG DB version → **v3 uses eggNOG 7**.
- A **MINOR** bump (`3.0` → `3.1`) means the DB structure or annotation sources
  changed — **you must re-download data**.
- **PATCH** releases (`3.0.0` → `3.0.1`) change only emapper code; the whole
  `3.0.x` series **reuses the same `emapper-3.0/` data**.

### Download

```bash
# Core databases: annotation DB + DIAMOND DB + taxonomy (+ prebuilt caches)
python download_eggnog_data.py -y --data_dir /path/to/data
```

`download_eggnog_data.py` derives the folder (`emapper-<MAJOR.MINOR>/`)
automatically from your emapper version. Flags:

| flag | effect |
|---|---|
| *(default)* | annotation DB, DIAMOND DB, taxonomy, GO OBO, prebuilt caches |
| `-D` | **skip** the DIAMOND DB (e.g. you'll only use MMseqs2) |
| `-M` | also install the **MMseqs2** DB (optional; large — see below) |
| `-P` | also install the **Pfam** DB (for `--pfam_realign`) |
| `-y` | assume "yes" |
| `-f` | force re-download |
| `--release 3.0` | pin a specific MAJOR.MINOR data release |
| `--data_dir DIR` | where to install |

**Approximate sizes (core):** annotation DB ~22 GB, DIAMOND DB ~23 GB,
taxonomy + GO OBO + caches < 0.4 GB → **~45 GB total**. Optional MMseqs2 DB
adds ~25 GB; Pfam is separate.

The download includes two prebuilt caches placed next to `eggnog.db`
(`eggnog.db.taxids.bin`, `eggnog.db.fieldpresence.bin`) so your **first run
starts instantly** instead of rebuilding them.

---

## Output files

Written as `<output_dir>/<prefix>.emapper.<suffix>`:

| file | contents |
|---|---|
| `.annotations` | the main functional-annotation table (TSV, 22 columns) |
| `.annotations.xlsx` | same, as Excel (only with `--excel`) |
| `.seed_orthologs` | the DIAMOND/MMseqs2 seed-ortholog hits |
| `.hits` | raw search hits |
| `.orthologs` | per-query ortholog list (only with `--report_orthologs`) |
| `.gff` | gene predictions / decorated GFF (genome/metagenome, or `--decorate_gff`) |
| `.dropped` | queries dropped + reason (only with `--report_dropped`) |

### The `.annotations` columns

| # | column | description |
|--:|---|---|
| 1 | `query` | your input sequence id |
| 2 | `seed_ortholog` | best eggNOG hit (the seed) |
| 3 | `evalue` | seed alignment e-value |
| 4 | `score` | seed alignment bit-score |
| 5 | `eggNOG_OGs` | orthologous groups the seed belongs to, at all levels |
| 6 | `tax_ceiling` | resolved taxonomic ceiling used for the transfer |
| 7 | `farthest_donor_lineage` | lineage of the most distant donor ortholog used |
| 8 | `COG_category` | functional category |
| 9 | `Preferred_name` | consensus gene name |
| 10 | `GOs` | Gene Ontology terms |
| 11 | `EC` | Enzyme Commission numbers |
| 12 | `KEGG_ko` | KEGG orthologs |
| 13 | `KEGG_Pathway` | KEGG pathways |
| 14 | `KEGG_Module` | KEGG modules |
| 15 | `KEGG_Reaction` | KEGG reactions |
| 16 | `KEGG_rclass` | KEGG reaction classes |
| 17 | `BRITE` | KEGG BRITE |
| 18 | `KEGG_TC` | transporter classification |
| 19 | `CAZy` | carbohydrate-active enzymes |
| 20 | `BiGG_Reaction` | BiGG model reactions |
| 21 | `PFAMs` | PFAM domains |
| 22 | `annotation_confidence` | per-field confidence (compact code — see below) |

With `--md5` a final `md5` column (md5 of the query sequence) is appended.

**`annotation_confidence` — compact encoding.** To avoid repeating field names
on every row, confidence is a **fixed positional string, one character per
functional column** (columns 9–21, in order), decoded via the legend written
into the output header:

```
## annotation_confidence: one char per annotation field
## confidence codes: h=high m=medium l=low -=not annotated
## confidence field order: Preferred_name GOs EC KEGG_ko KEGG_Pathway KEGG_Module KEGG_Reaction KEGG_rclass BRITE KEGG_TC CAZy BiGG_Reaction PFAMs
```

e.g. `hlhhhh--h---h` → `Preferred_name=high, GOs=low, EC=high, …, PFAMs=high`.
Decode with `zip(field_order, string)`. Confidence reflects the cascade tier the
value came from (closest donor = high).

> `--no_file_comments` suppresses all `##` header lines (including this legend);
> the field order is fixed and documented here.

---

## Verify your install (`--selftest`)

One command downloads a tiny reference dataset and checks that your build
reproduces the expected annotations for every input type:

```bash
emapper.py --selftest
# or, from the Apptainer image:
apptainer run eggnog-mapper-3.0.0.sif emapper.py --selftest
```

It exits non-zero if anything mismatches. If you already have the bundle
locally, point at it and skip the download: `emapper.py --selftest --data_dir /path/to/test_datasets`.

---

## Apptainer / Singularity

The image is fully self-contained: it bundles emapper **and** DIAMOND, MMseqs2,
HMMER and Prodigal, so nothing else needs installing. Its run-script **is**
`emapper.py`, so you pass emapper arguments straight through.

### Run a prebuilt image

```bash
# annotate (bind your data dir + working dir as needed)
apptainer run \
  --bind /path/to/data:/data \
  eggnog-mapper-3.0.0.sif \
  -i proteins.faa -o my_run --itype proteins \
  --data_dir /data --output_dir "$PWD" --cpu 8

# verify the image end-to-end
apptainer run eggnog-mapper-3.0.0.sif --selftest
```

Prebuilt images are published alongside the data — see
[Data & image downloads](#data--image-downloads).

### Build your own

```bash
cd apptainer
./build.sh                 # release image from committed HEAD  -> eggnog-mapper-<version>.sif
./build.sh --local         # dev image including uncommitted work -> eggnog-mapper-<version>-dev.sif
```
`build.sh` reads the version from `eggnogmapper/version.py` and can override tool
versions (`./build.sh <emapper> <diamond> <mmseqs> <hmmer> <prodigal>`).
Requires Apptainer ≥ 1.1.

> Apptainer needs no root and honours cluster filesystems — ideal for HPC where
> Docker isn't available. Bind your data dir and pass it with `--data_dir`.

---

## Advanced usage

### Search backends

```bash
# DIAMOND (default) — fast, low memory, recommended
emapper.py -i q.faa -m diamond -o out --data_dir DATA

# MMseqs2 — optional; needs the MMseqs2 DB (download with -M)
emapper.py -i q.faa -m mmseqs -o out --data_dir DATA

# No search — annotate an existing seed-orthologs table
emapper.py -m no_search --annotate_hits_table run.emapper.seed_orthologs \
           -o out --data_dir DATA
```

> **DIAMOND vs MMseqs2:** both produce the same annotations. DIAMOND is
> substantially **faster and lighter** at this DB scale and is the default;
> MMseqs2 is provided for users who prefer it and must be downloaded separately.

**Two-step (search once, annotate many ways):**

```bash
emapper.py -i q.faa -m diamond --no_annot -o step1 --data_dir DATA     # search only
emapper.py -m no_search --annotate_hits_table step1.emapper.seed_orthologs \
           -o step2 --data_dir DATA                                    # annotate
```

### Search sensitivity (DIAMOND)

DIAMOND runs an **iterative** search by default: easy queries are caught fast,
and only queries still unaligned escalate to higher sensitivity, up to a ceiling.

```bash
--dmnd_sensmode sensitive     # ceiling sensitivity (default); or fast / very-sensitive / …
--dmnd_iterate {yes,no}       # iterative escalation (default yes); no = single pass
--dmnd_algo ctg               # much faster for SMALL query sets (a few sequences)
```

### Gene prediction (genome / metagenome)

```bash
--genepred search             # default: genes inferred from blastx hits
--genepred prodigal           # predict with Prodigal (supports training)
--trans_table N               # genetic code
--training_genome FILE / --training_file FILE   # Prodigal training
```

### Controlling the annotation transfer

```bash
--tax_scope auto              # per-seed ceiling (default): Eukaryota / Prokaryota
--tax_scope Metazoa           # or a fixed clade name / NCBI taxid
--tax_scope auto-broad        # broad-first (tries Metazoa before Eukaryota for animals)
--donor_pool closest          # first non-empty cascade tier wins (default)
--donor_pool union            # union across tiers
--target_orthologs all        # one2one | many2one | one2many | many2many | all
--target_taxa 2,2157          # restrict donors to these taxa (+descendants)
--excluded_taxa 9606          # exclude donors from these taxa
--annot_evalue 0.001          # seed significance thresholds
--annot_score 60
```

### Extra outputs & QC

```bash
--report_orthologs            # write the per-query .orthologs file
--md5                         # add md5(query sequence) column
--excel                       # also write .xlsx (needs the xlsxwriter package)
--report_dropped              # log dropped queries + reason (.dropped)
--decorate_gff yes            # add hits/annotations to the GFF
--no_file_comments            # omit ## header/stat lines
```

### PFAM realignment (optional, needs the Pfam DB `-P`)

```bash
--pfam_realign realign        # realign queries to domains transferred from orthologs
--pfam_realign denovo         # realign against the whole Pfam DB
```
These spin up an internal `hmmpgmd` server (`--num_servers`, `--num_workers`,
`--port`, `--end_port`).

---

## Maximum performance on HPC clusters

eggNOG-mapper v3 is built for scale. A few facts drive the tuning below:

- **Memory does not scale linearly with `--cpu`.** The large lookup structures
  (the taxid array ~226 MB and the field-presence mask ~119 MB) are built once
  and **copy-on-write shared** across worker processes (fork). Adding CPUs adds
  little RAM.
- **Redundant inputs are deduplicated.** Seeds are sorted and identical seeds
  annotated once (result fanned out), removing ~3× redundant work on large
  proteome/UniProt-scale sets. This is the default.
- **Prebuilt caches ship with the DB** (`eggnog.db.taxids.bin`,
  `eggnog.db.fieldpresence.bin`) so you don't pay the one-time mask build
  (~50 min on the 59 M-protein DB) on first run.

### 1. Use all cores, tune DIAMOND to node RAM

```bash
emapper.py -i q.faa -o out --data_dir DATA --cpu 0 \
           --dmnd_block_size 8 --dmnd_index_chunks 1
```
`--dmnd_block_size`/`--dmnd_index_chunks` are **auto-picked from host RAM** when
unset (bigger block + fewer chunks = faster, more RAM). On a high-memory node,
forcing `--dmnd_block_size 8 --dmnd_index_chunks 1` keeps the whole DB resident.
Rough DIAMOND peak RAM ≈ `block_size × 6 + db_size/index_chunks + threads × 0.5 GB`.

### 2. Keep I/O off the network filesystem

```bash
--temp_dir  /local/scratch     # local SSD for temporary files
--scratch_dir /local/scratch   # compute on local disk, move results to the final dir at the end
```
On NFS/Lustre this matters a lot. If you run many tasks per node, also consider
copying the database to node-local scratch — the ~22 GB DB is memory-mapped and
random reads over a network FS can dominate wall-time.

### 3. Split huge inputs across array jobs

```bash
# split the FASTA, run one array task per chunk, then merge
seqkit split2 -p 100 proteins.faa -O chunks/          # or any splitter
# in each array task i:
emapper.py -i chunks/part_$i.faa -o chunk_$i --data_dir DATA --cpu $SLURM_CPUS_ON_NODE
# merge (keep a single header):
head -n 100 chunk_1.emapper.annotations | grep '^#' > all.emapper.annotations
grep -h -v '^#' chunk_*.emapper.annotations >> all.emapper.annotations
```
Alternatively, run the **search once** and fan out **annotation** with the
two-step `no_search` workflow.

### 4. Restart cheaply

`--resume` continues an interrupted run, skipping already-written results (it
refuses to resume across an output-schema change, telling you to `--override`).

### 5. Containerised, no root

Ship one `.sif` (all tools bundled), bind the data dir, and pass `--data_dir`.
Great for schedulers where you can't install system packages:

```bash
apptainer run --bind $DATA:/data eggnog-mapper-3.0.0.sif \
  -i q.faa -o out --data_dir /data --output_dir "$PWD" --cpu 0 \
  --temp_dir /local/scratch --scratch_dir /local/scratch
```

---

## Data & image downloads

Everything lives under one MAJOR.MINOR folder:

```
https://data.cgmlab.org/eggnog-mapper/emapper-3.0/
├── data/        # eggnog.db, eggnog_proteins.dmnd, taxonomy, go-basic.obo,
│                #   prebuilt caches (taxids.bin, fieldpresence.bin),
│                #   optional mmseqs.tar.gz
├── selftest/    # eggnog-mapper-selftest.tar.gz  (used by --selftest)
└── *.sif        # prebuilt Apptainer images (eggnog-mapper-3.0.x*.sif)
```

- Fetch databases with `download_eggnog_data.py` (it targets this folder
  automatically for your version).
- Fetch the image directly (e.g. `wget`/`curl`) and run it with `apptainer run`.
- The self-test bundle is fetched automatically by `emapper.py --selftest`
  (override the URL with `$EGGNOG_SELFTEST_URL`).

> All `3.0.x` releases share the `emapper-3.0/` data. You only re-download data
> when the **minor** version changes (e.g. `3.1`), which signals a database
> update.

---

## Citation

If you use eggNOG-mapper, please cite:

- Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J.
  **eggNOG-mapper v2: functional annotation, orthology assignments, and domain
  prediction at the metagenomic scale.** *Mol Biol Evol* (2021).
- Hernández-Plaza A, Deng Z, Robledo-Yagüe F, Szklarczyk D, von Mering C, Bork P,
  Huerta-Cepas J. **eggNOG v7: phylogeny-based orthology predictions and
  functional annotations.** *Nucleic Acids Res* (2026).
- Buchfink B, Reuter K, Drost HG. **Sensitive protein alignments at tree-of-life
  scale using DIAMOND.** *Nat Methods* (2021).

The end-of-run message prints the exact emapper + eggNOG DB versions to cite for
your run.

---

<sub>Publishing this page on **GitHub Pages**: enable Pages in the repository
settings (Source: your default branch, folder `/` or `/docs`). This file renders
as-is; to make it the landing page, copy it to `docs/index.md` (or symlink
`README`→`USAGE`). No Jekyll config is required for a single Markdown page.</sub>
