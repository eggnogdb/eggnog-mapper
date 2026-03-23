# eggnog-mapper Usage Guide

## Overview

eggnog-mapper is a tool for functional annotation of novel sequences (proteins,
CDS, genomes, or metagenomes) using precomputed eggNOG orthology assignments.
It maps query sequences to eggNOG orthologous groups (OGs) via fast sequence
searches, then transfers functional annotations — GO terms, KEGG pathways, Pfam
domains, EC numbers, and more — from known orthologs.

## Quick Start

### 1. Install databases

```bash
# Download the default eggNOG 5 databases (annotation DB + DIAMOND DB + taxonomy DB)
python download_eggnog_data.py -y
```

### 2. Annotate protein sequences

```bash
# Basic protein annotation using DIAMOND (default)
python emapper.py -i my_proteins.fasta -o my_output

# Use multiple CPUs
python emapper.py -i my_proteins.fasta -o my_output --cpu 8
```

### 3. Annotate a genome

```bash
# Gene prediction + annotation from a genome FASTA
python emapper.py -i genome.fna --itype genome -o my_output --cpu 8

# Using Prodigal for gene prediction instead of the default blastx-based method
python emapper.py -i genome.fna --itype genome --genepred prodigal -o my_output --cpu 8
```

### 4. Quick test run

```bash
# Run against the bundled test database (no downloads needed)
python emapper.py --db e5-test -i tests/fixtures/test_queries.fa -o test_run --output_dir /tmp
```

### Output files

| File                          | Contents                                          |
|-------------------------------|---------------------------------------------------|
| `*.emapper.seed_orthologs`    | Best hits (seed orthologs) from the search phase   |
| `*.emapper.annotations`       | Full functional annotations (TSV)                  |
| `*.emapper.orthologs`         | Orthologs per query (if `--report_orthologs`)      |
| `*.emapper.pfam`              | Pfam domain predictions (if `--pfam_realign`)      |
| `*.emapper.genepred.fasta`    | Predicted proteins (genome/metagenome input)       |
| `*.emapper.genepred.gff`      | Gene prediction GFF (genome/metagenome input)      |
| `*.emapper.decorated.gff`     | GFF with annotations (if `--decorate_gff`)         |

---

## Detailed / Advanced Usage

### Search Modes

eggnog-mapper supports several search engines, selected with `-m`:

```bash
# DIAMOND (default, fast and sensitive)
python emapper.py -m diamond -i seqs.fa -o out

# MMseqs2 (fast, good for very large datasets)
python emapper.py -m mmseqs -i seqs.fa -o out

# HMMER (profile-based, most accurate for specific taxonomic groups)
python emapper.py -m hmmer -i seqs.fa -d bact -o out

# Novel families (for uncharacterized protein families)
python emapper.py -m novel_fams -i seqs.fa -o out

# No search: re-annotate a previous seed_orthologs file
python emapper.py -m no_search --annotate_hits_table prev.emapper.seed_orthologs -o out

# Cache-based annotation
python emapper.py -m cache -i seqs.fa -c cache_file.txt -o out
```

### Database Backend Selection

The `--db` option selects which database backend to use. This determines the
data directory where all database files are looked up:

```bash
# Default: use the full eggNOG 5 databases in data/
python emapper.py --db eggnog5 -i seqs.fa -o out

# Use the bundled test database (tiny, for CI/testing)
python emapper.py --db e5-test -i tests/fixtures/test_queries.fa -o out --output_dir /tmp
```

The precedence for database path resolution is:

    --data_dir  >  EGGNOG_DATA_DIR env var  >  --db  >  default (eggnog5)

You can always override with `--data_dir` to point to a custom database
directory, regardless of `--db`:

```bash
python emapper.py --data_dir /path/to/my/custom/data -i seqs.fa -o out
```

### Input Types

```bash
# Protein sequences (default)
python emapper.py -i proteins.fa --itype proteins -o out

# CDS (coding DNA sequences)
python emapper.py -i cds.fna --itype CDS -o out

# CDS translated to proteins before search
python emapper.py -i cds.fna --itype CDS --translate -o out

# Genome
python emapper.py -i genome.fna --itype genome -o out

# Metagenome
python emapper.py -i metagenome.fna --itype metagenome -o out
```

### Gene Prediction (for genome/metagenome inputs)

When using `--itype genome` or `--itype metagenome`, genes must be predicted
before annotation:

```bash
# Blastx-based gene prediction (default): infers genes from DIAMOND/MMseqs2 hits
python emapper.py -i genome.fna --itype genome --genepred search -o out

# Prodigal-based gene prediction
python emapper.py -i genome.fna --itype genome --genepred prodigal -o out

# Prodigal with a training genome
python emapper.py -i genome.fna --itype genome --genepred prodigal \
    --training_genome reference.fna -o out

# Control overlap handling in blastx-based prediction
python emapper.py -i genome.fna --itype metagenome --genepred search \
    --allow_overlaps diff_frame --overlap_tol 0.5 -o out
```

### Taxonomic Scope and Ortholog Filtering

Taxonomic scope controls which clades are used for functional transfer:

```bash
# Restrict annotation to bacterial orthologs
python emapper.py -i seqs.fa --tax_scope bacteria -o out

# Use a specific list of taxonomic IDs
python emapper.py -i seqs.fa --tax_scope 2759,2157 -o out

# Control which clade level is used for annotation
python emapper.py -i seqs.fa --tax_scope bacteria --tax_scope_mode broadest -o out

# Only use one-to-one orthologs for annotation
python emapper.py -i seqs.fa --target_orthologs one2one -o out

# Restrict annotation to specific target taxa
python emapper.py -i seqs.fa --target_taxa 9606,10090 -o out

# Exclude specific taxa from annotation
python emapper.py -i seqs.fa --excluded_taxa 9606 -o out
```

### DIAMOND Tuning

```bash
# Ultra-sensitive search
python emapper.py -i seqs.fa --sensmode ultra-sensitive -o out

# Iterative search (starts fast, increases sensitivity)
python emapper.py -i seqs.fa --dmnd_iterate yes --sensmode very-sensitive -o out

# Small query optimization
python emapper.py -i small_set.fa --dmnd_algo ctg -o out

# Use a custom DIAMOND database
python emapper.py -i seqs.fa --dmnd_db /path/to/custom.dmnd -o out

# Enable frameshift-aware alignment (for error-prone sequences)
python emapper.py -i seqs.fa --dmnd_frameshift 15 -o out

# Tune memory/performance
python emapper.py -i seqs.fa --block_size 4.0 --index_chunks 1 -o out --cpu 16
```

### MMseqs2 Tuning

```bash
# Custom sensitivity ramp
python emapper.py -m mmseqs -i seqs.fa --start_sens 1 --sens_steps 5 --final_sens 9 -o out

# Custom database and substitution matrix
python emapper.py -m mmseqs -i seqs.fa --mmseqs_db /path/to/db --mmseqs_sub_mat blosum62.out -o out
```

### HMMER Usage

```bash
# Search against a specific taxonomic-level HMM database
python emapper.py -m hmmer -i seqs.fa -d bact -o out

# Load database into memory for faster searches
python emapper.py -m hmmer -i seqs.fa -d bact --usemem -o out

# Use a remote HMMER server
python emapper.py -m hmmer -i seqs.fa -d mydb:hostname:51700 -o out

# Multiple servers and workers
python emapper.py -m hmmer -i seqs.fa -d bact --usemem \
    --num_servers 2 --num_workers 4 --cpu 8 -o out
```

### Pfam Domain Annotation

```bash
# Realign queries to Pfam domains found via orthology transfer
python emapper.py -i seqs.fa --pfam_realign realign -o out

# De novo Pfam search against the entire Pfam database
python emapper.py -i seqs.fa --pfam_realign denovo -o out
```

### Search Filtering

```bash
# Filter by identity and coverage
python emapper.py -i seqs.fa --pident 40 --query_cover 60 --subject_cover 60 -o out

# Stricter e-value and score thresholds
python emapper.py -i seqs.fa --evalue 1e-10 --score 100 -o out

# Adjust seed ortholog thresholds for annotation
python emapper.py -i seqs.fa --seed_ortholog_evalue 1e-5 --seed_ortholog_score 60 -o out
```

### Output Options

```bash
# Also produce Excel output
python emapper.py -i seqs.fa -o out --excel

# Report orthologs used for annotation
python emapper.py -i seqs.fa -o out --report_orthologs

# Include MD5 hash of each query
python emapper.py -i seqs.fa -o out --md5

# No header/comment lines in output files
python emapper.py -i seqs.fa -o out --no_file_comments

# Use a scratch directory for faster I/O on network filesystems
python emapper.py -i seqs.fa -o out --scratch_dir /local/scratch --output_dir /nfs/results

# Decorate a GFF file with annotations
python emapper.py -i genome.fna --itype genome --decorate_gff yes -o out
```

### Resuming and Overriding

```bash
# Resume a previous run (skip already-computed results)
python emapper.py -i seqs.fa -o out --resume

# Force overwrite of existing output files
python emapper.py -i seqs.fa -o out --override
```

### Downloading Databases

```bash
# Download everything (annotation DB + DIAMOND + taxonomy)
python download_eggnog_data.py -y

# Skip DIAMOND database download
python download_eggnog_data.py -y -D

# Also install MMseqs2 database
python download_eggnog_data.py -y -M

# Also install Pfam database (required for --pfam_realign)
python download_eggnog_data.py -y -P

# Install novel families databases
python download_eggnog_data.py -y -F

# Install HMMER database for Bacteria
python download_eggnog_data.py -y -H -d 2 --dbname Bacteria

# Download to a custom directory
python download_eggnog_data.py -y --data_dir /path/to/data

# Simulate (show commands without downloading)
python download_eggnog_data.py -s
```

### Creating Custom Databases

Use `create_dbs.py` to build DIAMOND or MMseqs2 databases from a subset of
eggNOG proteins filtered by taxonomy:

```bash
# Create a DIAMOND database for Bacteria and Archaea
python create_dbs.py -y --dbname my_prok --taxa Bacteria,Archaea -m diamond

# Create an MMseqs2 database for specific tax IDs
python create_dbs.py -y --dbname my_custom --taxids 2,2157 -m mmseqs
```

---

## Reference: All Options

### General

| Option | Default | Description |
|--------|---------|-------------|
| `-v`, `--version` | | Show version and exit |
| `--list_taxa` | | List available taxonomic scopes and exit |

### Execution

| Option | Default | Description |
|--------|---------|-------------|
| `--cpu NUM` | `1` | Number of CPUs (0 = all available) |
| `--mp_start_method` | `spawn` | Python multiprocessing start method (`fork`, `spawn`, `forkserver`) |
| `--resume` | off | Resume a previous run, skipping existing results |
| `--override` | off | Overwrite existing output files |

### Input Data

| Option | Default | Description |
|--------|---------|-------------|
| `-i FILE` | | Input FASTA file (required unless `-m no_search`) |
| `--itype` | `proteins` | Input type: `proteins`, `CDS`, `genome`, `metagenome` |
| `--translate` | off | Translate CDS to proteins before search |
| `--annotate_hits_table FILE` | | Annotate an existing seed_orthologs file (requires `-m no_search`) |
| `-c`, `--cache FILE` | | Cache file for `-m cache` mode |
| `--db BACKEND` | `eggnog5` | Database backend (`eggnog5`, `e5-test`). Overridden by `--data_dir` |
| `--data_dir DIR` | | Custom path to database directory. Overrides `--db` and `EGGNOG_DATA_DIR` |

### Gene Prediction

| Option | Default | Description |
|--------|---------|-------------|
| `--genepred` | `search` | Gene prediction method: `search` (blastx-based) or `prodigal` |
| `--trans_table CODE` | program default | Genetic code table for Prodigal / DIAMOND / MMseqs2 |
| `--training_genome FILE` | | Genome for Prodigal training mode |
| `--training_file FILE` | | Pre-existing Prodigal training file |
| `--allow_overlaps` | `none` | Overlap handling: `none`, `strand`, `diff_frame`, `all` |
| `--overlap_tol FLOAT` | `0.0` | Overlap tolerance (0.0–1.0) |

### Search Mode

| Option | Default | Description |
|--------|---------|-------------|
| `-m MODE` | `diamond` | Search mode: `diamond`, `mmseqs`, `hmmer`, `no_search`, `cache`, `novel_fams` |

### Search Filtering (DIAMOND / MMseqs2)

| Option | Default | Description |
|--------|---------|-------------|
| `--pident FLOAT` | none | Min percent identity (0–100) |
| `--query_cover FLOAT` | none | Min query coverage (0–100) |
| `--subject_cover FLOAT` | none | Min subject coverage (0–100) |
| `--evalue FLOAT` | `0.001` | Max e-value threshold |
| `--score FLOAT` | none | Min bit score threshold |

### DIAMOND Options

| Option | Default | Description |
|--------|---------|-------------|
| `--dmnd_db FILE` | auto | Path to custom DIAMOND database |
| `--sensmode` | `sensitive` | Sensitivity: `fast`, `default`, `mid-sensitive`, `sensitive`, `more-sensitive`, `very-sensitive`, `ultra-sensitive` |
| `--dmnd_iterate` | `no` | Iterative search (`yes` / `no`) |
| `--dmnd_algo` | `auto` | Algorithm: `auto`, `0`, `1`, `ctg` |
| `--matrix` | none | Scoring matrix: `BLOSUM62`, `BLOSUM90`, `BLOSUM80`, `BLOSUM50`, `BLOSUM45`, `PAM250`, `PAM70`, `PAM30` |
| `--dmnd_frameshift INT` | none | Frameshift penalty (recommended: 15) |
| `--gapopen INT` | none | Gap open penalty |
| `--gapextend INT` | none | Gap extend penalty |
| `--block_size FLOAT` | diamond default | DIAMOND `-b` block size |
| `--index_chunks INT` | diamond default | DIAMOND `-c` index chunks |
| `--outfmt_short` | off | Minimal DIAMOND output (qseqid, sseqid, evalue, score only) |
| `--dmnd_ignore_warnings` | off | Ignore DIAMOND warnings |

### MMseqs2 Options

| Option | Default | Description |
|--------|---------|-------------|
| `--mmseqs_db FILE` | auto | Path to custom MMseqs2 database |
| `--start_sens FLOAT` | `3` | Starting sensitivity |
| `--sens_steps INT` | `3` | Number of sensitivity steps |
| `--final_sens FLOAT` | `7` | Final sensitivity |
| `--mmseqs_sub_mat` | mmseqs2 default | Substitution matrix |

### HMMER Options

| Option | Default | Description |
|--------|---------|-------------|
| `-d`, `--database DB` | | HMMER database: `euk`, `bact`, `arch`, or `db.hmm:host:port` |
| `--servers_list FILE` | | File with `host:port` entries for remote servers |
| `--qtype` | `seq` | Query type: `seq` or `hmm` |
| `--dbtype` | `hmmdb` | Database type: `hmmdb` or `seqdb` |
| `--usemem` | off | Load HMMER database into memory |
| `-p`, `--port PORT` | `51700` | Starting port for HMM server |
| `--end_port PORT` | `53200` | Last port for HMM server |
| `--num_servers INT` | `1` | Number of hmmpgmd servers |
| `--num_workers INT` | `1` | Workers per server |
| `--timeout_load_server INT` | `10` | Attempts to load a server on a port before trying the next |
| `--hmm_maxhits INT` | `1` | Max hits to report (0 = all) |
| `--report_no_hits` | off | Include queries without hits in output |
| `--hmm_maxseqlen INT` | `5000` | Skip queries longer than this |
| `--Z FLOAT` | `40000000` | Fixed database size for e-value calculation |
| `--cut_ga` | off | Use Pfam gathering thresholds |
| `--clean_overlaps` | none | Remove overlapping hits: `none`, `all`, `clans`, `hmmsearch_all`, `hmmsearch_clans` |

### Annotation Options

| Option | Default | Description |
|--------|---------|-------------|
| `--no_annot` | off | Skip annotation, report search hits only |
| `--dbmem` | off | Load `eggnog.db` into memory (~45 GB RAM) |
| `--seed_ortholog_evalue FLOAT` | `0.001` | Min e-value for seed ortholog acceptance |
| `--seed_ortholog_score FLOAT` | none | Min bit score for seed ortholog acceptance |
| `--tax_scope` | `auto` | Taxonomic scope: `auto`, `bacteria`, `eukaryota`, `archaea`, tax IDs, file path, or `none` |
| `--tax_scope_mode` | `inner_narrowest` | Scope mode: `broadest`, `inner_broadest`, `inner_narrowest`, `narrowest` |
| `--target_orthologs` | `all` | Ortholog types: `one2one`, `many2one`, `one2many`, `many2many`, `all` |
| `--target_taxa LIST` | none | Restrict to specific taxa (comma-separated tax IDs) |
| `--excluded_taxa LIST` | none | Exclude specific taxa (comma-separated tax IDs) |
| `--report_orthologs` | off | Output orthologs to a `.orthologs` file |
| `--go_evidence` | `non-electronic` | GO evidence filter: `experimental`, `non-electronic`, `all` |
| `--pfam_realign` | `none` | Pfam realignment: `none`, `realign`, `denovo` |
| `--md5` | off | Add MD5 hash of each query to output |

### Output Options

| Option | Default | Description |
|--------|---------|-------------|
| `-o`, `--output PREFIX` | | Base name for output files |
| `--output_dir DIR` | current dir | Output directory |
| `--scratch_dir DIR` | | Temporary scratch directory (for network filesystems) |
| `--temp_dir DIR` | current dir | Directory for temporary files |
| `--no_file_comments` | off | Omit header and stats from output files |
| `--decorate_gff` | `no` | Decorate GFF: `no`, `yes` (from gene prediction), or path to GFF file |
| `--decorate_gff_ID_field` | `ID` | GFF field used as feature ID |
| `--excel` | off | Also output annotations in `.xlsx` format |

### download_eggnog_data.py Options

| Option | Default | Description |
|--------|---------|-------------|
| `-D` | off | Skip DIAMOND database download |
| `-F` | off | Install novel families databases |
| `-P` | off | Install Pfam database |
| `-M` | off | Install MMseqs2 database |
| `-H` | off | Install HMMER database (requires `-d`) |
| `-d TAXID` | | Tax ID for HMMER database download |
| `--dbname NAME` | | Name for the HMMER database directory |
| `--db BACKEND` | `eggnog5` | Database backend |
| `--data_dir DIR` | | Custom data directory |
| `-y` | off | Assume "yes" to all prompts |
| `-f` | off | Force re-download |
| `-s` | off | Simulate (print commands, don't download) |
| `-q` | off | Quiet mode |

### create_dbs.py Options

| Option | Default | Description |
|--------|---------|-------------|
| `-m MODE` | `diamond` | Database type to create: `diamond` or `mmseqs` |
| `--dbname NAME` | | Prefix name for the output database (required) |
| `--taxids LIST` | | Comma-separated tax IDs to include |
| `--taxa LIST` | | Comma-separated taxon names to include |
| `-x` | off | Skip MMseqs2 index creation |
| `--db BACKEND` | `eggnog5` | Database backend |
| `--data_dir DIR` | | Custom data directory |
| `-y` | off | Assume "yes" to all prompts |
| `-s` | off | Simulate mode |

---

## Annotation Output Columns

The `*.emapper.annotations` file is a TSV with the following columns:

| # | Column | Description |
|---|--------|-------------|
| 1 | `query` | Query sequence name |
| 2 | `seed_ortholog` | Best matching eggNOG protein |
| 3 | `evalue` | E-value of the seed ortholog hit |
| 4 | `score` | Bit score of the seed ortholog hit |
| 5 | `eggNOG_OGs` | Orthologous groups at different taxonomic levels |
| 6 | `max_annot_lvl` | Broadest taxonomic level used for annotation |
| 7 | `COG_category` | COG functional category letter(s) |
| 8 | `Description` | Functional description from the best OG |
| 9 | `Preferred_name` | Short gene name |
| 10 | `GOs` | Gene Ontology terms (comma-separated) |
| 11 | `EC` | Enzyme Commission numbers |
| 12 | `KEGG_ko` | KEGG orthology identifiers |
| 13 | `KEGG_Pathway` | KEGG pathway maps |
| 14 | `KEGG_Module` | KEGG modules |
| 15 | `KEGG_Reaction` | KEGG reactions |
| 16 | `KEGG_rclass` | KEGG reaction classes |
| 17 | `BRITE` | KEGG BRITE hierarchies |
| 18 | `KEGG_TC` | Transporter Classification |
| 19 | `CAZy` | Carbohydrate-Active Enzymes |
| 20 | `BiGG_Reaction` | BiGG metabolic reactions |
| 21 | `PFAMs` | Pfam domain identifiers |
