# eggNOG-mapper v3 — offline self-test dataset

A small, self-consistent slice of the **final** eggNOG v7 mapper DB
(`data/e7/full/final/mapper/eggnog.db`) plus per-`--itype` query fixtures.
Running emapper against it reproduces fixed golden annotations, so the full
option surface can be regression-tested **offline, with no external DB**.

The heavy artefacts (`data/`, `fixtures/`, `golden/`, and the `.tar.gz` bundle)
are **not committed** to the repo — they are downloadable. Only the builder,
the golden generator, the pytest runner, and this README are tracked.

## Contents (downloadable bundle)

```
test_datasets/
├── make_test_dataset.py   # builds data/ + fixtures/ from the full DB (tracked)
├── gen_golden.py          # regenerates golden/ (tracked)
├── test_selftest.py       # pytest: live emapper vs golden (tracked)
├── data/                  # mini DB + DIAMOND db + taxonomy + OBO   (download)
│   ├── eggnog.db                     # 6 nif query seeds + annotation closure
│   ├── eggnog_proteins.dmnd          # query seeds only (deterministic self-hit)
│   ├── eggnog.taxa.db(+.traverse.pkl)
│   └── go-basic.obo
├── fixtures/              # one query file per input type            (download)
│   ├── proteins.faa       # --itype proteins
│   ├── cds.fna            # --itype CDS   (exact reverse-translation)
│   ├── genome.fna         # --itype genome     (1 ORF / contig, prodigal)
│   └── metagenome.fna     # --itype metagenome (multi-ORF, mixed strands)
└── golden/                # expected .annotations per itype          (download)
```

## The test genes

A genuine **nif operon** spanning Archaea + Bacteria, chosen for annotation
richness (GO / KEGG KO / EC / pathway / PFAM / BiGG) and taxonomic spread:

| query id        | gene | organism (taxid)                         | expected OG          |
|-----------------|------|------------------------------------------|----------------------|
| `605553_nifH`   | nifH | *Methanothermobacter marburgensis* (Archaea) | `Fer4_NifH`      |
| `604694_nifD`   | nifD | *M. marburgensis* (Archaea)              | `Oxidored_nitro`     |
| `604695_nifK`   | nifK | *M. marburgensis* (Archaea)              | `Oxidored_nitro`     |
| `105893_nifD`   | nifD | *Paenibacillus sabinae* (Bacteria)       | `Oxidored_nitro`     |
| `17828_nifJ`    | nifJ | Bacteria                                 | `EKR` / `POR`        |
| `17704775_nifJ` | nifJ | *Escherichia coli* IAI39 (Bacteria)      | `EKR` / `POR` (BiGG) |

## Running the self-test

### The easy way — `emapper.py --selftest`

Downloads this dataset (if not already local) and verifies the installation
reproduces the reference annotations for every input type, then exits non-zero
on any mismatch. It runs *inside* emapper, so from the Apptainer image the
bundled DIAMOND/prodigal are already on `PATH` — nothing else to install:

```bash
# from the built image (".sif"):
apptainer run eggnog-mapper-3.0.0-beta5.sif --selftest

# or a normal install / source checkout (needs DIAMOND on PATH):
emapper.py --selftest

# skip the download and use an already-downloaded bundle:
emapper.py --selftest --data_dir /path/to/test_datasets
```

The bundle is fetched from
`https://data.cgmlab.org/eggnog-mapper/emapper-<release>/selftest/eggnog-mapper-selftest.tar.gz`
(override with `$EGGNOG_SELFTEST_URL`).

### For maintainers — pytest

Runs the same four itype checks plus two nifH biological assertions against a
local `data/ fixtures/ golden/`, without downloading (DIAMOND must be on `PATH`;
the repo ships one under `/eggnog-eco/bin`):

```bash
pytest test_datasets/test_selftest.py -v
```

## Regenerating (maintainers)

```bash
# rebuild the mini DB + fixtures from the full final DB, then the goldens:
python test_datasets/make_test_dataset.py
python test_datasets/gen_golden.py
```

`make_test_dataset.py` documents the design: for each query seed it copies the
**annotation closure** (`event_index` → `sp_events` → decoded ortholog donors),
so the mini DB yields biologically meaningful annotations; `event_index` is
kept only for the query seeds (no dangling events), and the DIAMOND db holds
only the seeds so every fixture self-hits deterministically. Fixtures are a
fixed reverse-translation of each protein, so CDS/genome/metagenome inputs are
fully offline and reproducible.
