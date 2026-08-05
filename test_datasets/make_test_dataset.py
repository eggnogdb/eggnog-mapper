#!/usr/bin/env python3
"""Build the eggNOG-mapper v3 self-test dataset.

Produces a small, self-consistent ``test_datasets/data/`` directory (a mini
eggNOG DB + DIAMOND db + taxonomy + OBO) together with per-``itype`` query
fixtures (proteins / CDS / genome / metagenome). Running emapper against this
data reproduces fixed golden annotations, so the whole option surface can be
regression-tested offline.

Design
------
* **Self-consistent subsample.** For each query seed we walk its annotation
  closure: ``event_index`` -> ``sp_events`` -> decoded ``side1``/``side2``
  (the ortholog donors the cascade transfers from). Those donor proteins get
  their ``prots``/``protein_names``/``ogs`` rows copied, so the mini DB yields
  *biologically meaningful* annotations rather than empty ones.
* **event_index only for query seeds**, so every event referenced by the mini
  ``event_index`` is present in the mini ``sp_events`` (no dangling events).
* **DIAMOND db = query seeds only**, so each fixture self-hits its seed
  deterministically and prodigal-predicted ORFs still land on the right seed.
* **Deterministic fixtures.** CDS are a fixed reverse-translation of each query
  protein; genome/metagenome contigs embed those CDS in deterministic flanks.
  Fully offline and reproducible -- no network, no external nucleotide source.

This script is *not* shipped in the repo (``test_datasets/`` is gitignored);
it documents exactly how the downloadable self-test data was generated.
"""

from __future__ import annotations

import argparse
import array
import gzip
import os
import shutil
import sqlite3
import subprocess
import sys
from os.path import join as pjoin

# --- repo import path (for the delta-varint codec) -------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from eggnogmapper.annotator.codec import decode_intlist  # noqa: E402

# --- source data (the final, latest full DB) -------------------------------
FULL_DB = "/eggnog-eco/data/e7/full/final/mapper/eggnog.db"
PROTEOME_FA = "/eggnog-eco/data/e7/full/source/proteomes/proteins_intids.fa.gz"
TAXA_DB = "/eggnog-eco/data/e7/full/final/mapper/eggnog.taxa.db"
TAXA_TRAVERSE = TAXA_DB + ".traverse.pkl"
OBO = "/eggnog-eco/data/e7/full/source/reference/go-basic.obo"
DIAMOND = "/eggnog-eco/bin/diamond"

# --- query seeds: a genuine nif operon spanning Archaea + Bacteria ---------
# (id, note) -- chosen for annotation richness and taxonomic spread.
QUERY_SEEDS = [
    (605553, "nifH  Methanothermobacter marburgensis (Archaea) - OG Fer4_NifH"),
    (604694, "nifD  Methanothermobacter marburgensis (Archaea)"),
    (604695, "nifK  Methanothermobacter marburgensis (Archaea)"),
    (105893, "nifD  Paenibacillus sabinae (Bacteria)"),
    (17828,  "nifJ  Sulfurihydrogenibium / Chlorobi (Bacteria)"),
    (17704775, "nifJ  Escherichia coli IAI39 (Bacteria) - BiGG-annotated"),
]

# Standard genetic code, one deterministic codon per amino acid.
AA2CODON = {
    "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT", "Q": "CAA",
    "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT", "L": "CTT", "K": "AAA",
    "M": "ATG", "F": "TTT", "P": "CCT", "S": "TCT", "T": "ACT", "W": "TGG",
    "Y": "TAT", "V": "GTT", "U": "TGT", "B": "AAT", "Z": "CAA", "J": "CTT",
    "X": "NNN", "*": "",
}
_COMP = str.maketrans("ACGTN", "TGCAN")


def reverse_translate(prot: str) -> str:
    """Reverse-translate a protein to a deterministic CDS (no trailing stop)."""
    return "".join(AA2CODON.get(aa.upper(), "NNN") for aa in prot if aa != "*")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def det_flank(seed: int, n: int) -> str:
    """Deterministic pseudo-random nucleotide flank (LCG seeded by ``seed``)."""
    bases = "ACGT"
    x = (seed * 2654435761 + 12345) & 0xFFFFFFFF
    out = []
    for _ in range(n):
        x = (1103515245 * x + 12345) & 0x7FFFFFFF
        out.append(bases[x & 3])
    return "".join(out)


def collect_closure(src: sqlite3.Connection, seed_ids):
    """Return (event_ids, all_protein_ids, og_names) for the seeds' closure."""
    event_ids: set[int] = set()
    prot_ids: set[int] = set(seed_ids)
    for pid in seed_ids:
        row = src.execute(
            "SELECT events FROM event_index WHERE protein_id=?", (pid,)
        ).fetchone()
        if row is None:
            print(f"  WARNING: seed {pid} has no event_index row", file=sys.stderr)
            continue
        event_ids.update(decode_intlist(row[0]))

    for ev in event_ids:
        row = src.execute(
            "SELECT side1, side2 FROM sp_events WHERE i=?", (ev,)
        ).fetchone()
        if row is None:
            continue
        prot_ids.update(decode_intlist(row[0]))
        prot_ids.update(decode_intlist(row[1]))

    og_names: set[str] = set()
    # OGs referenced by the copied prots rows...
    q = ",".join("?" * len(prot_ids))
    for (ogs,) in src.execute(f"SELECT ogs FROM prots WHERE id IN ({q})", tuple(prot_ids)):
        if ogs:
            og_names.update(part.split("@")[0] + "@" + part.split("@")[1]
                            for part in ogs.split(",") if "@" in part)
    # ...and by the copied sp_events rows.
    if event_ids:
        qe = ",".join("?" * len(event_ids))
        for (og,) in src.execute(
            f"SELECT og FROM sp_events WHERE i IN ({qe})", tuple(event_ids)
        ):
            if og:
                og_names.add(og)
    return sorted(event_ids), sorted(prot_ids), og_names


def build_mini_db(src, dst_path, seed_ids, event_ids, prot_ids, og_names):
    if os.path.exists(dst_path):
        os.remove(dst_path)
    dst = sqlite3.connect(dst_path)
    # replicate every table's CREATE statement verbatim
    for (sql,) in src.execute(
        "SELECT sql FROM sqlite_master WHERE type='table' AND sql IS NOT NULL"
    ).fetchall():
        dst.execute(sql)

    def chunks(seq, n=900):
        seq = list(seq)
        for i in range(0, len(seq), n):
            yield seq[i:i + n]

    # prots + protein_names for the whole closure (annotation donors)
    for grp in chunks(prot_ids):
        q = ",".join("?" * len(grp))
        for row in src.execute(f"SELECT * FROM prots WHERE id IN ({q})", grp):
            dst.execute(f"INSERT INTO prots VALUES ({','.join('?'*len(row))})", row)
        for row in src.execute(f"SELECT * FROM protein_names WHERE id IN ({q})", grp):
            dst.execute(
                f"INSERT INTO protein_names VALUES ({','.join('?'*len(row))})", row
            )
    # event_index only for query seeds
    for grp in chunks(seed_ids):
        q = ",".join("?" * len(grp))
        for row in src.execute(
            f"SELECT * FROM event_index WHERE protein_id IN ({q})", grp
        ):
            dst.execute(
                f"INSERT INTO event_index VALUES ({','.join('?'*len(row))})", row
            )
    # sp_events for the seeds' events
    for grp in chunks(event_ids):
        q = ",".join("?" * len(grp))
        for row in src.execute(f"SELECT * FROM sp_events WHERE i IN ({q})", grp):
            dst.execute(
                f"INSERT INTO sp_events VALUES ({','.join('?'*len(row))})", row
            )
    # ogs referenced anywhere
    for grp in chunks(sorted(og_names)):
        q = ",".join("?" * len(grp))
        for row in src.execute(f"SELECT * FROM ogs WHERE name IN ({q})", grp):
            dst.execute(f"INSERT INTO ogs VALUES ({','.join('?'*len(row))})", row)
    # version
    for row in src.execute("SELECT * FROM version"):
        dst.execute(f"INSERT INTO version VALUES ({','.join('?'*len(row))})", row)

    dst.commit()
    n = {t: dst.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
         for t in ("prots", "protein_names", "event_index", "sp_events", "ogs")}
    dst.execute("VACUUM")
    dst.commit()
    dst.close()
    return n


def extract_sequences(seed_ids):
    """Scan proteins_intids.fa.gz once, returning {id: sequence} for seeds."""
    want = {str(i) for i in seed_ids}
    seqs: dict[int, str] = {}
    cur_id = None
    buf: list[str] = []
    with gzip.open(PROTEOME_FA, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(buf)
                hid = line[1:].strip().split()[0]
                cur_id = int(hid) if hid in want else None
                buf = []
                if len(seqs) == len(want):
                    break
            elif cur_id is not None:
                buf.append(line.strip())
        if cur_id is not None and cur_id not in seqs:
            seqs[cur_id] = "".join(buf)
    missing = set(seed_ids) - set(seqs)
    if missing:
        raise SystemExit(f"Sequences not found for seeds: {sorted(missing)}")
    return seqs


def label_for(seed_id):
    for sid, note in QUERY_SEEDS:
        if sid == seed_id:
            return note.split()[0]  # gene name, e.g. "nifH"
    return f"seed{seed_id}"


def write_fixtures(fx_dir, seed_ids, seqs, names):
    os.makedirs(fx_dir, exist_ok=True)
    # header id = "<intid>_<gene>" so golden outputs are human-readable & stable
    hid = {s: f"{s}_{label_for(s)}" for s in seed_ids}

    # 1) proteins
    with open(pjoin(fx_dir, "proteins.faa"), "w") as fh:
        for s in seed_ids:
            fh.write(f">{hid[s]}\n{seqs[s]}\n")

    # 2) CDS (exact reverse-translation -> auto-translate path)
    cds = {s: reverse_translate(seqs[s]) for s in seed_ids}
    with open(pjoin(fx_dir, "cds.fna"), "w") as fh:
        for s in seed_ids:
            fh.write(f">{hid[s]}\n{cds[s]}\n")

    # 3) genome: one contig per gene = flank + CDS + stop + flank (prodigal)
    with open(pjoin(fx_dir, "genome.fna"), "w") as fh:
        for s in seed_ids:
            contig = det_flank(s, 90) + cds[s] + "TAA" + det_flank(s + 7, 90)
            fh.write(f">contig_{hid[s]}\n{contig}\n")

    # 4) metagenome: few contigs, several ORFs each, mixed strands
    with open(pjoin(fx_dir, "metagenome.fna"), "w") as fh:
        per = 3
        for ci in range(0, len(seed_ids), per):
            group = seed_ids[ci:ci + per]
            parts = [det_flank(1000 + ci, 120)]
            for j, s in enumerate(group):
                orf = cds[s] + "TAA"
                parts.append(revcomp(orf) if j % 2 else orf)
                parts.append(det_flank(2000 + s, 120))
            fh.write(f">metacontig_{ci // per}\n{''.join(parts)}\n")

    return hid


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=pjoin(HERE, "data"),
                    help="output data dir (default: test_datasets/data)")
    ap.add_argument("--fixtures", default=pjoin(HERE, "fixtures"),
                    help="output fixtures dir (default: test_datasets/fixtures)")
    args = ap.parse_args()

    for p in (FULL_DB, PROTEOME_FA, TAXA_DB, OBO, DIAMOND):
        if not os.path.exists(p):
            raise SystemExit(f"Required source missing: {p}")

    os.makedirs(args.out, exist_ok=True)
    seed_ids = [s for s, _ in QUERY_SEEDS]

    print(f"[1/6] closure over {len(seed_ids)} query seeds ...")
    src = sqlite3.connect(FULL_DB)
    event_ids, prot_ids, og_names = collect_closure(src, seed_ids)
    print(f"      events={len(event_ids)} closure_prots={len(prot_ids)} ogs={len(og_names)}")

    print("[2/6] writing mini eggnog.db ...")
    counts = build_mini_db(src, pjoin(args.out, "eggnog.db"),
                           seed_ids, event_ids, prot_ids, og_names)
    src.close()
    print(f"      {counts}")

    print("[3/6] extracting query-seed sequences ...")
    seqs = extract_sequences(seed_ids)

    print("[4/6] writing fixtures ...")
    names = {}
    hid = write_fixtures(args.fixtures, seed_ids, seqs, names)

    print("[5/6] building DIAMOND db (query seeds only) ...")
    faa = pjoin(args.out, "_seeds.faa")
    with open(faa, "w") as fh:
        for s in seed_ids:
            fh.write(f">{s}\n{seqs[s]}\n")
    dmnd = pjoin(args.out, "eggnog_proteins.dmnd")
    subprocess.run([DIAMOND, "makedb", "--in", faa, "--db", dmnd],
                   check=True, capture_output=True)
    os.remove(faa)

    print("[6/6] copying taxonomy + OBO ...")
    shutil.copy2(TAXA_DB, pjoin(args.out, "eggnog.taxa.db"))
    shutil.copy2(TAXA_TRAVERSE, pjoin(args.out, "eggnog.taxa.db.traverse.pkl"))
    shutil.copy2(OBO, pjoin(args.out, "go-basic.obo"))

    print("\nDone. Data dir:", args.out)
    print("Query header ids:", hid)


if __name__ == "__main__":
    main()
