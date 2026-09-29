#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd

BASE_DIR = Path(".")
BLAST_DIR = BASE_DIR / "07_blast_results"

FILES = {
    "strong": BLAST_DIR / "strong_vs_refs.tsv",
    "moderate": BLAST_DIR / "moderate_vs_refs.tsv",
}

COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

for label, file in FILES.items():
    if not file.exists() or file.stat().st_size == 0:
        print(f"{label}: no BLAST hits file found or empty: {file}")
        continue

    df = pd.read_csv(file, sep="\t", header=None, names=COLS)

    if df.empty:
        print(f"{label}: no hits")
        continue

    # best hit per query by bitscore then evalue
    best = (
        df.sort_values(["qseqid", "bitscore", "evalue"], ascending=[True, False, True])
          .drop_duplicates(subset=["qseqid"], keep="first")
          .copy()
    )

    out_file = BLAST_DIR / f"{label}_best_hits.tsv"
    best.to_csv(out_file, sep="\t", index=False)

    print(f"{label}:")
    print(f"  total hits: {len(df)}")
    print(f"  unique query sequences with hits: {best['qseqid'].nunique()}")
    print(f"  best-hit file: {out_file}")
    print(f"  subject distribution:")
    print(best["sseqid"].value_counts().to_string())
    print()
