#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd

BASE_DIR = Path("<PROJECT_ROOT>")
BLAST_DIR = BASE_DIR / "07_blast_results"
OUT_DIR = BASE_DIR / "07_blast_results" / "filtered"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STRONG_FILE = BLAST_DIR / "strong_best_hits.tsv"
MODERATE_FILE = BLAST_DIR / "moderate_best_hits.tsv"

OUT_STRONG = OUT_DIR / "strong_best_hits_filtered.tsv"
OUT_MODERATE = OUT_DIR / "moderate_best_hits_filtered.tsv"
OUT_COMBINED = OUT_DIR / "priority_candidates_final.tsv"


def main():
    strong = pd.read_csv(STRONG_FILE, sep="\t")
    moderate = pd.read_csv(MODERATE_FILE, sep="\t")

    # Strong: keep all, since they are already very restricted
    strong_f = strong.copy()
    strong_f["priority_tier"] = "strong_blast_supported"

    # Moderate: apply more stringent biological filters
    moderate_f = moderate[
        (moderate["evalue"] <= 1e-8) &
        (moderate["bitscore"] >= 40) &
        (moderate["length"] >= 100)
    ].copy()
    moderate_f["priority_tier"] = "moderate_blast_supported"

    strong_f.to_csv(OUT_STRONG, sep="\t", index=False)
    moderate_f.to_csv(OUT_MODERATE, sep="\t", index=False)

    combined = pd.concat([strong_f, moderate_f], ignore_index=True)
    combined = combined.sort_values(
        by=["priority_tier", "bitscore", "evalue"],
        ascending=[True, False, True]
    )
    combined.to_csv(OUT_COMBINED, sep="\t", index=False)

    print("Filtered strong:", len(strong_f))
    print("Filtered moderate:", len(moderate_f))
    print("Combined priority candidates:", len(combined))
    print(f"\nFiles written:")
    print(f"  {OUT_STRONG}")
    print(f"  {OUT_MODERATE}")
    print(f"  {OUT_COMBINED}")


if __name__ == "__main__":
    main()
