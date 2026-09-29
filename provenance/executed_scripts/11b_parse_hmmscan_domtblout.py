#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd

BASE_DIR = Path("<PROJECT_ROOT>")
DOMTBLOUT = BASE_DIR / "13_hmmer_pfam" / "priority_candidates_vs_pfam.domtblout"
OUT_DIR = BASE_DIR / "13_hmmer_pfam"
OUT_TABLE = OUT_DIR / "priority_candidates_vs_pfam_parsed.tsv"
OUT_BEST = OUT_DIR / "priority_candidates_vs_pfam_best_hits.tsv"
OUT_SUMMARY = OUT_DIR / "priority_candidates_vs_pfam_summary.tsv"


def parse_domtblout(file_path):
    rows = []
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 23:
                continue

            rows.append({
                "target_name": parts[0],
                "target_accession": parts[1],
                "target_len": parts[2],
                "query_name": parts[3],
                "query_accession": parts[4],
                "query_len": parts[5],
                "full_evalue": parts[6],
                "full_score": parts[7],
                "full_bias": parts[8],
                "domain_num": parts[9],
                "domain_of": parts[10],
                "c_evalue": parts[11],
                "i_evalue": parts[12],
                "domain_score": parts[13],
                "domain_bias": parts[14],
                "hmm_from": parts[15],
                "hmm_to": parts[16],
                "ali_from": parts[17],
                "ali_to": parts[18],
                "env_from": parts[19],
                "env_to": parts[20],
                "acc": parts[21],
                "description": " ".join(parts[22:]),
            })
    return pd.DataFrame(rows)


def main():
    if not DOMTBLOUT.exists():
        raise FileNotFoundError(f"Missing domtblout: {DOMTBLOUT}")

    df = parse_domtblout(DOMTBLOUT)

    if df.empty:
        print("No Pfam hits parsed.")
        df.to_csv(OUT_TABLE, sep="\t", index=False)
        return

    # numeric conversion
    for col in ["full_evalue", "full_score", "i_evalue", "domain_score", "acc"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df.to_csv(OUT_TABLE, sep="\t", index=False)

    best = (
        df.sort_values(["query_name", "i_evalue", "domain_score"], ascending=[True, True, False])
          .drop_duplicates(subset=["query_name"], keep="first")
          .copy()
    )
    best.to_csv(OUT_BEST, sep="\t", index=False)

    summary = (
        best.groupby("target_name")
        .size()
        .reset_index(name="n_queries_best_assigned")
        .sort_values("n_queries_best_assigned", ascending=False)
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print(f"Parsed table written to: {OUT_TABLE}")
    print(f"Best-hit table written to: {OUT_BEST}")
    print(f"Summary written to: {OUT_SUMMARY}")
    print(f"Queries with best Pfam assignment: {len(best)}")


if __name__ == "__main__":
    main()
