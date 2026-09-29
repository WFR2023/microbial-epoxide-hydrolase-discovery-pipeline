#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd
from Bio import SeqIO


BASE_DIR = Path("<PROJECT_ROOT>")
HITS_TABLE = BASE_DIR / "06_candidate_mining" / "candidate_hits.tsv"
HITS_FASTA = BASE_DIR / "06_candidate_mining" / "candidate_hits.faa"

OUT_DIR = BASE_DIR / "06_candidate_mining" / "prioritized"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_TABLE = OUT_DIR / "candidate_hits_prioritized.tsv"
OUT_FASTA_STRONG = OUT_DIR / "candidate_hits_strong.faa"
OUT_FASTA_MODERATE = OUT_DIR / "candidate_hits_moderate.faa"
OUT_SUMMARY = OUT_DIR / "candidate_hits_prioritized_summary.tsv"


STRONG_TERMS = [
    "epoxide hydrolase",
    "epoxide hydrolase-like",
    "epoxide-hydrolase",
    "cif",
    "cif-like",
]

MODERATE_TERMS = [
    "alpha/beta hydrolase",
    "alpha beta hydrolase",
    "ab hydrolase",
    "lipid hydrolase",
    "xenobiotic hydrolase",
    "haloalkane",
]

WEAK_TERMS = [
    "hydrolase family protein",
]


def normalize_text(text):
    return str(text).lower().strip()


def classify_description(desc):
    desc_n = normalize_text(desc)

    strong_hits = [t for t in STRONG_TERMS if t in desc_n]
    moderate_hits = [t for t in MODERATE_TERMS if t in desc_n]
    weak_hits = [t for t in WEAK_TERMS if t in desc_n]

    if strong_hits:
        return "strong", "; ".join(strong_hits)
    if moderate_hits:
        return "moderate", "; ".join(moderate_hits)
    if weak_hits:
        return "weak", "; ".join(weak_hits)
    return "other", ""


def main():
    if not HITS_TABLE.exists():
        raise FileNotFoundError(f"Missing hits table: {HITS_TABLE}")
    if not HITS_FASTA.exists():
        raise FileNotFoundError(f"Missing hits FASTA: {HITS_FASTA}")

    hits = pd.read_csv(HITS_TABLE, sep="\t")
    hits["priority_class"] = ""
    hits["priority_match"] = ""

    for idx, row in hits.iterrows():
        priority, matched = classify_description(row["description"])
        hits.at[idx, "priority_class"] = priority
        hits.at[idx, "priority_match"] = matched

    hits.to_csv(OUT_TABLE, sep="\t", index=False)

    strong_ids = set(
        hits.loc[hits["priority_class"] == "strong", "assembly_accession"]
        + "|" +
        hits.loc[hits["priority_class"] == "strong", "protein_id"]
    )

    moderate_ids = set(
        hits.loc[hits["priority_class"] == "moderate", "assembly_accession"]
        + "|" +
        hits.loc[hits["priority_class"] == "moderate", "protein_id"]
    )

    strong_records = []
    moderate_records = []

    for record in SeqIO.parse(HITS_FASTA, "fasta"):
        rid = record.id
        if rid in strong_ids:
            strong_records.append(record)
        elif rid in moderate_ids:
            moderate_records.append(record)

    SeqIO.write(strong_records, OUT_FASTA_STRONG, "fasta")
    SeqIO.write(moderate_records, OUT_FASTA_MODERATE, "fasta")

    summary = (
        hits.groupby(["target_group", "target_species", "priority_class"])
        .agg(
            n_hits=("protein_id", "count"),
            n_assemblies=("assembly_accession", "nunique"),
        )
        .reset_index()
        .sort_values(["target_group", "target_species", "priority_class"])
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print(f"Prioritized table written to: {OUT_TABLE}")
    print(f"Strong FASTA written to: {OUT_FASTA_STRONG}")
    print(f"Moderate FASTA written to: {OUT_FASTA_MODERATE}")
    print(f"Summary written to: {OUT_SUMMARY}")
    print("\nPriority counts:")
    print(hits["priority_class"].value_counts().to_string())


if __name__ == "__main__":
    main()
