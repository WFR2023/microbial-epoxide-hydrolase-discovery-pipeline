#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd
from Bio import SeqIO


BASE_DIR = Path("<PROJECT_ROOT>")
PROTEIN_DIR = BASE_DIR / "04_derived" / "proteins_final"
SELECTED_FILE = BASE_DIR / "01_metadata" / "selected_assemblies.tsv"

OUT_DIR = BASE_DIR / "06_candidate_mining"
OUT_DIR.mkdir(parents=True, exist_ok=True)

HITS_TABLE = OUT_DIR / "candidate_hits.tsv"
HITS_FASTA = OUT_DIR / "candidate_hits.faa"
SUMMARY_FILE = OUT_DIR / "candidate_hits_summary.tsv"

KEYWORDS = [
    "epoxide hydrolase",
    "epoxide hydrolase-like",
    "epoxide-hydrolase",
    "alpha/beta hydrolase",
    "alpha beta hydrolase",
    "ab hydrolase",
    "hydrolase family protein",
    "cif",
    "cif-like",
    "haloalkane",
    "lipid hydrolase",
    "xenobiotic hydrolase",
]


def safe_text(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def normalize_text(text):
    return safe_text(text).lower()


def main():
    if not PROTEIN_DIR.exists():
        raise FileNotFoundError(f"Missing protein directory: {PROTEIN_DIR}")
    if not SELECTED_FILE.exists():
        raise FileNotFoundError(f"Missing selected assemblies file: {SELECTED_FILE}")

    selected = pd.read_csv(SELECTED_FILE, sep="\t")
    meta = selected.set_index("assembly_accession").to_dict(orient="index")

    hit_rows = []
    hit_records = []

    protein_files = sorted(PROTEIN_DIR.glob("*.faa"))

    for faa_file in protein_files:
        accession = faa_file.stem
        m = meta.get(accession, {})

        target_group = safe_text(m.get("target_group", ""))
        target_species = safe_text(m.get("target_species", ""))
        source_database = safe_text(m.get("source_database", ""))

        for record in SeqIO.parse(faa_file, "fasta"):
            desc = normalize_text(record.description)

            matched_keywords = [kw for kw in KEYWORDS if kw in desc]

            if matched_keywords:
                hit_rows.append({
                    "target_group": target_group,
                    "target_species": target_species,
                    "assembly_accession": accession,
                    "source_database": source_database,
                    "protein_id": record.id,
                    "description": record.description,
                    "sequence_length_aa": len(record.seq),
                    "matched_keywords": "; ".join(matched_keywords),
                    "protein_file": str(faa_file),
                })

                record.id = f"{accession}|{record.id}"
                record.description = f"{target_species} | {record.description}"
                hit_records.append(record)

    hits_df = pd.DataFrame(hit_rows)

    if hits_df.empty:
        hits_df = pd.DataFrame(columns=[
            "target_group",
            "target_species",
            "assembly_accession",
            "source_database",
            "protein_id",
            "description",
            "sequence_length_aa",
            "matched_keywords",
            "protein_file",
        ])
        hits_df.to_csv(HITS_TABLE, sep="\t", index=False)

        with open(HITS_FASTA, "w") as f:
            pass

        summary = pd.DataFrame(columns=[
            "target_group",
            "target_species",
            "n_candidate_hits",
            "n_unique_assemblies_with_hits",
        ])
        summary.to_csv(SUMMARY_FILE, sep="\t", index=False)

        print("No textual candidate hits found.")
        print(f"Empty table written to: {HITS_TABLE}")
        print(f"Empty FASTA written to: {HITS_FASTA}")
        print(f"Empty summary written to: {SUMMARY_FILE}")
        return

    hits_df.to_csv(HITS_TABLE, sep="\t", index=False)
    SeqIO.write(hit_records, HITS_FASTA, "fasta")

    summary = (
        hits_df.groupby(["target_group", "target_species"])
        .agg(
            n_candidate_hits=("protein_id", "count"),
            n_unique_assemblies_with_hits=("assembly_accession", "nunique"),
        )
        .reset_index()
        .sort_values(["target_group", "target_species"])
    )

    summary.to_csv(SUMMARY_FILE, sep="\t", index=False)

    print(f"Candidate hits written to: {HITS_TABLE}")
    print(f"Candidate FASTA written to: {HITS_FASTA}")
    print(f"Summary written to: {SUMMARY_FILE}")
    print(f"Total candidate proteins found: {len(hits_df)}")


if __name__ == "__main__":
    main()
