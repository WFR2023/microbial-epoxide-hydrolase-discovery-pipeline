#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd

BASE_DIR = Path(".")

DESC_FILE = BASE_DIR / "09_candidate_description" / "priority_candidates_descriptive_table.tsv"
MOTIF_FILE = BASE_DIR / "12_domain_motif_screen" / "priority_candidates_motif_scan.tsv"
OUT_DIR = BASE_DIR / "12_domain_motif_screen"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_FILE = OUT_DIR / "priority_candidates_integrated_evidence.tsv"

def main():
    desc = pd.read_csv(DESC_FILE, sep="\t")
    motif = pd.read_csv(MOTIF_FILE, sep="\t")

    merged = desc.merge(
        motif,
        left_on="qseqid",
        right_on="sequence_id",
        how="left"
    )

    def assign_priority(row):
        ref = str(row["best_reference_hit"])
        tier = str(row["priority_tier"])
        hcls = str(row["heuristic_class"])
        pident = float(row["pident"])

        if tier == "strong_blast_supported" and "pseudomonas_Cif" in ref:
            return "high_priority_cif_like"
        if tier == "strong_blast_supported" and "human_EPHX1" in ref:
            return "high_priority_ephx1_like"
        if tier == "moderate_blast_supported" and "pseudomonas_Cif" in ref and "hydrolase" in hcls:
            return "moderate_priority_cif_related"
        if tier == "moderate_blast_supported" and "human_EPHX2" in ref:
            return "moderate_priority_ephx2_related"
        return "exploratory"

    merged["integrated_priority"] = merged.apply(assign_priority, axis=1)
    merged.to_csv(OUT_FILE, sep="\t", index=False)

    print(f"Integrated evidence table written to: {OUT_FILE}")
    print("\nIntegrated priority counts:")
    print(merged["integrated_priority"].value_counts().to_string())

if __name__ == "__main__":
    main()
