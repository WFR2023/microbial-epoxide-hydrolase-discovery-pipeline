#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

BASE_DIR = Path(".")

IN_TABLE = BASE_DIR / "08_priority_candidates" / "priority_candidates_final_annotated.tsv"
IN_FASTA = BASE_DIR / "08_priority_candidates" / "priority_candidates_final.faa"
OUT_DIR = BASE_DIR / "09_candidate_description"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_TABLE = OUT_DIR / "priority_candidates_descriptive_table.tsv"
OUT_SUMMARY = OUT_DIR / "priority_candidates_descriptive_summary.tsv"


def safe_protein_metrics(seq):
    try:
        pa = ProteinAnalysis(str(seq))
        return {
            "molecular_weight": round(pa.molecular_weight(), 3),
            "isoelectric_point": round(pa.isoelectric_point(), 3),
            "aromaticity": round(pa.aromaticity(), 5),
            "instability_index": round(pa.instability_index(), 3),
            "gravy": round(pa.gravy(), 5),
        }
    except Exception:
        return {
            "molecular_weight": "",
            "isoelectric_point": "",
            "aromaticity": "",
            "instability_index": "",
            "gravy": "",
        }


def main():
    if not IN_TABLE.exists():
        raise FileNotFoundError(f"Missing annotated table: {IN_TABLE}")
    if not IN_FASTA.exists():
        raise FileNotFoundError(f"Missing candidate FASTA: {IN_FASTA}")

    ann = pd.read_csv(IN_TABLE, sep="\t")
    records = {rec.id: rec for rec in SeqIO.parse(IN_FASTA, "fasta")}

    rows = []

    for _, row in ann.iterrows():
        qseqid = row["qseqid"]
        rec = records.get(qseqid)

        if rec is None:
            continue

        metrics = safe_protein_metrics(rec.seq)

        rows.append({
            "assembly_accession": row["assembly_accession"],
            "protein_id": row["protein_id"],
            "qseqid": qseqid,
            "best_reference_hit": row["sseqid"],
            "pident": row["pident"],
            "blast_aln_length": row["length"],
            "evalue": row["evalue"],
            "bitscore": row["bitscore"],
            "priority_tier": row["priority_tier"],
            "sequence_length_aa": len(rec.seq),
            **metrics,
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT_TABLE, sep="\t", index=False)

    summary = (
        df.groupby(["best_reference_hit", "priority_tier"])
        .agg(
            n_candidates=("qseqid", "count"),
            mean_length_aa=("sequence_length_aa", "mean"),
            mean_pident=("pident", "mean"),
            mean_bitscore=("bitscore", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print(f"Descriptive table written to: {OUT_TABLE}")
    print(f"Summary written to: {OUT_SUMMARY}")
    print(f"Candidates described: {len(df)}")


if __name__ == "__main__":
    main()
