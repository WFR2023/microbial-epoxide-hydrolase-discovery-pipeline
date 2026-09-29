#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd
from Bio import SeqIO

BASE_DIR = Path(".")

PRIORITY_FILE = BASE_DIR / "07_blast_results" / "filtered" / "priority_candidates_final.tsv"
PROTEIN_DIR = BASE_DIR / "04_derived" / "proteins_final"
OUT_DIR = BASE_DIR / "08_priority_candidates"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_FASTA = OUT_DIR / "priority_candidates_final.faa"
OUT_TABLE = OUT_DIR / "priority_candidates_final_annotated.tsv"


def main():
    if not PRIORITY_FILE.exists():
        raise FileNotFoundError(f"Missing priority file: {PRIORITY_FILE}")

    df = pd.read_csv(PRIORITY_FILE, sep="\t")

    # Split qseqid into accession and protein id
    df[["assembly_accession", "protein_id"]] = df["qseqid"].str.split("|", n=1, expand=True)

    selected_records = []
    annotated_rows = []

    for _, row in df.iterrows():
        accession = row["assembly_accession"]
        protein_id = row["protein_id"]
        protein_file = PROTEIN_DIR / f"{accession}.faa"

        if not protein_file.exists():
            continue

        found = False
        for record in SeqIO.parse(protein_file, "fasta"):
            if record.id == protein_id:
                found = True
                record.id = f"{accession}|{protein_id}"
                record.description = (
                    f"subject={row['sseqid']} "
                    f"pident={row['pident']} "
                    f"length={row['length']} "
                    f"evalue={row['evalue']} "
                    f"bitscore={row['bitscore']} "
                    f"tier={row['priority_tier']}"
                )
                selected_records.append(record)

                annotated_rows.append({
                    "assembly_accession": accession,
                    "protein_id": protein_id,
                    "qseqid": row["qseqid"],
                    "sseqid": row["sseqid"],
                    "pident": row["pident"],
                    "length": row["length"],
                    "evalue": row["evalue"],
                    "bitscore": row["bitscore"],
                    "priority_tier": row["priority_tier"],
                    "sequence_length_aa": len(record.seq),
                })
                break

        if not found:
            annotated_rows.append({
                "assembly_accession": accession,
                "protein_id": protein_id,
                "qseqid": row["qseqid"],
                "sseqid": row["sseqid"],
                "pident": row["pident"],
                "length": row["length"],
                "evalue": row["evalue"],
                "bitscore": row["bitscore"],
                "priority_tier": row["priority_tier"],
                "sequence_length_aa": "",
            })

    SeqIO.write(selected_records, OUT_FASTA, "fasta")
    pd.DataFrame(annotated_rows).to_csv(OUT_TABLE, sep="\t", index=False)

    print(f"Priority candidate FASTA written to: {OUT_FASTA}")
    print(f"Priority candidate table written to: {OUT_TABLE}")
    print(f"Sequences extracted: {len(selected_records)}")


if __name__ == "__main__":
    main()
