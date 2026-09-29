#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
03_build_inventory.py

Purpose
-------
Inventory and QC all directly downloaded NCBI genome files organized under
03_downloads/direct_ncbi_files.

Outputs
-------
05_qc/file_inventory.tsv
05_qc/genome_capture_qc.tsv
05_qc/missing_expected_files.tsv
"""

from pathlib import Path
import pandas as pd


BASE_DIR = Path("<PROJECT_ROOT>")

SELECTED_FILE = BASE_DIR / "01_metadata" / "selected_assemblies.tsv"
DIRECT_ROOT = BASE_DIR / "03_downloads" / "direct_ncbi_files"
QC_DIR = BASE_DIR / "05_qc"

QC_DIR.mkdir(parents=True, exist_ok=True)


def find_first(paths):
    for p in paths:
        if p.exists() and p.is_file():
            return str(p)
    return ""


def main():
    if not SELECTED_FILE.exists():
        raise FileNotFoundError(f"Missing selected file: {SELECTED_FILE}")

    if not DIRECT_ROOT.exists():
        raise FileNotFoundError(f"Missing direct download root: {DIRECT_ROOT}")

    selected = pd.read_csv(SELECTED_FILE, sep="\t")

    inventory_rows = []

    for _, row in selected.iterrows():
        accession = str(row["assembly_accession"]).strip()
        species = str(row["target_species"]).strip()
        group = str(row["target_group"]).strip()
        source_database = str(row.get("source_database", "")).strip()
        assembly_name = str(row.get("assembly_name", "")).strip()

        species_safe = species.replace(" ", "_").replace("/", "_")
        accession_dir = DIRECT_ROOT / species_safe / accession

        genome_fna = find_first(list(accession_dir.glob("*_genomic.fna")))
        genome_fna_gz = find_first(list(accession_dir.glob("*_genomic.fna.gz")))

        protein_faa = find_first(list(accession_dir.glob("*_protein.faa")))
        protein_faa_gz = find_first(list(accession_dir.glob("*_protein.faa.gz")))

        cds_fna = find_first(list(accession_dir.glob("*_cds_from_genomic.fna")))
        cds_fna_gz = find_first(list(accession_dir.glob("*_cds_from_genomic.fna.gz")))

        gff = find_first(list(accession_dir.glob("*_genomic.gff")))
        gff_gz = find_first(list(accession_dir.glob("*_genomic.gff.gz")))

        gbff = find_first(list(accession_dir.glob("*_genomic.gbff")))
        gbff_gz = find_first(list(accession_dir.glob("*_genomic.gbff.gz")))

        assembly_report = find_first(list(accession_dir.glob("*_assembly_report.txt")))
        md5_file = find_first(list(accession_dir.glob("md5checksums.txt")))

        inventory_rows.append({
            "target_group": group,
            "target_species": species,
            "assembly_accession": accession,
            "assembly_name": assembly_name,
            "source_database": source_database,
            "accession_dir": str(accession_dir),

            "genome_fna": genome_fna,
            "genome_fna_gz": genome_fna_gz,
            "protein_faa": protein_faa,
            "protein_faa_gz": protein_faa_gz,
            "cds_fna": cds_fna,
            "cds_fna_gz": cds_fna_gz,
            "gff": gff,
            "gff_gz": gff_gz,
            "gbff": gbff,
            "gbff_gz": gbff_gz,
            "assembly_report": assembly_report,
            "md5checksums": md5_file,

            "has_accession_dir": accession_dir.exists(),
            "has_genome_fna": bool(genome_fna),
            "has_protein_faa": bool(protein_faa),
            "has_cds_fna": bool(cds_fna),
            "has_gff": bool(gff),
            "has_gbff": bool(gbff),
            "has_assembly_report": bool(assembly_report),
            "has_md5checksums": bool(md5_file),
        })

    inventory = pd.DataFrame(inventory_rows)

    inventory_file = QC_DIR / "file_inventory.tsv"
    inventory.to_csv(inventory_file, sep="\t", index=False)

    qc = (
        inventory.groupby(["target_group", "target_species"])
        .agg(
            n_assemblies=("assembly_accession", "nunique"),
            n_with_accession_dir=("has_accession_dir", "sum"),
            n_with_genome_fna=("has_genome_fna", "sum"),
            n_with_protein_faa=("has_protein_faa", "sum"),
            n_with_cds_fna=("has_cds_fna", "sum"),
            n_with_gff=("has_gff", "sum"),
            n_with_gbff=("has_gbff", "sum"),
            n_with_assembly_report=("has_assembly_report", "sum"),
            n_with_md5checksums=("has_md5checksums", "sum"),
        )
        .reset_index()
    )

    qc_file = QC_DIR / "genome_capture_qc.tsv"
    qc.to_csv(qc_file, sep="\t", index=False)

    missing_rows = []
    expected_cols = [
        ("has_accession_dir", "accession_dir"),
        ("has_genome_fna", "genome_fna"),
        ("has_protein_faa", "protein_faa"),
        ("has_cds_fna", "cds_fna"),
        ("has_gff", "gff"),
        ("has_gbff", "gbff"),
        ("has_assembly_report", "assembly_report"),
    ]

    for _, row in inventory.iterrows():
        for bool_col, file_label in expected_cols:
            if not bool(row[bool_col]):
                missing_rows.append({
                    "target_group": row["target_group"],
                    "target_species": row["target_species"],
                    "assembly_accession": row["assembly_accession"],
                    "missing_item": file_label,
                })

    missing = pd.DataFrame(missing_rows)
    missing_file = QC_DIR / "missing_expected_files.tsv"
    missing.to_csv(missing_file, sep="\t", index=False)

    print(f"Inventory written to: {inventory_file}")
    print(f"QC summary written to: {qc_file}")
    print(f"Missing-file report written to: {missing_file}")
    print(f"Assemblies inventoried: {len(inventory)}")


if __name__ == "__main__":
    main()
