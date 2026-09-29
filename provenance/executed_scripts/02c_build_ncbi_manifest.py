#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd

BASE_DIR = Path("<PROJECT_ROOT>")
SELECTED_FILE = BASE_DIR / "01_metadata" / "selected_assemblies.tsv"
OUT_FILE = BASE_DIR / "03_downloads" / "ncbi_download_manifest.tsv"

BASE_DIR.joinpath("03_downloads").mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SELECTED_FILE, sep="\t")

required = {
    "assembly_accession",
    "assembly_name",
    "source_database",
    "ftp_path_refseq",
    "ftp_path_genbank",
    "target_group",
    "target_species",
}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns in selected_assemblies.tsv: {missing}")

rows = []

for _, row in df.iterrows():
    accession = str(row["assembly_accession"]).strip()
    assembly_name = str(row["assembly_name"]).strip()
    source_database = str(row["source_database"]).strip()
    target_group = str(row["target_group"]).strip()
    target_species = str(row["target_species"]).strip()

    ftp_refseq = str(row.get("ftp_path_refseq", "")).strip()
    ftp_genbank = str(row.get("ftp_path_genbank", "")).strip()

    ftp_base = ftp_refseq if ftp_refseq and ftp_refseq != "nan" else ftp_genbank
    if not ftp_base or ftp_base == "nan":
        continue

    # Use HTTPS instead of FTP for better reliability
    https_base = ftp_base.replace("ftp://", "https://")

    base_name = ftp_base.rstrip("/").split("/")[-1]

    file_map = {
        "genome_fna_gz": f"{https_base}/{base_name}_genomic.fna.gz",
        "protein_faa_gz": f"{https_base}/{base_name}_protein.faa.gz",
        "cds_fna_gz": f"{https_base}/{base_name}_cds_from_genomic.fna.gz",
        "gff_gz": f"{https_base}/{base_name}_genomic.gff.gz",
        "gbff_gz": f"{https_base}/{base_name}_genomic.gbff.gz",
        "assembly_report_txt": f"{https_base}/{base_name}_assembly_report.txt",
        "md5checksums_txt": f"{https_base}/md5checksums.txt",
    }

    rows.append({
        "target_group": target_group,
        "target_species": target_species,
        "assembly_accession": accession,
        "assembly_name": assembly_name,
        "source_database": source_database,
        "ftp_base": ftp_base,
        "https_base": https_base,
        **file_map,
    })

manifest = pd.DataFrame(rows)
manifest.to_csv(OUT_FILE, sep="\t", index=False)

print(f"Manifest written to: {OUT_FILE}")
print(f"Assemblies in manifest: {len(manifest)}")
