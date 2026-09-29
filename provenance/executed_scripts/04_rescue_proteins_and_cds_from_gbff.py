#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


BASE_DIR = Path("<PROJECT_ROOT>")
INVENTORY_FILE = BASE_DIR / "05_qc" / "file_inventory.tsv"

DERIVED_DIR = BASE_DIR / "04_derived"
PROTEIN_DIR = DERIVED_DIR / "proteins_final"
CDS_DIR = DERIVED_DIR / "cds_final"
REPORT_DIR = DERIVED_DIR / "rescue_reports"

PROTEIN_DIR.mkdir(parents=True, exist_ok=True)
CDS_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)


def safe_text(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def extract_from_gbff(gbff_file, accession):
    proteins = []
    cds_records = []

    for record in SeqIO.parse(gbff_file, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue

            qualifiers = feature.qualifiers

            locus_tag = qualifiers.get("locus_tag", ["unknown_locus"])[0]
            protein_id = qualifiers.get("protein_id", [f"{accession}_{locus_tag}"])[0]
            product = qualifiers.get("product", ["hypothetical protein"])[0]
            gene = qualifiers.get("gene", [""])[0]

            # protein sequence
            translation = qualifiers.get("translation", [])
            if translation:
                aa_seq = translation[0]
                proteins.append(
                    SeqRecord(
                        Seq(aa_seq),
                        id=protein_id,
                        description=f"gene={gene} locus_tag={locus_tag} product={product} source={accession}"
                    )
                )

            # cds sequence
            try:
                cds_seq = feature.extract(record.seq)
                cds_records.append(
                    SeqRecord(
                        cds_seq,
                        id=protein_id,
                        description=f"gene={gene} locus_tag={locus_tag} product={product} source={accession}"
                    )
                )
            except Exception:
                pass

    return proteins, cds_records


def main():
    if not INVENTORY_FILE.exists():
        raise FileNotFoundError(f"Missing inventory file: {INVENTORY_FILE}")

    inv = pd.read_csv(INVENTORY_FILE, sep="\t")

    report_rows = []

    for _, row in inv.iterrows():
        accession = safe_text(row["assembly_accession"])
        species = safe_text(row["target_species"])
        group = safe_text(row["target_group"])

        protein_faa = safe_text(row.get("protein_faa", ""))
        cds_fna = safe_text(row.get("cds_fna", ""))
        gbff = safe_text(row.get("gbff", ""))

        out_protein = PROTEIN_DIR / f"{accession}.faa"
        out_cds = CDS_DIR / f"{accession}.fna"

        protein_source = ""
        cds_source = ""
        n_proteins = 0
        n_cds = 0

        # Proteins
        if protein_faa and Path(protein_faa).exists():
            records = list(SeqIO.parse(protein_faa, "fasta"))
            if records:
                SeqIO.write(records, out_protein, "fasta")
                protein_source = "native_protein_faa"
                n_proteins = len(records)
        elif gbff and Path(gbff).exists():
            proteins, cds_records_tmp = extract_from_gbff(gbff, accession)
            if proteins:
                SeqIO.write(proteins, out_protein, "fasta")
                protein_source = "rescued_from_gbff"
                n_proteins = len(proteins)

        # CDS
        if cds_fna and Path(cds_fna).exists():
            records = list(SeqIO.parse(cds_fna, "fasta"))
            if records:
                SeqIO.write(records, out_cds, "fasta")
                cds_source = "native_cds_fna"
                n_cds = len(records)
        elif gbff and Path(gbff).exists():
            proteins_tmp, cds_records = extract_from_gbff(gbff, accession)
            if cds_records:
                SeqIO.write(cds_records, out_cds, "fasta")
                cds_source = "rescued_from_gbff"
                n_cds = len(cds_records)

        report_rows.append({
            "target_group": group,
            "target_species": species,
            "assembly_accession": accession,
            "protein_source": protein_source if protein_source else "missing",
            "cds_source": cds_source if cds_source else "missing",
            "n_proteins": n_proteins,
            "n_cds": n_cds,
            "protein_output": str(out_protein) if out_protein.exists() else "",
            "cds_output": str(out_cds) if out_cds.exists() else "",
        })

    report = pd.DataFrame(report_rows)
    report_file = REPORT_DIR / "protein_cds_rescue_report.tsv"
    report.to_csv(report_file, sep="\t", index=False)

    summary = (
        report.groupby(["target_group", "target_species"])
        .agg(
            n_assemblies=("assembly_accession", "nunique"),
            n_with_proteins=( "protein_output", lambda x: sum(bool(v) for v in x) ),
            n_with_cds=( "cds_output", lambda x: sum(bool(v) for v in x) ),
        )
        .reset_index()
    )

    summary_file = REPORT_DIR / "protein_cds_rescue_summary.tsv"
    summary.to_csv(summary_file, sep="\t", index=False)

    print(f"Rescue report written to: {report_file}")
    print(f"Summary written to: {summary_file}")
    print(f"Assemblies processed: {len(report)}")


if __name__ == "__main__":
    main()
