#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
01_select_refseq_genomes.py

Pipeline step:
Selection of genome assemblies and metadata from NCBI Assembly.

Purpose:
Select up to N genome assemblies per target species, prioritizing RefSeq/GCF
assemblies and using GenBank/GCA assemblies as a documented fallback when
RefSeq availability is insufficient.

Main principles:
1. Use NCBI Assembly through Entrez.
2. Search by taxon and latest assemblies.
3. Avoid the refseq[filter] term because it may cause HTTP 400 errors.
4. Prefer RefSeq/GCF assemblies.
5. Use GenBank/GCA only as fallback to reach target_n when RefSeq is limited.
6. Prioritize non-anomalous assemblies, FTP availability, reference or
   representative status, higher assembly level, and recent updates.
7. Produce traceable metadata tables and accession list for NCBI Datasets CLI.

Input:
00_config/target_species.tsv

Required columns:
group
species
taxon_id
target_n
priority_note

Outputs:
01_metadata/all_candidate_assemblies.tsv
01_metadata/selected_assemblies.tsv
01_metadata/selection_log.tsv
01_metadata/run_summary.txt
02_accessions/selected_accessions.txt
"""

import os
import time
from pathlib import Path
from datetime import datetime

import pandas as pd
from tqdm import tqdm
from Bio import Entrez


# =============================================================================
# Paths
# =============================================================================

BASE_DIR = Path("<PROJECT_ROOT>")

CONFIG_FILE = BASE_DIR / "00_config" / "target_species.tsv"
METADATA_DIR = BASE_DIR / "01_metadata"
ACCESSION_DIR = BASE_DIR / "02_accessions"
LOG_DIR = BASE_DIR / "logs"

METADATA_DIR.mkdir(parents=True, exist_ok=True)
ACCESSION_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# NCBI Entrez configuration
# =============================================================================

NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "<SET_NCBI_EMAIL>")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")

Entrez.email = NCBI_EMAIL
Entrez.tool = "seh_epoxide_hydrolase_genome_capture"

if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

# Conservative delay:
# - without API key: NCBI allows 3 requests/s
# - with API key: NCBI allows higher throughput, but we keep it gentle
REQUEST_DELAY = 0.15 if NCBI_API_KEY else 0.45


# =============================================================================
# Ranking schemes
# =============================================================================

ASSEMBLY_LEVEL_RANK = {
    "Complete Genome": 4,
    "Chromosome": 3,
    "Scaffold": 2,
    "Contig": 1,
}

REFSEQ_CATEGORY_RANK = {
    "reference genome": 3,
    "representative genome": 2,
    "na": 1,
    "": 1,
}


# =============================================================================
# Helper functions
# =============================================================================

def clean_text(value):
    """Return clean string from Entrez fields."""
    if value is None:
        return ""
    return str(value).strip()


def get_source_database(accession):
    """Classify assembly accession source."""
    accession = clean_text(accession)
    if accession.startswith("GCF_"):
        return "RefSeq"
    if accession.startswith("GCA_"):
        return "GenBank"
    return "Other"


def get_source_rank(source_database):
    """RefSeq is preferred over GenBank."""
    if source_database == "RefSeq":
        return 2
    if source_database == "GenBank":
        return 1
    return 0


def entrez_esearch_assembly(taxon_id, retmax=3000):
    """
    Search NCBI Assembly by taxon using latest assemblies.

    Important:
    We do not use refseq[filter] because this can cause HTTP 400 errors in
    some Entrez Assembly queries. RefSeq/GenBank filtering is done after
    retrieving Assembly summaries.
    """
    query = f"txid{taxon_id}[Organism:exp] AND latest[filter]"

    handle = Entrez.esearch(
        db="assembly",
        term=query,
        retmax=retmax,
        sort="relevance",
    )
    result = Entrez.read(handle)
    handle.close()

    time.sleep(REQUEST_DELAY)

    return result.get("IdList", []), query, result.get("Count", "0")


def entrez_esummary_assembly(id_list, batch_size=100):
    """Fetch NCBI Assembly summaries in batches."""
    records = []

    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i + batch_size]

        handle = Entrez.esummary(
            db="assembly",
            id=",".join(batch),
            report="full",
        )
        result = Entrez.read(handle)
        handle.close()

        time.sleep(REQUEST_DELAY)

        docs = result.get("DocumentSummarySet", {}).get("DocumentSummary", [])
        records.extend(docs)

    return records


def normalize_assembly_record(doc, species, group, taxon_id, target_n, priority_note):
    """Normalize one Entrez Assembly record into a flat metadata row."""

    assembly_accession = clean_text(doc.get("AssemblyAccession"))
    source_database = get_source_database(assembly_accession)

    assembly_name = clean_text(doc.get("AssemblyName"))
    organism = clean_text(doc.get("Organism"))

    taxid = clean_text(doc.get("Taxid"))
    species_taxid = clean_text(doc.get("SpeciesTaxid"))

    assembly_level = clean_text(doc.get("AssemblyStatus"))
    refseq_category = clean_text(doc.get("RefSeq_category"))

    submitter = clean_text(doc.get("SubmitterOrganization"))
    release_date = clean_text(doc.get("SeqReleaseDate"))
    update_date = clean_text(doc.get("AsmReleaseDate")) or release_date

    ftp_path_refseq = clean_text(doc.get("FtpPath_RefSeq"))
    ftp_path_genbank = clean_text(doc.get("FtpPath_GenBank"))

    biosample = clean_text(doc.get("BioSampleAccn"))
    bioproject = clean_text(doc.get("BioProjectAccn"))

    assembly_type = clean_text(doc.get("AssemblyType"))

    excluded_from_refseq = clean_text(doc.get("ExclFromRefSeq"))
    anomalous = clean_text(doc.get("AnomalousList"))

    has_any_ftp = bool(ftp_path_refseq) or bool(ftp_path_genbank)
    anomaly_flag = bool(excluded_from_refseq) or bool(anomalous)

    assembly_level_rank = ASSEMBLY_LEVEL_RANK.get(assembly_level, 0)
    refseq_category_rank = REFSEQ_CATEGORY_RANK.get(refseq_category.lower(), 1)

    return {
        "target_group": group,
        "target_species": species,
        "target_taxon_id": taxon_id,
        "target_n": target_n,
        "priority_note": priority_note,

        "organism": organism,
        "taxid": taxid,
        "species_taxid": species_taxid,

        "assembly_accession": assembly_accession,
        "source_database": source_database,
        "source_database_rank": get_source_rank(source_database),

        "assembly_name": assembly_name,
        "assembly_level": assembly_level,
        "assembly_level_rank": assembly_level_rank,

        "refseq_category": refseq_category,
        "refseq_category_rank": refseq_category_rank,

        "assembly_type": assembly_type,
        "biosample": biosample,
        "bioproject": bioproject,
        "submitter": submitter,

        "release_date": release_date,
        "update_date": update_date,

        "ftp_path_refseq": ftp_path_refseq,
        "ftp_path_genbank": ftp_path_genbank,
        "has_any_ftp": has_any_ftp,

        "excluded_from_refseq": excluded_from_refseq,
        "anomalous": anomalous,
        "anomaly_flag": anomaly_flag,
    }


def score_assemblies(df_species):
    """Assign a reproducible quality/preference score to candidate assemblies."""

    df = df_species.copy()

    df["update_date_parsed"] = pd.to_datetime(df["update_date"], errors="coerce")
    df["release_date_parsed"] = pd.to_datetime(df["release_date"], errors="coerce")

    df["has_any_ftp_score"] = df["has_any_ftp"].astype(int)
    df["non_anomalous_score"] = (~df["anomaly_flag"]).astype(int)

    # Composite score:
    # RefSeq/GCF strongly preferred.
    # GenBank/GCA is allowed as fallback.
    df["selection_score"] = (
        df["source_database_rank"] * 1000
        + df["non_anomalous_score"] * 500
        + df["has_any_ftp_score"] * 250
        + df["refseq_category_rank"] * 100
        + df["assembly_level_rank"] * 50
    )

    return df


def select_assemblies(df_species, target_n):
    """
    Select up to target_n assemblies using quality-aware ranking.

    The sorting keeps RefSeq as priority, but allows GenBank fallback.
    """
    df = score_assemblies(df_species)

    # Prefer non-anomalous assemblies.
    non_anomalous = df[~df["anomaly_flag"]].copy()
    if not non_anomalous.empty:
        df_for_selection = non_anomalous
    else:
        df_for_selection = df.copy()

    df_for_selection = df_for_selection.sort_values(
        by=[
            "selection_score",
            "source_database_rank",
            "refseq_category_rank",
            "assembly_level_rank",
            "update_date_parsed",
            "release_date_parsed",
            "assembly_accession",
        ],
        ascending=[False, False, False, False, False, False, True],
    )

    selected = df_for_selection.head(target_n).copy()
    selected["selected"] = True
    selected["selection_rank_within_species"] = range(1, len(selected) + 1)

    df["selected"] = df["assembly_accession"].isin(selected["assembly_accession"])
    df["selection_rank_within_species"] = ""

    rank_map = dict(
        zip(
            selected["assembly_accession"],
            selected["selection_rank_within_species"],
        )
    )

    df.loc[df["selected"], "selection_rank_within_species"] = df.loc[
        df["selected"], "assembly_accession"
    ].map(rank_map)

    n_refseq_selected = int((selected["source_database"] == "RefSeq").sum())
    n_genbank_selected = int((selected["source_database"] == "GenBank").sum())

    return df, selected, n_refseq_selected, n_genbank_selected


def build_empty_all_candidates():
    return pd.DataFrame(columns=[
        "target_group",
        "target_species",
        "target_taxon_id",
        "target_n",
        "priority_note",
        "organism",
        "taxid",
        "species_taxid",
        "assembly_accession",
        "source_database",
        "source_database_rank",
        "assembly_name",
        "assembly_level",
        "assembly_level_rank",
        "refseq_category",
        "refseq_category_rank",
        "assembly_type",
        "biosample",
        "bioproject",
        "submitter",
        "release_date",
        "update_date",
        "ftp_path_refseq",
        "ftp_path_genbank",
        "has_any_ftp",
        "excluded_from_refseq",
        "anomalous",
        "anomaly_flag",
        "selection_score",
        "selected",
        "selection_rank_within_species",
    ])


def build_empty_selected():
    return pd.DataFrame(columns=[
        "target_group",
        "target_species",
        "target_taxon_id",
        "target_n",
        "priority_note",
        "organism",
        "taxid",
        "species_taxid",
        "assembly_accession",
        "source_database",
        "source_database_rank",
        "assembly_name",
        "assembly_level",
        "assembly_level_rank",
        "refseq_category",
        "refseq_category_rank",
        "assembly_type",
        "biosample",
        "bioproject",
        "submitter",
        "release_date",
        "update_date",
        "ftp_path_refseq",
        "ftp_path_genbank",
        "has_any_ftp",
        "excluded_from_refseq",
        "anomalous",
        "anomaly_flag",
        "selection_score",
        "selected",
        "selection_rank_within_species",
    ])


def main():
    started = datetime.now()

    if not CONFIG_FILE.exists():
        raise FileNotFoundError(f"Config file not found: {CONFIG_FILE}")

    targets = pd.read_csv(CONFIG_FILE, sep="\t")

    required_cols = {"group", "species", "taxon_id", "target_n"}
    missing = required_cols - set(targets.columns)
    if missing:
        raise ValueError(f"Missing required columns in target_species.tsv: {missing}")

    if "priority_note" not in targets.columns:
        targets["priority_note"] = ""

    all_candidate_rows = []
    selected_rows = []
    log_rows = []

    for _, row in tqdm(targets.iterrows(), total=len(targets), desc="Selecting genomes"):
        group = str(row["group"])
        species = str(row["species"])
        taxon_id = str(row["taxon_id"])
        target_n = int(row["target_n"])
        priority_note = str(row.get("priority_note", ""))

        try:
            ids, query, ncbi_count = entrez_esearch_assembly(taxon_id)

            if not ids:
                log_rows.append({
                    "species": species,
                    "taxon_id": taxon_id,
                    "query": query,
                    "ncbi_count": ncbi_count,
                    "n_assembly_ids_retrieved": 0,
                    "n_candidates_after_summary": 0,
                    "n_refseq_candidates": 0,
                    "n_genbank_candidates": 0,
                    "n_selected": 0,
                    "n_refseq_selected": 0,
                    "n_genbank_selected": 0,
                    "status": "no_assemblies_found",
                    "message": "No assembly IDs returned by Entrez.",
                })
                continue

            docs = entrez_esummary_assembly(ids)

            normalized = [
                normalize_assembly_record(
                    doc=doc,
                    species=species,
                    group=group,
                    taxon_id=taxon_id,
                    target_n=target_n,
                    priority_note=priority_note,
                )
                for doc in docs
            ]

            df_species = pd.DataFrame(normalized)

            if df_species.empty:
                log_rows.append({
                    "species": species,
                    "taxon_id": taxon_id,
                    "query": query,
                    "ncbi_count": ncbi_count,
                    "n_assembly_ids_retrieved": len(ids),
                    "n_candidates_after_summary": 0,
                    "n_refseq_candidates": 0,
                    "n_genbank_candidates": 0,
                    "n_selected": 0,
                    "n_refseq_selected": 0,
                    "n_genbank_selected": 0,
                    "status": "empty_summary",
                    "message": "Entrez returned IDs but no assembly summaries.",
                })
                continue

            # Keep only GCF_ or GCA_ assemblies.
            df_species = df_species[
                df_species["assembly_accession"].str.startswith(("GCF_", "GCA_"), na=False)
            ].copy()

            if df_species.empty:
                log_rows.append({
                    "species": species,
                    "taxon_id": taxon_id,
                    "query": query,
                    "ncbi_count": ncbi_count,
                    "n_assembly_ids_retrieved": len(ids),
                    "n_candidates_after_summary": 0,
                    "n_refseq_candidates": 0,
                    "n_genbank_candidates": 0,
                    "n_selected": 0,
                    "n_refseq_selected": 0,
                    "n_genbank_selected": 0,
                    "status": "no_gcf_or_gca",
                    "message": "No GCF_ or GCA_ assemblies retained.",
                })
                continue

            n_refseq_candidates = int((df_species["source_database"] == "RefSeq").sum())
            n_genbank_candidates = int((df_species["source_database"] == "GenBank").sum())

            scored, selected, n_refseq_selected, n_genbank_selected = select_assemblies(
                df_species=df_species,
                target_n=target_n,
            )

            all_candidate_rows.extend(scored.to_dict("records"))
            selected_rows.extend(selected.to_dict("records"))

            if len(selected) >= target_n:
                status = "selection_ok"
                message = (
                    f"Selected {len(selected)}/{target_n}: "
                    f"{n_refseq_selected} RefSeq and {n_genbank_selected} GenBank."
                )
            elif len(selected) >= 5:
                status = "selection_partial_but_acceptable"
                message = (
                    f"Selected {len(selected)}/{target_n}; at least 5 available: "
                    f"{n_refseq_selected} RefSeq and {n_genbank_selected} GenBank."
                )
            else:
                status = "selection_limited"
                message = (
                    f"Selected only {len(selected)}/{target_n}: "
                    f"{n_refseq_selected} RefSeq and {n_genbank_selected} GenBank. "
                    "This species may have limited public assemblies."
                )

            log_rows.append({
                "species": species,
                "taxon_id": taxon_id,
                "query": query,
                "ncbi_count": ncbi_count,
                "n_assembly_ids_retrieved": len(ids),
                "n_candidates_after_summary": len(df_species),
                "n_refseq_candidates": n_refseq_candidates,
                "n_genbank_candidates": n_genbank_candidates,
                "n_selected": len(selected),
                "n_refseq_selected": n_refseq_selected,
                "n_genbank_selected": n_genbank_selected,
                "status": status,
                "message": message,
            })

        except Exception as e:
            log_rows.append({
                "species": species,
                "taxon_id": taxon_id,
                "query": "",
                "ncbi_count": "",
                "n_assembly_ids_retrieved": "",
                "n_candidates_after_summary": "",
                "n_refseq_candidates": "",
                "n_genbank_candidates": "",
                "n_selected": 0,
                "n_refseq_selected": 0,
                "n_genbank_selected": 0,
                "status": "error",
                "message": repr(e),
            })

    all_candidates = pd.DataFrame(all_candidate_rows)
    selected_df = pd.DataFrame(selected_rows)
    log_df = pd.DataFrame(log_rows)

    if all_candidates.empty:
        all_candidates = build_empty_all_candidates()

    if selected_df.empty:
        selected_df = build_empty_selected()

    all_candidates_file = METADATA_DIR / "all_candidate_assemblies.tsv"
    selected_file = METADATA_DIR / "selected_assemblies.tsv"
    log_file = METADATA_DIR / "selection_log.tsv"
    accessions_file = ACCESSION_DIR / "selected_accessions.txt"
    summary_file = METADATA_DIR / "run_summary.txt"

    all_candidates.to_csv(all_candidates_file, sep="\t", index=False)
    selected_df.to_csv(selected_file, sep="\t", index=False)

    selected_df["assembly_accession"].dropna().drop_duplicates().to_csv(
        accessions_file,
        index=False,
        header=False,
    )

    log_df.to_csv(log_file, sep="\t", index=False)

    finished = datetime.now()
    elapsed = finished - started

    with open(summary_file, "w", encoding="utf-8") as f:
        f.write("Genome capture selection summary\n")
        f.write("================================\n\n")
        f.write(f"Started: {started}\n")
        f.write(f"Finished: {finished}\n")
        f.write(f"Elapsed: {elapsed}\n")
        f.write(f"NCBI email configured: {NCBI_EMAIL}\n")
        f.write(f"NCBI API key configured: {'yes' if bool(NCBI_API_KEY) else 'no'}\n\n")

        f.write("Selection strategy:\n")
        f.write("- Latest NCBI Assembly records were queried by taxon.\n")
        f.write("- RefSeq/GCF assemblies were prioritized.\n")
        f.write("- GenBank/GCA assemblies were used as documented fallback when RefSeq availability was insufficient.\n")
        f.write("- Non-anomalous assemblies, FTP availability, reference/representative category, assembly level, and update date were prioritized.\n\n")

        f.write(f"Target species: {len(targets)}\n")
        f.write(f"Candidate assemblies retained: {len(all_candidates)}\n")
        f.write(f"Selected assemblies: {len(selected_df)}\n\n")

        if not selected_df.empty:
            f.write("Selected assemblies by species and source database:\n")
            counts = (
                selected_df
                .groupby(["target_group", "target_species", "source_database"])
                .size()
                .reset_index(name="n_selected")
            )
            f.write(counts.to_string(index=False))
            f.write("\n\n")

            f.write("Total selected by source database:\n")
            source_counts = (
                selected_df
                .groupby("source_database")
                .size()
                .reset_index(name="n_selected")
            )
            f.write(source_counts.to_string(index=False))
            f.write("\n\n")

        f.write("Output files:\n")
        f.write(f"- {all_candidates_file}\n")
        f.write(f"- {selected_file}\n")
        f.write(f"- {log_file}\n")
        f.write(f"- {accessions_file}\n")

    print("\nDone.")
    print(f"All candidates: {all_candidates_file}")
    print(f"Selected assemblies: {selected_file}")
    print(f"Accessions: {accessions_file}")
    print(f"Log: {log_file}")
    print(f"Summary: {summary_file}")


if __name__ == "__main__":
    main()
