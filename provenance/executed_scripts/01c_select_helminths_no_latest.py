#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
01c_select_helminths_no_latest.py

Purpose
-------
Rescue genome selection for the originally selected helminths
(Schistosoma mansoni and Strongyloides stercoralis) without using
latest[filter], because that filter may be too restrictive for some
eukaryotic parasites in NCBI Assembly.

Selection strategy
------------------
1. Search NCBI Assembly by taxon WITHOUT latest[filter].
2. Retrieve full assembly summaries.
3. Keep GCF_ and GCA_ assemblies.
4. Prioritize:
   - RefSeq over GenBank
   - non-anomalous assemblies
   - assemblies with FTP available
   - reference/representative genomes
   - assembly level (Complete Genome > Chromosome > Scaffold > Contig)
   - recent update/release dates
5. Select up to target_n assemblies per species.

Outputs
-------
01_metadata/helminths_no_latest_candidates.tsv
01_metadata/helminths_no_latest_selected.tsv
01_metadata/helminths_no_latest_selection_log.tsv
"""

import os
import time
from pathlib import Path
from datetime import datetime

import pandas as pd
from tqdm import tqdm
from Bio import Entrez


BASE_DIR = Path("<PROJECT_ROOT>")
CONFIG_FILE = BASE_DIR / "00_config" / "target_species_helminths_only.tsv"
METADATA_DIR = BASE_DIR / "01_metadata"
LOG_DIR = BASE_DIR / "logs"

METADATA_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "<SET_NCBI_EMAIL>")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY", "")

Entrez.email = NCBI_EMAIL
Entrez.tool = "seh_epoxide_hydrolase_helminth_rescue_no_latest"

if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

REQUEST_DELAY = 0.15 if NCBI_API_KEY else 0.45

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


def clean_text(value):
    if value is None:
        return ""
    return str(value).strip()


def get_source_database(accession):
    accession = clean_text(accession)
    if accession.startswith("GCF_"):
        return "RefSeq"
    if accession.startswith("GCA_"):
        return "GenBank"
    return "Other"


def get_source_rank(source_database):
    if source_database == "RefSeq":
        return 2
    if source_database == "GenBank":
        return 1
    return 0


def entrez_esearch_assembly_no_latest(taxon_id, retmax=5000):
    """
    Search NCBI Assembly by taxon WITHOUT latest[filter].
    This is intentionally more permissive for helminths.
    """
    query = f"txid{taxon_id}[Organism:exp]"

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


def score_assemblies(df):
    df = df.copy()

    df["update_date_parsed"] = pd.to_datetime(df["update_date"], errors="coerce")
    df["release_date_parsed"] = pd.to_datetime(df["release_date"], errors="coerce")

    df["has_any_ftp_score"] = df["has_any_ftp"].astype(int)
    df["non_anomalous_score"] = (~df["anomaly_flag"]).astype(int)

    df["selection_score"] = (
        df["source_database_rank"] * 1000
        + df["non_anomalous_score"] * 500
        + df["has_any_ftp_score"] * 250
        + df["refseq_category_rank"] * 100
        + df["assembly_level_rank"] * 50
    )

    return df


def select_assemblies(df_species, target_n):
    df = score_assemblies(df_species)

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

    rank_map = dict(zip(selected["assembly_accession"], selected["selection_rank_within_species"]))

    df.loc[df["selected"], "selection_rank_within_species"] = df.loc[
        df["selected"], "assembly_accession"
    ].map(rank_map)

    return df, selected


def main():
    started = datetime.now()

    if not CONFIG_FILE.exists():
        raise FileNotFoundError(f"Config file not found: {CONFIG_FILE}")

    targets = pd.read_csv(CONFIG_FILE, sep="\t")

    all_candidate_rows = []
    selected_rows = []
    log_rows = []

    for _, row in tqdm(targets.iterrows(), total=len(targets), desc="Helminth no-latest rescue"):
        group = str(row["group"])
        species = str(row["species"])
        taxon_id = str(row["taxon_id"])
        target_n = int(row["target_n"])
        priority_note = str(row.get("priority_note", ""))

        try:
            ids, query, ncbi_count = entrez_esearch_assembly_no_latest(taxon_id)

            if not ids:
                log_rows.append({
                    "species": species,
                    "taxon_id": taxon_id,
                    "query": query,
                    "ncbi_count": ncbi_count,
                    "n_selected": 0,
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
                    "n_selected": 0,
                    "status": "empty_summary",
                    "message": "No assembly summaries returned.",
                })
                continue

            df_species = df_species[
                df_species["assembly_accession"].str.startswith(("GCF_", "GCA_"), na=False)
            ].copy()

            if df_species.empty:
                log_rows.append({
                    "species": species,
                    "taxon_id": taxon_id,
                    "query": query,
                    "ncbi_count": ncbi_count,
                    "n_selected": 0,
                    "status": "no_gcf_or_gca",
                    "message": "No GCF/GCA accessions after filtering.",
                })
                continue

            scored, selected = select_assemblies(df_species, target_n)

            all_candidate_rows.extend(scored.to_dict("records"))
            selected_rows.extend(selected.to_dict("records"))

            n_refseq = int((selected["source_database"] == "RefSeq").sum())
            n_genbank = int((selected["source_database"] == "GenBank").sum())

            log_rows.append({
                "species": species,
                "taxon_id": taxon_id,
                "query": query,
                "ncbi_count": ncbi_count,
                "n_selected": len(selected),
                "n_refseq_selected": n_refseq,
                "n_genbank_selected": n_genbank,
                "status": "selection_ok" if len(selected) > 0 else "selection_limited",
                "message": f"Selected {len(selected)} assemblies: {n_refseq} RefSeq, {n_genbank} GenBank.",
            })

        except Exception as e:
            log_rows.append({
                "species": species,
                "taxon_id": taxon_id,
                "query": "",
                "ncbi_count": "",
                "n_selected": 0,
                "status": "error",
                "message": repr(e),
            })

    candidates = pd.DataFrame(all_candidate_rows)
    selected = pd.DataFrame(selected_rows)
    log = pd.DataFrame(log_rows)

    candidates_file = METADATA_DIR / "helminths_no_latest_candidates.tsv"
    selected_file = METADATA_DIR / "helminths_no_latest_selected.tsv"
    log_file = METADATA_DIR / "helminths_no_latest_selection_log.tsv"

    candidates.to_csv(candidates_file, sep="\t", index=False)
    selected.to_csv(selected_file, sep="\t", index=False)
    log.to_csv(log_file, sep="\t", index=False)

    print("\nHelminth no-latest rescue completed.")
    print(f"Started: {started}")
    print(f"Finished: {datetime.now()}")
    print(f"Candidates: {candidates_file}")
    print(f"Selected: {selected_file}")
    print(f"Log: {log_file}")


if __name__ == "__main__":
    main()
