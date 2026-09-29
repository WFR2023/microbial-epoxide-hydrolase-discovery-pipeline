#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO

BASE_DIR = Path("<PROJECT_ROOT>")

# Inputs already produced in the project
SELECTED_ASM = BASE_DIR / "01_metadata" / "selected_assemblies.tsv"
PRIORITY_ANN = BASE_DIR / "08_priority_candidates" / "priority_candidates_final_annotated.tsv"
PRIORITY_FASTA = BASE_DIR / "08_priority_candidates" / "priority_candidates_final.faa"
INTEGRATED_EVIDENCE = BASE_DIR / "12_domain_motif_screen" / "priority_candidates_integrated_evidence.tsv"

# Main source of downloaded assemblies
DIRECT_FILES = BASE_DIR / "03_downloads" / "direct_ncbi_files"

# Derived protein files
PROTEIN_DIR = BASE_DIR / "04_derived" / "proteins_final"

OUT_DIR = BASE_DIR / "15_remaining_context_analysis"
NEIGH_DIR = OUT_DIR / "neighborhoods"
PRES_DIR = OUT_DIR / "presence_absence"
LOC_DIR = OUT_DIR / "localization"
FIG_DIR = OUT_DIR / "figures"

for d in [OUT_DIR, NEIGH_DIR, PRES_DIR, LOC_DIR, FIG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

MASTER_TABLE = OUT_DIR / "candidate_master_table.tsv"
NEIGH_TABLE = NEIGH_DIR / "candidate_neighborhood_table.tsv"
NEIGH_SUMMARY = NEIGH_DIR / "candidate_neighborhood_summary.tsv"
PRES_LONG = PRES_DIR / "within_species_presence_long.tsv"
PRES_MATRIX = PRES_DIR / "within_species_presence_matrix.tsv"
LOC_STATUS = LOC_DIR / "localization_status.txt"
PRES_FIG = FIG_DIR / "Figure_remaining_presence_counts.png"
PRES_FIG_PDF = FIG_DIR / "Figure_remaining_presence_counts.pdf"


# ------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------
def load_data():
    selected = pd.read_csv(SELECTED_ASM, sep="\t")
    priority = pd.read_csv(PRIORITY_ANN, sep="\t")
    integ = pd.read_csv(INTEGRATED_EVIDENCE, sep="\t")

    # Add assembly and protein ids
    if "assembly_accession" not in priority.columns or "protein_id" not in priority.columns:
        priority[["assembly_accession", "protein_id"]] = priority["qseqid"].str.split("|", n=1, expand=True)

    selected_small = selected[[
        "assembly_accession", "target_group", "target_species", "source_database",
        "assembly_name", "assembly_level", "refseq_category"
    ]].drop_duplicates()

    integ_small = integ[[
        "qseqid", "integrated_priority", "heuristic_class"
    ]].drop_duplicates()

    master = (
        priority.merge(selected_small, on="assembly_accession", how="left")
                .merge(integ_small, on="qseqid", how="left")
    )

    master.to_csv(MASTER_TABLE, sep="\t", index=False)
    return selected, master


def find_candidate_seq_records():
    return {rec.id: rec for rec in SeqIO.parse(PRIORITY_FASTA, "fasta")}


def parse_gff_attributes(attr_text):
    attrs = {}
    for item in attr_text.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def find_annotation_files(assembly_accession):
    """
    Search recursively inside direct_ncbi_files for GFF/GBFF corresponding to assembly accession.
    """
    gff_matches = list(DIRECT_FILES.rglob(f"{assembly_accession}*genomic.gff")) + \
                  list(DIRECT_FILES.rglob(f"{assembly_accession}*genomic.gff.gz"))
    gbff_matches = list(DIRECT_FILES.rglob(f"{assembly_accession}*genomic.gbff")) + \
                   list(DIRECT_FILES.rglob(f"{assembly_accession}*genomic.gbff.gz"))

    # Prefer uncompressed if present
    gff = next((p for p in gff_matches if p.suffix != ".gz"), gff_matches[0] if gff_matches else None)
    gbff = next((p for p in gbff_matches if p.suffix != ".gz"), gbff_matches[0] if gbff_matches else None)
    return gff, gbff


def parse_gff_cds_features(gff_path):
    features = []
    if gff_path is None or not gff_path.exists():
        return features

    opener = open
    if gff_path.suffix == ".gz":
        import gzip
        opener = gzip.open

    with opener(gff_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype != "CDS":
                continue
            ad = parse_gff_attributes(attrs)
            features.append({
                "seqid": seqid,
                "feature_type": ftype,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "protein_id": ad.get("protein_id", ""),
                "locus_tag": ad.get("locus_tag", ""),
                "gene": ad.get("gene", ""),
                "product": ad.get("product", ""),
                "raw_attributes": attrs
            })
    features.sort(key=lambda x: (x["seqid"], x["start"], x["end"]))
    return features


def find_focal_feature_index(features, protein_id):
    for i, feat in enumerate(features):
        if feat["protein_id"] == protein_id:
            return i
    # fallback looser search
    for i, feat in enumerate(features):
        if protein_id in feat["protein_id"]:
            return i
    return None


def build_neighborhood_rows(master):
    rows = []
    summary = []

    for _, row in master.iterrows():
        asm = row["assembly_accession"]
        pid = row["protein_id"]
        qseqid = row["qseqid"]

        gff_path, gbff_path = find_annotation_files(asm)
        features = parse_gff_cds_features(gff_path)

        if not features:
            summary.append({
                "qseqid": qseqid,
                "assembly_accession": asm,
                "protein_id": pid,
                "target_species": row.get("target_species", ""),
                "status": "no_gff_features",
                "gff_path": str(gff_path) if gff_path else "",
                "gbff_path": str(gbff_path) if gbff_path else ""
            })
            continue

        idx = find_focal_feature_index(features, pid)
        if idx is None:
            summary.append({
                "qseqid": qseqid,
                "assembly_accession": asm,
                "protein_id": pid,
                "target_species": row.get("target_species", ""),
                "status": "protein_not_found_in_gff",
                "gff_path": str(gff_path) if gff_path else "",
                "gbff_path": str(gbff_path) if gbff_path else ""
            })
            continue

        focal = features[idx]
        seqid = focal["seqid"]
        seq_features = [f for f in features if f["seqid"] == seqid]
        seq_positions = [i for i, f in enumerate(seq_features) if f["protein_id"] == focal["protein_id"]]
        if not seq_positions:
            summary.append({
                "qseqid": qseqid,
                "assembly_accession": asm,
                "protein_id": pid,
                "target_species": row.get("target_species", ""),
                "status": "protein_not_found_on_seqid",
                "gff_path": str(gff_path) if gff_path else "",
                "gbff_path": str(gbff_path) if gbff_path else ""
            })
            continue

        focal_seq_idx = seq_positions[0]
        left = max(0, focal_seq_idx - 5)
        right = min(len(seq_features), focal_seq_idx + 6)
        window = seq_features[left:right]

        summary.append({
            "qseqid": qseqid,
            "assembly_accession": asm,
            "protein_id": pid,
            "target_species": row.get("target_species", ""),
            "status": "ok",
            "gff_path": str(gff_path) if gff_path else "",
            "gbff_path": str(gbff_path) if gbff_path else "",
            "seqid": seqid,
            "focal_start": focal["start"],
            "focal_end": focal["end"],
            "n_window_features": len(window)
        })

        for j, feat in enumerate(window):
            rows.append({
                "qseqid": qseqid,
                "assembly_accession": asm,
                "protein_id": pid,
                "target_species": row.get("target_species", ""),
                "best_reference_hit": row.get("sseqid", ""),
                "integrated_priority": row.get("integrated_priority", ""),
                "seqid": seqid,
                "window_rank": j - (focal_seq_idx - left),
                "is_focal": int(feat["protein_id"] == pid),
                "neighbor_protein_id": feat["protein_id"],
                "neighbor_locus_tag": feat["locus_tag"],
                "neighbor_gene": feat["gene"],
                "neighbor_product": feat["product"],
                "start": feat["start"],
                "end": feat["end"],
                "strand": feat["strand"],
                "raw_attributes": feat["raw_attributes"]
            })

    neigh_df = pd.DataFrame(rows)
    summ_df = pd.DataFrame(summary)
    neigh_df.to_csv(NEIGH_TABLE, sep="\t", index=False)
    summ_df.to_csv(NEIGH_SUMMARY, sep="\t", index=False)
    return neigh_df, summ_df


def plot_neighborhoods(neigh_df):
    if neigh_df.empty:
        return

    for qseqid, sub in neigh_df.groupby("qseqid"):
        sub = sub.sort_values("start")
        fig, ax = plt.subplots(figsize=(12, 2.6))

        ymin, height = 0.35, 0.3
        xmin = sub["start"].min()
        xmax = sub["end"].max()

        for _, r in sub.iterrows():
            color = "crimson" if r["is_focal"] == 1 else "steelblue"
            width = r["end"] - r["start"] + 1
            ax.add_patch(plt.Rectangle((r["start"], ymin), width, height))
            ax.text(
                (r["start"] + r["end"]) / 2,
                ymin + height + 0.06,
                r["neighbor_protein_id"] if r["neighbor_protein_id"] else r["neighbor_locus_tag"],
                ha="center",
                va="bottom",
                fontsize=6,
                rotation=45
            )

        ax.set_xlim(xmin - 100, xmax + 100)
        ax.set_ylim(0, 1.2)
        ax.set_yticks([])
        ax.set_xlabel("Genomic coordinate")
        ax.set_title(f"Neighborhood: {qseqid}")
        for spine in ["left", "right", "top"]:
            ax.spines[spine].set_visible(False)

        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", qseqid)
        fig.tight_layout()
        fig.savefig(NEIGH_DIR / f"{safe_name}_neighborhood.png", dpi=300, bbox_inches="tight")
        fig.savefig(NEIGH_DIR / f"{safe_name}_neighborhood.pdf", bbox_inches="tight")
        plt.close(fig)


def run_local_blast_presence(master, seq_records):
    """
    For each priority candidate, search same-species assemblies using local blastp -subject.
    Presence rule:
      evalue <= 1e-10
      pident >= 30
      qcov >= 70
    """
    if shutil.which("blastp") is None:
        raise RuntimeError("blastp not found in PATH")

    long_rows = []

    # Build species -> assemblies mapping from selected metadata
    selected = pd.read_csv(SELECTED_ASM, sep="\t")
    species_map = (
        selected[["assembly_accession", "target_species"]]
        .drop_duplicates()
        .groupby("target_species")["assembly_accession"].apply(list).to_dict()
    )

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        for _, row in master.iterrows():
            qseqid = row["qseqid"]
            asm = row["assembly_accession"]
            species = row["target_species"]
            rec = seq_records.get(qseqid)
            if rec is None:
                continue

            query_fa = td / "query.faa"
            SeqIO.write([rec], query_fa, "fasta")

            for target_asm in species_map.get(species, []):
                subj_fa = PROTEIN_DIR / f"{target_asm}.faa"
                if not subj_fa.exists():
                    long_rows.append({
                        "query_qseqid": qseqid,
                        "query_species": species,
                        "query_source_assembly": asm,
                        "target_assembly": target_asm,
                        "status": "missing_subject_proteins",
                        "best_hit_subject_id": "",
                        "pident": "",
                        "length": "",
                        "evalue": "",
                        "bitscore": "",
                        "qlen": "",
                        "qcov": "",
                        "present_by_threshold": 0
                    })
                    continue

                out_tsv = td / "blast.tsv"
                cmd = [
                    "blastp",
                    "-query", str(query_fa),
                    "-subject", str(subj_fa),
                    "-evalue", "1e-5",
                    "-max_target_seqs", "1",
                    "-outfmt",
                    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
                ]

                with open(out_tsv, "w") as fout:
                    subprocess.run(cmd, check=False, stdout=fout, stderr=subprocess.DEVNULL)

                if out_tsv.stat().st_size == 0:
                    long_rows.append({
                        "query_qseqid": qseqid,
                        "query_species": species,
                        "query_source_assembly": asm,
                        "target_assembly": target_asm,
                        "status": "no_hit",
                        "best_hit_subject_id": "",
                        "pident": "",
                        "length": "",
                        "evalue": "",
                        "bitscore": "",
                        "qlen": "",
                        "qcov": "",
                        "present_by_threshold": 0
                    })
                    continue

                cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]
                blast_df = pd.read_csv(out_tsv, sep="\t", header=None, names=cols)
                best = blast_df.iloc[0]
                qcov = (float(best["length"]) / float(best["qlen"])) * 100.0 if float(best["qlen"]) > 0 else 0.0
                present = int(
                    (float(best["evalue"]) <= 1e-10) and
                    (float(best["pident"]) >= 30.0) and
                    (qcov >= 70.0)
                )

                long_rows.append({
                    "query_qseqid": qseqid,
                    "query_species": species,
                    "query_source_assembly": asm,
                    "target_assembly": target_asm,
                    "status": "hit",
                    "best_hit_subject_id": best["sseqid"],
                    "pident": best["pident"],
                    "length": best["length"],
                    "evalue": best["evalue"],
                    "bitscore": best["bitscore"],
                    "qlen": best["qlen"],
                    "qcov": round(qcov, 3),
                    "present_by_threshold": present
                })

    long_df = pd.DataFrame(long_rows)
    long_df.to_csv(PRES_LONG, sep="\t", index=False)

    if long_df.empty:
        pd.DataFrame().to_csv(PRES_MATRIX, sep="\t", index=False)
        return long_df, pd.DataFrame()

    matrix = long_df.pivot_table(
        index="query_qseqid",
        columns="target_assembly",
        values="present_by_threshold",
        aggfunc="max",
        fill_value=0
    ).reset_index()
    matrix.to_csv(PRES_MATRIX, sep="\t", index=False)

    return long_df, matrix


def plot_presence_summary(long_df, master):
    if long_df.empty:
        return

    counts = (
        long_df.groupby("query_qseqid")["present_by_threshold"]
        .sum()
        .reset_index(name="n_same_species_assemblies_present")
    )
    counts = counts.merge(master[["qseqid", "integrated_priority"]].drop_duplicates(),
                          left_on="query_qseqid", right_on="qseqid", how="left")

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.bar(range(len(counts)), counts["n_same_species_assemblies_present"])
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels([q.split("|", 1)[1] if "|" in q else q for q in counts["query_qseqid"]],
                       rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Assemblies with within-species homolog above threshold")
    ax.set_title("Within-species presence of prioritized candidates")
    fig.tight_layout()
    fig.savefig(PRES_FIG, dpi=300, bbox_inches="tight")
    fig.savefig(PRES_FIG_PDF, bbox_inches="tight")
    plt.close(fig)


def run_optional_localization():
    lines = []
    candidates = []

    # Check tools
    signalp_bin = shutil.which("signalp6") or shutil.which("signalp")
    tm_bin = shutil.which("deeptmhmm") or shutil.which("tmhmm") or shutil.which("DeepTMHMM")
    psortb_bin = shutil.which("psortb")
    targetp_bin = shutil.which("targetp")
    wolfpsort_bin = shutil.which("wolfpsort")

    lines.append("Optional localization scan status")
    lines.append("================================")
    lines.append(f"signalp: {signalp_bin if signalp_bin else 'NOT_FOUND'}")
    lines.append(f"tm predictor: {tm_bin if tm_bin else 'NOT_FOUND'}")
    lines.append(f"psortb: {psortb_bin if psortb_bin else 'NOT_FOUND'}")
    lines.append(f"targetp: {targetp_bin if targetp_bin else 'NOT_FOUND'}")
    lines.append(f"wolfpsort: {wolfpsort_bin if wolfpsort_bin else 'NOT_FOUND'}")
    lines.append("")
    lines.append("No automatic localization results were generated unless the required external tools were available in PATH.")
    lines.append("This status file documents whether the optional localization layer can be executed on the current system.")

    LOC_STATUS.write_text("\n".join(lines), encoding="utf-8")


def main():
    print("Loading metadata and candidate tables...")
    selected, master = load_data()
    seq_records = find_candidate_seq_records()

    print("Building genomic neighborhood tables...")
    neigh_df, neigh_summary = build_neighborhood_rows(master)
    plot_neighborhoods(neigh_df)

    print("Running within-species presence/absence analysis...")
    long_df, matrix_df = run_local_blast_presence(master, seq_records)
    plot_presence_summary(long_df, master)

    print("Checking optional localization tools...")
    run_optional_localization()

    print("\nRemaining context analysis completed.")
    print(f"Master table: {MASTER_TABLE}")
    print(f"Neighborhood table: {NEIGH_TABLE}")
    print(f"Neighborhood summary: {NEIGH_SUMMARY}")
    print(f"Presence long table: {PRES_LONG}")
    print(f"Presence matrix: {PRES_MATRIX}")
    print(f"Localization status: {LOC_STATUS}")
    print(f"Figures: {FIG_DIR}")


if __name__ == "__main__":
    main()
