#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict, Counter

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = Path("<PROJECT_ROOT>")

MASTER = BASE_DIR / "15_remaining_context_analysis" / "candidate_master_table.tsv"
NEIGH = BASE_DIR / "15_remaining_context_analysis" / "neighborhoods" / "candidate_neighborhood_table.tsv"
NEIGH_SUM = BASE_DIR / "15_remaining_context_analysis" / "neighborhoods" / "candidate_neighborhood_summary.tsv"
PRES_LONG = BASE_DIR / "15_remaining_context_analysis" / "presence_absence" / "within_species_presence_long.tsv"
FINAL_CAND = BASE_DIR / "14_outputs" / "tables" / "final_priority_candidates.tsv"
PRIORITY_FASTA = BASE_DIR / "08_priority_candidates" / "priority_candidates_final.faa"

OUT_DIR = BASE_DIR / "16_functional_context"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"
OPT_DIR = OUT_DIR / "optional_tools"

for d in [OUT_DIR, FIG_DIR, TAB_DIR, OPT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

KEYWORD_CATEGORIES = {
    "transport": [
        "transporter", "permease", "abc", "efflux", "porin", "channel", "import", "export"
    ],
    "lipid_metabolism": [
        "lipid", "fatty acid", "acyl", "phospholipid", "membrane lipid", "esterase", "lipase"
    ],
    "oxidative_stress": [
        "oxidative", "peroxidase", "superoxide", "catalase", "glutathione", "thioredoxin", "redox"
    ],
    "detoxification_xenobiotic": [
        "xenobiotic", "detox", "dehalogenase", "drug", "multidrug", "toxin", "resistance"
    ],
    "regulation": [
        "regulator", "transcriptional", "sigma factor", "two-component", "repressor", "activator"
    ],
    "virulence_host_interaction": [
        "virulence", "pathogenic", "adhesin", "invasion", "secretion", "toxin", "host"
    ],
    "mobile_element": [
        "transposase", "integrase", "recombinase", "phage", "insertion sequence", "mobile element"
    ],
    "metabolism_general": [
        "dehydrogenase", "synthetase", "transferase", "isomerase", "kinase", "hydrolase", "metabolic"
    ]
}

def categorize_product(text: str):
    text = str(text).lower()
    hits = []
    for cat, kws in KEYWORD_CATEGORIES.items():
        for kw in kws:
            if kw in text:
                hits.append(cat)
                break
    if not hits:
        hits = ["unclassified"]
    return hits

def make_candidate_function_table():
    final_df = pd.read_csv(FINAL_CAND, sep="\t")
    neigh_df = pd.read_csv(NEIGH, sep="\t")
    pres_df = pd.read_csv(PRES_LONG, sep="\t")

    # neighborhood categorization
    neigh_df["neighbor_product"] = neigh_df["neighbor_product"].fillna("")
    neigh_df["functional_categories"] = neigh_df["neighbor_product"].apply(categorize_product)
    neigh_exploded = neigh_df.explode("functional_categories")

    # summarize neighborhood categories per candidate
    neigh_cat_counts = (
        neigh_exploded.groupby(["qseqid", "functional_categories"])
        .size()
        .reset_index(name="n_neighbors")
    )

    neigh_summary_wide = (
        neigh_cat_counts.pivot_table(
            index="qseqid",
            columns="functional_categories",
            values="n_neighbors",
            fill_value=0
        )
        .reset_index()
    )

    # focal vs neighbors summary
    neigh_stats = (
        neigh_df.groupby("qseqid")
        .agg(
            n_neighbor_features=("neighbor_protein_id", "count"),
            n_focal_records=("is_focal", "sum")
        )
        .reset_index()
    )

    # within-species presence summary
    pres_ok = pres_df.copy()
    if "present_by_threshold" in pres_ok.columns:
        pres_summary = (
            pres_ok.groupby("query_qseqid")
            .agg(
                n_target_assemblies=("target_assembly", "nunique"),
                n_present_same_species=("present_by_threshold", "sum")
            )
            .reset_index()
            .rename(columns={"query_qseqid": "qseqid"})
        )
    else:
        pres_summary = pd.DataFrame(columns=["qseqid", "n_target_assemblies", "n_present_same_species"])

    merged = (
        final_df.merge(neigh_stats, on="qseqid", how="left")
                .merge(neigh_summary_wide, on="qseqid", how="left")
                .merge(pres_summary, on="qseqid", how="left")
    )

    merged = merged.fillna(0)

    # heuristic contextual score
    def score_row(r):
        score = 0
        if str(r.get("integrated_priority", "")).startswith("high_priority"):
            score += 3
        elif str(r.get("integrated_priority", "")).startswith("moderate_priority"):
            score += 2
        else:
            score += 1

        if r.get("n_present_same_species", 0) >= 3:
            score += 2
        elif r.get("n_present_same_species", 0) >= 1:
            score += 1

        for c in ["transport", "lipid_metabolism", "oxidative_stress",
                  "detoxification_xenobiotic", "virulence_host_interaction", "regulation"]:
            if r.get(c, 0) > 0:
                score += 1
        return score

    merged["functional_context_score"] = merged.apply(score_row, axis=1)

    def score_class(x):
        if x >= 8:
            return "high_context_support"
        elif x >= 5:
            return "moderate_context_support"
        else:
            return "limited_context_support"

    merged["functional_context_class"] = merged["functional_context_score"].apply(score_class)

    merged.to_csv(TAB_DIR / "candidate_functional_context_master.tsv", sep="\t", index=False)
    neigh_cat_counts.to_csv(TAB_DIR / "candidate_neighborhood_category_counts_long.tsv", sep="\t", index=False)
    neigh_exploded.to_csv(TAB_DIR / "candidate_neighborhood_categorized_long.tsv", sep="\t", index=False)

    return merged, neigh_cat_counts, neigh_exploded

def run_optional_eggnog():
    emapper = shutil.which("emapper.py") or shutil.which("emapper")
    status_lines = []
    status_lines.append("Optional eggNOG layer")
    status_lines.append("====================")

    if not emapper:
        status_lines.append("emapper.py/emapper not found in PATH")
        (OPT_DIR / "eggnog_status.txt").write_text("\n".join(status_lines), encoding="utf-8")
        return

    status_lines.append(f"Found: {emapper}")
    out_prefix = OPT_DIR / "eggnog_priority"

    cmd = [
        emapper,
        "-i", str(PRIORITY_FASTA),
        "--output", str(out_prefix.name),
        "--output_dir", str(OPT_DIR),
        "--itype", "proteins",
        "--cpu", "2"
    ]

    try:
        subprocess.run(cmd, check=True)
        status_lines.append("eggNOG-mapper completed successfully.")
    except Exception as e:
        status_lines.append(f"eggNOG-mapper failed: {e}")

    (OPT_DIR / "eggnog_status.txt").write_text("\n".join(status_lines), encoding="utf-8")

def plot_functional_context(master_df, neigh_cat_counts):
    # Figure 1: context score
    plot_df = master_df.sort_values(["functional_context_score", "bitscore"], ascending=[False, False]).copy()
    labels = [x.split("|", 1)[1] if "|" in x else x for x in plot_df["qseqid"]]

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(range(len(plot_df)), plot_df["functional_context_score"])
    ax.set_xticks(range(len(plot_df)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Functional context score")
    ax.set_title("Functional-context support of prioritized candidates")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "Figure_5_functional_context_score.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "Figure_5_functional_context_score.pdf", bbox_inches="tight")
    plt.close(fig)

    # Figure 2: neighborhood category heatmap-like matrix
    wide = (
        neigh_cat_counts.pivot_table(
            index="qseqid",
            columns="functional_categories",
            values="n_neighbors",
            fill_value=0
        )
    )
    if not wide.empty:
        wide = wide.loc[plot_df["qseqid"].tolist()]
        fig, ax = plt.subplots(figsize=(10, max(4, 0.45 * len(wide))))
        im = ax.imshow(wide.values, aspect="auto")
        ax.set_xticks(range(len(wide.columns)))
        ax.set_xticklabels(wide.columns, rotation=45, ha="right")
        ax.set_yticks(range(len(wide.index)))
        ax.set_yticklabels([x.split("|", 1)[1] if "|" in x else x for x in wide.index], fontsize=8)
        ax.set_title("Genomic neighborhood functional categories")
        plt.colorbar(im, ax=ax, label="Neighbor count")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "Figure_6_neighborhood_function_heatmap.png", dpi=300, bbox_inches="tight")
        fig.savefig(FIG_DIR / "Figure_6_neighborhood_function_heatmap.pdf", bbox_inches="tight")
        plt.close(fig)

def export_xlsx():
    xlsx = TAB_DIR / "functional_context_tables.xlsx"
    master = pd.read_csv(TAB_DIR / "candidate_functional_context_master.tsv", sep="\t")
    long1 = pd.read_csv(TAB_DIR / "candidate_neighborhood_category_counts_long.tsv", sep="\t")
    long2 = pd.read_csv(TAB_DIR / "candidate_neighborhood_categorized_long.tsv", sep="\t")

    for engine in ["xlsxwriter", "openpyxl"]:
        try:
            with pd.ExcelWriter(xlsx, engine=engine) as writer:
                master.to_excel(writer, sheet_name="Functional_Context_Master", index=False)
                long1.to_excel(writer, sheet_name="Neighborhood_Cat_Counts", index=False)
                long2.to_excel(writer, sheet_name="Neighborhood_Categorized", index=False)
            return xlsx, engine
        except Exception:
            continue
    return None, None

def main():
    print("Building functional-context master tables...")
    master_df, neigh_cat_counts, neigh_exploded = make_candidate_function_table()

    print("Running optional eggNOG layer if available...")
    run_optional_eggnog()

    print("Generating figures...")
    plot_functional_context(master_df, neigh_cat_counts)

    print("Exporting workbook...")
    xlsx, engine = export_xlsx()

    summary = [
        "Functional context layer completed.",
        f"Master table: {TAB_DIR / 'candidate_functional_context_master.tsv'}",
        f"Neighborhood counts: {TAB_DIR / 'candidate_neighborhood_category_counts_long.tsv'}",
        f"Categorized neighborhood table: {TAB_DIR / 'candidate_neighborhood_categorized_long.tsv'}",
        f"Figures directory: {FIG_DIR}",
        f"Excel workbook: {xlsx if xlsx else 'not created'}",
        f"Excel engine: {engine if engine else 'none'}"
    ]
    (OUT_DIR / "functional_context_summary.txt").write_text("\n".join(summary), encoding="utf-8")
    print("\n".join(summary))

if __name__ == "__main__":
    main()
