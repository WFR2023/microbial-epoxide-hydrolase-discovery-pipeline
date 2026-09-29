#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import math
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo

BASE_DIR = Path("<PROJECT_ROOT>")
OUT_DIR = BASE_DIR / "14_outputs"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

FILES = {
    "priority_final": BASE_DIR / "07_blast_results" / "filtered" / "priority_candidates_final.tsv",
    "integrated": BASE_DIR / "12_domain_motif_screen" / "priority_candidates_integrated_evidence.tsv",
    "desc": BASE_DIR / "09_candidate_description" / "priority_candidates_descriptive_table.tsv",
    "pfam_best": BASE_DIR / "13_hmmer_pfam" / "priority_candidates_vs_pfam_best_hits.tsv",
    "pfam_summary": BASE_DIR / "13_hmmer_pfam" / "priority_candidates_vs_pfam_summary.tsv",
    "tree_contree": BASE_DIR / "11_phylogeny" / "priority_candidates_tree.contree",
    "tree_iqtree": BASE_DIR / "11_phylogeny" / "priority_candidates_tree.iqtree",
}

for name, path in FILES.items():
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")

priority_df = pd.read_csv(FILES["priority_final"], sep="\t")
integrated_df = pd.read_csv(FILES["integrated"], sep="\t")
desc_df = pd.read_csv(FILES["desc"], sep="\t")
pfam_best_df = pd.read_csv(FILES["pfam_best"], sep="\t")
pfam_summary_df = pd.read_csv(FILES["pfam_summary"], sep="\t")

# ----------------------------
# Helpers
# ----------------------------
def simplify_label(label: str) -> str:
    if label.startswith("human_EPHX2"):
        return "Human EPHX2"
    if label.startswith("human_EPHX1"):
        return "Human EPHX1"
    if label.startswith("pseudomonas_Cif"):
        return "Pseudomonas Cif"
    if "|" in label:
        acc, prot = label.split("|", 1)
        return f"{prot}\n({acc})"
    return label

def save_plot(fig, path_base: Path):
    fig.tight_layout()
    fig.savefig(path_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

# ----------------------------
# Build final integrated table
# ----------------------------
keep_cols_desc = [
    "assembly_accession", "protein_id", "qseqid", "best_reference_hit", "pident",
    "blast_aln_length", "evalue", "bitscore", "priority_tier",
    "sequence_length_aa", "molecular_weight", "isoelectric_point",
    "aromaticity", "instability_index", "gravy"
]
desc_small = desc_df[keep_cols_desc].copy()

keep_cols_int = [
    "qseqid", "integrated_priority", "heuristic_class"
]
int_small = integrated_df[keep_cols_int].copy()

keep_cols_pfam = [
    "query_name", "target_name", "target_accession", "full_evalue",
    "full_score", "description"
]
pfam_small = pfam_best_df[keep_cols_pfam].copy().rename(columns={
    "query_name": "qseqid",
    "target_name": "pfam_best_domain",
    "target_accession": "pfam_accession",
    "full_evalue": "pfam_evalue",
    "full_score": "pfam_score",
    "description": "pfam_description",
})

final_df = (
    desc_small
    .merge(int_small, on="qseqid", how="left")
    .merge(pfam_small, on="qseqid", how="left")
)

priority_order = {
    "high_priority_cif_like": 1,
    "high_priority_ephx1_like": 2,
    "moderate_priority_cif_related": 3,
    "moderate_priority_ephx2_related": 4,
    "exploratory": 5
}
final_df["priority_rank"] = final_df["integrated_priority"].map(priority_order).fillna(99)
final_df = final_df.sort_values(
    by=["priority_rank", "bitscore", "evalue"],
    ascending=[True, False, True]
).drop(columns=["priority_rank"])

final_tsv = TAB_DIR / "final_priority_candidates.tsv"
final_csv = TAB_DIR / "final_priority_candidates.csv"
final_df.to_csv(final_tsv, sep="\t", index=False)
final_df.to_csv(final_csv, index=False)

# ----------------------------
# Figure 1: phylogenetic tree
# ----------------------------
tree = Phylo.read(str(FILES["tree_contree"]), "newick")
for clade in tree.find_clades():
    if clade.name:
        clade.name = simplify_label(clade.name)

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=ax, do_show=False, show_confidence=True)
ax.set_title("Maximum-likelihood phylogeny of prioritized candidates and curated references")
save_plot(fig, FIG_DIR / "Figure_1_phylogenetic_tree")

# ----------------------------
# Figure 2: integrated priority counts
# ----------------------------
priority_counts = (
    final_df["integrated_priority"]
    .value_counts()
    .rename_axis("integrated_priority")
    .reset_index(name="n_candidates")
)

fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(priority_counts["integrated_priority"], priority_counts["n_candidates"])
ax.set_title("Integrated priority classes of final candidates")
ax.set_ylabel("Number of candidates")
ax.set_xlabel("Integrated priority")
ax.tick_params(axis="x", rotation=30)
save_plot(fig, FIG_DIR / "Figure_2_integrated_priority_counts")

# ----------------------------
# Figure 3: BLAST support scatter
# ----------------------------
fig, ax = plt.subplots(figsize=(8, 6))
for label, subdf in final_df.groupby("integrated_priority"):
    ax.scatter(
        subdf["pident"],
        subdf["bitscore"],
        label=label,
        s=60
    )
for _, row in final_df.iterrows():
    ax.annotate(
        row["protein_id"],
        (row["pident"], row["bitscore"]),
        fontsize=7,
        alpha=0.85
    )
ax.set_title("BLAST support of prioritized candidates")
ax.set_xlabel("Percent identity")
ax.set_ylabel("Bitscore")
ax.legend(fontsize=8)
save_plot(fig, FIG_DIR / "Figure_3_blast_support_scatter")

# ----------------------------
# Figure 4: Pfam best-hit summary
# ----------------------------
fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(pfam_summary_df["target_name"], pfam_summary_df["n_queries_best_assigned"])
ax.set_title("Best Pfam domain assignments among prioritized sequences")
ax.set_xlabel("Pfam domain")
ax.set_ylabel("Number of sequences")
ax.tick_params(axis="x", rotation=25)
save_plot(fig, FIG_DIR / "Figure_4_pfam_best_hits")

# ----------------------------
# Methods/phylogeny metadata
# ----------------------------
model_line = ""
loglik_line = ""
with open(FILES["tree_iqtree"], "r", encoding="utf-8") as f:
    for line in f:
        if "Best-fit model according to BIC:" in line:
            model_line = line.strip()
        if "Log-likelihood of the tree:" in line:
            loglik_line = line.strip()

meta_df = pd.DataFrame([
    {"parameter": "best_fit_model", "value": model_line.replace("Best-fit model according to BIC:", "").strip()},
    {"parameter": "log_likelihood", "value": loglik_line.replace("Log-likelihood of the tree:", "").strip()},
])

meta_df.to_csv(TAB_DIR / "phylogeny_metadata.tsv", sep="\t", index=False)

# ----------------------------
# Excel workbook
# ----------------------------
xlsx_path = TAB_DIR / "paper_tables.xlsx"
xlsx_written = False
xlsx_error = ""

for engine in ["xlsxwriter", "openpyxl"]:
    try:
        with pd.ExcelWriter(xlsx_path, engine=engine) as writer:
            final_df.to_excel(writer, sheet_name="Final_Candidates", index=False)
            desc_df.to_excel(writer, sheet_name="Descriptive_Table", index=False)
            integrated_df.to_excel(writer, sheet_name="Integrated_Evidence", index=False)
            priority_df.to_excel(writer, sheet_name="BLAST_Filtered", index=False)
            pfam_best_df.to_excel(writer, sheet_name="Pfam_Best_Hits", index=False)
            pfam_summary_df.to_excel(writer, sheet_name="Pfam_Summary", index=False)
            meta_df.to_excel(writer, sheet_name="Phylogeny_Metadata", index=False)
        xlsx_written = True
        break
    except Exception as e:
        xlsx_error = f"{engine}: {e}"

summary_lines = [
    f"Final candidate table: {final_tsv}",
    f"Figures directory: {FIG_DIR}",
    f"Tables directory: {TAB_DIR}",
    f"Excel workbook written: {xlsx_written}",
]
if not xlsx_written:
    summary_lines.append(f"Excel error: {xlsx_error}")
    summary_lines.append("CSV/TSV tables were still exported successfully.")

summary_file = OUT_DIR / "artifact_generation_summary.txt"
summary_file.write_text("\n".join(summary_lines), encoding="utf-8")

print("\n".join(summary_lines))
