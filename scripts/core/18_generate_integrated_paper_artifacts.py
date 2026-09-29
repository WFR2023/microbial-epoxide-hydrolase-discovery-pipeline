#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

BASE_DIR = Path(".")

FINAL_CAND = BASE_DIR / "14_outputs" / "tables" / "final_priority_candidates.tsv"
PHYLO_TREE = BASE_DIR / "11_phylogeny" / "priority_candidates_tree.contree"
PRES_LONG = BASE_DIR / "15_remaining_context_analysis" / "presence_absence" / "within_species_presence_long.tsv"
NEIGH_COUNTS = BASE_DIR / "16_functional_context" / "tables" / "candidate_neighborhood_category_counts_long.tsv"
FUNC_MASTER = BASE_DIR / "16_functional_context" / "tables" / "candidate_functional_context_master.tsv"
ADV_MASTER = BASE_DIR / "17_optional_advanced_layers" / "tables" / "advanced_layers_master_status.tsv"

OUT_DIR = BASE_DIR / "18_integrated_artifacts"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"

for d in [OUT_DIR, FIG_DIR, TAB_DIR]:
    d.mkdir(parents=True, exist_ok=True)

def simplify_label(label: str) -> str:
    if str(label).startswith("human_EPHX2"):
        return "Human EPHX2"
    if str(label).startswith("human_EPHX1"):
        return "Human EPHX1"
    if str(label).startswith("pseudomonas_Cif"):
        return "Pseudomonas Cif"
    if "|" in str(label):
        acc, prot = str(label).split("|", 1)
        return f"{prot}\n({acc})"
    return str(label)

def save_plot(fig, outbase: Path):
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

# ------------------------------------------------------------
# Load input tables
# ------------------------------------------------------------
final_df = pd.read_csv(FINAL_CAND, sep="\t")
func_df = pd.read_csv(FUNC_MASTER, sep="\t")
pres_long_df = pd.read_csv(PRES_LONG, sep="\t")
neigh_counts_df = pd.read_csv(NEIGH_COUNTS, sep="\t")
adv_df = pd.read_csv(ADV_MASTER, sep="\t")

# ------------------------------------------------------------
# Build integrated master table robustly
# ------------------------------------------------------------
master = final_df.copy()

# Add functional context columns if present
func_keep = [c for c in [
    "qseqid", "functional_context_score", "functional_context_class"
] if c in func_df.columns]
if "qseqid" in func_keep:
    func_small = func_df[func_keep].drop_duplicates()
    master = master.merge(func_small, on="qseqid", how="left")

# Add advanced-layer status columns if present
adv_keep = [c for c in [
    "qseqid",
    "localization_layer",
    "formal_pathway_layer",
    "structural_modeling_layer",
    "docking_layer"
] if c in adv_df.columns]
if "qseqid" in adv_keep:
    adv_small = adv_df[adv_keep].drop_duplicates()
    master = master.merge(adv_small, on="qseqid", how="left")

# Presence summary
if {"query_qseqid", "target_assembly", "present_by_threshold"}.issubset(set(pres_long_df.columns)):
    pres_summary = (
        pres_long_df.groupby("query_qseqid")
        .agg(
            n_target_assemblies=("target_assembly", "nunique"),
            n_present_same_species=("present_by_threshold", "sum")
        )
        .reset_index()
        .rename(columns={"query_qseqid": "qseqid"})
    )
    master = master.merge(pres_summary, on="qseqid", how="left")
else:
    master["n_target_assemblies"] = 0
    master["n_present_same_species"] = 0

# Ensure required columns exist
for col in [
    "functional_context_score",
    "functional_context_class",
    "localization_layer",
    "formal_pathway_layer",
    "structural_modeling_layer",
    "docking_layer",
    "n_target_assemblies",
    "n_present_same_species"
]:
    if col not in master.columns:
        master[col] = 0 if col in ["functional_context_score", "n_target_assemblies", "n_present_same_species"] else "not_available"

master["n_target_assemblies"] = pd.to_numeric(master["n_target_assemblies"], errors="coerce").fillna(0)
master["n_present_same_species"] = pd.to_numeric(master["n_present_same_species"], errors="coerce").fillna(0)
master["functional_context_score"] = pd.to_numeric(master["functional_context_score"], errors="coerce").fillna(0)

# Save integrated master
master.to_csv(TAB_DIR / "integrated_candidate_master_table.tsv", sep="\t", index=False)
master.to_csv(TAB_DIR / "integrated_candidate_master_table.csv", index=False)

# ------------------------------------------------------------
# Figure 1: updated phylogeny
# ------------------------------------------------------------
tree = Phylo.read(str(PHYLO_TREE), "newick")
for clade in tree.find_clades():
    if clade.name:
        clade.name = simplify_label(clade.name)

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=ax, do_show=False, show_confidence=True)
ax.set_title("Figure 1. Maximum-likelihood phylogeny of prioritized candidates and references")
save_plot(fig, FIG_DIR / "Figure_1_updated_phylogeny")

# ------------------------------------------------------------
# Figure 2: integrated priority counts
# ------------------------------------------------------------
priority_counts = (
    master["integrated_priority"].astype(str).value_counts()
    .rename_axis("integrated_priority")
    .reset_index(name="n_candidates")
)

fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(priority_counts["integrated_priority"], priority_counts["n_candidates"])
ax.set_title("Figure 2. Integrated priority classes of prioritized candidates")
ax.set_ylabel("Number of candidates")
ax.set_xlabel("Integrated priority")
ax.tick_params(axis="x", rotation=30)
save_plot(fig, FIG_DIR / "Figure_2_integrated_priority_counts")

# ------------------------------------------------------------
# Figure 3: within-species presence
# ------------------------------------------------------------
plot_df = master.copy().sort_values(
    ["n_present_same_species", "bitscore"],
    ascending=[False, False]
)
labels = [x.split("|", 1)[1] if "|" in str(x) else str(x) for x in plot_df["qseqid"]]

fig, ax = plt.subplots(figsize=(12, 5))
ax.bar(range(len(plot_df)), plot_df["n_present_same_species"])
ax.set_xticks(range(len(plot_df)))
ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Assemblies with same-species homolog")
ax.set_xlabel("Prioritized candidate")
ax.set_title("Figure 3. Within-species distribution of prioritized candidates")
save_plot(fig, FIG_DIR / "Figure_3_within_species_presence")

# ------------------------------------------------------------
# Figure 4: neighborhood functional categories heatmap
# ------------------------------------------------------------
if {"qseqid", "functional_categories", "n_neighbors"}.issubset(set(neigh_counts_df.columns)):
    wide = (
        neigh_counts_df.pivot_table(
            index="qseqid",
            columns="functional_categories",
            values="n_neighbors",
            fill_value=0
        )
    )

    if not wide.empty:
        order = [q for q in master["qseqid"].tolist() if q in wide.index]
        wide = wide.loc[order]

        fig, ax = plt.subplots(figsize=(10, max(4, 0.45 * len(wide))))
        im = ax.imshow(wide.values, aspect="auto")
        ax.set_xticks(range(len(wide.columns)))
        ax.set_xticklabels(wide.columns, rotation=45, ha="right")
        ax.set_yticks(range(len(wide.index)))
        ax.set_yticklabels([x.split("|", 1)[1] if "|" in str(x) else str(x) for x in wide.index], fontsize=8)
        ax.set_title("Figure 4. Functional categories in candidate genomic neighborhoods")
        plt.colorbar(im, ax=ax, label="Neighbor count")
        save_plot(fig, FIG_DIR / "Figure_4_neighborhood_function_heatmap")

# ------------------------------------------------------------
# Figure 5: functional context score
# ------------------------------------------------------------
plot_df2 = master.copy().sort_values(
    ["functional_context_score", "bitscore"],
    ascending=[False, False]
)
labels2 = [x.split("|", 1)[1] if "|" in str(x) else str(x) for x in plot_df2["qseqid"]]

fig, ax = plt.subplots(figsize=(12, 5))
ax.bar(range(len(plot_df2)), plot_df2["functional_context_score"])
ax.set_xticks(range(len(plot_df2)))
ax.set_xticklabels(labels2, rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Functional context score")
ax.set_xlabel("Prioritized candidate")
ax.set_title("Figure 5. Functional-context support of prioritized candidates")
save_plot(fig, FIG_DIR / "Figure_5_functional_context_score")

# ------------------------------------------------------------
# Figure 6: advanced layer status summary
# ------------------------------------------------------------
adv_cols = [
    "localization_layer",
    "formal_pathway_layer",
    "structural_modeling_layer",
    "docking_layer"
]
adv_counts = []
for col in adv_cols:
    if col not in master.columns:
        master[col] = "not_executed_or_not_available"
    n_done = (master[col].fillna("not_executed_or_not_available") != "not_executed_or_not_available").sum()
    n_not = (master[col].fillna("not_executed_or_not_available") == "not_executed_or_not_available").sum()
    adv_counts.append({
        "layer": col,
        "executed_or_available": int(n_done),
        "not_executed_or_not_available": int(n_not)
    })

adv_counts_df = pd.DataFrame(adv_counts)

fig, ax = plt.subplots(figsize=(9, 5))
x = range(len(adv_counts_df))
ax.bar(x, adv_counts_df["not_executed_or_not_available"], label="Not executed / not available")
ax.bar(
    x,
    adv_counts_df["executed_or_available"],
    bottom=adv_counts_df["not_executed_or_not_available"],
    label="Executed / available"
)
ax.set_xticks(list(x))
ax.set_xticklabels(adv_counts_df["layer"], rotation=20, ha="right")
ax.set_ylabel("Count of prioritized candidates")
ax.set_title("Figure 6. Status of advanced optional analytical layers")
ax.legend()
save_plot(fig, FIG_DIR / "Figure_6_advanced_layer_status")

# ------------------------------------------------------------
# Excel workbook
# ------------------------------------------------------------
xlsx_path = TAB_DIR / "integrated_paper_tables.xlsx"
used_engine = None

for engine in ["xlsxwriter", "openpyxl"]:
    try:
        with pd.ExcelWriter(xlsx_path, engine=engine) as writer:
            master.to_excel(writer, sheet_name="Integrated_Master", index=False)
            final_df.to_excel(writer, sheet_name="Final_Candidates", index=False)
            func_df.to_excel(writer, sheet_name="Functional_Context", index=False)
            pres_long_df.to_excel(writer, sheet_name="Presence_Long", index=False)
            neigh_counts_df.to_excel(writer, sheet_name="Neighborhood_Categories", index=False)
            adv_df.to_excel(writer, sheet_name="Advanced_Layers", index=False)
            adv_counts_df.to_excel(writer, sheet_name="Advanced_Layer_Summary", index=False)
        used_engine = engine
        break
    except Exception:
        continue

summary = [
    "Integrated paper artifacts generated.",
    f"Figures directory: {FIG_DIR}",
    f"Tables directory: {TAB_DIR}",
    f"Workbook: {xlsx_path if xlsx_path.exists() else 'not created'}",
    f"Workbook engine: {used_engine if used_engine else 'none'}",
]
(TAB_DIR / "integrated_artifacts_summary.txt").write_text("\n".join(summary), encoding="utf-8")
print("\n".join(summary))
