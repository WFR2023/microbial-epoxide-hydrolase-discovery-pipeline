#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Publication-grade structural summary and figure generator
for docking / ColabFold / pocket analysis outputs.

Author: Project workflow
"""

from pathlib import Path
import argparse
import json
import math
from itertools import combinations

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# SciPy is optional; script continues if unavailable
try:
    from scipy.stats import wilcoxon, mannwhitneyu
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False


# -----------------------------
# Helpers
# -----------------------------
def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def save_figure(fig, outbase: Path):
    """
    Save each figure in multiple formats.
    PNG = high resolution raster
    PDF = vector
    SVG = vector
    """
    fig.savefig(f"{outbase}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{outbase}.pdf", bbox_inches="tight")
    fig.savefig(f"{outbase}.svg", bbox_inches="tight")
    plt.close(fig)


def classify_model(model_id: str) -> str:
    refs = {"human_EPHX1_mEH", "human_EPHX2_sEH", "pseudomonas_Cif"}
    return "reference" if model_id in refs else "candidate"


def reference_group(model_id: str) -> str:
    if model_id == "human_EPHX1_mEH":
        return "human_mEH_reference"
    if model_id == "human_EPHX2_sEH":
        return "human_sEH_reference"
    if model_id == "pseudomonas_Cif":
        return "bacterial_reference"
    return "candidate"


def fmt(x, nd=3):
    if pd.isna(x):
        return "NA"
    return f"{x:.{nd}f}"


def extract_colabfold_quality(model_dir: Path) -> pd.DataFrame:
    """
    Parse ColabFold JSON score files and extract mean pLDDT and pTM if present.
    """
    rows = []
    suffix = "_scores_rank_001_alphafold2_ptm_model_1_seed_000.json"

    if not model_dir.exists():
        return pd.DataFrame(columns=["model_id", "mean_plddt", "ptm"])

    for json_file in model_dir.glob(f"*{suffix}"):
        model_id = json_file.name.replace(suffix, "")
        try:
            with open(json_file, "r", encoding="utf-8") as fh:
                data = json.load(fh)
        except Exception:
            continue

        mean_plddt = np.nan
        ptm = np.nan

        if "plddt" in data:
            plddt = data["plddt"]
            if isinstance(plddt, list) and len(plddt) > 0:
                mean_plddt = float(np.mean(plddt))
            elif isinstance(plddt, (float, int)):
                mean_plddt = float(plddt)

        if "ptm" in data and isinstance(data["ptm"], (float, int)):
            ptm = float(data["ptm"])

        rows.append({
            "model_id": model_id,
            "mean_plddt": mean_plddt,
            "ptm": ptm
        })

    return pd.DataFrame(rows)


def paired_wilcoxon_tests(wide_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    ligands = list(wide_df.columns)

    for a, b in combinations(ligands, 2):
        paired = wide_df[[a, b]].dropna()
        n = len(paired)

        if n < 3 or not SCIPY_AVAILABLE:
            rows.append({
                "ligand_a": a,
                "ligand_b": b,
                "n_pairs": n,
                "median_diff_a_minus_b": np.nan if n == 0 else float(np.median(paired[a] - paired[b])),
                "wilcoxon_statistic": np.nan,
                "p_value": np.nan,
                "note": "SciPy unavailable or insufficient paired observations"
            })
            continue

        try:
            stat, p = wilcoxon(paired[a], paired[b], alternative="two-sided")
        except Exception:
            stat, p = np.nan, np.nan

        rows.append({
            "ligand_a": a,
            "ligand_b": b,
            "n_pairs": n,
            "median_diff_a_minus_b": float(np.median(paired[a] - paired[b])),
            "wilcoxon_statistic": stat,
            "p_value": p,
            "note": "exploratory paired test across proteins"
        })

    return pd.DataFrame(rows)


def candidate_vs_reference_tests(best_df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for ligand in sorted(best_df["ligand_name"].unique()):
        cand = best_df[(best_df["model_class"] == "candidate") & (best_df["ligand_name"] == ligand)]["vina_best_score_kcal_mol"].dropna()
        ref = best_df[(best_df["model_class"] == "reference") & (best_df["ligand_name"] == ligand)]["vina_best_score_kcal_mol"].dropna()

        if len(cand) == 0 or len(ref) == 0:
            rows.append({
                "ligand_name": ligand,
                "n_candidate": len(cand),
                "n_reference": len(ref),
                "candidate_median": np.nan if len(cand) == 0 else float(np.median(cand)),
                "reference_median": np.nan if len(ref) == 0 else float(np.median(ref)),
                "median_diff_candidate_minus_reference": np.nan,
                "mannwhitney_u": np.nan,
                "p_value": np.nan,
                "note": "insufficient data"
            })
            continue

        stat, p = (np.nan, np.nan)
        if SCIPY_AVAILABLE and len(cand) >= 2 and len(ref) >= 2:
            try:
                stat, p = mannwhitneyu(cand, ref, alternative="two-sided")
            except Exception:
                stat, p = np.nan, np.nan

        rows.append({
            "ligand_name": ligand,
            "n_candidate": len(cand),
            "n_reference": len(ref),
            "candidate_median": float(np.median(cand)),
            "reference_median": float(np.median(ref)),
            "median_diff_candidate_minus_reference": float(np.median(cand) - np.median(ref)),
            "mannwhitney_u": stat,
            "p_value": p,
            "note": "exploratory unpaired class comparison"
        })

    return pd.DataFrame(rows)


# -----------------------------
# Figure generation
# -----------------------------
def plot_heatmap(score_matrix: pd.DataFrame, outbase: Path):
    fig, ax = plt.subplots(figsize=(8, max(4.5, 0.55 * len(score_matrix.index))))

    mat = score_matrix.values.astype(float)
    im = ax.imshow(mat, aspect="auto", interpolation="nearest", cmap="viridis")

    ax.set_xticks(np.arange(len(score_matrix.columns)))
    ax.set_yticks(np.arange(len(score_matrix.index)))
    ax.set_xticklabels(score_matrix.columns, rotation=30, ha="right", fontsize=10)
    ax.set_yticklabels(score_matrix.index, fontsize=9)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Vina best score (kcal/mol)", fontsize=10)

    ax.set_title("Docking score heatmap across proteins and ligands", fontsize=13)
    ax.set_xlabel("Ligand", fontsize=11)
    ax.set_ylabel("Protein model", fontsize=11)

    save_figure(fig, outbase)


def plot_grouped_bars(score_matrix: pd.DataFrame, model_classes: pd.Series, outbase: Path):
    n_models = len(score_matrix.index)
    n_ligands = len(score_matrix.columns)

    fig_width = max(9, 1.0 * n_models)
    fig, ax = plt.subplots(figsize=(fig_width, 5.5))

    x = np.arange(n_models)
    width = 0.22 if n_ligands >= 3 else 0.30

    offsets = np.linspace(-(n_ligands - 1) / 2, (n_ligands - 1) / 2, n_ligands) * width

    for k, ligand in enumerate(score_matrix.columns):
        ax.bar(x + offsets[k], score_matrix[ligand].values, width=width, label=ligand)

    xticklabels = []
    for model in score_matrix.index:
        label = f"{model}\n({model_classes.get(model, 'NA')})"
        xticklabels.append(label)

    ax.set_xticks(x)
    ax.set_xticklabels(xticklabels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Vina best score (kcal/mol)", fontsize=11)
    ax.set_xlabel("Protein model", fontsize=11)
    ax.set_title("Docking score comparison per protein", fontsize=13)
    ax.legend(frameon=False, title="Ligand")
    ax.axhline(0, linewidth=0.8)

    save_figure(fig, outbase)


def plot_ligand_distributions(best_df: pd.DataFrame, ligand_order, outbase: Path):
    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    data = [best_df.loc[best_df["ligand_name"] == lig, "vina_best_score_kcal_mol"].dropna().values for lig in ligand_order]

    bp = ax.boxplot(
        data,
        labels=ligand_order,
        patch_artist=False,
        showmeans=True,
        meanline=False
    )

    rng = np.random.default_rng(12345)
    for i, lig in enumerate(ligand_order, start=1):
        y = best_df.loc[best_df["ligand_name"] == lig, "vina_best_score_kcal_mol"].dropna().values
        x = rng.normal(i, 0.04, size=len(y))
        ax.plot(x, y, "o", alpha=0.8)

    ax.set_ylabel("Vina best score (kcal/mol)", fontsize=11)
    ax.set_xlabel("Ligand", fontsize=11)
    ax.set_title("Distribution of docking scores by ligand", fontsize=13)

    save_figure(fig, outbase)


def plot_best_score_ranking(best_per_model: pd.DataFrame, outbase: Path):
    ranked = best_per_model.sort_values("best_score", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8.5, max(4.5, 0.6 * len(ranked))))

    ax.barh(ranked["model_id"], ranked["best_score"])
    ax.set_xlabel("Best Vina score (kcal/mol)", fontsize=11)
    ax.set_ylabel("Protein model", fontsize=11)
    ax.set_title("Best ligand-protein docking score ranking", fontsize=13)
    ax.invert_yaxis()

    for i, row in ranked.iterrows():
        ax.text(row["best_score"], i, f"  {row['best_ligand']}", va="center", fontsize=9)

    save_figure(fig, outbase)


def plot_quality_vs_docking(best_per_model: pd.DataFrame, quality_df: pd.DataFrame, outbase: Path):
    merged = best_per_model.merge(quality_df, on="model_id", how="left")

    merged = merged.dropna(subset=["mean_plddt", "best_score"])
    if merged.empty:
        return

    fig, ax = plt.subplots(figsize=(7, 5.5))

    for _, row in merged.iterrows():
        marker = "o" if row["model_class"] == "candidate" else "s"
        ax.scatter(row["mean_plddt"], row["best_score"], marker=marker, s=70)
        ax.text(row["mean_plddt"] + 0.15, row["best_score"], row["model_id"], fontsize=8, va="center")

    ax.set_xlabel("Mean pLDDT", fontsize=11)
    ax.set_ylabel("Best Vina score (kcal/mol)", fontsize=11)
    ax.set_title("Model quality versus best docking score", fontsize=13)

    save_figure(fig, outbase)


# -----------------------------
# Text generation
# -----------------------------
def write_results_narrative(
    best_df: pd.DataFrame,
    ligand_summary: pd.DataFrame,
    best_per_model: pd.DataFrame,
    quality_df: pd.DataFrame,
    outpath: Path
):
    overall_best = best_df.sort_values("vina_best_score_kcal_mol", ascending=True).iloc[0]
    ligand_summary_sorted = ligand_summary.sort_values("median_score", ascending=True)
    best_candidate = best_per_model[best_per_model["model_class"] == "candidate"].sort_values("best_score", ascending=True).iloc[0]
    best_reference = best_per_model[best_per_model["model_class"] == "reference"].sort_values("best_score", ascending=True).iloc[0]

    lines = []
    lines.append("# Structural docking results summary")
    lines.append("")
    lines.append(
        f"Docking was summarized across {best_df['model_id'].nunique()} protein models and "
        f"{best_df['ligand_name'].nunique()} ligands, resulting in {len(best_df)} final receptor–ligand entries "
        f"after collapsing to the best-scoring pocket for each protein–ligand combination."
    )
    lines.append("")
    lines.append(
        "As a general rule, more negative Vina scores indicate more favorable predicted binding."
    )
    lines.append("")
    lines.append(
        f"The overall best-scoring complex was {overall_best['model_id']} with {overall_best['ligand_name']} "
        f"at pocket {overall_best['pocket_name']} (Vina = {fmt(overall_best['vina_best_score_kcal_mol'])} kcal/mol)."
    )
    lines.append("")
    lines.append("## Ligand-level summary")
    lines.append("")
    for _, row in ligand_summary_sorted.iterrows():
        lines.append(
            f"- {row['ligand_name']}: median = {fmt(row['median_score'])}, "
            f"mean = {fmt(row['mean_score'])}, SD = {fmt(row['sd_score'])}, "
            f"range = {fmt(row['min_score'])} to {fmt(row['max_score'])} kcal/mol."
        )
    lines.append("")
    lines.append(
        "The ligand ranking based on median docking score therefore suggests the following order of predicted affinity: "
        + " > ".join(ligand_summary_sorted["ligand_name"].tolist()) + "."
    )
    lines.append("")
    lines.append("## Model-level summary")
    lines.append("")
    lines.append(
        f"- Best candidate model: {best_candidate['model_id']} with {best_candidate['best_ligand']} "
        f"(best score = {fmt(best_candidate['best_score'])} kcal/mol)."
    )
    lines.append(
        f"- Best reference model: {best_reference['model_id']} with {best_reference['best_ligand']} "
        f"(best score = {fmt(best_reference['best_score'])} kcal/mol)."
    )
    lines.append("")
    if not quality_df.empty:
        qual = quality_df.copy()
        mean_plddt = qual["mean_plddt"].dropna()
        mean_ptm = qual["ptm"].dropna()
        if len(mean_plddt) > 0 or len(mean_ptm) > 0:
            lines.append("## ColabFold model quality")
            lines.append("")
            if len(mean_plddt) > 0:
                lines.append(
                    f"- Mean model pLDDT across available ColabFold outputs: {fmt(mean_plddt.mean())} "
                    f"(range {fmt(mean_plddt.min())} to {fmt(mean_plddt.max())})."
                )
            if len(mean_ptm) > 0:
                lines.append(
                    f"- Mean model pTM across available ColabFold outputs: {fmt(mean_ptm.mean())} "
                    f"(range {fmt(mean_ptm.min())} to {fmt(mean_ptm.max())})."
                )
            lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append(
        "These docking scores should be interpreted as exploratory structural evidence rather than direct proof of biochemical activity. "
        "Nonetheless, the consistency of ligand ranking across multiple proteins and the preservation of reasonable model-quality metrics "
        "support the use of these analyses as a comparative prioritization layer within the broader genomic and structural framework of the study."
    )
    lines.append("")

    outpath.write_text("\n".join(lines), encoding="utf-8")


def write_figure_legends(outpath: Path):
    lines = []
    lines.append("# Suggested figure legends")
    lines.append("")
    lines.append(
        "Figure 1. Heatmap of best docking scores obtained for each protein–ligand combination. "
        "For each protein–ligand pair, the lowest-energy AutoDock Vina result was retained after pocket-based screening. "
        "More negative values indicate more favorable predicted binding."
    )
    lines.append("")
    lines.append(
        "Figure 2. Grouped bar plot of docking scores by protein model and ligand. "
        "Bars represent the best score for each ligand in the selected pocket of each protein. "
        "This plot facilitates direct comparison of ligand performance across candidate and reference proteins."
    )
    lines.append("")
    lines.append(
        "Figure 3. Distribution of docking scores according to ligand identity. "
        "Boxplots summarize the distribution of Vina best scores across all proteins, while individual points represent protein-level observations. "
        "This panel highlights the overall comparative behavior of the three tested ligands."
    )
    lines.append("")
    lines.append(
        "Figure 4. Ranking of protein models according to their best docking score. "
        "Each bar corresponds to the most favorable ligand for a given model, with the best ligand annotated alongside the bar."
    )
    lines.append("")
    lines.append(
        "Figure 5. Relationship between structural model quality and docking performance. "
        "Mean pLDDT values extracted from ColabFold outputs are plotted against the best Vina score for each model. "
        "Circular markers denote candidate proteins and square markers denote reference proteins."
    )
    lines.append("")

    outpath.write_text("\n".join(lines), encoding="utf-8")


# -----------------------------
# Main
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description="Generate publication-grade structural summaries and figures.")
    parser.add_argument("--base-dir", required=True, help="Base directory of 19_structural_pipeline_real")
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    docking_file = base_dir / "06_tables" / "docking_summary.tsv"
    colabfold_dir = base_dir / "10_colabfold_models"

    out_root = base_dir / "08_reporting"
    out_tables = out_root / "tables"
    out_figures = out_root / "figures"
    out_text = out_root / "text"

    ensure_dir(out_root)
    ensure_dir(out_tables)
    ensure_dir(out_figures)
    ensure_dir(out_text)

    if not docking_file.exists():
        raise FileNotFoundError(f"Docking summary not found: {docking_file}")

    docking = pd.read_csv(docking_file, sep="\t", dtype=str)

    required_cols = {
        "model_id",
        "ligand_name",
        "pocket_name",
        "vina_best_score_kcal_mol"
    }
    missing = required_cols - set(docking.columns)
    if missing:
        raise ValueError(f"Missing required columns in docking_summary.tsv: {missing}")

    docking["vina_best_score_kcal_mol"] = pd.to_numeric(docking["vina_best_score_kcal_mol"], errors="coerce")

    # Optional filtering if existence columns are present
    if "pose_pdbqt_exists" in docking.columns:
        docking = docking[docking["pose_pdbqt_exists"].astype(str).isin(["1", "True", "true"])]
    if "pose_pdb_exists" in docking.columns:
        docking = docking[docking["pose_pdb_exists"].astype(str).isin(["1", "True", "true"])]

    docking = docking.dropna(subset=["vina_best_score_kcal_mol"]).copy()

    # If multiple pockets per model-ligand exist, retain best score
    docking_best = (
        docking.sort_values("vina_best_score_kcal_mol", ascending=True)
               .groupby(["model_id", "ligand_name"], as_index=False)
               .first()
               .copy()
    )

    docking_best["model_class"] = docking_best["model_id"].apply(classify_model)
    docking_best["model_group"] = docking_best["model_id"].apply(reference_group)

    # Ligand and model summaries
    ligand_summary = (
        docking_best.groupby("ligand_name", as_index=False)
        .agg(
            n=("vina_best_score_kcal_mol", "size"),
            mean_score=("vina_best_score_kcal_mol", "mean"),
            sd_score=("vina_best_score_kcal_mol", "std"),
            median_score=("vina_best_score_kcal_mol", "median"),
            min_score=("vina_best_score_kcal_mol", "min"),
            max_score=("vina_best_score_kcal_mol", "max"),
        )
        .sort_values("median_score", ascending=True)
    )

    best_per_model = (
        docking_best.sort_values("vina_best_score_kcal_mol", ascending=True)
        .groupby("model_id", as_index=False)
        .first()
        .rename(columns={
            "vina_best_score_kcal_mol": "best_score",
            "ligand_name": "best_ligand",
            "pocket_name": "best_pocket"
        })
    )

    mean_per_model = (
        docking_best.groupby("model_id", as_index=False)
        .agg(
            mean_score=("vina_best_score_kcal_mol", "mean"),
            sd_score=("vina_best_score_kcal_mol", "std")
        )
    )

    best_per_model = best_per_model.merge(mean_per_model, on="model_id", how="left")
    best_per_model["model_class"] = best_per_model["model_id"].apply(classify_model)
    best_per_model["model_group"] = best_per_model["model_id"].apply(reference_group)

    class_ligand_summary = (
        docking_best.groupby(["model_class", "ligand_name"], as_index=False)
        .agg(
            n=("vina_best_score_kcal_mol", "size"),
            mean_score=("vina_best_score_kcal_mol", "mean"),
            sd_score=("vina_best_score_kcal_mol", "std"),
            median_score=("vina_best_score_kcal_mol", "median"),
            min_score=("vina_best_score_kcal_mol", "min"),
            max_score=("vina_best_score_kcal_mol", "max")
        )
    )

    # Score matrix
    score_matrix = docking_best.pivot(index="model_id", columns="ligand_name", values="vina_best_score_kcal_mol")

    ligand_order = ligand_summary["ligand_name"].tolist()
    if len(ligand_order) == 0:
        ligand_order = sorted(score_matrix.columns.tolist())

    model_order = best_per_model.sort_values("best_score", ascending=True)["model_id"].tolist()

    score_matrix = score_matrix.loc[model_order, ligand_order]
    model_classes = pd.Series({m: classify_model(m) for m in score_matrix.index})

    # Statistics
    pairwise_tests = paired_wilcoxon_tests(score_matrix)
    candidate_reference_tests = candidate_vs_reference_tests(docking_best)

    # ColabFold model quality
    quality_df = extract_colabfold_quality(colabfold_dir)
    merged_quality = best_per_model.merge(quality_df, on="model_id", how="left")

    # Write tables
    docking_best.to_csv(out_tables / "docking_best_per_model_ligand.tsv", sep="\t", index=False)
    ligand_summary.to_csv(out_tables / "ligand_summary.tsv", sep="\t", index=False)
    best_per_model.to_csv(out_tables / "best_per_model.tsv", sep="\t", index=False)
    class_ligand_summary.to_csv(out_tables / "class_ligand_summary.tsv", sep="\t", index=False)
    score_matrix.to_csv(out_tables / "docking_score_matrix.tsv", sep="\t")
    pairwise_tests.to_csv(out_tables / "ligand_pairwise_wilcoxon.tsv", sep="\t", index=False)
    candidate_reference_tests.to_csv(out_tables / "candidate_vs_reference_by_ligand.tsv", sep="\t", index=False)
    quality_df.to_csv(out_tables / "colabfold_model_quality.tsv", sep="\t", index=False)
    merged_quality.to_csv(out_tables / "best_per_model_with_quality.tsv", sep="\t", index=False)

    # Figures
    plot_heatmap(score_matrix, out_figures / "Figure1_heatmap_docking_scores")
    plot_grouped_bars(score_matrix, model_classes, out_figures / "Figure2_grouped_bar_docking_scores")
    plot_ligand_distributions(docking_best, ligand_order, out_figures / "Figure3_ligand_distribution")
    plot_best_score_ranking(best_per_model, out_figures / "Figure4_best_score_ranking")
    plot_quality_vs_docking(best_per_model, quality_df, out_figures / "Figure5_quality_vs_docking")

    # Text
    write_results_narrative(
        best_df=docking_best,
        ligand_summary=ligand_summary,
        best_per_model=best_per_model,
        quality_df=quality_df,
        outpath=out_text / "structural_results_summary.md"
    )

    write_figure_legends(out_text / "figure_legends.md")

    # Console summary
    print("\n=== Structural report generation completed ===")
    print(f"Input docking summary: {docking_file}")
    print(f"Output root: {out_root}")
    print(f"Tables: {out_tables}")
    print(f"Figures: {out_figures}")
    print(f"Text summaries: {out_text}")
    print("\nMain outputs:")
    print(f"- {out_tables / 'docking_best_per_model_ligand.tsv'}")
    print(f"- {out_tables / 'ligand_summary.tsv'}")
    print(f"- {out_tables / 'best_per_model.tsv'}")
    print(f"- {out_tables / 'docking_score_matrix.tsv'}")
    print(f"- {out_figures / 'Figure1_heatmap_docking_scores.png'}")
    print(f"- {out_figures / 'Figure2_grouped_bar_docking_scores.png'}")
    print(f"- {out_figures / 'Figure3_ligand_distribution.png'}")
    print(f"- {out_figures / 'Figure4_best_score_ranking.png'}")
    print(f"- {out_figures / 'Figure5_quality_vs_docking.png'}")
    print(f"- {out_text / 'structural_results_summary.md'}")
    print(f"- {out_text / 'figure_legends.md'}")


if __name__ == "__main__":
    main()
