#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile
import pandas as pd
from Bio import SeqIO

BASE_DIR = Path("<PROJECT_ROOT>")

PRIORITY_FASTA = BASE_DIR / "08_priority_candidates" / "priority_candidates_final.faa"
FINAL_CAND = BASE_DIR / "14_outputs" / "tables" / "final_priority_candidates.tsv"
FUNC_MASTER = BASE_DIR / "16_functional_context" / "tables" / "candidate_functional_context_master.tsv"

OUT_DIR = BASE_DIR / "17_optional_advanced_layers"
RAW_DIR = OUT_DIR / "raw_outputs"
TAB_DIR = OUT_DIR / "tables"
STATUS_DIR = OUT_DIR / "status"

for d in [OUT_DIR, RAW_DIR, TAB_DIR, STATUS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

STATUS_FILE = STATUS_DIR / "advanced_layers_status.txt"
SUMMARY_FILE = OUT_DIR / "advanced_layers_summary.txt"

# ------------------------------------------------------------
# Helper utilities
# ------------------------------------------------------------
def run_cmd(cmd, cwd=None, stdout_path=None, stderr_path=None):
    try:
        with open(stdout_path, "w") if stdout_path else subprocess.DEVNULL as out_f, \
             open(stderr_path, "w") if stderr_path else subprocess.DEVNULL as err_f:
            proc = subprocess.run(
                cmd,
                cwd=cwd,
                stdout=out_f if stdout_path else subprocess.DEVNULL,
                stderr=err_f if stderr_path else subprocess.DEVNULL,
                check=False,
                text=True
            )
        return proc.returncode
    except Exception:
        return 999

def write_status(lines):
    STATUS_FILE.write_text("\n".join(lines), encoding="utf-8")

def fasta_ids(path):
    return [rec.id for rec in SeqIO.parse(path, "fasta")]

# ------------------------------------------------------------
# 1. Localization layer
# ------------------------------------------------------------
def run_localization_layer():
    out_sub = RAW_DIR / "localization"
    out_sub.mkdir(parents=True, exist_ok=True)

    status = []
    results_rows = []

    signalp = shutil.which("signalp6") or shutil.which("signalp")
    deeptmhmm = shutil.which("deeptmhmm") or shutil.which("DeepTMHMM")
    tmhmm = shutil.which("tmhmm")
    psortb = shutil.which("psortb")
    targetp = shutil.which("targetp")
    wolfpsort = shutil.which("wolfpsort")

    status.append("[Localization layer]")
    status.append(f"signalp: {signalp if signalp else 'NOT_FOUND'}")
    status.append(f"deeptmhmm: {deeptmhmm if deeptmhmm else 'NOT_FOUND'}")
    status.append(f"tmhmm: {tmhmm if tmhmm else 'NOT_FOUND'}")
    status.append(f"psortb: {psortb if psortb else 'NOT_FOUND'}")
    status.append(f"targetp: {targetp if targetp else 'NOT_FOUND'}")
    status.append(f"wolfpsort: {wolfpsort if wolfpsort else 'NOT_FOUND'}")

    ids = fasta_ids(PRIORITY_FASTA)
    for qid in ids:
        results_rows.append({
            "sequence_id": qid,
            "signalp_status": "not_run",
            "tm_status": "not_run",
            "psort_status": "not_run",
            "targetp_status": "not_run",
            "wolfpsort_status": "not_run"
        })

    # Minimal raw execution only if tools exist
    if signalp:
        status.append("signalp detected; raw execution attempted.")
        sig_out = out_sub / "signalp_stdout.txt"
        sig_err = out_sub / "signalp_stderr.txt"
        # Different SignalP versions differ a lot; keep raw output only
        rc = run_cmd([signalp, "-h"], stdout_path=sig_out, stderr_path=sig_err)
        status.append(f"signalp probe return code: {rc}")
    else:
        status.append("signalp not available; no predictions generated.")

    if deeptmhmm:
        status.append("DeepTMHMM detected; raw execution support available.")
    elif tmhmm:
        status.append("TMHMM detected; raw execution support available.")
    else:
        status.append("No TM predictor available.")

    if psortb:
        status.append("psortb detected.")
    else:
        status.append("psortb not available.")

    if targetp:
        status.append("targetp detected.")
    else:
        status.append("targetp not available.")

    if wolfpsort:
        status.append("wolfpsort detected.")
    else:
        status.append("wolfpsort not available.")

    loc_df = pd.DataFrame(results_rows)
    loc_df.to_csv(TAB_DIR / "optional_localization_status_table.tsv", sep="\t", index=False)
    return status

# ------------------------------------------------------------
# 2. eggNOG / KO / pathway layer
# ------------------------------------------------------------
def run_annotation_layer():
    out_sub = RAW_DIR / "annotation"
    out_sub.mkdir(parents=True, exist_ok=True)

    status = []
    emapper = shutil.which("emapper.py") or shutil.which("emapper")
    kofam = shutil.which("exec_annotation") or shutil.which("kofam_scan") or shutil.which("kofamscan")

    status.append("[Functional annotation layer]")
    status.append(f"eggNOG-mapper: {emapper if emapper else 'NOT_FOUND'}")
    status.append(f"KOfam/KO scanner: {kofam if kofam else 'NOT_FOUND'}")

    # eggNOG layer
    if emapper:
        out_prefix = out_sub / "eggnog_priority"
        stdout = out_sub / "eggnog_stdout.txt"
        stderr = out_sub / "eggnog_stderr.txt"
        cmd = [
            emapper,
            "-i", str(PRIORITY_FASTA),
            "--itype", "proteins",
            "--output", out_prefix.name,
            "--output_dir", str(out_sub),
            "--cpu", "2"
        ]
        rc = run_cmd(cmd, stdout_path=stdout, stderr_path=stderr)
        status.append(f"eggNOG return code: {rc}")

        annot_candidates = list(out_sub.glob("eggnog_priority*.annotations"))
        if annot_candidates:
            ann = annot_candidates[0]
            try:
                df = pd.read_csv(ann, sep="\t", comment="#", header=None)
                df.to_csv(TAB_DIR / "eggnog_annotations_raw.tsv", sep="\t", index=False)
                status.append(f"eggNOG annotation file parsed: {ann.name}")
            except Exception as e:
                status.append(f"eggNOG annotation parse failed: {e}")
        else:
            status.append("eggNOG annotations file not found.")
    else:
        status.append("eggNOG-mapper not available; no eggNOG/GO/KEGG-like functional assignment generated.")

    # KOfam / KO layer
    if kofam:
        stdout = out_sub / "kofam_stdout.txt"
        stderr = out_sub / "kofam_stderr.txt"
        # We do not assume installed database paths; just probe command availability.
        rc = run_cmd([kofam, "-h"], stdout_path=stdout, stderr_path=stderr)
        status.append(f"KOfam probe return code: {rc}")
        status.append("KOfam command detected, but full KO assignment requires a properly configured local database.")
    else:
        status.append("KOfam/KO scanner not available; no formal KO-based pathway mapping generated.")

    return status

# ------------------------------------------------------------
# 3. Structural layer
# ------------------------------------------------------------
def run_structure_layer():
    out_sub = RAW_DIR / "structure"
    out_sub.mkdir(parents=True, exist_ok=True)

    status = []
    colabfold = shutil.which("colabfold_batch")
    hhsearch = shutil.which("hhsearch")
    foldseek = shutil.which("foldseek")
    mmseqs = shutil.which("mmseqs")
    chimerax = shutil.which("chimerax")

    status.append("[Structural layer]")
    status.append(f"colabfold_batch: {colabfold if colabfold else 'NOT_FOUND'}")
    status.append(f"hhsearch: {hhsearch if hhsearch else 'NOT_FOUND'}")
    status.append(f"foldseek: {foldseek if foldseek else 'NOT_FOUND'}")
    status.append(f"mmseqs: {mmseqs if mmseqs else 'NOT_FOUND'}")
    status.append(f"chimerax: {chimerax if chimerax else 'NOT_FOUND'}")

    # No heavy prediction automatically run unless ColabFold is installed
    if colabfold:
        status.append("ColabFold detected. Structure prediction could be launched on prioritized candidates.")
        status.append("Automatic prediction was not forced by default to avoid heavy runtime and model downloads.")
    else:
        status.append("No local structure-prediction engine detected; structural modeling not executed.")

    # create a candidate list for future structural work
    final_df = pd.read_csv(FINAL_CAND, sep="\t")
    struct_df = final_df[
        final_df["integrated_priority"].isin(["high_priority_cif_like", "high_priority_ephx1_like"])
    ].copy()
    struct_df.to_csv(TAB_DIR / "structure_priority_candidates.tsv", sep="\t", index=False)
    status.append("A structure-priority candidate table was exported for future modeling.")

    return status

# ------------------------------------------------------------
# 4. Docking layer
# ------------------------------------------------------------
def run_docking_layer():
    out_sub = RAW_DIR / "docking"
    out_sub.mkdir(parents=True, exist_ok=True)

    status = []
    vina = shutil.which("vina")
    smina = shutil.which("smina")
    gnina = shutil.which("gnina")
    obabel = shutil.which("obabel")

    status.append("[Docking layer]")
    status.append(f"vina: {vina if vina else 'NOT_FOUND'}")
    status.append(f"smina: {smina if smina else 'NOT_FOUND'}")
    status.append(f"gnina: {gnina if gnina else 'NOT_FOUND'}")
    status.append(f"obabel: {obabel if obabel else 'NOT_FOUND'}")

    if vina or smina or gnina:
        status.append("A docking engine is available, but docking was not automatically executed because prepared receptor structures and curated ligand files were not yet available.")
    else:
        status.append("No docking engine detected; docking layer not executed.")

    status.append("Docking should remain exploratory and should only be performed after structural models and binding-site candidates are curated.")
    return status

# ------------------------------------------------------------
# 5. Final integration
# ------------------------------------------------------------
def build_master_status_table():
    final_df = pd.read_csv(FINAL_CAND, sep="\t")
    func_df = pd.read_csv(FUNC_MASTER, sep="\t")

    merged = final_df.merge(
        func_df[["qseqid", "functional_context_score", "functional_context_class"]],
        on="qseqid",
        how="left"
    )

    # placeholders for advanced layers
    merged["localization_layer"] = "not_executed_or_not_available"
    merged["formal_pathway_layer"] = "not_executed_or_not_available"
    merged["structural_modeling_layer"] = "not_executed_or_not_available"
    merged["docking_layer"] = "not_executed_or_not_available"

    merged.to_csv(TAB_DIR / "advanced_layers_master_status.tsv", sep="\t", index=False)
    return merged

def main():
    all_status = []

    all_status.extend(run_localization_layer())
    all_status.append("")

    all_status.extend(run_annotation_layer())
    all_status.append("")

    all_status.extend(run_structure_layer())
    all_status.append("")

    all_status.extend(run_docking_layer())
    all_status.append("")

    master = build_master_status_table()

    write_status(all_status)

    summary_lines = [
        "Optional advanced layers executed.",
        f"Status file: {STATUS_FILE}",
        f"Master status table: {TAB_DIR / 'advanced_layers_master_status.tsv'}",
        f"Localization status table: {TAB_DIR / 'optional_localization_status_table.tsv'}",
        f"Structure-priority candidate table: {TAB_DIR / 'structure_priority_candidates.tsv'}",
        f"Raw outputs directory: {RAW_DIR}",
        f"Rows in master status table: {len(master)}"
    ]
    SUMMARY_FILE.write_text("\n".join(summary_lines), encoding="utf-8")
    print("\n".join(summary_lines))

if __name__ == "__main__":
    main()
