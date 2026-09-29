#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="<PROJECT_ROOT>"
OUT_DIR="${BASE_DIR}/19_structural_pipeline"

INPUT_DIR="${OUT_DIR}/01_inputs"
MODEL_DIR="${OUT_DIR}/02_models"
MODEL_PDB_DIR="${MODEL_DIR}/pdb"
MANUAL_MODEL_DIR="${MODEL_DIR}/manual_models"
POCKET_DIR="${OUT_DIR}/03_pockets"
DOCK_DIR="${OUT_DIR}/04_docking"
VIS_DIR="${OUT_DIR}/05_visualization"
TAB_DIR="${OUT_DIR}/06_tables"
LOG_DIR="${OUT_DIR}/07_logs"
STATUS_DIR="${OUT_DIR}/08_status"

mkdir -p "${INPUT_DIR}" "${MODEL_DIR}" "${MODEL_PDB_DIR}" "${MANUAL_MODEL_DIR}" \
         "${POCKET_DIR}" "${DOCK_DIR}" "${VIS_DIR}" "${TAB_DIR}" "${LOG_DIR}" "${STATUS_DIR}"

MAX_CANDIDATES="${MAX_CANDIDATES:-4}"
BOX_SIZE="${BOX_SIZE:-22}"
CPU_THREADS="${CPU_THREADS:-2}"

STATUS_FILE="${STATUS_DIR}/structural_pipeline_status.txt"
SUMMARY_FILE="${OUT_DIR}/structural_pipeline_summary.txt"

# ------------------------------------------------------------
# Tool detection
# ------------------------------------------------------------
COLABFOLD="$(command -v colabfold_batch || true)"
FPOCKET="$(command -v fpocket || true)"
OBABEL="$(command -v obabel || true)"
VINA="$(command -v vina || true)"
PYMOL="$(command -v pymol || true)"
CHIMERAX="$(command -v chimerax || true)"

{
  echo "Structural pipeline status"
  echo "========================="
  echo "Base directory: ${BASE_DIR}"
  echo "Output directory: ${OUT_DIR}"
  echo "MAX_CANDIDATES: ${MAX_CANDIDATES}"
  echo "BOX_SIZE: ${BOX_SIZE}"
  echo "CPU_THREADS: ${CPU_THREADS}"
  echo
  echo "Detected tools:"
  echo "colabfold_batch: ${COLABFOLD:-NOT_FOUND}"
  echo "fpocket: ${FPOCKET:-NOT_FOUND}"
  echo "obabel: ${OBABEL:-NOT_FOUND}"
  echo "vina: ${VINA:-NOT_FOUND}"
  echo "pymol: ${PYMOL:-NOT_FOUND}"
  echo "chimerax: ${CHIMERAX:-NOT_FOUND}"
  echo
} > "${STATUS_FILE}"

# ------------------------------------------------------------
# 1. Prepare candidate and reference inputs
# ------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import os
import re
import pandas as pd
from Bio import SeqIO

BASE = Path("<PROJECT_ROOT>")
OUT = BASE / "19_structural_pipeline"
INPUT = OUT / "01_inputs"
INPUT.mkdir(parents=True, exist_ok=True)

MAX_CANDIDATES = int(os.environ.get("MAX_CANDIDATES", "4"))

candidate_table_candidates = [
    BASE / "18_integrated_artifacts" / "tables" / "integrated_candidate_master_table.tsv",
    BASE / "16_functional_context" / "tables" / "candidate_functional_context_master.tsv",
    BASE / "14_outputs" / "tables" / "final_priority_candidates.tsv",
]
candidate_table = None
for p in candidate_table_candidates:
    if p.exists():
        candidate_table = p
        break
if candidate_table is None:
    raise FileNotFoundError("No candidate master/final table found.")

priority_fasta = BASE / "08_priority_candidates" / "priority_candidates_final.faa"
ref_fasta = BASE / "07_blast_reference" / "epoxide_hydrolase_references.faa"

df = pd.read_csv(candidate_table, sep="\t")

# Ensure columns
if "qseqid" not in df.columns:
    raise ValueError("Candidate table lacks qseqid column.")

if "bitscore" not in df.columns:
    df["bitscore"] = 0

if "integrated_priority" not in df.columns:
    df["integrated_priority"] = "exploratory"

if "functional_context_score" not in df.columns:
    df["functional_context_score"] = 0

priority_rank = {
    "high_priority_cif_like": 1,
    "high_priority_ephx1_like": 2,
    "moderate_priority_cif_related": 3,
    "moderate_priority_ephx2_related": 4,
    "exploratory": 5
}
df["priority_rank"] = df["integrated_priority"].map(priority_rank).fillna(99)
df["functional_context_score"] = pd.to_numeric(df["functional_context_score"], errors="coerce").fillna(0)
df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce").fillna(0)

selected = (
    df.sort_values(["priority_rank", "functional_context_score", "bitscore"], ascending=[True, False, False])
      .drop_duplicates("qseqid")
      .head(MAX_CANDIDATES)
      .copy()
)

def safe_name(x: str):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))

seqs = {rec.id: rec for rec in SeqIO.parse(priority_fasta, "fasta")}
refs = {rec.id: rec for rec in SeqIO.parse(ref_fasta, "fasta")}

# Export selected candidates
manifest_rows = []
combined_records = []

for _, row in selected.iterrows():
    qid = row["qseqid"]
    if qid not in seqs:
        continue
    rec = seqs[qid]
    safe = safe_name(qid)
    out_fa = INPUT / f"{safe}.faa"
    SeqIO.write([rec], out_fa, "fasta")
    combined_records.append(rec)

    manifest_rows.append({
        "record_type": "candidate",
        "qseqid": qid,
        "safe_id": safe,
        "integrated_priority": row.get("integrated_priority", ""),
        "bitscore": row.get("bitscore", ""),
        "functional_context_score": row.get("functional_context_score", "")
    })

# Export core references
for ref_id in ["human_EPHX2_sEH", "human_EPHX1_mEH", "pseudomonas_Cif"]:
    if ref_id in refs:
        rec = refs[ref_id]
        safe = safe_name(ref_id)
        out_fa = INPUT / f"{safe}.faa"
        SeqIO.write([rec], out_fa, "fasta")
        combined_records.append(rec)

        manifest_rows.append({
            "record_type": "reference",
            "qseqid": ref_id,
            "safe_id": safe,
            "integrated_priority": "reference",
            "bitscore": "",
            "functional_context_score": ""
        })

combined_fa = INPUT / "structural_priority_plus_references.faa"
SeqIO.write(combined_records, combined_fa, "fasta")

manifest = pd.DataFrame(manifest_rows)
manifest.to_csv(INPUT / "structural_manifest.tsv", sep="\t", index=False)

# Ligand template
ligand_template = INPUT / "ligands.tsv"
if not ligand_template.exists():
    ligand_template.write_text(
        "ligand_name\tsmiles\tenabled\tcomment\n"
        "ethylene_oxide\tC1CO1\t1\tminimal epoxide exploratory ligand\n"
        "propylene_oxide\tCC1CO1\t1\tminimal epoxide exploratory ligand\n"
        "glycidol\tOCC1CO1\t1\tminimal epoxide exploratory ligand\n"
        "manual_add_TPPU\t\t0\tfill the SMILES manually if you want to test TPPU\n"
        "manual_add_AUDA\t\t0\tfill the SMILES manually if you want to test AUDA\n",
        encoding="utf-8"
    )

print(f"Prepared structural manifest: {INPUT / 'structural_manifest.tsv'}")
print(f"Prepared combined FASTA: {combined_fa}")
print(f"Prepared ligand template: {ligand_template}")
PY

# ------------------------------------------------------------
# 2. Optional structure modeling with ColabFold
# ------------------------------------------------------------
if [[ -n "${COLABFOLD}" ]]; then
  echo "[2/8] ColabFold detected. Running structure prediction..." | tee -a "${STATUS_FILE}"
  CF_OUT="${MODEL_DIR}/colabfold_raw"
  mkdir -p "${CF_OUT}"

  "${COLABFOLD}" \
    "${INPUT_DIR}/structural_priority_plus_references.faa" \
    "${CF_OUT}" \
    > "${LOG_DIR}/colabfold_stdout.log" 2> "${LOG_DIR}/colabfold_stderr.log" || true

  echo "ColabFold run attempted." >> "${STATUS_FILE}"
else
  echo "[2/8] ColabFold not detected. Skipping automatic structure prediction." | tee -a "${STATUS_FILE}"
fi

# ------------------------------------------------------------
# 3. Collect PDB models from ColabFold or manual directory
# ------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import re
import shutil

BASE = Path("<PROJECT_ROOT>")
OUT = BASE / "19_structural_pipeline"
MANIFEST = OUT / "01_inputs" / "structural_manifest.tsv"
MODEL_PDB_DIR = OUT / "02_models" / "pdb"
MANUAL_MODEL_DIR = OUT / "02_models" / "manual_models"
CF_OUT = OUT / "02_models" / "colabfold_raw"

import pandas as pd
manifest = pd.read_csv(MANIFEST, sep="\t")

MODEL_PDB_DIR.mkdir(parents=True, exist_ok=True)

def safe_name(x: str):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))

# First copy manual models if present
for pdb in MANUAL_MODEL_DIR.glob("*.pdb"):
    target = MODEL_PDB_DIR / pdb.name
    if not target.exists():
        shutil.copy2(pdb, target)

# Then search ColabFold outputs and map by safe_id
for _, row in manifest.iterrows():
    safe_id = row["safe_id"]
    target = MODEL_PDB_DIR / f"{safe_id}.pdb"
    if target.exists():
        continue
    candidates = list(CF_OUT.rglob(f"{safe_id}*rank_001*.pdb")) + list(CF_OUT.rglob(f"{safe_id}*.pdb"))
    if candidates:
        shutil.copy2(candidates[0], target)

print(f"Collected PDB models into: {MODEL_PDB_DIR}")
PY

# ------------------------------------------------------------
# 4. Pocket detection with fpocket
# ------------------------------------------------------------
if [[ -n "${FPOCKET}" ]]; then
  echo "[4/8] Running fpocket on available PDB models..." | tee -a "${STATUS_FILE}"
  for pdb in "${MODEL_PDB_DIR}"/*.pdb; do
    [[ -e "$pdb" ]] || continue
    bn="$(basename "$pdb" .pdb)"
    workdir="${POCKET_DIR}/${bn}"
    mkdir -p "${workdir}"
    cp "$pdb" "${workdir}/${bn}.pdb"
    (
      cd "${workdir}"
      "${FPOCKET}" -f "${bn}.pdb" > "${LOG_DIR}/${bn}_fpocket_stdout.log" 2> "${LOG_DIR}/${bn}_fpocket_stderr.log" || true
    )
  done
else
  echo "[4/8] fpocket not detected. Pocket analysis skipped." | tee -a "${STATUS_FILE}"
fi

# ------------------------------------------------------------
# 5. Parse pocket centroids and create docking job table
# ------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import pandas as pd

BASE = Path("<PROJECT_ROOT>")
OUT = BASE / "19_structural_pipeline"
POCKET_DIR = OUT / "03_pockets"
TAB_DIR = OUT / "06_tables"
INPUT_DIR = OUT / "01_inputs"

TAB_DIR.mkdir(parents=True, exist_ok=True)

rows = []
for model_dir in POCKET_DIR.iterdir():
    if not model_dir.is_dir():
        continue
    out_dirs = [d for d in model_dir.iterdir() if d.is_dir() and d.name.endswith("_out")]
    for od in out_dirs:
        pockets_dir = od / "pockets"
        if not pockets_dir.exists():
            continue
        for pocket_file in sorted(pockets_dir.glob("pocket*_atm.pdb")):
            pocket_name = pocket_file.stem
            xs, ys, zs = [], [], []
            with open(pocket_file, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            xs.append(x); ys.append(y); zs.append(z)
                        except Exception:
                            pass
            if xs:
                rows.append({
                    "model_id": model_dir.name,
                    "pocket_name": pocket_name,
                    "center_x": sum(xs)/len(xs),
                    "center_y": sum(ys)/len(ys),
                    "center_z": sum(zs)/len(zs),
                    "n_atoms": len(xs),
                    "pocket_atoms_file": str(pocket_file)
                })

pocket_df = pd.DataFrame(rows)
pocket_df.to_csv(TAB_DIR / "pocket_centroids.tsv", sep="\t", index=False)

# Create docking jobs against the top pocket of each model
ligands = pd.read_csv(INPUT_DIR / "ligands.tsv", sep="\t")
ligands = ligands[(ligands["enabled"] == 1) & ligands["smiles"].fillna("").ne("")].copy()

jobs = []
if not pocket_df.empty and not ligands.empty:
    top_pockets = pocket_df.sort_values(["model_id", "n_atoms"], ascending=[True, False]).drop_duplicates("model_id")
    for _, prow in top_pockets.iterrows():
        for _, lrow in ligands.iterrows():
            jobs.append({
                "model_id": prow["model_id"],
                "pocket_name": prow["pocket_name"],
                "center_x": prow["center_x"],
                "center_y": prow["center_y"],
                "center_z": prow["center_z"],
                "ligand_name": lrow["ligand_name"],
                "smiles": lrow["smiles"]
            })

jobs_df = pd.DataFrame(jobs)
jobs_df.to_csv(TAB_DIR / "docking_jobs.tsv", sep="\t", index=False)

print(f"Pocket centroids: {TAB_DIR / 'pocket_centroids.tsv'}")
print(f"Docking jobs: {TAB_DIR / 'docking_jobs.tsv'}")
PY

# ------------------------------------------------------------
# 6. Prepare ligands and receptors with Open Babel
# ------------------------------------------------------------
if [[ -n "${OBABEL}" ]]; then
  echo "[6/8] Preparing ligands and receptors with Open Babel..." | tee -a "${STATUS_FILE}"
  mkdir -p "${DOCK_DIR}/ligands" "${DOCK_DIR}/receptors"

  # Ligands from SMILES
  python - <<'PY'
from pathlib import Path
import pandas as pd

BASE = Path("<PROJECT_ROOT>")
OUT = BASE / "19_structural_pipeline"
INPUT = OUT / "01_inputs"
TAB = OUT / "06_tables"

ligands = pd.read_csv(INPUT / "ligands.tsv", sep="\t")
ligands = ligands[(ligands["enabled"] == 1) & ligands["smiles"].fillna("").ne("")].copy()
ligands.to_csv(TAB / "enabled_ligands.tsv", sep="\t", index=False)
print(f"Enabled ligands: {len(ligands)}")
PY

  if [[ -f "${TAB_DIR}/enabled_ligands.tsv" ]]; then
    tail -n +2 "${TAB_DIR}/enabled_ligands.tsv" | while IFS=$'\t' read -r ligand_name smiles enabled comment; do
      [[ -n "${ligand_name}" ]] || continue
      "${OBABEL}" -:"${smiles}" -osdf -O "${DOCK_DIR}/ligands/${ligand_name}.sdf" --gen3d \
        > "${LOG_DIR}/${ligand_name}_obabel_stdout.log" 2> "${LOG_DIR}/${ligand_name}_obabel_stderr.log" || true

      "${OBABEL}" "${DOCK_DIR}/ligands/${ligand_name}.sdf" -opdbqt -O "${DOCK_DIR}/ligands/${ligand_name}.pdbqt" \
        >> "${LOG_DIR}/${ligand_name}_obabel_stdout.log" 2>> "${LOG_DIR}/${ligand_name}_obabel_stderr.log" || true

      "${OBABEL}" "${DOCK_DIR}/ligands/${ligand_name}.sdf" -opdb -O "${DOCK_DIR}/ligands/${ligand_name}.pdb" \
        >> "${LOG_DIR}/${ligand_name}_obabel_stdout.log" 2>> "${LOG_DIR}/${ligand_name}_obabel_stderr.log" || true
    done
  fi

  # Receptors from PDB
  for pdb in "${MODEL_PDB_DIR}"/*.pdb; do
    [[ -e "$pdb" ]] || continue
    bn="$(basename "$pdb" .pdb)"
    "${OBABEL}" "$pdb" -opdbqt -O "${DOCK_DIR}/receptors/${bn}.pdbqt" \
      > "${LOG_DIR}/${bn}_receptor_prep_stdout.log" 2> "${LOG_DIR}/${bn}_receptor_prep_stderr.log" || true
  done
else
  echo "[6/8] Open Babel not detected. Ligand/receptor preparation skipped." | tee -a "${STATUS_FILE}"
fi

# ------------------------------------------------------------
# 7. Docking with Vina
# ------------------------------------------------------------
if [[ -n "${VINA}" && -n "${OBABEL}" && -f "${TAB_DIR}/docking_jobs.tsv" ]]; then
  echo "[7/8] Running exploratory docking with Vina..." | tee -a "${STATUS_FILE}"
  tail -n +2 "${TAB_DIR}/docking_jobs.tsv" | while IFS=$'\t' read -r model_id pocket_name center_x center_y center_z ligand_name smiles; do
    [[ -n "${model_id}" ]] || continue
    receptor="${DOCK_DIR}/receptors/${model_id}.pdbqt"
    ligand="${DOCK_DIR}/ligands/${ligand_name}.pdbqt"
    out_pose="${DOCK_DIR}/${model_id}__${ligand_name}__${pocket_name}.pdbqt"
    out_log="${DOCK_DIR}/${model_id}__${ligand_name}__${pocket_name}.log"

    if [[ -f "${receptor}" && -f "${ligand}" ]]; then
      "${VINA}" \
        --receptor "${receptor}" \
        --ligand "${ligand}" \
        --center_x "${center_x}" \
        --center_y "${center_y}" \
        --center_z "${center_z}" \
        --size_x "${BOX_SIZE}" \
        --size_y "${BOX_SIZE}" \
        --size_z "${BOX_SIZE}" \
        --cpu "${CPU_THREADS}" \
        --out "${out_pose}" \
        --log "${out_log}" \
        > "${LOG_DIR}/${model_id}__${ligand_name}__${pocket_name}_vina_stdout.log" \
        2> "${LOG_DIR}/${model_id}__${ligand_name}__${pocket_name}_vina_stderr.log" || true

      # convert pose to PDB for visualization
      if [[ -f "${out_pose}" ]]; then
        "${OBABEL}" "${out_pose}" -opdb -O "${DOCK_DIR}/${model_id}__${ligand_name}__${pocket_name}.pdb" \
          >> "${LOG_DIR}/${model_id}__${ligand_name}__${pocket_name}_vina_stdout.log" \
          2>> "${LOG_DIR}/${model_id}__${ligand_name}__${pocket_name}_vina_stderr.log" || true
      fi
    fi
  done
else
  echo "[7/8] Vina and/or Open Babel not detected, or docking_jobs.tsv missing. Docking skipped." | tee -a "${STATUS_FILE}"
fi

# ------------------------------------------------------------
# 8. Collect results and generate PyMOL / ChimeraX scripts
# ------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import re
import pandas as pd

BASE = Path("<PROJECT_ROOT>")
OUT = BASE / "19_structural_pipeline"
INPUT = OUT / "01_inputs"
MODEL_PDB = OUT / "02_models" / "pdb"
DOCK = OUT / "04_docking"
VIS = OUT / "05_visualization"
TAB = OUT / "06_tables"

VIS.mkdir(parents=True, exist_ok=True)

manifest = pd.read_csv(INPUT / "structural_manifest.tsv", sep="\t")

# Model status
model_rows = []
for _, row in manifest.iterrows():
    safe_id = row["safe_id"]
    pdb = MODEL_PDB / f"{safe_id}.pdb"
    model_rows.append({
        "qseqid": row["qseqid"],
        "safe_id": safe_id,
        "record_type": row["record_type"],
        "model_pdb_exists": int(pdb.exists()),
        "model_pdb_path": str(pdb) if pdb.exists() else ""
    })
model_df = pd.DataFrame(model_rows)
model_df.to_csv(TAB / "structural_model_status.tsv", sep="\t", index=False)

# Docking summary
dock_rows = []
for logf in DOCK.glob("*.log"):
    name = logf.stem
    parts = name.split("__")
    if len(parts) < 3:
        continue
    model_id, ligand_name, pocket_name = parts[0], parts[1], parts[2]
    best_score = None
    with open(logf, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if re.match(r"^\s*1\s+", line):
                toks = line.split()
                if len(toks) >= 2:
                    try:
                        best_score = float(toks[1])
                    except Exception:
                        pass
                break

    pose_pdb = DOCK / f"{model_id}__{ligand_name}__{pocket_name}.pdb"
    pose_pdbqt = DOCK / f"{model_id}__{ligand_name}__{pocket_name}.pdbqt"
    receptor_pdb = MODEL_PDB / f"{model_id}.pdb"

    dock_rows.append({
        "model_id": model_id,
        "ligand_name": ligand_name,
        "pocket_name": pocket_name,
        "vina_best_score_kcal_mol": best_score,
        "pose_pdb_exists": int(pose_pdb.exists()),
        "pose_pdbqt_exists": int(pose_pdbqt.exists()),
        "pose_pdb": str(pose_pdb) if pose_pdb.exists() else "",
        "receptor_pdb": str(receptor_pdb) if receptor_pdb.exists() else ""
    })

dock_df = pd.DataFrame(dock_rows)
dock_df.to_csv(TAB / "docking_summary.tsv", sep="\t", index=False)

# Generate visualization scripts
for _, row in dock_df.iterrows():
    if row["pose_pdb_exists"] != 1:
        continue
    model_id = row["model_id"]
    ligand_name = row["ligand_name"]
    pocket_name = row["pocket_name"]
    receptor_pdb = Path(row["receptor_pdb"])
    pose_pdb = Path(row["pose_pdb"])

    safe_vis = f"{model_id}__{ligand_name}__{pocket_name}"

    pml = VIS / f"{safe_vis}.pml"
    pml.write_text(
        f"reinitialize\n"
        f"load {receptor_pdb}, receptor\n"
        f"load {pose_pdb}, ligand\n"
        f"hide everything\n"
        f"show cartoon, receptor\n"
        f"color gray80, receptor\n"
        f"show sticks, ligand\n"
        f"color yellow, ligand\n"
        f"select pocketres, byres (receptor within 4 of ligand)\n"
        f"show sticks, pocketres\n"
        f"color cyan, pocketres\n"
        f"zoom pocketres, 6\n"
        f"bg_color white\n"
        f"ray 1800,1400\n"
        f"png {VIS / (safe_vis + '.png')}, dpi=300\n",
        encoding="utf-8"
    )

    cxc = VIS / f"{safe_vis}.cxc"
    cxc.write_text(
        f"open {receptor_pdb}\n"
        f"open {pose_pdb}\n"
        f"cartoon\n"
        f"style stick ligand\n"
        f"color gray target #1\n"
        f"color yellow target #2\n"
        f"select zone #2 range 4\n"
        f"style stick sel\n"
        f"color cyan sel\n"
        f"view\n"
        f"save {VIS / (safe_vis + '.cxs')}\n",
        encoding="utf-8"
    )

print(f"Structural model status: {TAB / 'structural_model_status.tsv'}")
print(f"Docking summary: {TAB / 'docking_summary.tsv'}")
print(f"Visualization scripts: {VIS}")
PY

# Optional: render PyMOL PNGs if pymol is available
if [[ -n "${PYMOL}" ]]; then
  echo "[8/8] Rendering PyMOL PNGs if docking poses exist..." | tee -a "${STATUS_FILE}"
  for pml in "${VIS_DIR}"/*.pml; do
    [[ -e "$pml" ]] || continue
    "${PYMOL}" -cq "$pml" > "${LOG_DIR}/$(basename "$pml" .pml)_pymol_stdout.log" 2> "${LOG_DIR}/$(basename "$pml" .pml)_pymol_stderr.log" || true
  done
else
  echo "[8/8] PyMOL not detected. Visualization scripts were generated but not rendered." | tee -a "${STATUS_FILE}"
fi

# Final summary
{
  echo "Structural pipeline completed."
  echo "Input directory: ${INPUT_DIR}"
  echo "Model PDB directory: ${MODEL_PDB_DIR}"
  echo "Pocket directory: ${POCKET_DIR}"
  echo "Docking directory: ${DOCK_DIR}"
  echo "Visualization directory: ${VIS_DIR}"
  echo "Tables directory: ${TAB_DIR}"
  echo "Status file: ${STATUS_FILE}"
  echo
  echo "Key outputs:"
  echo "  - ${INPUT_DIR}/structural_manifest.tsv"
  echo "  - ${TAB_DIR}/pocket_centroids.tsv"
  echo "  - ${TAB_DIR}/docking_jobs.tsv"
  echo "  - ${TAB_DIR}/structural_model_status.tsv"
  echo "  - ${TAB_DIR}/docking_summary.tsv"
} > "${SUMMARY_FILE}"

cat "${SUMMARY_FILE}"
