#!/usr/bin/env bash
set -euo pipefail

BASE="<PROJECT_ROOT>/19_structural_pipeline_real"
PDB_DIR="${BASE}/02_models/pdb"
POCKET_DIR="${BASE}/03_pockets"
DOCK_DIR="${BASE}/04_docking"
VIS_DIR="${BASE}/05_visualization"
TAB_DIR="${BASE}/06_tables"
LOG_DIR="${BASE}/07_logs"

BOX_SIZE="${BOX_SIZE:-22}"
CPU_THREADS="${CPU_THREADS:-2}"

mkdir -p "${DOCK_DIR}/ligands" "${DOCK_DIR}/receptors" "${VIS_DIR}" "${TAB_DIR}" "${LOG_DIR}"

echo "Using Vina without --log flag"

# 1) prepare ligands
tail -n +2 "${TAB_DIR}/enabled_ligands.tsv" | while IFS=$'\t' read -r ligand_name smiles enabled comment; do
  [[ -n "${ligand_name}" ]] || continue
  obabel -:"${smiles}" -osdf -O "${DOCK_DIR}/ligands/${ligand_name}.sdf" --gen3d \
    > "${LOG_DIR}/${ligand_name}_lig_stdout.log" 2> "${LOG_DIR}/${ligand_name}_lig_stderr.log"
  obabel "${DOCK_DIR}/ligands/${ligand_name}.sdf" -opdbqt -O "${DOCK_DIR}/ligands/${ligand_name}.pdbqt" \
    >> "${LOG_DIR}/${ligand_name}_lig_stdout.log" 2>> "${LOG_DIR}/${ligand_name}_lig_stderr.log"
  obabel "${DOCK_DIR}/ligands/${ligand_name}.sdf" -opdb -O "${DOCK_DIR}/ligands/${ligand_name}.pdb" \
    >> "${LOG_DIR}/${ligand_name}_lig_stdout.log" 2>> "${LOG_DIR}/${ligand_name}_lig_stderr.log"
done

# 2) prepare receptors
for pdb in "${PDB_DIR}"/*.pdb; do
  [[ -e "$pdb" ]] || continue
  bn="$(basename "$pdb" .pdb)"
  obabel "$pdb" -opdbqt -O "${DOCK_DIR}/receptors/${bn}.pdbqt" \
    > "${LOG_DIR}/${bn}_rec_stdout.log" 2> "${LOG_DIR}/${bn}_rec_stderr.log"
done

# 3) run vina
tail -n +2 "${TAB_DIR}/docking_jobs.tsv" | while IFS=$'\t' read -r model_id pocket_name center_x center_y center_z ligand_name smiles; do
  [[ -n "${model_id}" ]] || continue

  receptor="${DOCK_DIR}/receptors/${model_id}.pdbqt"
  ligand="${DOCK_DIR}/ligands/${ligand_name}.pdbqt"
  out_pose="${DOCK_DIR}/${model_id}__${ligand_name}__${pocket_name}.pdbqt"
  out_pdb="${DOCK_DIR}/${model_id}__${ligand_name}__${pocket_name}.pdb"
  out_log="${DOCK_DIR}/${model_id}__${ligand_name}__${pocket_name}.log"

  vina \
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
    > "${out_log}" 2> "${LOG_DIR}/${model_id}__${ligand_name}_vina_stderr.log"

  obabel "${out_pose}" -opdb -O "${out_pdb}" \
    > "${LOG_DIR}/${model_id}__${ligand_name}_poseconv_stdout.log" \
    2> "${LOG_DIR}/${model_id}__${ligand_name}_poseconv_stderr.log"
done

# 4) summarize
python - <<'PY'
from pathlib import Path
import re
import pandas as pd

BASE = Path("<PROJECT_ROOT>/19_structural_pipeline_real")
DOCK = BASE / "04_docking"
TAB = BASE / "06_tables"
VIS = BASE / "05_visualization"
PDB = BASE / "02_models" / "pdb"
VIS.mkdir(parents=True, exist_ok=True)

rows = []
for logf in DOCK.glob("*.log"):
    parts = logf.stem.split("__")
    if len(parts) < 3:
        continue
    model_id, ligand_name, pocket_name = parts[0], parts[1], parts[2]
    score = None
    with open(logf, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if re.match(r"^\s*1\s+", s):
                toks = s.split()
                if len(toks) >= 2:
                    try:
                        score = float(toks[1])
                    except Exception:
                        pass
                break
    pose_pdb = DOCK / f"{model_id}__{ligand_name}__{pocket_name}.pdb"
    pose_pdbqt = DOCK / f"{model_id}__{ligand_name}__{pocket_name}.pdbqt"
    receptor_pdb = PDB / f"{model_id}.pdb"
    rows.append({
        "model_id": model_id,
        "ligand_name": ligand_name,
        "pocket_name": pocket_name,
        "vina_best_score_kcal_mol": score,
        "pose_pdbqt_exists": int(pose_pdbqt.exists()),
        "pose_pdb_exists": int(pose_pdb.exists()),
        "pose_pdbqt": str(pose_pdbqt),
        "pose_pdb": str(pose_pdb),
        "receptor_pdb": str(receptor_pdb),
    })

df = pd.DataFrame(rows)
df.to_csv(TAB / "docking_summary.tsv", sep="\t", index=False)

for _, row in df.iterrows():
    if row["pose_pdb_exists"] != 1:
        continue
    model_id = row["model_id"]
    ligand_name = row["ligand_name"]
    pocket_name = row["pocket_name"]
    receptor_pdb = row["receptor_pdb"]
    pose_pdb = row["pose_pdb"]
    safe = f"{model_id}__{ligand_name}__{pocket_name}"
    pml = VIS / f"{safe}.pml"
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
        f"png {VIS / (safe + '.png')}, dpi=300\n",
        encoding="utf-8"
    )

print(f"Docking summary written to: {TAB / 'docking_summary.tsv'}")
print(f"Visualization scripts written to: {VIS}")
PY

echo "Done."
