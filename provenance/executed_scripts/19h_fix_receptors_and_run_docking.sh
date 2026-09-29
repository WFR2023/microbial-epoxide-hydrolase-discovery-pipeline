#!/usr/bin/env bash
set -euo pipefail

BASE="<PROJECT_ROOT>/19_structural_pipeline_real"
DOCK_DIR="${BASE}/04_docking"
VIS_DIR="${BASE}/05_visualization"
TAB_DIR="${BASE}/06_tables"
LOG_DIR="${BASE}/07_logs"
PDB_DIR="${BASE}/02_models/pdb"

mkdir -p "${DOCK_DIR}/receptors_clean" "${VIS_DIR}" "${TAB_DIR}" "${LOG_DIR}"

BOX_SIZE="${BOX_SIZE:-22}"
CPU_THREADS="${CPU_THREADS:-2}"

echo "Step 1: cleaning receptor PDBQT files for Vina compatibility..."

python - <<'PY'
from pathlib import Path

base = Path("<PROJECT_ROOT>/19_structural_pipeline_real")
src_dir = base / "04_docking" / "receptors"
dst_dir = base / "04_docking" / "receptors_clean"
dst_dir.mkdir(parents=True, exist_ok=True)

remove_prefixes = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")

for src in sorted(src_dir.glob("*.pdbqt")):
    dst = dst_dir / src.name
    kept = []
    with open(src, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith(remove_prefixes):
                continue
            kept.append(line)
    with open(dst, "w", encoding="utf-8") as out:
        out.writelines(kept)

print(f"Clean receptors written to: {dst_dir}")
PY

echo "Step 2: running docking with cleaned receptors..."

tail -n +2 "${TAB_DIR}/docking_jobs.tsv" | while IFS=$'\t' read -r model_id pocket_name center_x center_y center_z ligand_name smiles; do
  [[ -n "${model_id}" ]] || continue

  receptor="${DOCK_DIR}/receptors_clean/${model_id}.pdbqt"
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

echo "Step 3: summarizing docking results..."

python - <<'PY'
from pathlib import Path
import re
import pandas as pd

base = Path("<PROJECT_ROOT>/19_structural_pipeline_real")
dock = base / "04_docking"
tab = base / "06_tables"
vis = base / "05_visualization"
pdb = base / "02_models" / "pdb"

vis.mkdir(parents=True, exist_ok=True)

rows = []
for logf in sorted(dock.glob("*.log")):
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

    pose_pdbqt = dock / f"{model_id}__{ligand_name}__{pocket_name}.pdbqt"
    pose_pdb = dock / f"{model_id}__{ligand_name}__{pocket_name}.pdb"
    receptor_pdb = pdb / f"{model_id}.pdb"

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
df.to_csv(tab / "docking_summary.tsv", sep="\t", index=False)

for _, row in df.iterrows():
    if row["pose_pdb_exists"] != 1:
        continue
    model_id = row["model_id"]
    ligand_name = row["ligand_name"]
    pocket_name = row["pocket_name"]
    receptor_pdb = row["receptor_pdb"]
    pose_pdb = row["pose_pdb"]
    safe = f"{model_id}__{ligand_name}__{pocket_name}"

    pml = vis / f"{safe}.pml"
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
        f"bg_color white\n",
        encoding="utf-8"
    )

print(f"Docking summary written to: {tab / 'docking_summary.tsv'}")
print(f"Visualization scripts written to: {vis}")
PY

echo "Done."
