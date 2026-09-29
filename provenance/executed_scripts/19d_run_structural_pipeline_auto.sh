#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="<PROJECT_ROOT>"
OUT_DIR="${BASE_DIR}/19_structural_pipeline_auto"

INPUT_DIR="${OUT_DIR}/01_inputs"
MODEL_DIR="${OUT_DIR}/02_models"
MANUAL_MODEL_DIR="${MODEL_DIR}/manual_models"
AUTO_MODEL_DIR="${MODEL_DIR}/auto_models"
PDB_DIR="${MODEL_DIR}/pdb"
POCKET_DIR="${OUT_DIR}/03_pockets"
DOCK_DIR="${OUT_DIR}/04_docking"
VIS_DIR="${OUT_DIR}/05_visualization"
TAB_DIR="${OUT_DIR}/06_tables"
LOG_DIR="${OUT_DIR}/07_logs"
STATUS_DIR="${OUT_DIR}/08_status"

mkdir -p "${INPUT_DIR}" "${MANUAL_MODEL_DIR}" "${AUTO_MODEL_DIR}" "${PDB_DIR}" \
         "${POCKET_DIR}" "${DOCK_DIR}" "${VIS_DIR}" "${TAB_DIR}" "${LOG_DIR}" "${STATUS_DIR}"

MAX_CANDIDATES="${MAX_CANDIDATES:-4}"
BOX_SIZE="${BOX_SIZE:-22}"
CPU_THREADS="${CPU_THREADS:-2}"

STATUS_FILE="${STATUS_DIR}/pipeline_status.txt"
SUMMARY_FILE="${OUT_DIR}/pipeline_summary.txt"

COLABFOLD="$(command -v colabfold_batch || true)"
OBABEL="$(command -v obabel || true)"
VINA="$(command -v vina || true)"
FPOCKET="$(command -v fpocket || true)"
PYMOL="$(command -v pymol || true)"
PYTHON_BIN="$(command -v python || true)"

{
  echo "Automatic structural pipeline"
  echo "============================"
  echo "Base dir: ${BASE_DIR}"
  echo "Output dir: ${OUT_DIR}"
  echo "MAX_CANDIDATES=${MAX_CANDIDATES}"
  echo "BOX_SIZE=${BOX_SIZE}"
  echo "CPU_THREADS=${CPU_THREADS}"
  echo
  echo "Detected tools:"
  echo "python: ${PYTHON_BIN:-NOT_FOUND}"
  echo "colabfold_batch: ${COLABFOLD:-NOT_FOUND}"
  echo "fpocket: ${FPOCKET:-NOT_FOUND}"
  echo "obabel: ${OBABEL:-NOT_FOUND}"
  echo "vina: ${VINA:-NOT_FOUND}"
  echo "pymol: ${PYMOL:-NOT_FOUND}"
  echo
} > "${STATUS_FILE}"

# ------------------------------------------------------------
# 1. Build manifest and FASTAs
# ------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import os
import re
import pandas as pd
from Bio import SeqIO

BASE = Path("<PROJECT_ROOT>")
OUT = BASE / "19_structural_pipeline_auto"
INPUT = OUT / "01_inputs"
INPUT.mkdir(parents=True, exist_ok=True)

MAX_CANDIDATES = int(os.environ.get("MAX_CANDIDATES", "4"))

table_candidates = [
    BASE / "18_integrated_artifacts" / "tables" / "integrated_candidate_master_table.tsv",
    BASE / "16_functional_context" / "tables" / "candidate_functional_context_master.tsv",
    BASE / "14_outputs" / "tables" / "final_priority_candidates.tsv",
]
candidate_table = None
for p in table_candidates:
    if p.exists():
        candidate_table = p
        break
if candidate_table is None:
    raise FileNotFoundError("No candidate table found.")

priority_fasta = BASE / "08_priority_candidates" / "priority_candidates_final.faa"
ref_fasta = BASE / "07_blast_reference" / "epoxide_hydrolase_references.faa"

df = pd.read_csv(candidate_table, sep="\t")
if "integrated_priority" not in df.columns:
    df["integrated_priority"] = "exploratory"
if "functional_context_score" not in df.columns:
    df["functional_context_score"] = 0
if "bitscore" not in df.columns:
    df["bitscore"] = 0

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

rows = []
combined = []

for _, row in selected.iterrows():
    qid = row["qseqid"]
    if qid not in seqs:
        continue
    safe = safe_name(qid)
    rec = seqs[qid]
    SeqIO.write([rec], INPUT / f"{safe}.faa", "fasta")
    combined.append(rec)
    rows.append({
        "record_type": "candidate",
        "qseqid": qid,
        "safe_id": safe,
        "sequence_length": len(rec.seq)
    })

for ref_id in ["human_EPHX2_sEH", "human_EPHX1_mEH", "pseudomonas_Cif"]:
    if ref_id in refs:
        safe = safe_name(ref_id)
        rec = refs[ref_id]
        SeqIO.write([rec], INPUT / f"{safe}.faa", "fasta")
        combined.append(rec)
        rows.append({
            "record_type": "reference",
            "qseqid": ref_id,
            "safe_id": safe,
            "sequence_length": len(rec.seq)
        })

SeqIO.write(combined, INPUT / "structural_priority_plus_references.faa", "fasta")
pd.DataFrame(rows).to_csv(INPUT / "structural_manifest.tsv", sep="\t", index=False)

ligands = INPUT / "ligands.tsv"
if not ligands.exists():
    ligands.write_text(
        "ligand_name\tsmiles\tenabled\tcomment\n"
        "ethylene_oxide\tC1CO1\t1\tminimal epoxide ligand\n"
        "propylene_oxide\tCC1CO1\t1\tminimal epoxide ligand\n"
        "glycidol\tOCC1CO1\t1\tepoxide-alcohol exploratory ligand\n",
        encoding="utf-8"
    )

print(f"Manifest written: {INPUT / 'structural_manifest.tsv'}")
print(f"FASTA written: {INPUT / 'structural_priority_plus_references.faa'}")
PY

# ------------------------------------------------------------
# 2. Try to reuse manual models
# ------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import pandas as pd
import shutil

BASE = Path("<PROJECT_ROOT>")
OUT = BASE / "19_structural_pipeline_auto"
manifest = pd.read_csv(OUT / "01_inputs" / "structural_manifest.tsv", sep="\t")
manual = OUT / "02_models" / "manual_models"
pdb_dir = OUT / "02_models" / "pdb"
pdb_dir.mkdir(parents=True, exist_ok=True)

rows = []
for _, row in manifest.iterrows():
    safe = row["safe_id"]
    src = manual / f"{safe}.pdb"
    dst = pdb_dir / f"{safe}.pdb"
    if src.exists() and src.stat().st_size > 0:
        shutil.copy2(src, dst)
        rows.append({"safe_id": safe, "source": "manual_models", "pdb_exists": 1, "pdb_path": str(dst)})
    else:
        rows.append({"safe_id": safe, "source": "manual_models", "pdb_exists": 0, "pdb_path": ""})

pd.DataFrame(rows).to_csv(OUT / "06_tables" / "initial_model_status.tsv", sep="\t", index=False)
print("Initial manual-model status saved.")
PY

# ------------------------------------------------------------
# 3. Try automatic modeling with ColabFold if available
# ------------------------------------------------------------
if [[ -n "${COLABFOLD}" ]]; then
  echo "[3] ColabFold detected, attempting automatic modeling..." | tee -a "${STATUS_FILE}"
  "${COLABFOLD}" \
    "${INPUT_DIR}/structural_priority_plus_references.faa" \
    "${AUTO_MODEL_DIR}" \
    > "${LOG_DIR}/colabfold_stdout.log" 2> "${LOG_DIR}/colabfold_stderr.log" || true
else
  echo "[3] ColabFold not detected." | tee -a "${STATUS_FILE}"
fi

# ------------------------------------------------------------
# 4. Try automatic modeling with local ESMFold if available
# ------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import importlib.util
import subprocess
import sys
import pandas as pd
from Bio import SeqIO

BASE = Path("<PROJECT_ROOT>")
OUT = BASE / "19_structural_pipeline_auto"
INPUT = OUT / "01_inputs"
AUTO = OUT / "02_models" / "auto_models"
PDB_DIR = OUT / "02_models" / "pdb"
LOG = OUT / "07_logs"
STATUS = OUT / "08_status" / "esmfold_status.txt"

AUTO.mkdir(parents=True, exist_ok=True)
PDB_DIR.mkdir(parents=True, exist_ok=True)

have_torch = importlib.util.find_spec("torch") is not None
have_esm = importlib.util.find_spec("esm") is not None

lines = []
lines.append(f"torch_available={have_torch}")
lines.append(f"esm_available={have_esm}")

if not (have_torch and have_esm):
    STATUS.write_text("\n".join(lines), encoding="utf-8")
    print("ESMFold local modeling not available.")
    sys.exit(0)

# Local esmfold inference
try:
    import torch
    import esm
    model = esm.pretrained.esmfold_v1()
    model = model.eval()
    if torch.cuda.is_available():
        model = model.cuda()
        lines.append("device=cuda")
    else:
        lines.append("device=cpu")

    manifest = pd.read_csv(INPUT / "structural_manifest.tsv", sep="\t")
    for _, row in manifest.iterrows():
        safe = row["safe_id"]
        dst = PDB_DIR / f"{safe}.pdb"
        if dst.exists() and dst.stat().st_size > 0:
            continue
        fa = INPUT / f"{safe}.faa"
        if not fa.exists():
            continue
        recs = list(SeqIO.parse(fa, "fasta"))
        if not recs:
            continue
        seq = str(recs[0].seq)
        with torch.no_grad():
            pdb_string = model.infer_pdb(seq)
        outp = AUTO / f"{safe}.pdb"
        outp.write_text(pdb_string, encoding="utf-8")
        dst.write_text(pdb_string, encoding="utf-8")
        lines.append(f"modeled={safe}")

except Exception as e:
    lines.append(f"esmfold_error={e}")

STATUS.write_text("\n".join(lines), encoding="utf-8")
print("ESMFold automatic modeling step finished.")
PY

# ------------------------------------------------------------
# 5. Collect ColabFold outputs if created
# ------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import pandas as pd
import shutil

BASE = Path("<PROJECT_ROOT>")
OUT = BASE / "19_structural_pipeline_auto"
manifest = pd.read_csv(OUT / "01_inputs" / "structural_manifest.tsv", sep="\t")
auto = OUT / "02_models" / "auto_models"
pdb_dir = OUT / "02_models" / "pdb"
pdb_dir.mkdir(parents=True, exist_ok=True)

for _, row in manifest.iterrows():
    safe = row["safe_id"]
    dst = pdb_dir / f"{safe}.pdb"
    if dst.exists() and dst.stat().st_size > 0:
        continue
    cands = list(auto.rglob(f"{safe}*rank_001*.pdb")) + list(auto.rglob(f"{safe}.pdb")) + list(auto.rglob(f"{safe}*.pdb"))
    if cands:
        shutil.copy2(cands[0], dst)

rows = []
for _, row in manifest.iterrows():
    safe = row["safe_id"]
    p = pdb_dir / f"{safe}.pdb"
    rows.append({
        "qseqid": row["qseqid"],
        "safe_id": safe,
        "record_type": row["record_type"],
        "pdb_exists": int(p.exists() and p.stat().st_size > 0),
        "pdb_path": str(p) if p.exists() and p.stat().st_size > 0 else ""
    })

pd.DataFrame(rows).to_csv(OUT / "06_tables" / "final_model_status.tsv", sep="\t", index=False)
print("Final model-status table written.")
PY

# ------------------------------------------------------------
# 6. Require at least all candidates modeled
# ------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import pandas as pd
import sys

BASE = Path("<PROJECT_ROOT>")
OUT = BASE / "19_structural_pipeline_auto"
df = pd.read_csv(OUT / "06_tables" / "final_model_status.tsv", sep="\t")

cand = df[df["record_type"] == "candidate"].copy()
missing = cand[cand["pdb_exists"] != 1]

if not missing.empty:
    print("MISSING_CANDIDATE_MODELS_AFTER_AUTO")
    print(missing[["safe_id"]].to_string(index=False))
    sys.exit(2)
else:
    print("ALL_CANDIDATE_MODELS_AVAILABLE")
PY
RC=$?

if [[ "${RC}" -eq 2 ]]; then
  echo "ERROR: automatic modeling did not recover all candidate PDBs." | tee -a "${STATUS_FILE}"
  echo "Place missing models manually into:" | tee -a "${STATUS_FILE}"
  echo "  ${MANUAL_MODEL_DIR}" | tee -a "${STATUS_FILE}"
  echo "Then rerun this same script." | tee -a "${STATUS_FILE}"
  exit 2
fi

# ------------------------------------------------------------
# 7. Require structural tools for downstream steps
# ------------------------------------------------------------
for cmd in fpocket obabel vina; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: downstream structural command not found: ${cmd}" | tee -a "${STATUS_FILE}"
    exit 3
  fi
done

# ------------------------------------------------------------
# 8. Pocket analysis
# ------------------------------------------------------------
echo "[8] Running fpocket..." | tee -a "${STATUS_FILE}"
for pdb in "${PDB_DIR}"/*.pdb; do
  [[ -e "$pdb" ]] || continue
  bn="$(basename "$pdb" .pdb)"
  workdir="${POCKET_DIR}/${bn}"
  mkdir -p "${workdir}"
  cp "$pdb" "${workdir}/${bn}.pdb"
  (
    cd "${workdir}"
    fpocket -f "${bn}.pdb" > "${LOG_DIR}/${bn}_fpocket_stdout.log" 2> "${LOG_DIR}/${bn}_fpocket_stderr.log"
  )
done

# ------------------------------------------------------------
# 9. Prepare pocket table and docking jobs
# ------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import pandas as pd

BASE = Path("<PROJECT_ROOT>")
OUT = BASE / "19_structural_pipeline_auto"
POCKET = OUT / "03_pockets"
TAB = OUT / "06_tables"
INPUT = OUT / "01_inputs"

rows = []
for model_dir in POCKET.iterdir():
    if not model_dir.is_dir():
        continue
    outdirs = [d for d in model_dir.iterdir() if d.is_dir() and d.name.endswith("_out")]
    for od in outdirs:
        pdir = od / "pockets"
        if not pdir.exists():
            continue
        for pf in sorted(pdir.glob("pocket*_atm.pdb")):
            xs, ys, zs = [], [], []
            with open(pf, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        try:
                            xs.append(float(line[30:38]))
                            ys.append(float(line[38:46]))
                            zs.append(float(line[46:54]))
                        except Exception:
                            pass
            if xs:
                rows.append({
                    "model_id": model_dir.name,
                    "pocket_name": pf.stem,
                    "center_x": sum(xs)/len(xs),
                    "center_y": sum(ys)/len(ys),
                    "center_z": sum(zs)/len(zs),
                    "n_atoms": len(xs),
                    "pocket_file": str(pf)
                })

pocket_df = pd.DataFrame(rows)
pocket_df.to_csv(TAB / "pocket_centroids.tsv", sep="\t", index=False)

lig = pd.read_csv(INPUT / "ligands.tsv", sep="\t")
lig = lig[(lig["enabled"] == 1) & lig["smiles"].fillna("").ne("")].copy()
lig.to_csv(TAB / "enabled_ligands.tsv", sep="\t", index=False)

jobs = []
if not pocket_df.empty and not lig.empty:
    top = pocket_df.sort_values(["model_id", "n_atoms"], ascending=[True, False]).drop_duplicates("model_id")
    for _, prow in top.iterrows():
        for _, lrow in lig.iterrows():
            jobs.append({
                "model_id": prow["model_id"],
                "pocket_name": prow["pocket_name"],
                "center_x": prow["center_x"],
                "center_y": prow["center_y"],
                "center_z": prow["center_z"],
                "ligand_name": lrow["ligand_name"],
                "smiles": lrow["smiles"]
            })

pd.DataFrame(jobs).to_csv(TAB / "docking_jobs.tsv", sep="\t", index=False)
print("Pocket and docking-job tables written.")
PY

# ------------------------------------------------------------
# 10. Prepare ligands and receptors
# ------------------------------------------------------------
mkdir -p "${DOCK_DIR}/ligands" "${DOCK_DIR}/receptors"

tail -n +2 "${TAB_DIR}/enabled_ligands.tsv" | while IFS=$'\t' read -r ligand_name smiles enabled comment; do
  [[ -n "${ligand_name}" ]] || continue
  obabel -:"${smiles}" -osdf -O "${DOCK_DIR}/ligands/${ligand_name}.sdf" --gen3d \
    > "${LOG_DIR}/${ligand_name}_lig_stdout.log" 2> "${LOG_DIR}/${ligand_name}_lig_stderr.log"
  obabel "${DOCK_DIR}/ligands/${ligand_name}.sdf" -opdbqt -O "${DOCK_DIR}/ligands/${ligand_name}.pdbqt" \
    >> "${LOG_DIR}/${ligand_name}_lig_stdout.log" 2>> "${LOG_DIR}/${ligand_name}_lig_stderr.log"
  obabel "${DOCK_DIR}/ligands/${ligand_name}.sdf" -opdb -O "${DOCK_DIR}/ligands/${ligand_name}.pdb" \
    >> "${LOG_DIR}/${ligand_name}_lig_stdout.log" 2>> "${LOG_DIR}/${ligand_name}_lig_stderr.log"
done

for pdb in "${PDB_DIR}"/*.pdb; do
  [[ -e "$pdb" ]] || continue
  bn="$(basename "$pdb" .pdb)"
  obabel "$pdb" -opdbqt -O "${DOCK_DIR}/receptors/${bn}.pdbqt" \
    > "${LOG_DIR}/${bn}_rec_stdout.log" 2> "${LOG_DIR}/${bn}_rec_stderr.log"
done

# ------------------------------------------------------------
# 11. Run docking
# ------------------------------------------------------------
echo "[11] Running AutoDock Vina..." | tee -a "${STATUS_FILE}"
tail -n +2 "${TAB_DIR}/docking_jobs.tsv" | while IFS=$'\t' read -r model_id pocket_name center_x center_y center_z ligand_name smiles; do
  [[ -n "${model_id}" ]] || continue
  receptor="${DOCK_DIR}/receptors/${model_id}.pdbqt"
  ligand="${DOCK_DIR}/ligands/${ligand_name}.pdbqt"
  out_pose="${DOCK_DIR}/${model_id}__${ligand_name}__${pocket_name}.pdbqt"
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
    --log "${out_log}" \
    > "${LOG_DIR}/${model_id}__${ligand_name}_vina_stdout.log" \
    2> "${LOG_DIR}/${model_id}__${ligand_name}_vina_stderr.log"

  obabel "${out_pose}" -opdb -O "${DOCK_DIR}/${model_id}__${ligand_name}__${pocket_name}.pdb" \
    >> "${LOG_DIR}/${model_id}__${ligand_name}_vina_stdout.log" \
    2>> "${LOG_DIR}/${model_id}__${ligand_name}_vina_stderr.log"
done

# ------------------------------------------------------------
# 12. Summarize docking and generate visualization scripts
# ------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import re
import pandas as pd

BASE = Path("<PROJECT_ROOT>")
OUT = BASE / "19_structural_pipeline_auto"
DOCK = OUT / "04_docking"
PDB = OUT / "02_models" / "pdb"
VIS = OUT / "05_visualization"
TAB = OUT / "06_tables"

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
            line = line.strip()
            if re.match(r"^\s*1\s+", line):
                toks = line.split()
                if len(toks) >= 2:
                    try:
                        score = float(toks[1])
                    except Exception:
                        pass
                break

    pose_pdb = DOCK / f"{model_id}__{ligand_name}__{pocket_name}.pdb"
    receptor_pdb = PDB / f"{model_id}.pdb"

    rows.append({
        "model_id": model_id,
        "ligand_name": ligand_name,
        "pocket_name": pocket_name,
        "vina_best_score_kcal_mol": score,
        "pose_exists": int(pose_pdb.exists()),
        "pose_pdb": str(pose_pdb),
        "receptor_pdb": str(receptor_pdb)
    })

dock_df = pd.DataFrame(rows)
dock_df.to_csv(TAB / "docking_summary.tsv", sep="\t", index=False)

for _, row in dock_df.iterrows():
    if row["pose_exists"] != 1:
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

print("Docking summary and PyMOL scripts generated.")
PY

# Optional render
if [[ -n "${PYMOL}" ]]; then
  for pml in "${VIS_DIR}"/*.pml; do
    [[ -e "$pml" ]] || continue
    pymol -cq "$pml" > "${LOG_DIR}/$(basename "$pml" .pml)_pymol_stdout.log" 2> "${LOG_DIR}/$(basename "$pml" .pml)_pymol_stderr.log" || true
  done
fi

{
  echo "Automatic structural pipeline completed."
  echo "Key outputs:"
  echo "  ${TAB_DIR}/final_model_status.tsv"
  echo "  ${TAB_DIR}/pocket_centroids.tsv"
  echo "  ${TAB_DIR}/docking_jobs.tsv"
  echo "  ${TAB_DIR}/docking_summary.tsv"
  echo "  ${VIS_DIR}"
} | tee "${SUMMARY_FILE}"
