#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="."
PIPE_DIR="${BASE_DIR}/19_structural_pipeline_real"
INPUT_DIR="${PIPE_DIR}/01_inputs"
MODEL_DIR="${PIPE_DIR}/02_models"
MANUAL_DIR="${MODEL_DIR}/manual_models"
PREP_DIR="${PIPE_DIR}/09_model_prep"
FASTA_DIR="${PREP_DIR}/candidate_fastas"
REF_DIR="${PREP_DIR}/reference_downloads"
mkdir -p "${FASTA_DIR}" "${REF_DIR}" "${MANUAL_DIR}"

MANIFEST="${INPUT_DIR}/structural_manifest.tsv"
COMBINED_FASTA="${INPUT_DIR}/structural_priority_plus_references.faa"
STATUS="${PREP_DIR}/model_prep_status.txt"
CHECKLIST="${PREP_DIR}/missing_model_checklist.tsv"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "ERROR: manifest not found: ${MANIFEST}"
  exit 1
fi

if [[ ! -f "${COMBINED_FASTA}" ]]; then
  echo "ERROR: combined FASTA not found: ${COMBINED_FASTA}"
  exit 1
fi

python - <<'PY'
from pathlib import Path
import pandas as pd
from Bio import SeqIO

BASE = Path(".")
PIPE = BASE / "19_structural_pipeline_real"
INPUT = PIPE / "01_inputs"
MANUAL = PIPE / "02_models" / "manual_models"
PREP = PIPE / "09_model_prep"
FASTA_DIR = PREP / "candidate_fastas"
CHECKLIST = PREP / "missing_model_checklist.tsv"

manifest = pd.read_csv(INPUT / "structural_manifest.tsv", sep="\t")
seqs = {rec.id: rec for rec in SeqIO.parse(INPUT / "structural_priority_plus_references.faa", "fasta")}

rows = []

for _, row in manifest.iterrows():
    qseqid = row["qseqid"]
    safe_id = row["safe_id"]
    record_type = row["record_type"]
    expected_pdb = MANUAL / f"{safe_id}.pdb"
    exists = expected_pdb.exists()

    if qseqid in seqs:
        SeqIO.write([seqs[qseqid]], FASTA_DIR / f"{safe_id}.faa", "fasta")

    rows.append({
        "record_type": record_type,
        "qseqid": qseqid,
        "safe_id": safe_id,
        "expected_pdb": str(expected_pdb),
        "pdb_exists": int(exists),
        "candidate_fasta": str(FASTA_DIR / f"{safe_id}.faa")
    })

pd.DataFrame(rows).to_csv(CHECKLIST, sep="\t", index=False)
print(f"Checklist written to: {CHECKLIST}")
print(f"Candidate FASTAs written to: {FASTA_DIR}")
PY

# Download reference structures when possible
# Human EPHX2 AlphaFold
wget -O "${REF_DIR}/human_EPHX2_sEH_alphafold.pdb" \
  "https://alphafold.ebi.ac.uk/files/AF-P34913-F1-model_v4.pdb" \
  > "${REF_DIR}/wget_ephx2_stdout.log" 2> "${REF_DIR}/wget_ephx2_stderr.log" || true

# Human EPHX1 AlphaFold
wget -O "${REF_DIR}/human_EPHX1_mEH_alphafold.pdb" \
  "https://alphafold.ebi.ac.uk/files/AF-P07099-F1-model_v4.pdb" \
  > "${REF_DIR}/wget_ephx1_stdout.log" 2> "${REF_DIR}/wget_ephx1_stderr.log" || true

# Pseudomonas Cif experimental structure if available from RCSB direct download
# If this specific ID fails, you can replace manually later.
wget -O "${REF_DIR}/pseudomonas_Cif_rcsb.pdb" \
  "https://files.rcsb.org/download/3KD2.pdb" \
  > "${REF_DIR}/wget_cif_stdout.log" 2> "${REF_DIR}/wget_cif_stderr.log" || true

# Copy downloaded refs into manual_models if non-empty
for src in \
  "${REF_DIR}/human_EPHX2_sEH_alphafold.pdb" \
  "${REF_DIR}/human_EPHX1_mEH_alphafold.pdb" \
  "${REF_DIR}/pseudomonas_Cif_rcsb.pdb"
do
  if [[ -s "$src" ]]; then
    case "$(basename "$src")" in
      human_EPHX2_sEH_alphafold.pdb) cp "$src" "${MANUAL_DIR}/human_EPHX2_sEH.pdb" ;;
      human_EPHX1_mEH_alphafold.pdb) cp "$src" "${MANUAL_DIR}/human_EPHX1_mEH.pdb" ;;
      pseudomonas_Cif_rcsb.pdb) cp "$src" "${MANUAL_DIR}/pseudomonas_Cif.pdb" ;;
    esac
  fi
done

{
  echo "Model preparation completed."
  echo "Checklist: ${CHECKLIST}"
  echo "Candidate FASTAs: ${FASTA_DIR}"
  echo "Reference downloads: ${REF_DIR}"
  echo "Manual models dir: ${MANUAL_DIR}"
  echo
  echo "Next action:"
  echo "1. Check which candidate PDBs are still missing in ${CHECKLIST}"
  echo "2. Generate those models externally (e.g., ColabFold web)"
  echo "3. Save them into ${MANUAL_DIR} with the exact safe_id.pdb filename"
  echo "4. Re-run 19b_run_structural_pipeline_real.sh"
} > "${STATUS}"

cat "${STATUS}"
