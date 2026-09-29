#!/usr/bin/env bash
set -euo pipefail

BASE="<PROJECT_ROOT>/19_structural_pipeline_real"
SRC="${BASE}/10_colabfold_models"
DST="${BASE}/02_models/manual_models"
MANIFEST="${BASE}/01_inputs/structural_manifest.tsv"
OUTTAB="${BASE}/06_tables/colabfold_model_collection_status.tsv"

mkdir -p "${DST}"
mkdir -p "$(dirname "${OUTTAB}")"

python - <<'PY'
from pathlib import Path
import pandas as pd
import shutil

BASE = Path("<PROJECT_ROOT>/19_structural_pipeline_real")
SRC = BASE / "10_colabfold_models"
DST = BASE / "02_models" / "manual_models"
MANIFEST = BASE / "01_inputs" / "structural_manifest.tsv"
OUTTAB = BASE / "06_tables" / "colabfold_model_collection_status.tsv"

manifest = pd.read_csv(MANIFEST, sep="\t")
rows = []

for _, row in manifest.iterrows():
    safe = row["safe_id"]
    target = DST / f"{safe}.pdb"

    if target.exists() and target.stat().st_size > 0:
        rows.append({"safe_id": safe, "status": "already_present", "source_file": str(target)})
        continue

    candidates = (
        list(SRC.rglob(f"{safe}*relaxed_rank_001*.pdb")) +
        list(SRC.rglob(f"{safe}*unrelaxed_rank_001*.pdb")) +
        list(SRC.rglob(f"{safe}*rank_001*.pdb")) +
        list(SRC.rglob(f"{safe}*.pdb"))
    )

    if candidates:
        best = candidates[0]
        shutil.copy2(best, target)
        rows.append({"safe_id": safe, "status": "copied", "source_file": str(best)})
    else:
        rows.append({"safe_id": safe, "status": "missing", "source_file": ""})

out = pd.DataFrame(rows)
out.to_csv(OUTTAB, sep="\t", index=False)
print(out.to_string(index=False))
print(f"\nStatus table written to: {OUTTAB}")
PY
