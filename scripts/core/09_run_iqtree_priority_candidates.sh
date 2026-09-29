#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="."
ALIGNMENT="${BASE_DIR}/10_alignment/priority_candidates_plus_references_mafft.faa"
OUT_DIR="${BASE_DIR}/11_phylogeny"

mkdir -p "${OUT_DIR}"

if [[ ! -s "${ALIGNMENT}" ]]; then
    echo "ERROR: alignment file missing:"
    echo "  ${ALIGNMENT}"
    exit 1
fi

if command -v iqtree2 >/dev/null 2>&1; then
    IQTREE_BIN="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
    IQTREE_BIN="iqtree"
else
    echo "ERROR: neither iqtree2 nor iqtree was found in PATH."
    exit 1
fi

cd "${OUT_DIR}"

echo "Running ${IQTREE_BIN}..."
"${IQTREE_BIN}" \
  -s "${ALIGNMENT}" \
  -m MFP \
  -bb 1000 \
  -alrt 1000 \
  -nt AUTO \
  --prefix priority_candidates_tree

echo "IQ-TREE completed."
echo "Output directory: ${OUT_DIR}"
