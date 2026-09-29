#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="."
IN_FASTA="${BASE_DIR}/08_priority_candidates/priority_candidates_plus_references.faa"
OUT_DIR="${BASE_DIR}/10_alignment"
OUT_FASTA="${OUT_DIR}/priority_candidates_plus_references_mafft.faa"

mkdir -p "${OUT_DIR}"

if [[ ! -s "${IN_FASTA}" ]]; then
    echo "ERROR: input FASTA missing:"
    echo "  ${IN_FASTA}"
    exit 1
fi

if ! command -v mafft >/dev/null 2>&1; then
    echo "ERROR: mafft not found in PATH."
    exit 1
fi

echo "Running MAFFT alignment..."
mafft --auto "${IN_FASTA}" > "${OUT_FASTA}"

echo "Alignment written to: ${OUT_FASTA}"
