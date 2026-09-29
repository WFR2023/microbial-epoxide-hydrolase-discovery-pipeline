#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="<PROJECT_ROOT>"
PFAM_DB="<PFAM_DB>/Pfam-A.hmm"
QUERY_FASTA="${BASE_DIR}/08_priority_candidates/priority_candidates_plus_references.faa"
OUT_DIR="${BASE_DIR}/13_hmmer_pfam"

mkdir -p "${OUT_DIR}"

DOMTBLOUT="${OUT_DIR}/priority_candidates_vs_pfam.domtblout"
TBLOUT="${OUT_DIR}/priority_candidates_vs_pfam.tblout"
LOGFILE="${OUT_DIR}/priority_candidates_vs_pfam.log"

if [[ ! -s "${PFAM_DB}" ]]; then
    echo "ERROR: Pfam database not found:"
    echo "  ${PFAM_DB}"
    exit 1
fi

if [[ ! -s "${QUERY_FASTA}" ]]; then
    echo "ERROR: Query FASTA not found:"
    echo "  ${QUERY_FASTA}"
    exit 1
fi

if ! command -v hmmscan >/dev/null 2>&1; then
    echo "ERROR: hmmscan not found in PATH."
    exit 1
fi

echo "Running hmmscan against Pfam..."
hmmscan \
  --cpu 2 \
  --tblout "${TBLOUT}" \
  --domtblout "${DOMTBLOUT}" \
  "${PFAM_DB}" \
  "${QUERY_FASTA}" \
  > "${LOGFILE}"

echo "hmmscan completed."
echo "Outputs:"
echo "  ${TBLOUT}"
echo "  ${DOMTBLOUT}"
echo "  ${LOGFILE}"
