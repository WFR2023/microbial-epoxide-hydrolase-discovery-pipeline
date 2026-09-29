#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="."
REF_FASTA="${BASE_DIR}/07_blast_reference/epoxide_hydrolase_references.faa"
QUERY_STRONG="${BASE_DIR}/06_candidate_mining/deduplicated/strong_deduplicated.faa"
QUERY_MODERATE="${BASE_DIR}/06_candidate_mining/deduplicated/moderate_deduplicated.faa"
OUT_DIR="${BASE_DIR}/07_blast_results"
DB_PREFIX="${OUT_DIR}/ref_db"

mkdir -p "${OUT_DIR}"

if [[ ! -s "${REF_FASTA}" ]]; then
    echo "ERROR: Reference FASTA missing or empty:"
    echo "  ${REF_FASTA}"
    exit 1
fi

if [[ ! -s "${QUERY_STRONG}" ]]; then
    echo "ERROR: Strong query FASTA missing:"
    echo "  ${QUERY_STRONG}"
    exit 1
fi

if [[ ! -s "${QUERY_MODERATE}" ]]; then
    echo "ERROR: Moderate query FASTA missing:"
    echo "  ${QUERY_MODERATE}"
    exit 1
fi

echo "Building BLAST database..."
makeblastdb -in "${REF_FASTA}" -dbtype prot -out "${DB_PREFIX}"

echo "Running BLASTp for strong candidates..."
blastp \
  -query "${QUERY_STRONG}" \
  -db "${DB_PREFIX}" \
  -out "${OUT_DIR}/strong_vs_refs.tsv" \
  -evalue 1e-5 \
  -max_target_seqs 10 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

echo "Running BLASTp for moderate candidates..."
blastp \
  -query "${QUERY_MODERATE}" \
  -db "${DB_PREFIX}" \
  -out "${OUT_DIR}/moderate_vs_refs.tsv" \
  -evalue 1e-5 \
  -max_target_seqs 10 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

echo "BLASTp completed."
echo "Outputs:"
echo "  ${OUT_DIR}/strong_vs_refs.tsv"
echo "  ${OUT_DIR}/moderate_vs_refs.tsv"
