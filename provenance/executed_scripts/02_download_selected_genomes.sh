#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="<PROJECT_ROOT>"
ACCESSIONS_FILE="${BASE_DIR}/02_accessions/selected_accessions.txt"
DOWNLOAD_DIR="${BASE_DIR}/03_downloads"
LOG_DIR="${BASE_DIR}/logs"

ZIP_FILE="${DOWNLOAD_DIR}/ncbi_dataset.zip"
EXTRACT_DIR="${DOWNLOAD_DIR}/ncbi_dataset"

mkdir -p "${DOWNLOAD_DIR}" "${LOG_DIR}"

echo "============================================================"
echo "NCBI genome download pipeline"
echo "============================================================"
echo "Base directory   : ${BASE_DIR}"
echo "Accessions file  : ${ACCESSIONS_FILE}"
echo "Download dir     : ${DOWNLOAD_DIR}"
echo "ZIP output       : ${ZIP_FILE}"
echo "Extract dir      : ${EXTRACT_DIR}"
echo "============================================================"

if [[ ! -s "${ACCESSIONS_FILE}" ]]; then
    echo "ERROR: Accessions file not found or empty:"
    echo "  ${ACCESSIONS_FILE}"
    exit 1
fi

if ! command -v datasets >/dev/null 2>&1; then
    echo "ERROR: NCBI Datasets CLI is not installed or not in PATH."
    echo "Install it first, for example:"
    echo "  conda install -c conda-forge ncbi-datasets-cli -y"
    exit 1
fi

if ! command -v unzip >/dev/null 2>&1; then
    echo "ERROR: unzip is not installed."
    echo "Install it with:"
    echo "  sudo apt-get update && sudo apt-get install -y unzip"
    exit 1
fi

echo "Selected accession count:"
wc -l "${ACCESSIONS_FILE}"

echo ""
echo "Preview of the first selected accessions:"
head -n 10 "${ACCESSIONS_FILE}" || true

cd "${DOWNLOAD_DIR}"

# Clean previous package to avoid confusion
if [[ -f "${ZIP_FILE}" ]]; then
    echo ""
    echo "Removing previous ZIP package..."
    rm -f "${ZIP_FILE}"
fi

# Clean previous extracted folder to ensure reproducibility
if [[ -d "${EXTRACT_DIR}" ]]; then
    echo "Removing previous extracted dataset folder..."
    rm -rf "${EXTRACT_DIR}"
fi

echo ""
echo "Starting NCBI Datasets download..."
echo "This may take some time depending on the number and size of assemblies."

datasets download genome accession \
    --inputfile "${ACCESSIONS_FILE}" \
    --include genome,protein,cds,gff3,gbff,seq-report \
    --filename "$(basename "${ZIP_FILE}")" \
    --no-progressbar \
    2>&1 | tee "${LOG_DIR}/02_datasets_download.log"

if [[ ! -f "${ZIP_FILE}" ]]; then
    echo "ERROR: Download did not produce ${ZIP_FILE}"
    exit 1
fi

echo ""
echo "Download completed. Unzipping package..."
unzip -q -o "${ZIP_FILE}" -d "${DOWNLOAD_DIR}"

if [[ ! -d "${EXTRACT_DIR}" ]]; then
    echo "ERROR: Extraction folder not found after unzip:"
    echo "  ${EXTRACT_DIR}"
    exit 1
fi

echo ""
echo "Download and extraction completed successfully."
echo "ZIP file        : ${ZIP_FILE}"
echo "Extracted folder: ${EXTRACT_DIR}"

echo ""
echo "Top-level extracted contents:"
find "${EXTRACT_DIR}" -maxdepth 2 -type f | head -n 20 || true

echo ""
echo "Done."
