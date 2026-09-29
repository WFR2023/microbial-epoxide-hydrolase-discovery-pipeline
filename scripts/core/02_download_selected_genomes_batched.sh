#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="."
ACCESSIONS_FILE="${BASE_DIR}/02_accessions/selected_accessions.txt"
DOWNLOAD_DIR="${BASE_DIR}/03_downloads"
BATCH_DIR="${DOWNLOAD_DIR}/batches"
EXTRACT_DIR="${DOWNLOAD_DIR}/extracted_batches"
LOG_DIR="${BASE_DIR}/logs"

BATCH_SIZE=10
MAX_RETRIES=3

mkdir -p "${DOWNLOAD_DIR}" "${BATCH_DIR}" "${EXTRACT_DIR}" "${LOG_DIR}"

echo "============================================================"
echo "NCBI batched genome download pipeline"
echo "============================================================"
echo "Base directory   : ${BASE_DIR}"
echo "Accessions file  : ${ACCESSIONS_FILE}"
echo "Batch size       : ${BATCH_SIZE}"
echo "Max retries      : ${MAX_RETRIES}"
echo "Batch dir        : ${BATCH_DIR}"
echo "Extract dir      : ${EXTRACT_DIR}"
echo "============================================================"

if [[ ! -s "${ACCESSIONS_FILE}" ]]; then
    echo "ERROR: Accessions file not found or empty:"
    echo "  ${ACCESSIONS_FILE}"
    exit 1
fi

if ! command -v datasets >/dev/null 2>&1; then
    echo "ERROR: NCBI Datasets CLI not found in PATH."
    exit 1
fi

if ! command -v unzip >/dev/null 2>&1; then
    echo "ERROR: unzip not found."
    exit 1
fi

echo "Total accession count:"
wc -l "${ACCESSIONS_FILE}"

echo
echo "Cleaning previous batch files..."
rm -f "${BATCH_DIR}"/batch_*.txt
rm -f "${BATCH_DIR}"/batch_*.zip
rm -rf "${EXTRACT_DIR}"/batch_*

echo
echo "Splitting accession list into batches of ${BATCH_SIZE}..."
split -d -l "${BATCH_SIZE}" "${ACCESSIONS_FILE}" "${BATCH_DIR}/batch_"

BATCH_LIST=$(ls "${BATCH_DIR}"/batch_* | grep -v '\.zip$' || true)

if [[ -z "${BATCH_LIST}" ]]; then
    echo "ERROR: No batch files were created."
    exit 1
fi

SUCCESS_LOG="${LOG_DIR}/02_batched_download_success.tsv"
FAIL_LOG="${LOG_DIR}/02_batched_download_fail.tsv"

echo -e "batch_file\tzip_file\tn_accessions\tstatus\tattempts" > "${SUCCESS_LOG}"
echo -e "batch_file\tzip_file\tn_accessions\tstatus\tattempts" > "${FAIL_LOG}"

download_batch () {
    local batch_txt="$1"
    local zip_file="$2"
    local attempt="$3"

    datasets download genome accession \
        --inputfile "${batch_txt}" \
        --include genome,protein,cds,gff3,gbff,seq-report \
        --filename "${zip_file}" \
        --no-progressbar
}

echo
echo "Starting batched downloads..."

for batch_txt in ${BATCH_LIST}; do
    batch_name="$(basename "${batch_txt}")"
    zip_file="${BATCH_DIR}/${batch_name}.zip"
    out_dir="${EXTRACT_DIR}/${batch_name}"
    n_accessions=$(wc -l < "${batch_txt}")

    echo
    echo "------------------------------------------------------------"
    echo "Batch: ${batch_name}"
    echo "Accessions: ${n_accessions}"
    echo "TXT: ${batch_txt}"
    echo "ZIP: ${zip_file}"
    echo "OUT: ${out_dir}"
    echo "------------------------------------------------------------"

    success=0

    for attempt in $(seq 1 "${MAX_RETRIES}"); do
        echo "Attempt ${attempt}/${MAX_RETRIES} for ${batch_name}..."

        rm -f "${zip_file}"
        rm -rf "${out_dir}"

        if download_batch "${batch_txt}" "${zip_file}" "${attempt}" \
            > "${LOG_DIR}/${batch_name}.attempt${attempt}.log" 2>&1; then

            if [[ -s "${zip_file}" ]]; then
                mkdir -p "${out_dir}"

                if unzip -q -o "${zip_file}" -d "${out_dir}" \
                    >> "${LOG_DIR}/${batch_name}.attempt${attempt}.log" 2>&1; then
                    echo "Batch ${batch_name} completed successfully."
                    echo -e "${batch_txt}\t${zip_file}\t${n_accessions}\tsuccess\t${attempt}" >> "${SUCCESS_LOG}"
                    success=1
                    break
                else
                    echo "Unzip failed for ${batch_name} on attempt ${attempt}."
                fi
            else
                echo "ZIP file missing or empty for ${batch_name} on attempt ${attempt}."
            fi
        else
            echo "Download failed for ${batch_name} on attempt ${attempt}."
        fi

        echo "Sleeping before retry..."
        sleep 5
    done

    if [[ "${success}" -eq 0 ]]; then
        echo "Batch ${batch_name} failed after ${MAX_RETRIES} attempts."
        echo -e "${batch_txt}\t${zip_file}\t${n_accessions}\tfailed\t${MAX_RETRIES}" >> "${FAIL_LOG}"
    fi
done

echo
echo "============================================================"
echo "Batched download completed"
echo "============================================================"

echo
echo "Successful batches:"
cat "${SUCCESS_LOG}"

echo
echo "Failed batches:"
cat "${FAIL_LOG}"

echo
echo "Extracted batch folders:"
find "${EXTRACT_DIR}" -maxdepth 1 -type d | sort
