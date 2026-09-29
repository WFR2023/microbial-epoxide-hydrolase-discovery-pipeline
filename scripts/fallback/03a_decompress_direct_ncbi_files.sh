#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="."
IN_ROOT="${BASE_DIR}/03_downloads/direct_ncbi_files"
LOG_DIR="${BASE_DIR}/logs"

mkdir -p "${LOG_DIR}"

if [[ ! -d "${IN_ROOT}" ]]; then
    echo "ERROR: Input directory not found: ${IN_ROOT}"
    exit 1
fi

echo "============================================================"
echo "Decompressing direct NCBI files"
echo "============================================================"
echo "Input root: ${IN_ROOT}"
echo "============================================================"

find "${IN_ROOT}" -type f \( -name "*.gz" \) | while read -r gzfile; do
    out="${gzfile%.gz}"

    if [[ -f "${out}" && -s "${out}" ]]; then
        echo "Skipping already decompressed file: ${out}"
        continue
    fi

    echo "Decompressing: ${gzfile}"
    gunzip -c "${gzfile}" > "${out}"
done

echo
echo "Decompression completed."
