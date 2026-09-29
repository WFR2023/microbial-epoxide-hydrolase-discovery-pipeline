#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="<PROJECT_ROOT>"
MANIFEST="${BASE_DIR}/03_downloads/ncbi_download_manifest.tsv"
OUT_ROOT="${BASE_DIR}/03_downloads/direct_ncbi_files"
LOG_DIR="${BASE_DIR}/logs"

mkdir -p "${OUT_ROOT}" "${LOG_DIR}"

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Manifest not found or empty: ${MANIFEST}"
    exit 1
fi

if ! command -v wget >/dev/null 2>&1; then
    echo "ERROR: wget not found."
    echo "Install with: sudo apt-get update && sudo apt-get install -y wget"
    exit 1
fi

echo "============================================================"
echo "Direct NCBI download from manifest"
echo "============================================================"
echo "Manifest : ${MANIFEST}"
echo "Output   : ${OUT_ROOT}"
echo "============================================================"

python - <<'PY'
import pandas as pd
df = pd.read_csv("<PROJECT_ROOT>/03_downloads/ncbi_download_manifest.tsv", sep="\t")
print(f"Assemblies to download: {len(df)}")
print(df[["assembly_accession","target_species","source_database"]].head(10).to_string(index=False))
PY

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r \
    target_group target_species assembly_accession assembly_name source_database ftp_base https_base \
    genome_fna_gz protein_faa_gz cds_fna_gz gff_gz gbff_gz assembly_report_txt md5checksums_txt
do
    species_safe=$(echo "${target_species}" | tr ' /' '__')
    accession_safe=$(echo "${assembly_accession}" | tr ' /' '__')
    asm_dir="${OUT_ROOT}/${species_safe}/${accession_safe}"

    mkdir -p "${asm_dir}"

    echo
    echo "------------------------------------------------------------"
    echo "Species    : ${target_species}"
    echo "Accession  : ${assembly_accession}"
    echo "Source     : ${source_database}"
    echo "Directory  : ${asm_dir}"
    echo "------------------------------------------------------------"

    for url in \
        "${genome_fna_gz}" \
        "${protein_faa_gz}" \
        "${cds_fna_gz}" \
        "${gff_gz}" \
        "${gbff_gz}" \
        "${assembly_report_txt}" \
        "${md5checksums_txt}"
    do
        filename=$(basename "${url}")
        outfile="${asm_dir}/${filename}"

        echo "Downloading: ${filename}"

        if ! wget -q -c --tries=3 --timeout=60 -O "${outfile}" "${url}"; then
            echo "WARNING: failed to download ${url}" | tee -a "${LOG_DIR}/02d_download_warnings.log"
            rm -f "${outfile}"
        fi
    done
done

echo
echo "============================================================"
echo "Direct download completed"
echo "============================================================"
echo "Files root: ${OUT_ROOT}"
