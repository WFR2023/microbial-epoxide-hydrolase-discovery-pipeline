#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-seh_struct}"
BASE_DIR="."
OUT_DIR="${BASE_DIR}/19_structural_pipeline_setup"
mkdir -p "${OUT_DIR}"

REPORT="${OUT_DIR}/install_report.txt"
VERSIONS="${OUT_DIR}/installed_tool_versions.tsv"

echo "Creating/using conda environment: ${ENV_NAME}" | tee "${REPORT}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found in PATH." | tee -a "${REPORT}"
  exit 1
fi

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "Conda environment already exists: ${ENV_NAME}" | tee -a "${REPORT}"
else
  conda create -n "${ENV_NAME}" -y python=3.11 | tee -a "${REPORT}"
fi

conda activate "${ENV_NAME}"

conda install -y -c conda-forge -c bioconda \
  pandas biopython matplotlib openpyxl requests \
  openbabel vina fpocket pymol-open-source \
  > "${OUT_DIR}/conda_install_stdout.log" 2> "${OUT_DIR}/conda_install_stderr.log" || {
    echo "ERROR: package installation failed. See logs in ${OUT_DIR}" | tee -a "${REPORT}"
    exit 1
  }

{
  echo -e "tool\tversion_or_status"
  for cmd in python obabel vina fpocket pymol; do
    if command -v "${cmd}" >/dev/null 2>&1; then
      ver="$("${cmd}" -h 2>&1 | head -n 1 | tr '\t' ' ' | sed 's/[[:space:]]\+/ /g')"
      echo -e "${cmd}\t${ver}"
    else
      echo -e "${cmd}\tNOT_FOUND"
    fi
  done
  python - <<'PY'
import sys
import pandas
import Bio
import matplotlib
import openpyxl
print(f"python_pkg\t{sys.version.split()[0]}")
print(f"pandas_pkg\t{pandas.__version__}")
print(f"biopython_pkg\t{Bio.__version__}")
print(f"matplotlib_pkg\t{matplotlib.__version__}")
print(f"openpyxl_pkg\t{openpyxl.__version__}")
PY
} > "${VERSIONS}"

echo "" | tee -a "${REPORT}"
echo "Installation finished." | tee -a "${REPORT}"
echo "Environment: ${ENV_NAME}" | tee -a "${REPORT}"
echo "Versions table: ${VERSIONS}" | tee -a "${REPORT}"
echo "Activate later with:" | tee -a "${REPORT}"
echo "  conda activate ${ENV_NAME}" | tee -a "${REPORT}"
