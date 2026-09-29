#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-colabfold}"
BASE_DIR="<PROJECT_ROOT>"
OUT_DIR="${BASE_DIR}/19_colabfold_install"
mkdir -p "${OUT_DIR}"

REPORT="${OUT_DIR}/install_report.txt"
VERSIONS="${OUT_DIR}/tool_versions.tsv"

echo "Starting ColabFold installation into env: ${ENV_NAME}" | tee "${REPORT}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found in PATH." | tee -a "${REPORT}"
  exit 1
fi

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "Conda environment already exists: ${ENV_NAME}" | tee -a "${REPORT}"
else
  conda create -n "${ENV_NAME}" -y -c conda-forge -c bioconda \
    python=3.13 kalign2=2.04 hhsuite=3.3.0 mmseqs2=18.8cc5c \
    > "${OUT_DIR}/conda_create_stdout.log" 2> "${OUT_DIR}/conda_create_stderr.log"
fi

conda activate "${ENV_NAME}"

# CPU-only official install path from ColabFold README
python -m pip install --upgrade pip \
  > "${OUT_DIR}/pip_upgrade_stdout.log" 2> "${OUT_DIR}/pip_upgrade_stderr.log"

python -m pip install "colabfold[alphafold,openmm]" \
  > "${OUT_DIR}/pip_install_stdout.log" 2> "${OUT_DIR}/pip_install_stderr.log"

{
  echo -e "tool\tstatus_or_version"
  for cmd in python pip colabfold_batch mmseqs hhsearch kalign; do
    if command -v "${cmd}" >/dev/null 2>&1; then
      ver="$("${cmd}" -h 2>&1 | head -n 1 | tr '\t' ' ' | sed 's/[[:space:]]\+/ /g')"
      echo -e "${cmd}\t${ver}"
    else
      echo -e "${cmd}\tNOT_FOUND"
    fi
  done
} > "${VERSIONS}"

echo "" | tee -a "${REPORT}"
echo "Installation complete." | tee -a "${REPORT}"
echo "Environment: ${ENV_NAME}" | tee -a "${REPORT}"
echo "Version table: ${VERSIONS}" | tee -a "${REPORT}"
echo "" | tee -a "${REPORT}"
echo "Next commands:" | tee -a "${REPORT}"
echo "  source \"\$(conda info --base)/etc/profile.d/conda.sh\"" | tee -a "${REPORT}"
echo "  conda activate ${ENV_NAME}" | tee -a "${REPORT}"
echo "  colabfold_batch --help" | tee -a "${REPORT}"
