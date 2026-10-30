#!/bin/bash
set -euo pipefail


if [ $# -lt 1 ]; then
  echo "Usage: $0 <SESSION_DIR>"
  exit 1
fi

SESSION_DIR="$1"   # es: /media/.../voxelwise_hurst/FIN/session_0
OUTDIR="/path/to/external/data/results/residuals/one_samp_per_session"
MASK="group_param_increase_mask.nii.gz"
N_PERM=5000
TFCE_FLAG="-T"

mkdir -p "${OUTDIR}"

SESSION_NAME="$(basename "$(dirname "$SESSION_DIR")")_$(basename "$SESSION_DIR")"

echo ">>> Session directory: ${SESSION_DIR}"
echo ">>> Creating 4D file ..."

FOURD="${OUTDIR}/${SESSION_NAME}_4D.nii.gz"
echo ">>> 4D output: ${FOURD}"

fslmerge -t "${FOURD}" "${SESSION_DIR}"/*.nii.gz

echo ">>> Running randomise ..."
randomise -i "${FOURD}" \
          -o "${OUTDIR}/masked/OneSampT_${SESSION_NAME}" \
          -1 \
          -m "${MASK}" \
          -n ${N_PERM} \
          ${TFCE_FLAG}

echo ">>> Done."
