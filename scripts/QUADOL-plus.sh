#!/bin/bash
# Usage: bash QUADOL-plus.sh

wdir=$(cd "$(dirname "$0")" && pwd)

SCRIPT_DIR=${wdir}

ALSRAC_SCRIPT="${SCRIPT_DIR}/ALSRAC.sh"
QUADOL_SCRIPT="${SCRIPT_DIR}/QUADOL.sh"

MODE="mred"

ERROR_TH=0.001

# Run ALSRAC
bash "${ALSRAC_SCRIPT}" "${MODE}" "${ERROR_TH}"

# Run QUADOL
bash "${QUADOL_SCRIPT}" "${MODE}" "${ERROR_TH}"

wait
