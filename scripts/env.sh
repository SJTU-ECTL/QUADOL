#!/bin/bash
# Usage: bash env.sh [Error Metric] [Error Threshold]

wdir=$(cd "$(dirname "$0")" && pwd)

ENV_DIR="${wdir}/.."
SRC_DIR="${ENV_DIR}/src"
UTILS_DIR="${ENV_DIR}/utils"
INPUT_DIR="${ENV_DIR}/inputs"
RESULTS_DIR="${ENV_DIR}/results"

# Location for ALSRAC project
ALSRAC_DIR="${ENV_DIR}/../ALSRAC"

MCNC_LIB="${INPUT_DIR}/mcnc.genlib"

ALSRAC_BIN="${ALSRAC_DIR}/main"
SHAREDLUT_BIN="${SRC_DIR}/main.py"
AF_ALSRAC_BIN="${SRC_DIR}/dataProcess.py"
SIMULATOR_BIN="${UTILS_DIR}/simulator.out"
APPDECOMP_BIN="${UTILS_DIR}/AppDecompose.out"

BENCHMARK_DIR="${INPUT_DIR}/Arith"

ALSRAC_RST_DIR="${RESULTS_DIR}/ALSRAC"
QUADOL_RST_DIR="${RESULTS_DIR}/QUADOL"

mkdir -p ${RESULTS_DIR}

mkdir -p ${ALSRAC_RST_DIR}
mkdir -p ${QUADOL_RST_DIR}

# Check if the input argument is provided
if [ -z "${1}" ]; then
  echo "ERROR: Error Metric not provided."
  exit
else
  MODE="${1}"
fi

# Check if the error threshold is provided
if [ -z "${2}" ]; then
  echo "ERROR: Error Threshold not provided."
  exit
else
  E_TH="${2}"
fi
