#!/bin/bash
# Usage: bash QUADOL.sh [Error Metric] [Error Threshold]

wdir=$(cd "$(dirname "$0")" && pwd)
source ${wdir}/env.sh ${1} ${2}

echo "TestCase,BestALSRAC,BestQUADOL"

# LUT Sharing Flow
for EXACT_BLIF in "${BENCHMARK_DIR}"/*
do
  BLIF_FILENAME=$(basename "${EXACT_BLIF}")
  PREFIX=${BLIF_FILENAME%%_*}

  INPUTS_DIR="${ALSRAC_RST_DIR}/${PREFIX}"
  OUTPUT_DIR="${QUADOL_RST_DIR}/${PREFIX}"

  RESULTS_CSV="${RESULTS_DIR}/${PREFIX}.csv"

  mkdir -p ${OUTPUT_DIR}

  start=$(date +%s)
  python ${SHAREDLUT_BIN} -a ${EXACT_BLIF} -o ${OUTPUT_DIR} -s ${SIMULATOR_BIN} -e ${E_TH} -m ${MODE} >> ${RESULTS_CSV}
  wait
  for BLIF in "${INPUTS_DIR}"/*.blif
  do
    python ${SHAREDLUT_BIN} -i ${EXACT_BLIF} -a ${BLIF} -o ${OUTPUT_DIR} -s ${SIMULATOR_BIN} -e ${E_TH} -m ${MODE} >> ${RESULTS_CSV}
  done

  wait

  end=$(date +%s)
  runtime=$((end - start))
  # Record the runtime
  echo "${PREFIX},${runtime}" >> ${QUADOL_RST_DIR}/runtime.log

  # Find the smallest number in the first column, including the header
  ALSRAC_MIN=$(cut -d',' -f1 "${RESULTS_CSV}" | sort -n | head -n1)

  # Find the smallest number in the third column, including the header
  QUADOL_MIN=$(cut -d',' -f3 "${RESULTS_CSV}" | sort -n | head -n1)

  # Output the results
  echo "${PREFIX},${ALSRAC_MIN},${QUADOL_MIN}"

done

echo "QUADOL are done."
