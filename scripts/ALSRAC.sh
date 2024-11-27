#!/bin/bash
# Usage: bash ALSRAC.sh [Error Metric] [Error Threshold]

wdir=$(cd "$(dirname "$0")" && pwd)
source ${wdir}/env.sh ${1} ${2}

NUM_TRIAL=10
MAP_TYPE=1

mkdir -p ${ALSRAC_RST_DIR}

# Run ALSRAC Flow for ${NUM_TRIAL} times
for BLIF in "${BENCHMARK_DIR}"/*
do
  BLIF_FILENAME=$(basename "${BLIF}")
  PREFIX=${BLIF_FILENAME%%_*}
  start=$(date +%s)
  for i in $(seq 1 ${NUM_TRIAL})
  do
    OUTPUT_DIR="${ALSRAC_RST_DIR}/${PREFIX}_${i}"
    mkdir -p ${OUTPUT_DIR}
    rm -f ${OUTPUT_DIR}/*

    ${ALSRAC_BIN} -i ${BLIF} -l ${MCNC_LIB} -m ${MODE} -o ${OUTPUT_DIR} -t ${MAP_TYPE} -f 64 -b ${E_TH} >> ${ALSRAC_RST_DIR}/${PREFIX}_${i}.log &
  done
  wait
  end=$(date +%s)
  runtime=$((end - start))
  echo "${PREFIX},${runtime}" >> ${ALSRAC_RST_DIR}/runtime.log
done

wait

python ${AF_ALSRAC_BIN} ${2}

echo "ALSRAC are done."
