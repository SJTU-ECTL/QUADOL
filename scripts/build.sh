#!/bin/bash
# Usage: bash build.sh

wdir=$(cd "$(dirname "$0")" && pwd)
WORKDING_DIR="${wdir}/.."

UTILS_DIR="${wdir}/utils"

UTILS_BUILD="${UTILS_DIR}/build"

mkdir -p ${UTILS_BUILD}
cmake -DCMAKE_BUILD_TYPE=Release -B "${UTILS_BUILD}" "${UTILS_DIR}"
make -C "${UTILS_BUILD}" -j$(nproc)
