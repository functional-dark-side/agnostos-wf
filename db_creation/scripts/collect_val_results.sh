#!/bin/bash


DIR="${1}"
STATS="${2}"
DIR2=$(dirname "${DIR}")
DIR_REJ="${DIR2}"/rejected
REJ="${3}"
PARALLEL="${4}"
NSLOTS="${5}"

# Collect compositional validation stats
find "${DIR}" -name '*_SSN_filt_stats.tsv' | "${PARALLEL}" -j "${NSLOTS}" cat {} > "${STATS}"

# Collect all ORFs identified as bad-aligned, "rejected"

find "${DIR_REJ}" -name '*_rejected.txt' | "${PARALLEL}" -j "${NSLOTS}" awk \'{print FILENAME,\$0}\' {} > "${REJ}"

sed -i -e 's|.*rejected/\(.*\)_rejected.txt|\1|' "${REJ}"

sed -i 's/ /\t/g' "${REJ}"
