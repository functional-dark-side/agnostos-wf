#!/bin/bash

set -x
set -e

HHBLITS=${1}
DB=${2}
DIR=${3}
mkdir -p "${DIR}"

perl -ne 'print $_' | "${HHBLITS}" -i stdin -n 2 -v 0 -d "${DB}" -cpu 2 -Z 10000000 -B 10000000 -e 1 -o "${DIR}"/"${MMSEQS_ENTRY_NAME}".hhr

QLEN=$(grep 'Match_' "${DIR}"/${MMSEQS_ENTRY_NAME}.hhr | awk '{print $2}')

sed -n '/No Hit/,/No 1/p' "${DIR}"/"${MMSEQS_ENTRY_NAME}".hhr \
  | grep -v 'No' \
  | sed '/^\s*$/d' | sed 's/[()]/ /g' | tr -s ' ' \
  | awk -v ql="${QLEN}" -v q="${MMSEQS_ENTRY_NAME}" \
  '{split($9,a,"-");split($10,b,"-"); print q"\t"$2"\t"$3"\t"$4"\t"$6"\t"$8"\t"a[1]"\t"a[2]"\t"b[1]"\t"b[2]"\t"ql"\t"$11"\t"(a[2]-a[1]+1)/ql"\t"(b[2]-b[1]+1)/$11}'
