#!/bin/bash

REMOV="${1}"
INDEX="${2}"
DIR=$(dirname "${REMOV}")
SEQS=$(perl -ne 'print $_')

echo "${SEQS}" | awk '/^>/{split($1,a,"- OS"); print a[1]; next}1' >"${DIR}"/"${MMSEQS_ENTRY_NAME}"_fa

# MMseqs2-index - Cluster-name correspondance
join -12 -22 <(sort -k2,2 "${REMOV}") <(sort -k2,2 "${INDEX}") >"${REMOV}".index

awk -v N=$MMSEQS_ENTRY_NAME '$3==N{print $2}' "${REMOV}".index >"${DIR}"/"${MMSEQS_ENTRY_NAME}"_rem.txt

if [[ -s "${DIR}"/"${MMSEQS_ENTRY_NAME}"_rem.txt ]]; then
  seqkit fx2tab "${DIR}"/"${MMSEQS_ENTRY_NAME}"_fa | LC_ALL=C grep -v -f "${DIR}"/"${MMSEQS_ENTRY_NAME}"_rem.txt | seqkit tab2fx
  rm "${DIR}"/"${MMSEQS_ENTRY_NAME}"_rem.txt
  rm "${DIR}"/"${MMSEQS_ENTRY_NAME}"_fa
else
  awk '{print $0}' "${DIR}"/"${MMSEQS_ENTRY_NAME}"_fa
  rm "${DIR}"/"${MMSEQS_ENTRY_NAME}"_rem.txt
  rm "${DIR}"/"${MMSEQS_ENTRY_NAME}"_fa
fi
