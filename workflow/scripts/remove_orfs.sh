#!/bin/bash

REMOV="${1}"
DIR=$(dirname "${REMOV}")
#SEQS=$(perl -ne 'print $_') #tr -d '\0' |



cat /dev/stdin | awk '/^>/{split($1,a,"- OS"); print a[1]; next}1' | seqkit replace -s -p "\*" -r "" >"${DIR}"/"${MMSEQS_ENTRY_NAME}"_fa

# get the eventual sequences to be removed from the cluster
awk -v N=$MMSEQS_ENTRY_NAME '$2==N{print $1}' "${REMOV}" > "${DIR}"/"${MMSEQS_ENTRY_NAME}"_rem.txt

if [[ -s "${DIR}"/"${MMSEQS_ENTRY_NAME}"_rem.txt ]]; then
  seqkit fx2tab "${DIR}"/"${MMSEQS_ENTRY_NAME}"_fa | LC_ALL=C grep -v -f "${DIR}"/"${MMSEQS_ENTRY_NAME}"_rem.txt | seqkit tab2fx
  rm "${DIR}"/"${MMSEQS_ENTRY_NAME}"_rem.txt
  rm "${DIR}"/"${MMSEQS_ENTRY_NAME}"_fa
else
  awk '{print $0}' "${DIR}"/"${MMSEQS_ENTRY_NAME}"_fa
  rm "${DIR}"/"${MMSEQS_ENTRY_NAME}"_rem.txt
  rm "${DIR}"/"${MMSEQS_ENTRY_NAME}"_fa
fi
