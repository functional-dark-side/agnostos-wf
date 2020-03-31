#!/bin/bash

set -x

#export PATH=$PATH:$HHLIB/bin:$HHLIB/scripts

#module load gcc

DB=${1}
DIR=${PWD}/${2}
HHBLITS=${3}
THREADS=${4}

perl -ne 'print $_' | "${HHBLITS}" -i stdin -n 2 -v 0 -d "${DB}" -cpu "${THREADS}" -Z 10000000 -B 10000000 -e 1 -o "${DIR}"/"${MMSEQS_ENTRY_NAME}".hhr

QLEN=$(grep 'Match_' "${DIR}"/"${MMSEQS_ENTRY_NAME}".hhr | awk '{print $2}')

sed -n '/No Hit/,/No 1/p' "${DIR}"/"${MMSEQS_ENTRY_NAME}".hhr |
  grep -v 'No' |
  sed -e 's/^\s*//' | sed 's/[()]/ /g' | tr -s ' ' |
  perl -e 'while(<>){@matched = $_ =~ m/^(\d+)\s([\S\s]{30})\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/; print "$_\t" foreach @matched; print "\n"}' |
  sed '/^\s*$/d' |
  awk -v ql="${QLEN}" -v q="${MMSEQS_ENTRY_NAME}" \
    'BEGIN{FS="\t"}{split($2, c, " "); split($9, a, "-"); split($10, b, "-"); print q"\t"c[1]"\t"$3"\t"$4"\t"$6"\t"$8"\t"a[1]"\t"a[2]"\t"b[1]"\t"b[2]"\t"ql"\t"$11"\t"(a[2]-a[1]+1)/ql"\t"(b[2]-b[1]+1)/$11}'
