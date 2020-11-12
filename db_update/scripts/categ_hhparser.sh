#!/bin/bash

set -x
set -e


FILE=$(perl -ne 'print $_')

QLEN=$(grep 'Match_' <(echo "${FILE}") | awk '{print $2}')

sed -n '/No Hit/,/No 1/p' <(echo "${FILE}") \
  | grep -v 'No' \
  | sed '/^\s*$/d' | sed 's/[()]/ /g' | tr -s ' ' \
  | awk -v ql="${QLEN}" -v q="${MMSEQS_ENTRY_NAME}" \
  '{split($9,a,"-");split($10,b,"-"); print q"\t"$2"\t"$3"\t"$4"\t"$6"\t"$8"\t"a[1]"\t"a[2]"\t"b[1]"\t"b[2]"\t"ql"\t"$11"\t"(a[2]-a[1]+1)/ql"\t"(b[2]-b[1]+1)/$11}'
