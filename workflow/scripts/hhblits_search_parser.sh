#!/bin/bash

set -x

DIR=${1}
THREADS=${2}
PROB=${3}

mkdir -p "${DIR}"
FILE=$(perl -ne 'print $_')
QLEN=$(grep 'Match_' <(echo "${FILE}") | awk '{print $2}')
TMP="${DIR}"/"${MMSEQS_ENTRY_NAME}"


sed -n '/No Hit/,/No 1/p' <(echo "${FILE}") | \
  grep -v 'No' | \
  sed -e 's/^\s*//' | sed 's/[()]/ /g' | tr -s ' ' | \
  perl -e 'while(<>){@matched = $_ =~ m/^(\d+)\s([\S\s]{30})\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/; print "$_\t" foreach @matched; print "\n"}' | \
  sed '/^\s*$/d' | sed '/^$/d' | \
  awk -vql="${QLEN}" -vq="${MMSEQS_ENTRY_NAME}" \
  'BEGIN{FS="\t"}{split($2, c, " "); split($9, a, "-"); split($10, b, "-"); print q"\t"c[1]"\t"$3"\t"$4"\t"$6"\t"$8"\t"a[1]"\t"a[2]"\t"b[1]"\t"b[2]"\t"ql"\t"$11"\t"(a[2]-a[1]+1)/ql"\t"(b[2]-b[1]+1)/$11}' | \
  awk -vP="${PROB}" '$3 >= P' >"${TMP}".tmp1

grep '^>' <(echo "${FILE}") | perl -e 'while(<>){
                                    if (@matched = $_ =~ m/^>(\S+) (\S.+) OS=\S.*/){
                                            print "$_\t" foreach @matched; print "\n";
                                        }elsif(@matched = $_ =~ m/^>(\S.+) ; (\S.+) ; (.*)/) {
                                            print "$_\t" foreach @matched; print "\n";
                                        }else{
                                            @matched = $_ =~ m/^>(\.*)/;
                                            print "$_\t" foreach @matched; print "\n";
                                        }
}' | \
  sed '/^\s*$/d' >"${TMP}".tmp2

if [[ ! -s "${TMP}".tmp2 ]]; then
    grep '^>' <(echo "${FILE}") | sed 's/^>//' | \
        sed 's/ /:/' | sed 's/ /_/g' | sed 's/:/\t/g' > "${TMP}".tmp2
fi

if [[ -s "${TMP}".tmp1 ]]; then
    join -t $'\t' -1 2 -2 1 <(sort -k2,2 --parallel "${THREADS}" -S25% "${TMP}".tmp1) <(sort -k1,1 --parallel "${THREADS}" -S25% "${TMP}".tmp2)
    rm "${TMP}".tmp1 "${TMP}".tmp2
else
    rm "${TMP}".tmp1 "${TMP}".tmp2
fi
