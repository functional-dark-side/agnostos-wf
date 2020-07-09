#!/bin/bash

set -x

DB=${1}
DIR=${2}
HHBLITS=${3}
THREADS=${4}
PROB=${5}

HHR="${DIR}"/"${MMSEQS_ENTRY_NAME}".hhr

perl -ne 'print $_' | "${HHBLITS}" -i stdin -n 2 -v 0 -d "${DB}" -cpu "${THREADS}" -Z 10000000 -B 10000000 -e 1 -o "${HHR}"

QLEN=$(grep 'Match_' "${HHR}" | awk '{print $2}')

sed -n '/No Hit/,/No 1/p' "${HHR}" | \
  grep -v 'No' | \
  sed -e 's/^\s*//' | sed 's/[()]/ /g' | tr -s ' ' | \
  perl -e 'while(<>){@matched = $_ =~ m/^(\d+)\s([\S\s]{30})\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/; print "$_\t" foreach @matched; print "\n"}' | \
  sed '/^\s*$/d' | \
  awk -vql="${QLEN}" -vq="${MMSEQS_ENTRY_NAME}" \
  'BEGIN{FS="\t"}{split($2, c, " "); split($9, a, "-"); split($10, b, "-"); print q"\t"c[1]"\t"$3"\t"$4"\t"$6"\t"$8"\t"a[1]"\t"a[2]"\t"b[1]"\t"b[2]"\t"ql"\t"$11"\t"(a[2]-a[1]+1)/ql"\t"(b[2]-b[1]+1)/$11}' | \
  awk -vP="${PROB}" '$3 >= P' >"${HHR}".tmp1

grep '^>' "${HHR}" | perl -e 'while(<>){
                                    if (@matched = $_ =~ m/^>(\S+) (\S.+) OS=\S.*/){
                                            print "$_\t" foreach @matched; print "\n";
                                        }elsif(@matched = $_ =~ m/^>(\S.+) ; (\S.+) ; (.*)/) {
                                            print "$_\t" foreach @matched; print "\n";
                                        }else{
                                            @matched = $_ =~ m/^>(\.*)/;
                                            print "$_\t" foreach @matched; print "\n";
                                        }
}' | \
  sed '/^\s*$/d' >"${HHR}".tmp2

join -t $'\t' -1 2 -2 1 <(sort -k2,2 --parallel "${THREADS}" -S25% "${HHR}".tmp1) <(sort -k1,1 --parallel "${THREADS}" -S25% "${HHR}".tmp2)

#rm "${HHR}".tmp1 "${HHR}".tmp2
