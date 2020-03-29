#!/bin/bash

seqDB="${1}"
CLSTR="${2}"
CLSTR_INFO="${3}"
MMSEQS="${4}"
NSLOTS="${5}"

set -x
set -e

# Cluster naming
# The official cluster names are going to be based on the line number of the wide formatted file
# We are also going to produce a correspondence file to access the clusters in the MMseqs2 indices

# 1. Get representative MMseqs2 DB (same index as cluDB)
"${MMSEQS}" result2repseq "${seqDB}" "${CLSTR}" "${CLSTR}"_rep --threads "${NSLOTS}"
"${MMSEQS}" result2flat "${seqDB}" "${seqDB}" "${CLSTR}"_rep "${CLSTR}"_rep.fasta --use-fasta-header
grep '^>' "${CLSTR}"_rep.fasta | sed 's/^>//' | cut -d' ' -f1 >"${CLSTR}"_rep.tsv

# 2. Combine MMseqs2 index with representative header
paste <(awk '{print $1}' "${CLSTR}"_rep.index) \
  "${CLSTR}"_rep.tsv > mmseqs_reps.tsv

# 3. Combine our cluster names with the original MMSEQS2 names
join -12 -22 <(sort -k2,2 --parallel="${NSLOTS}" -S20% mmseqs_reps.tsv) \
  <(sort -k2,2 --parallel="${NSLOTS}" -S20% <(awk '{print $1"\t"$2}' "${CLSTR_INFO}")) |
sort -k2,2n --parallel="${NSLOTS}" -S20% > mmseqs_reps_clname.tsv

# 4. Create index file with <MMseqs-index-num> <Cluster-name>
awk '{print $2"\t"$3}' mmseqs_reps_clname.tsv | sort -k2,2 --parallel="${NSLOTS}" -S20% > "${CLSTR}"_name_index.txt

# Clean-up
rm mmseqs_reps.tsv mmseqs_reps_clname.tsv "${CLSTR}"_rep "${CLSTR}"_rep.index
