#!/bin/bash

set -e
set -x

usage() {
  echo -e "Usage: $(basename "$0") [OPTIONS...]
Options:
--search             program to perform the search (MMSEQS) \
--mpi_runner         mpi runner \
--ltmp               local template folder \
--cons               query consensus sequences fasta format \
--db_fasta           database sequences fasta format \
--db_info            database sequence additional information \
--evalue_filter      script to filter hitsbased on evalue \
--evalue_threshold   evalue threshold \
--hypo_threshold     threshold for the number of hypothetical proteins per cluster \
--hypo_patterns      patterns to retrieve hypothetical proteins \
--grep               parallel grep program \
--output             path to the output file \
--threads            number of threads to use" 1>&2
  exit 1
}

OPT_LIST="search:,mpi_runner:,ltmp:,cons:,db_fasta:,db_info:,evalue_filter:,evalue_threshold:,hypo_threshold:,hypo_patterns:,grep:,output:,threads:,outdir:"

eval set -- "$(getopt -o '' --long "${OPT_LIST}" -- "$@")"
while true; do
  case "$1" in
  --search)
    MMSEQS_BIN=$2
    shift 2
    ;;
  --mpi_runner)
    MPI=$2
    shift 2
    ;;
  --ltmp)
    LTMP=$2
    shift 2
    ;;
  --cons)
    CONS=$2
    shift 2
    ;;
  --db_fasta)
    DB_FA=$2
    shift 2
    ;;
  --db_info)
    DB_PROT=$2
    shift 2
    ;;
  --evalue_filter)
    EFILTER=$2
    shift 2
    ;;
  --evalue_threshold)
    E=$2
    shift 2
    ;;
  --hypo_threshold)
    HYPO=$2
    shift 2
    ;;
  --hypo_patterns)
    PATTERNS=$2
    shift 2
    ;;
  --grep)
    PGREP=$2
    shift 2
    ;;
  --output)
    RES=$2
    shift 2
    ;;
  --outdir)
    DIR=$2
    shift 2
    ;;
  --threads)
    NSLOTS=$2
    shift 2
    ;;
  --)
    shift
    break
    ;;
  *)
    usage
    ;;
  esac
done

if [ -z "${MMSEQS_BIN}" ] || [ -z "${MPI}" ] || [ -z "${LTMP}" ] || [ -z "${CONS}" ] ||
  [ -z "${DB_FA}" ] || [ -z "${DB_PROT}" ] || [ -z "${EFILTER}" ] || [ -z "${DIR}" ] ||
  [ -z "${PATTERNS}" ] || [ -z "${PGREP}" ] || [ -z "${RES}" ]; then
  usage
fi

queryDB="${DIR}"/cons_db
targetDB="${DIR}"/$(md5sum <(echo "${DB_FA}") | awk '{print $1}')

# Create query and target DBs
"${MMSEQS_BIN}" createdb "${CONS}" "${queryDB}"

if [ ! -f "${targetDB}" ]; then
  "${MMSEQS_BIN}" createdb "${DB_FA}" "${targetDB}"
fi


res1="${DIR}"/round1_res


top1="${DIR}"/top1

aligned="${DIR}"/aligned

aligned_db="${DIR}"/aligned_db

res2="${DIR}"/round2_res

merged="${DIR}"/merged

aln_2b="${DIR}"/aln_2b
# The search is performed following a similar strategy to the 2bLCA protocol (dual BLAST based last common ancestor)
# 1. search query with e-vlaue < 1e-5
# 2. search again with aligned regions of best hits and same or lower e-evalue
# 3. merge results
# First search [round 1]
mkdir -p "${DIR}"/tmp_hsp1

"${MMSEQS_BIN}" search "${queryDB}" "${targetDB}" "${res1}" "${DIR}"/tmp_hsp1 \
  --local-tmp "${LTMP}" \
  --max-seqs 300 --threads "${NSLOTS}" -a -e 1e-5 \
  --cov-mode 2 -c 0.6 \
  --mpi-runner "${MPI}"

# Filter search results:
## 2bLCA mode
"${MMSEQS_BIN}" filterdb "${res1}" "${top1}" --extract-lines 1
## Top-hit mode
"${MMSEQS_BIN}" filterdb "${res1}" "${top1}" --beats-first \
  --filter-column 4 --comparison-operator le

"${MMSEQS_BIN}" createsubdb "${top1}".index "${queryDB}"_h "${top1}"_h
# Extract aligned regions from target DB
"${MMSEQS_BIN}" extractalignedregion "${queryDB}" "${targetDB}" \
  "${top1}" "${aligned}" --extract-mode 2 --threads "${NSLOTS}"

# Second search [round 2]
mkdir -p "${DIR}"/tmp_hsp2

"${MMSEQS_BIN}" createsubdb "${aligned}" "${queryDB}" "${aligned_db}"

#lower e-value?
"${MMSEQS_BIN}" search "${aligned_db}" "${targetDB}" "${res2}" "${DIR}"/tmp_hsp2 \
  --local-tmp "${LTMP}" \
  --max-seqs 300 --threads "${NSLOTS}" -a -e 1e-5 \
  --cov-mode 2 -c 0.6 \
  --mpi-runner "${MPI}"

# Concat top hit from first search with all the results from second search
# We can use filterdb --beats-first to filter out all entries from second search that do not reach the evalue of the top1 hit
"${MMSEQS_BIN}" mergedbs "${top1}" "${merged}" "${top1}" "${res2}"

"${MMSEQS_BIN}" filterdb "${merged}" "${aln_2b}" --beats-first --filter-column 4 --comparison-operator le

#The results can be converted into flat format both with createtsv or convertalis
"${MMSEQS_BIN}" convertalis "${queryDB}" "${targetDB}" "${aln_2b}" "${RES}" \
  --threads "${NSLOTS}" \
  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov'

# Filter e-value inside the 60% of the best evalue
OUT=$(basename "${RES}" '.tsv')
F="${DIR}"/"${OUT}"_E"${E}"
awk -v P="${E}" -f "${EFILTER}" <(sort -k1,1 -k11,11g "${RES}" | awk '!seen[$1,$2]++') >"${F}"

# Retrive database protein information
FP="${F}"_prot.tsv
join -12 -21 <(sort -k2,2 "${F}") <(zcat "${DB_PROT}" | sort -k1,1 --parallel="${NSLOTS}" -S25% -T "${LTMP}") > "${FP}"

FOUT="${F}"_hypo_char
LC_ALL=C "${PGREP}" -j 4 -i -f "${PATTERNS}" "${FP}" | awk '{print $0"\thypo"}' >"${FOUT}"
LC_ALL=C "${PGREP}" -j 4 -v -i -f "${PATTERNS}" "${FP}" | awk '{print $0"\tchar"}' >>"${FOUT}"
sed -i 's/ /\t/g' "${FOUT}"
RES_HYPO="${DIR}"/"${OUT}"_hypo.txt
awk -v P="${HYPO}" 'BEGIN{FS="\t"}{a[$2][$16]++}END{for (i in a) {N=a[i]["hypo"]/(a[i]["hypo"]+a[i]["char"]); if (N >= P){print i}}}' "${FOUT}" >"${RES_HYPO}"

awk '{print $2"\t"$1"\t"$11"\t"$15}' "${FOUT}" >"${FP}"

rm -rf "${F}" "${FOUT}" "${queryDB}"* "${targetDB}"* "${merged}"* "${aln_2b}"* "${top1}"* "${res1}"* "${res2}"* "${aligned}"* "${aligned_db}"*
rm -rf "${DIR}"/tmp_hsp1 "${DIR}"/tmp_hsp2
