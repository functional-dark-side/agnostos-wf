#!/bin/bash

set -e
set -x

usage() { echo -e "Usage: $(basename $0) [OPTIONS...]
Options:
--query            genes fasta sequences \
--info             genes additional information \
--mmseqs           location of mmseqs binary \
  --threads          number of threads to use" 1>&2; exit 1; }

OPT_LIST="query:,info:,mmseqs:,threads:"

eval set -- $(getopt -o '' --long "${OPT_LIST}" -- "$@")
while true; do
  case "$1" in
    --query)
      QUERY=$2; shift 2 ;;
    --info)
      INFO=$2; shift 2 ;;
    --mmseqs)
      MMSEQS_BIN=$2; shift 2 ;;
    --threads)
      NSLOTS=$2; shift 2 ;;
    --)
      shift; break ;;
    *)
      usage ;;
  esac
done

if [ -z "${QUERY}" ] ; then
  usage
fi

MMSEQS_BIN=${MMSEQS_BIN:=~/opt/MMseqs2/bin/mmseqs}

# Fixed variables
CLHMM="${PWD}"/mg_gtdb_hmm_db #/bioinf/projects/megx/UNKNOWNS/2017_11/mg_gtdb_hmm_db
CL_CATEG="${PWD}"/mg_gtdb_kept_cluster_categories_class.tsv.gz
MPI="srun --mpi=pmi2"
EFILTER="${PWD}"/scripts/evalue_filter.awk
MVOTE="${PWD}"/scripts/majority_vote_categ.R
LTMP=/vol/scratch/tmp
OUTDIR="${PWD}"/results

mkdir -p "${OUTDIR}"

# Create query DB
# Create sequence databases
NAME=$(basename "${QUERY}" .fasta)
QUERY_DB="${OUTDIR}"/"${NAME}"_db
"${MMSEQS_BIN}" createdb "${QUERY}" "${QUERY_DB}"

# Decide to use mpi search based on the number of query sequences:
NSEQS=$(grep -c '^>' "${QUERY}")

if [[ ${NSEQS} -ge 1000000 ]]; then
  # Sequernce-profile search against the Cluster HMM profiles
  "${MMSEQS_BIN}" search "${QUERY_DB}" "${CLHMM}" "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_db "${OUTDIR}"/tmp \
    --local-tmp "${LTMP}" \
    --threads "${NSLOTS}" -e 1e-20 --cov-mode 2 -c 0.6 \
    --mpi-runner "${MPI}"
else
  # Sequernce-profile search against the Cluster HMM profiles
  "${MMSEQS_BIN}" search "${QUERY_DB}" "${CLHMM}" "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_db "${OUTDIR}"/tmp\
    --threads "${NSLOTS}" -e 1e-20 --cov-mode 2 -c 0.6
fi


# Convert result database in tsv-file
"${MMSEQS_BIN}" convertalis "${QUERY_DB}" "${CLHMM}" "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_db \
  "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_hits.tsv \
  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov' \
  --threads "${NSLOTS}"
# Filter hits within 90% of the Log(best(e-value))
awk -v P=0.9 -f "${EFILTER}" "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_hits.tsv > "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_res_qcov06_e90.tsv
# Add functional category information
join -12 -21 <(sort -k2,2 "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_res_qcov06_e90.tsv ) \
  <(sort -k1,1  <(zcat "${CL_CATEG}")) \
  > "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_search_res.tsv

gzip "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_search_res.tsv

# Parse results to get best-hits and category proportions per sample (or genome or contig..)
# NB: the info file should have the following format: <gene> <sample> (or <genome>, or <contig>)
Rscript --vanilla "${MVOTE}" "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_search_res.tsv.gz "${INFO}"

rm "${MAG}"_genes* "${MAG}"_mg_gtdb_db* "${MAG}"_mg_gtdb_qcov06_e90.tsv
