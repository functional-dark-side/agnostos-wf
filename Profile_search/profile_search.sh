#!/bin/bash

set -e
set -x

usage() { echo -e "Usage: $(basename $0) [OPTIONS...]
Options:
--query            genes fasta sequences \
--clu_hmm          cluster HMMs db \
--clu_cat        cluster categories table \
--info             genes additional information \
--mmseqs           location of mmseqs binary \
--mpi              logic, default is FALSE \
--threads          number of threads to use" 1>&2; exit 1; }

OPT_LIST="query:,clu_hmm:,clu_cat:,info:,mmseqs:,mpi:,threads:"

eval set -- $(getopt -o '' --long "${OPT_LIST}" -- "$@")
while true; do
  case "$1" in
    --query)
      QUERY=$2; shift 2 ;;
    --clu_hmm)
      CLHMM=$2; shift 2 ;;
    --clu_cat)
      CLCAT=$2; shift 2 ;;
    --info)
      INFO=$2; shift 2 ;;
    --mmseqs)
      MMSEQS_BIN=$2; shift 2 ;;
    --mpi)
      MPI=$2; shift 2 ;;
    --threads)
      NSLOTS=$2; shift 2 ;;
    --)
      shift; break ;;
    *)
      usage ;;
  esac
done

if [ -z "${QUERY}" ] || [ -z "${CLHMM}" ] || [ -z "${CLCAT}" ]; then
  usage
fi

MMSEQS_BIN=${MMSEQS_BIN:=mmseqs}

MPI=${MPI:=FALSE}

INFO=${INFO:="none"}
# Fixed variables
MPI="srun --mpi=pmi2"
EFILTER="${PWD}"/Profile_search/evalue_filter.awk
MVOTE="${PWD}"/Profile_search/majority_vote_categ.R
LTMP=/vol/scratch/tmp
OUTDIR="${PWD}"/profile_search_res

mkdir -p "${OUTDIR}"

# Create query DB
# Create sequence databases
NAME=$(basename "${QUERY}" .fasta)
QUERY_DB="${OUTDIR}"/"${NAME}"_db
"${MMSEQS_BIN}" createdb "${QUERY}" "${QUERY_DB}"

# Decide to use mpi search based on the number of query sequences:
NSEQS=$(grep -c '^>' "${QUERY}")

if [[ ${NSEQS} -ge 100000 ]]; then
  if [[ ${MPI}=="mpi" ]]; then
    # Sequernce-profile search against the Cluster HMM profiles
    sbatch --ntasks-per-node 1 --wait \
      --nodes 9 --cpus-per-task 28 \
      --job-name profs --partition nomaster \
      --wrap " ${MMSEQS_BIN} search ${QUERY_DB} ${CLHMM} ${OUTDIR}/${NAME}_vs_mg_gtdb_hmm_db ${OUTDIR}/tmp \
      --local-tmp ${LTMP} \
      --threads ${NSLOTS} -e 1e-20 --cov-mode 2 -c 0.6 \
      --split-mode 0 --split-memory-limit 150G \
      --mpi-runner \"${MPI}\" "
  else
    "${MMSEQS_BIN}" search "${QUERY_DB}" "${CLHMM}" "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_db "${OUTDIR}"/tmp \
      --threads "${NSLOTS}" -e 1e-20 --cov-mode 2 -c 0.6
  fi
else
  # Sequernce-profile search against the Cluster HMM profiles
  "${MMSEQS_BIN}" search "${QUERY_DB}" "${CLHMM}" "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_db "${OUTDIR}"/tmp \
   --threads "${NSLOTS}" -e 1e-20 --cov-mode 2 -c 0.6
fi

wait

# Convert result database in tsv-file
"${MMSEQS_BIN}" convertalis "${QUERY_DB}" "${CLHMM}" "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_db \
  "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_hits.tsv \
  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov' \
  --threads "${NSLOTS}"
# Filter hits within 90% of the Log(best(e-value))
awk -v P=0.9 -f "${EFILTER}" "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_hits.tsv > "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_res_qcov06_e90.tsv
# Add functional category information
join -12 -21 <(sort -k2,2 "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_res_qcov06_e90.tsv ) \
  <(sort -k1,1   "${CLCAT}") \
  > "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_search_res.tsv

gzip "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_search_res.tsv

# Parse results to get best-hits and category proportions per sample (or genome or contig..)
# NB: the info file should have the following format: <gene> <sample> (or <genome>, or <contig>)

"${MVOTE}"  --res "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_search_res.tsv.gz \
            --info "${INFO}"

#Rscript --vanilla "${MVOTE}" "${OUTDIR}"/"${NAME}"_vs_mg_gtdb_hmm_search_res.tsv.gz "${INFO}"
