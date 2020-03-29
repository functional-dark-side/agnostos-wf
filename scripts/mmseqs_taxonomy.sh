#!/bin/bash

set -e
set -x

usage() {
  echo -e "Usage: $(basename $0) [OPTIONS...]
Options:
--search             program to perform the search (MMSEQS)
--input              query sequences (refined db) \
--db_fasta           database sequences fasta format with unirpot ids \
--cl_info            file containing the corresp. orf-cl_name-category \
--output             path to the output file \
    --threads            number of threads to use" 1>&2
  exit 1
}

OPT_LIST="search:,input:,taxdb:,output:,threads:,cl_info:,mpi_runner:"

eval set -- $(getopt -o '' --long "${OPT_LIST}" -- "$@")
while true; do
  case "$1" in
    --search)
      MMSEQS_BIN=$2
      shift 2
      ;;
    --mpi_runner)
      RUNNER=$2
      shift 2
      ;;
    --input)
      QUERY=$2
      shift 2
      ;;
    --taxdb)
      TAXDB=$2
      shift 2
      ;;
    --cl_info)
      CL_INFO=$2
      shift 2
      ;;
    --output)
      RES=$2
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

if [ -z "${MMSEQS_BIN}" ] ||
[ -z "${QUERY}" ] ||
[ -z "${RUNNER}" ] ||
[ -z "${TAXDB}" ] ||
[ -z "${CL_INFO}" ] ||
[ -z "${RES}" ] ||
[ -z "${NSLOTS}" ]; then
  usage
fi

DIR=$(dirname "${RES}")
queryDB="${DIR}"/query_db
targetDB="${DIR}"/taxa_db
taxresDB="${DIR}"/taxres_db
taxres="${DIR}"/taxres.tsv
report="${DIR}"/taxres_report

sed -e 's/\x0//g' "${QUERY}" >"${DIR}"/query.fasta

# Create query and target DBs
"${MMSEQS_BIN}" createdb "${DIR}"/query.fasta "${queryDB}"

if [[ ! -e "${targetDB}" ]]; then
  "${MMSEQS_BIN}" databases "${TAXDB}" "${targetDB}" "${DIR}"/tmp
  #"${MMSEQS_BIN}" createdb "${DB_FA}" "${targetDB}"
  # Create a taxadb, the target DB should contain UniProt IDS (This scripts is using the UniProtKB)
  # The next module will download the Uniprot idmapping and ncbi-taxdump and map the identifier of the seqTaxDB to NCBI taxonomic identifier.
  # By default, createtaxdb downloads the Uniprot id mapping file (idmapping.dat.gz), and thus only supports Uniprot identifiers.
  "${MMSEQS_BIN}" createtaxdb "${targetDB}" "${DIR}"/tmp
fi

# Run the taxonomy search
# --lca-mode 4 (top-hit, default) -> assigns taxonomic labels based on the lowest common ancestor of all equal scoring top hits
rm -rf "${DIR}"/taxres*

"${MMSEQS_BIN}" taxonomy "${queryDB}" "${targetDB}" "${taxresDB}" "${DIR}"/tmp --threads "${NSLOTS}" -e 1e-05 --cov-mode 0 -c 0.6  --lca-mode 2 --tax-lineage --mpi-runner "${RUNNER}"

# Parse the results
"${MMSEQS_BIN}" createtsv "${queryDB}" "${taxresDB}" "${taxres}" --threads "${NSLOTS}"
sed -i 's/ /__/g' "${taxres}"

# Add cluster information
join -13 -21 <(sort -k3,3 "${CL_INFO}") \
  <(sort -k1,1 "${taxres}") >"${RES}"

sed -i 's/ /\t/g' "${RES}"

# Optional: create a taxonomy Kraken-style report that can be visualized using the interactive metagenomics data explorer Pavian
#"${MMSEQS_BIN}" taxonomyreport "${targetDB}" "${taxresDB}" "${report}" --threads "${NSLOTS}"
