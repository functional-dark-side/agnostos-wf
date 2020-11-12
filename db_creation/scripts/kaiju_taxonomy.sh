#!/bin/bash

set -e
set -x

usage() {
  echo -e "Usage: $(basename $0) [OPTIONS...]
Options:
--search             program to perform the search (MMSEQS)
--input              query sequences (refined db) \
--db_fasta           database sequences fasta format with unirpot ids \
--parsing            file containing the corresp. orf-cl_name-category \
--output             path to the output file \
--tmpl               template folder \
    --threads            number of threads to use" 1>&2
  exit 1
}

OPT_LIST="search:,input:,taxdb:,output:,threads:,parsing:,tmpl:"

eval set -- $(getopt -o '' --long "${OPT_LIST}" -- "$@")
while true; do
  case "$1" in
    --search)
      KAIJU_BIN=$2
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
    --parsing)
      PARSER=$2
      shift 2
      ;;
    --output)
      RES=$2
      shift 2
      ;;
    --tmpl)
      TMP=$2
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

if [ -z "${KAIJU_BIN}" ] ||
[ -z "${QUERY}" ] ||
[ -z "${TAXDB}" ] ||
[ -z "${PARSER}" ] ||
[ -z "${RES}" ] ||
[ -z "${TMP}" ] ||
[ -z "${NSLOTS}" ]; then
  usage
fi

DIR_DB=$(dirname "${TAXDB}")
DIR=$(dirname "${RES}")
TMP="${DIR}"/tmp
taxresDB="${DIR}"/kaiju_taxres.out

sed -e 's/\x0//g' "${QUERY}" > "${DIR}"/query.fasta

"${KAIJU_BIN}" -f "${TAXDB}" \
  -z  "${NSLOTS}" \
  -t "${DIR_DB}"/nodes.dmp \
  -o "${taxresDB}" \
  -i "${DIR}"/query.fasta \
  -p -v


"${PARSER}" -r "${taxresDB}" \
  -o "${taxresDB}".parsed \
  --names "${DIR_DB}"/names.dmp \
  --nodes "${DIR_DB}"/nodes.dmp \
  --tmp "${TMP}"

sed 's/ /_/g' "${taxresDB}".parsed | \
  awk '{print $2"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' > "${RES}"

rm -rf "${DIR}"/query.fasta
