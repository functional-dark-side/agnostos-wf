#!/bin/bash

set -e
set -x

usage() {
    echo -e "Usage: $(basename "$0") [OPTIONS...]
Options:
--cluseq_db  path to cluster sequence db \
--step       name of the refinement step \
--mpi_runner mpi-mode runner \
--mmseqs      program to extract subdb \
--aln        msa program \
--hhsuite    path to the hh-suite location \
--reformat   script to reformat msa \
--consensus  script to extract consensus sequences \
--hhmake     script to build the HMM profiles \
--output     path to the output folder \
--idir       path to the input folder \
--clu_hhm    path to cluster hmm db \
  --threads    number of threads to use" 1>&2
    exit 1
}

OPT_LIST="cluseq_db:,step:,mpi_runner:,mmseqs:,aln:,hhsuite:,reformat:,consensus:,hhmake:,outdir:,idir:,clu_hhm:,threads:"

eval set -- "$(getopt -o '' --long "${OPT_LIST}" -- "$@")"

while true; do
    case "$1" in
    --cluseq_db)
        CLU=$2
        shift 2
        ;;
    --step)
        STEP=$2
        shift 2
        ;;
    --mpi_runner)
        RUNNER=$2
        shift 2
        ;;
    --mmseqs)
        MMSEQS_BIN=$2
        shift 2
        ;;
    --aln)
        ALN=$2
        shift 2
        ;;
    --hhsuite)
        HHSUITE=$2
        shift 2
        ;;
    --reformat)
        REFORM=$2
        shift 2
        ;;
    --consensus)
        CONS=$2
        shift 2
        ;;
    --hhmake)
        HHMAKE=$2
        shift 2
        ;;
    --outdir)
        RES=$2
        shift 2
        ;;
    --idir)
        IDIR=$2
        shift 2
        ;;
    --clu_hhm)
        CLU_HHM=$2
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

if [ -z "${CLU}" ] || [ -z "${STEP}" ] || [ -z "${RUNNER}" ] ||
    [ -z "${MMSEQS_BIN}" ] || [ -z "${ALN}" ] || [ -z "${HHSUITE}" ] || [ -z "${IDIR}" ] ||
    [ -z "${REFORM}" ] || [ -z "${CONS}" ] || [ -z "${HHMAKE}" ] || [ -z "${RES}" ] || [ -z "${NSLOTS}" ]; then
    usage
fi

if [[ ${STEP} = "known refinement" ]]; then
    CATEG=$(echo -e "kwp")
elif [[ ${STEP} = "unknown refinement" ]]; then
    CATEG=$(echo -e "eu")
elif [[ ! -s "${IDIR}"/eu_ids.txt ]]; then
    CATEG=$(echo -e "gu\nkwp\nk")
elif [[ ! -s "${IDIR}"/gu_ids.txt ]]; then
    CATEG=$(echo -e "eu\nkwp\nk")
elif [[ ! -s "${IDIR}"/kwp_ids.txt ]]; then
    CATEG=$(echo -e "eu\ngu\nk")
else
    CATEG=$(echo -e "eu\ngu\nkwp\nk")
fi

#Create subDB for each cluster categ to retrieve the MSAs, the consensus and the profiles; using the cluster ids/names
#Ex: Genomic unknowns = gu

for categ in $CATEG; do

    "${MMSEQS_BIN}" createsubdb "${IDIR}"/"${categ}"_ids.txt \
        "${CLU}" "${RES}"/"${categ}"_clseqdb

    if [[ ${STEP} = "category database" ]]; then
        # Retrieve set of ORFs for each category
        sed -e 's/\x0//g' "${RES}"/"${categ}"_clseqdb >"${RES}"/"${categ}"_clu_orfs.fasta

        grep '^>' "${RES}"/"${categ}"_clu_orfs.fasta | sed 's/^>//' >"${RES}"/"${categ}"_orfs.txt
    fi

    # Retrieve alignments, consensus sequences and HMMs
    ${RUNNER} "${MMSEQS_BIN}" apply "${RES}"/"${categ}"_clseqdb "${RES}"/"${categ}"_aln \
        --threads 1 \
        -- "${ALN}" STDIN STDOUT 2>/dev/null

    ${RUNNER} "${MMSEQS_BIN}" apply "${RES}"/"${categ}"_aln "${RES}"/"${categ}"_a3m_db \
        --threads 1 \
        -- "${REFORM}" "${HHSUITE}"

    ln -sf "${RES}"/"${categ}"_aln "${RES}"/"${categ}"_aln.ffdata
    ln -sf "${RES}"/"${categ}"_aln.index "${RES}"/"${categ}"_aln.ffindex
    ln -sf "${RES}"/"${categ}"_a3m_db "${RES}"/"${categ}"_a3m.ffdata
    ln -sf "${RES}"/"${categ}"_a3m_db.index "${RES}"/"${categ}"_a3m.ffindex

    if [[ ${STEP} = "category database" ]]; then
        ${RUNNER} "${MMSEQS_BIN}" apply "${RES}"/"${categ}"_aln "${RES}"/"${categ}"_cons \
                                        --threads 1 \
                                        -- "${CONS}" "${HHSUITE}"/bin/hhconsensus

         NCL=$(cat "${RES}"/"${categ}"_aln.index | wc -l)
         if [[ "${NCL}" -le 10000 ]]; then
            "${HHSUITE}"/bin/cstranslate -i "${RES}"/"${categ}"_a3m -o "${RES}"/"${categ}"_cs219 \
             -A "${HHSUITE}"/data/cs219.lib -D "${HHSUITE}"/data/context_data.lib \
             -x 0.3 -c 4 -b -f -I a3m &>/dev/null
        else
             ${RUNNER} "${HHSUITE}"/bin/cstranslate_mpi -i "${RES}"/"${categ}"_a3m -o "${RES}"/"${categ}"_cs219 \
              -A "${HHSUITE}"/data/cs219.lib -D "${HHSUITE}"/data/context_data.lib \
              -x 0.3 -c 4 -b -f -I a3m &>/dev/null
         fi
         ln -sf "${RES}"/"${categ}"_cs219.ffdata "${RES}"/"${categ}"_cs219
         ln -sf "${RES}"/"${categ}"_cs219.ffindex "${RES}"/"${categ}"_cs219.index
    fi

    ${RUNNER} "${MMSEQS_BIN}" apply "${RES}"/"${categ}"_a3m_db "${RES}"/"${categ}"_hhm_db \
        --threads 1 \
        -- "${HHMAKE}" "${HHSUITE}"/bin/hhmake

    ln -sf "${RES}"/"${categ}"_hhm_db "${RES}"/"${categ}"_hhm_db.ffdata
    ln -sf "${RES}"/"${categ}"_hhm_db.index "${RES}"/"${categ}"_hhm_db.ffindex
    ln -sf "${RES}"/"${categ}"_hhm_db "${RES}"/"${categ}"_hhm.ffdata
    ln -sf "${RES}"/"${categ}"_hhm_db.index "${RES}"/"${categ}"_hhm.ffindex
done

echo "Done building hh-suite databases for ${categ}"

if [[ ${STEP} = "known refinement" ]]; then
    rm -rf "${RES}"/"${categ}"_aln* "${RES}"/"${categ}"_a3m* "${RES}"/"${categ}"_clseqdb*
elif [[ ${STEP} = "unknown refinement" ]]; then
    rm -rf "${RES}"/"${categ}"_aln* "${RES}"/"${categ}"_a3m* "${RES}"/"${categ}"_clseqdb*
fi

if [[ ${STEP} = "category database" ]]; then

  for categ in $CATEG; do
    #mv "${RES}"/"${categ}"_aln "${RES}"/"${categ}"_a3m_db
    #mv "${RES}"/"${categ}"_aln.index "${RES}"/"${categ}"_a3m_db.index
    #mv "${RES}"/"${categ}"_aln.dbtype "${RES}"/"${categ}"_a3m_db.dbtype
    #ln -sf "${RES}"/"${categ}"_a3m_db "${RES}"/"${categ}"_a3m.ffdata
    #ln -sf "${RES}"/"${categ}"_a3m_db.index "${RES}"/"${categ}"_a3m.ffindex
    #rm "${RES}"/"${categ}"_aln.ff*
    if [[ -f "${RES}"/"${categ}"_cs219.log* ]]; then
      rm "${RES}"/"${categ}"_cs219.log*
    fi
  done

  if [[ ! -s "${IDIR}"/eu_ids.txt ]]; then
    "${MMSEQS_BIN}" concatdbs "${RES}"/k_hhm_db "${RES}"/kwp_hhm_db "${RES}"/tmp_hhm_db --threads 1 --preserve-keys

    "${MMSEQS_BIN}" concatdbs "${RES}"/gu_hhm_db "${RES}"/tmp_hhm_db "${RES}"/clu_hhm_db --threads 1 --preserve-keys

    rm "${RES}"/tmp_hhm_db*
  else
    "${MMSEQS_BIN}" concatdbs "${RES}"/eu_hhm_db "${RES}"/gu_hhm_db "${RES}"/tmp_hhm_db --threads 1 --preserve-keys

    "${MMSEQS_BIN}" concatdbs "${RES}"/k_hhm_db "${RES}"/kwp_hhm_db "${RES}"/tmp1_hhm_db --threads 1 --preserve-keys

    "${MMSEQS_BIN}" concatdbs "${RES}"/tmp_hhm_db "${RES}"/tmp1_hhm_db "${RES}"/clu_hhm_db --threads 1 --preserve-keys

    rm "${RES}"/tmp_hhm_db* "${RES}"/tmp1_hhm_db*
  fi
  echo "Done building clusters hhm databases"
fi
