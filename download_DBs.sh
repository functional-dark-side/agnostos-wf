#!/bin/#!/usr/bin/env bash

MMSEQS="${PWD}"/bin/mmseqs

mkdir -p databases

cd databases


# Pfam database
if [ ! -s Pfam-A.hmm ]; then
  echo "Dowloading Pfam-A database"
  wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
  gunzip Pfam-A.hmm.gz
fi
# Pfam clans
if [ ! -s Pfam-A.clans.tsv.gz ]; then
  echo "Dowloading Pfam-A clan information"
  wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
fi
# Pfam HHBLITS DB
if [ ! -s pfam_hhm.ffdata ]; then
  echo "Dowloading Pfam hh-suite database"
  wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_31.0.tgz
  tar xvfz pfamA_31.0.tgz
fi

# Pfam list common domain terms
if [ ! -s Pfam-31_names_mod_01122019.tsv ]; then
  echo "Dowloading Pfam list of shared domain names"
  wget https://ndownloader.figshare.com/files/23756204 -O Pfam-31_names_mod_01122019.tsv
fi
# Antifam databases
if [ ! -s AntiFam.hmm ]; then
  echo "Dowloading AntiFam database"
  wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz
  tar xvfz Antifam.tar.gz
  ../bin/hmmpress AntiFam.hmm
fi
# Uniref90
if [ ! -s uniref90.proteins.tsv.gz ]; then
  echo "Dowloading UniRef90 DB"
  aria2c --file-allocation=none -c -x 10 -s 10 ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
  echo "Creating the MMseqs2 database"
  "${MMSEQS}" createdb uniref90.fasta.gz uniref90.db --write-lookup 0 -v 0
  # Create the protein description file:
  echo "extracting protein information"
  zcat uniref90.fasta.gz | grep '^>' | sed 's/^>//' | sed 's/ /\t/' | sed 's/ /_/g' | gzip > uniref90.proteins.tsv.gz
  rm uniref90.fasta.gz
fi

# NCBI nr
if [ ! -s nr.proteins.tsv.gz ]; then
  echo "Dowloading NCBI nr DB"
  aria2c --file-allocation=none -c -x 10 -s 10 ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
  echo "Creating the MMseqs2 database"
  "${MMSEQS}" createdb nr.gz nr.db --write-lookup 0 -v 0
  # Create the protein description file:
  echo "extracting protein information"
  zcat nr.gz | grep '^>' | sed 's/^>//' | sed 's/ /\t/' | sed 's/ /_/g' | gzip > nr.proteins.tsv.gz
  rm nr.gz
fi
# Uniclust HHBLITS DB
if [ ! -s uniclust30_2018_08/uniclust30_2018_08_a3m.ffdata ]; then
  echo "Dowloading Uniclust30 hh-suite database"
  aria2c --file-allocation=none -c -x 10 -s 10 http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz
  tar xvfz uniclust30_2018_08_hhsuite.tar.gz
fi

# DPD and info
if [ ! -s dpd_uniprot_sprot.fasta.gz ]; then
  echo "Dowloading DPD"
  wget https://ndownloader.figshare.com/files/23756312 -O dpd_uniprot_sprot.fasta.gz
  wget https://ndownloader.figshare.com/files/23756306 -O dpd_ids_all_info.tsv.gz
fi

# GTDB r89 taxonomy DB
if [ ! -s gtdb-r89_54k/gtdb-r89_54k.fmi ]; then
  echo "Dowloading GTDB-r89 kaiju DB"
  wget https://ndownloader.figshare.com/files/24745184 -O gtdb-r89_54k.tar.gz
  tar xzvf gtdb-r89_54k.tar.gz
fi

# Uniprot KB for MMseqs2 taxonomy
if [ ! -s uniprotKB ]; then
  echo "Dowloading UniProtKB mmseqs taxonomy DB"
  "${MMSEQS}" databases "UniProtKB" uniprotKB tmp --remove-tmp-files 1 -v 0
  rm uniprotKB.lookout
fi


cd ..
