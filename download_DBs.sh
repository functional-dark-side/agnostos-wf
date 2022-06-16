#!/bin/#!/usr/bin/env bash

MMSEQS="${PWD}"/bin/mmseqs

mkdir -p databases

cd databases

# Pfam database
if [ ! -s Pfam-A.hmm ]; then
  echo "Dowloading Pfam-A database"
  wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz
  gunzip Pfam-A.hmm.gz
fi
# Pfam clans
if [ ! -s Pfam-A.clans.tsv.gz ]; then
  echo "Dowloading Pfam-A clan information"
  wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.clans.tsv.gz
fi
# Pfam HHBLITS DB
if [ ! -s pfam_hhm.ffdata ]; then
  echo "Dowloading Pfam hh-suite database"
  wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_34.0.tar.gz
  tar xvfz pfamA_34.0.tar.gz
fi

# Pfam list common domain terms
if [ ! -s Pfam-34_names_mod_20102021.tsv ]; then
  echo "Dowloading Pfam list of shared domain names"
  wget https://figshare.com/ndownloader/files/31127782 -O Pfam-34_names_mod_20102021.tsv
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
  echo "Dowloading UniRef90 DB using mmseqs"
  "${MMSEQS}" databases UniRef90 uniref90.db tmp --threads 28 --remove-tmp-files
  # Create the protein description file:
  echo "extracting protein information"
  sed 's/\x0//g' uniref90.db_h | sed 's/ /__/g' | sed 's/__/\t/' | sed 's/__/_/g' | gzip > uniref90.proteins.tsv.gz
fi

# NCBI nr
if [ ! -s nr.proteins.tsv.gz ]; then
  echo "Dowloading NR DB using mmseqs"
  "${MMSEQS}" databases NR nr.db tmp --threads 28 --remove-tmp-files
  # Create the protein description file:
  echo "extracting protein information"
  sed 's/\x0//g' nr.db_h | sed 's/ /__/g' | sed 's/__/\t/' | sed 's/__/_/g' | gzip > nr.proteins.tsv.gz
fi

# Uniclust HHBLITS DB
if [ ! -s UniRef30_2021_03_a3m.ffdata ]; then
  echo "Dowloading Uniclust hh-suite database"
  aria2c --file-allocation=none -c -x 10 -s 10 http://wwwuser.gwdg.de/~compbiol/uniclust/2021_03/UniRef30_2021_03.tar.gz
  tar xvfz UniRef30_2021_03.tar.gz
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

# Mutant phenotypes (Price et al. 2018)
if [ ! -s aaseqs ]; then
  ## Amino acid sequences
  wget https://fit.genomics.lbl.gov/cgi_data/aaseqs
fi
if [ ! -s feba.db ]; then
  ## Contextual data
  wget https://fit.genomics.lbl.gov/cgi_data/feba.db
fi

# Uniprot KB for MMseqs2 taxonomy
#if [ ! -s uniprotKB ]; then
#  echo "Dowloading UniProtKB mmseqs taxonomy DB"
#  "${MMSEQS}" databases "UniProtKB" uniprotKB tmp --remove-tmp-files 1 -v 0
#  rm uniprotKB.lookout
#fi


cd ..
