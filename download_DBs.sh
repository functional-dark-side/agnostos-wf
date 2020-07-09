#!/bin/#!/usr/bin/env bash

mkdir -p databases

cd databases

# Pfam database
if [ ! -s Pfam-A.hmm ]; then
  wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
  gunzip Pfam-A.hmm.gz
fi
# Pfam clans
if [ ! -s Pfam-A.clans.tsv.gz ]; then
  wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
fi
# Pfam HHBLITS DB
if [ ! -s pfam_hhm.ffdata ]; then
  wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_31.0.tgz
  tar xvfz pfamA_31.0.tgz
fi

# Pfam list common domain terms
if [ ! -s Pfam-31_names_mod_01122019.tsv ]; then
  wget https://ndownloader.figshare.com/files/23756204 -O Pfam-31_names_mod_01122019.tsv
fi
# Antifam databases
if [ ! -s AntiFam.hmm ]; then
  wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz
  tar xvfz Antifam.tar.gz
  ../bin/hmmpress AntiFam.hmm
fi
# Uniref90
if [ ! -s uniref90.fasta.gz ]; then
  aria2c --file-allocation=none -c -x 10 -s 10 ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
  # Create the protein description file:
  zcat uniref90.fasta.gz | grep '^>' | sed 's/^>//' | sed 's/ /\t/' | sed 's/ /_/g' | gzip > uniref90.proteins.tsv.gz
fi

# NCBI nr
if [ ! -s nr.fasta.gz ]; then
  aria2c --file-allocation=none -c -x 10 -s 10 ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
  mv nr.gz nr.fasta.gz
  # Create the protein description file:
  zcat nr.fasta.gz | grep '^>' | sed 's/^>//' | sed 's/ /\t/' | sed 's/ /_/g' | gzip > nr.proteins.tsv.gz
fi
# Uniclust HHBLITS DB
if [ ! -s uniclust30_2018_08/uniclust30_2018_08_a3m.ffdata ]; then
  aria2c --file-allocation=none -c -x 10 -s 10 http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08.tar.gz
  tar xvfz uniclust30_2018_08.tar.gz
fi

# DPD and info
if [ ! -s dpd_uniprot_sprot.fasta.gz ]; then
  wget https://ndownloader.figshare.com/files/23756312 -O dpd_uniprot_sprot.fasta.gz
  wget https://ndownloader.figshare.com/files/23756306 -O dpd_ids_all_info.tsv.gz
fi
# Uniprot KB
# wget ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
# wget wget ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
# zcat uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz | gzip > uniprotKB.fasta.gz
# rm uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz

# NCBI Taxonomy
# mkdir ncbi-taxdump && cd ncbi-taxdump
# wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
# tar xzvf taxdump.tar.gz
# cd ..

cd ..
