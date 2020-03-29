#!/bin/bash

# Cluster naming
# The official cluster names are going to be based on the line number of the wide formatted file
# We are also going to produce a correspondence file to access the clusters in the MMseqs2 indices

SEQS=$(perl -ne 'print $_')

NAME=$MMSEQS_ENTRY_NAME

echo "${SEQS}" | awk '/^>/{split($1,a,"- OS"); print a[1]; next}1' | grep '^>' | sed "s/^>/$NAME /"
