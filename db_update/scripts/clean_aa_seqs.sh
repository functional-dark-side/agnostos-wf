#!/bin/bash

set -e

SEQS=$(perl -ne 'print $_')

echo "${SEQS}" | awk '/^>/{split($1,a,"- OS"); print a[1]; next}1' \
| seqkit fx2tab | awk '!seen[$0]++' | seqkit tab2fx | seqkit replace -s -p "\*" -r ""
