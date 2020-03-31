#!/bin/bash

SEQS=$(perl -ne 'print $_')
TMP=$(mktemp -q)
HHSUITE="${1}"

"${HHSUITE}"/scripts/reformat.pl -v 0 -M 50 fas a3m <(echo "${SEQS}") "${TMP}"

cat "${TMP}"

rm "${TMP}"
