#!/bin/bash

HHCONS="${1}"

"${HHCONS}" -M 50 -maxres 65535 -i stdin -s stdout -v 0 \
  | awk -v name="${MMSEQS_ENTRY_NAME}" 'NR==1{print ">"name; next}{print $0}'
