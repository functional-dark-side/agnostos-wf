#!/bin/bash

HHMAKE="${1}"

"${HHMAKE}" -nocontxt -diff 1000 -add_cons -i stdin -o stdout -v 0 -maxres 65535 -name "${MMSEQS_ENTRY_NAME}"

# Add -M50 if using the *_aln files
