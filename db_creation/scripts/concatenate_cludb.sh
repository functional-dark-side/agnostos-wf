#!/bin/bash

cludb="${1}"
threads="${2}"

N=$(("${threads}" - 1))

cat "${cludb}".{0..$N} > "${cludb}"

rm "${cludb}".{0..$N}
