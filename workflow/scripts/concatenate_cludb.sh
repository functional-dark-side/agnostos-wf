#!/bin/bash

set -x
set -e

cludb="${1}"
threads="${2}"

N=$(($threads-1))

for i in $(seq 0 $N); do cat "${cludb}".$i >> "${cludb}"; done

for i in $(seq 0 $N); do rm "${cludb}".$i; done
