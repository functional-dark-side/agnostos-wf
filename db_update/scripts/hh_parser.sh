#!/bin/bash
FILE=$(perl -ne 'print $_')

python "${1}" <(echo "$FILE") | awk '$2 >= 90' | awk -F"OS=" '{$0=$1}1'
