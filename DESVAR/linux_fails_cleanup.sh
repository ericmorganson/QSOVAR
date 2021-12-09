#!/bin/bash

# Run this script for finding fits row runs  that didn't end successfully.
# {DIR} should be any of the prefixes for the fits files, like x3, c1, etc.
# {DIR}_fail_lines.txt is used in python_cleanup.sl to create a script to
# finish processing these unfinished fits. 
DIR=$1  #"e2"

cd /home/thrush2/QSOVAR/DESVAR/${DIR}_runlogs/.
find . -type f -empty -delete
find . -type f -name "*.txt" -exec grep -HL 'Successfully finished!' '{}' ';' > ../${DIR}_fails.txt

while read p; do
    grep "Running object" "$p" | tail -n 1
    echo "$p"
done < ../${DIR}_fails.txt > ../${DIR}_fail_lines.txt
