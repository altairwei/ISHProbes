#!/usr/bin/env bash

source activate bioinfo

CSVFILE=$1

# Filter table with minimum MFE
cat ${CSVFILE} \
    | csvtk filter -H -f "2>-50" \
    | csvtk pretty -H
