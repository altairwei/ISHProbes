#!/usr/bin/env bash

source activate bioinfo

gfftools seq --type gene \
    -g data/Triticum_aestivum.IWGSC.dna.toplevel.fa \
    --line-length 60 \
    --fasta-header "attributes['ID']" \
    data/Triticum_aestivum.IWGSC.48.gff3 > data/genes.fa