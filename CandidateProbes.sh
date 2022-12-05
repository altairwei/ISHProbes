#!/usr/bin/env bash

source activate probeMining

python OligoMiner/blockParse.py \
    -l 150 -L 200 -g 50 -G 60 \
    -n RNA_NN3 --verbose \
    -f data/Triticum_aestivum.IWGSC.dna.toplevel.fa