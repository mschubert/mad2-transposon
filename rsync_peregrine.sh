#!/bin/bash

rsync -auvr --include='*.pdf' \
    --exclude '.snakemake' \
#    --include 'data/singlecell/BB150521_I_merged.bam{,.bai}' \
#    --exclude 'data/*' \
    --include '*.RData' \
    --include='*/' \
    --include='docs/*' \
    --exclude='*' \
    rug:/data/p282396/mad2-transposon/* $(dirname $0)

#rsync -avur figure/* yoda:/nfs/research2/saezrodriguez/mike/speed2/figure/
