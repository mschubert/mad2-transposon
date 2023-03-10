#!/bin/bash

rsync -auvr $@ --include='*.pdf' \
    --exclude '.snakemake' \
    --include '*.html' \
    --include '*.rds' \
    --include '*.xls' \
    --include '*.xlsx' \
    --include '*.tsv' \
    --include='*/' \
    --include='docs/*' \
    --exclude='*' \
    harris:/home/m.schubert/mad2-transposon/* $(dirname $0)
