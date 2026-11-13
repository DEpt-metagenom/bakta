#!/usr/bin/env bash

for db in $(ls bakta_db_light/*.dmnd); do
    base=$(basename $db .dmnd)
    echo "Converting $db to FASTA format..."
    diamond getseq --db $db --out bakta_db_light/$base.fasta
done