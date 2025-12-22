#!/bin/sh

# Following directions from CLEANSER
# Paper: https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00022-9 
# Github: https://github.com/Gersbachlab-Bioinformatics/CLEANSER?tab=readme-ov-file

# File directory for 10X output
tenX_out="/mnt/mass_storage1/mass_storage_projects/compass/COM-iN-none-aggr/count/batch1/"
echo $tenX_out

# Converting to a suitable matrix format for guide information only from CellRanger output
cr2cleanser --matrix-market ${tenX_out}matrix.mtx.gz \
    --features ${tenX_out}features.tsv.gz \
    --output CLEANSER/matrix.mtx

# Running cleanser - this outputs two tsv files. You can use "read.csv" with the "sep = '\t'" option in R to import them!
cleanser --input CLEANSER/matrix.mtx \
    --posteriors-output CLEANSER/posteriors-output \
    --samples-output CLEANSER/sample-output \
    --seed 67 \
    --parallel-runs 10 \
    --direct-capture


