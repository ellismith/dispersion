#!/bin/bash
# Applies filter_percent_animals.R (cutoffs 0.5, 1, 2) to every
# *_pseudobulk.csv file in a directory -- i.e. every louvain-subcluster x
# region combo produced by pseudobulk.py.
#
# Usage: ./filter_all_pseudobulk.sh <pseudobulk_dir> <outdir>

set -e
PB_DIR=$1
OUTDIR=$2
mkdir -p "$OUTDIR"

for f in "$PB_DIR"/*_pseudobulk.csv; do
    label=$(basename "$f" | sed 's/_pseudobulk.csv$//')
    echo "=== $label ==="
    Rscript /scratch/easmit31/dispersion/dglm/scripts/filter_percent_animals.R "$f" "$OUTDIR" "$label"
done
