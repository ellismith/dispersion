#!/bin/bash
# Stages one subcluster's cutoff-filtered pseudobulk files into the exact
# naming convention dglm_model.R expects ({label}_{region}_pseudobulk.csv,
# {label}_metadata.csv, {label}_gene_names.csv), without touching
# dglm_model.R itself.
#
# Usage: ./stage_filtered_for_dglm.sh <raw_dir> <filtered_dir> <label> <cutoff> <staging_outdir>

set -e
RAW_DIR=$1
FILTERED_DIR=$2
LABEL=$3
CUTOFF=$4
STAGE_DIR=$5

mkdir -p "$STAGE_DIR"

cp "$RAW_DIR/${LABEL}_metadata.csv" "$STAGE_DIR/"
cp "$RAW_DIR/${LABEL}_gene_names.csv" "$STAGE_DIR/"

n_staged=0
for f in "$FILTERED_DIR/${LABEL}_"*"_filtered_cutoff${CUTOFF}.csv"; do
    [ -e "$f" ] || continue
    # strip "{LABEL}_" prefix and "_filtered_cutoff{CUTOFF}.csv" suffix to get the region
    base=$(basename "$f")
    region=${base#${LABEL}_}
    region=${region%_filtered_cutoff${CUTOFF}.csv}
    cp "$f" "$STAGE_DIR/${LABEL}_${region}_pseudobulk.csv"
    n_staged=$((n_staged + 1))
done

echo "Staged $n_staged region file(s) for $LABEL (cutoff $CUTOFF) into $STAGE_DIR"
