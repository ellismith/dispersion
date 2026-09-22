#!/bin/bash
# For EVERY subcluster of one cell type (natural numeric order): stages its
# cutoff-filtered pseudobulk files, runs dglm_model.R (with the given
# --dispersion mode), and copies the resulting *_dglm_results.rds into one
# SHARED checkpoints directory for that cell type + dispersion mode.
#
# IDEMPOTENT: skips any subcluster whose .rds already exists in the
# checkpoints dir.
#
# Usage: ./run_dglm_for_celltype.sh <raw_dir> <filtered_dir> <cutoff> <stage_base_dir> <shared_checkpoints_dir> [dispersion_mode]
#   dispersion_mode: 'age' (default, V2) or 'age_sex' (V3)

set -e
RAW_DIR=$1
FILTERED_DIR=$2
CUTOFF=$3
STAGE_BASE=$4
CHECKPOINTS=$5
DISPERSION=${6:-age}

SCRIPTS=/scratch/easmit31/dispersion/dglm/scripts
mkdir -p "$CHECKPOINTS"

for meta in $(ls "$RAW_DIR"/*_metadata.csv 2>/dev/null | sort -V); do
    label=$(basename "$meta" | sed 's/_metadata.csv$//')

    if [ -f "$CHECKPOINTS/${label}_dglm_results.rds" ]; then
        echo "=== $label: already in checkpoints, skipping ==="
        continue
    fi

    stage_dir="$STAGE_BASE/$label"
    echo "=== $label (cutoff $CUTOFF, dispersion $DISPERSION) ==="

    if [ -f "$stage_dir/${label}_dglm_results.rds" ]; then
        echo "  dglm already ran (found in $stage_dir), just copying to checkpoints"
    else
        "$SCRIPTS/stage_filtered_for_dglm.sh" "$RAW_DIR" "$FILTERED_DIR" "$label" "$CUTOFF" "$stage_dir"
        conda run -n mashr_env Rscript "$SCRIPTS/dglm_model.R" \
            --cell_type "$label" \
            --outdir "$stage_dir" \
            --dispersion "$DISPERSION"
    fi

    cp "$stage_dir/${label}_dglm_results.rds" "$CHECKPOINTS/"
done

echo ""
echo "All subclusters done. Shared checkpoints dir:"
ls "$CHECKPOINTS"
