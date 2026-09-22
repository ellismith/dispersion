#!/bin/bash
# Runs pseudobulk.py (all louvain subclusters x all regions) + the CPM
# filter (cutoffs 0.5/1/2, >= "at least" logic, audit trail) for every
# remaining cell type. astrocytes and microglia already done separately.
#
# Usage: ./run_all_celltypes.sh

set -e

CELL_TYPES=(
  basket_cells
  cerebellar_neurons
  ependymal_cells
  GABAergic_neurons
  glutamatergic_neurons
  medium_spiny_neurons
  midbrain_neurons
  oligodendrocytes
  opc
  vascular_cells
)

BASE_OUT=/scratch/easmit31/dispersion/dglm/disp_age__raw_counts
SCRIPTS=/scratch/easmit31/dispersion/dglm/scripts

for ct in "${CELL_TYPES[@]}"; do
  echo "=================================================="
  echo "=== $ct : pseudobulk ==="
  echo "=================================================="
  mkdir -p "$BASE_OUT/$ct"
  conda run -n latent_analysis python "$SCRIPTS/pseudobulk.py" \
    --cell_type "$ct" \
    --outdir "$BASE_OUT/$ct" \
    --subcluster_col ct_louvain

  echo "=================================================="
  echo "=== $ct : filtering ==="
  echo "=================================================="
  mkdir -p "$BASE_OUT/$ct/filtered"
  "$SCRIPTS/filter_all_pseudobulk.sh" \
    "$BASE_OUT/$ct" \
    "$BASE_OUT/$ct/filtered"

  echo "=== $ct done ==="
  echo ""
done

echo "ALL CELL TYPES COMPLETE"
