#!/bin/bash
# V3: louvain resolution, dispersion ~ age + sex. Reuses the SAME existing
# cutoff-0.5 pseudobulk/filtered CSVs already on disk for each cell type
# (filtering doesn't depend on the dispersion formula, so no regeneration
# needed) -- just runs dglm with --dispersion age_sex into separate
# *_agesex checkpoints dirs so V2 (age-only) results are never touched.

SCRIPTS=/scratch/easmit31/dispersion/dglm/scripts
BASE=/scratch/easmit31/dispersion/dglm/disp_age__raw_counts
LOGDIR=/scratch/easmit31/dispersion/dglm/logs
RSCRIPT=/home/easmit31/.conda/envs/mashr_env/bin/Rscript
CUTOFF=0.5

mkdir -p "$LOGDIR"

declare -A MEM_MAP
MEM_MAP[glutamatergic_neurons]=256G
MEM_MAP[GABAergic_neurons]=128G
MEM_MAP[cerebellar_neurons]=128G

CELL_TYPES=(astrocytes basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons glutamatergic_neurons medium_spiny_neurons microglia midbrain_neurons oligodendrocytes opc vascular_cells)

for CT in "${CELL_TYPES[@]}"; do
    MEM=${MEM_MAP[$CT]:-64G}
    RAW_DIR="$BASE/$CT"
    FILTERED_DIR="$BASE/$CT/filtered"
    STAGE_BASE="$BASE/$CT/dglm_run_cutoff${CUTOFF}_agesex"
    CHECKPOINTS="$BASE/$CT/dglm_checkpoints_cutoff${CUTOFF}_agesex"

    sbatch \
        --job-name=dglm_v3_${CT} \
        --output=${LOGDIR}/dglm_v3_${CT}_%j.out \
        --error=${LOGDIR}/dglm_v3_${CT}_%j.err \
        --mem=${MEM} --time=03:59:00 -p htc \
        --cpus-per-task=16 \
        --wrap="\
set -e; \
if [ ! -d '${RAW_DIR}' ] || [ -z \"\$(ls -A ${RAW_DIR}/*_metadata.csv 2>/dev/null)\" ]; then \
    echo 'ERROR: no pseudobulk metadata found in ${RAW_DIR}' >&2; \
    exit 1; \
fi; \
${SCRIPTS}/run_dglm_for_celltype.sh '${RAW_DIR}' '${FILTERED_DIR}' '${CUTOFF}' '${STAGE_BASE}' '${CHECKPOINTS}' age_sex; \
${RSCRIPT} ${SCRIPTS}/dglm_mashr.R --mode combined --checkpoints '${CHECKPOINTS}' --shat_mode sqrt --fast --random_subset_n 3000 --gridmult 2.5; \
${RSCRIPT} ${SCRIPTS}/dglm_fdr_combined.R --checkpoints '${CHECKPOINTS}'; \
echo 'DONE: ${CT}'"
done

echo "Submitted ${#CELL_TYPES[@]} V3 (age+sex) jobs. Check with: squeue -u \$USER"
