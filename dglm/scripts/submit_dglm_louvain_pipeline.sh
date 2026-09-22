#!/bin/bash
# Submits one Slurm job per remaining cell type, each running the FULL
# chain: dglm (every louvain subcluster) -> combined mashr (one fit across
# that cell type's subcluster x region conditions) -> master TSV.
# astrocytes and microglia excluded -- already done/in progress separately.
#
# Assumes pseudobulk + CPM filtering (cutoff 0.5) already ran for each cell
# type via run_all_celltypes.sh / submit_remaining_pseudobulk.sh -- fails
# fast per-job if that's missing rather than silently proceeding on
# incomplete data.

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

CELL_TYPES=(basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons glutamatergic_neurons medium_spiny_neurons midbrain_neurons oligodendrocytes opc vascular_cells)

for CT in "${CELL_TYPES[@]}"; do
    MEM=${MEM_MAP[$CT]:-64G}
    RAW_DIR="$BASE/$CT"
    FILTERED_DIR="$BASE/$CT/filtered"
    STAGE_BASE="$BASE/$CT/dglm_run_cutoff${CUTOFF}"
    CHECKPOINTS="$BASE/$CT/dglm_checkpoints_cutoff${CUTOFF}"

    sbatch \
        --job-name=dglm_pipeline_${CT} \
        --output=${LOGDIR}/dglm_pipeline_${CT}_%j.out \
        --error=${LOGDIR}/dglm_pipeline_${CT}_%j.err \
        --mem=${MEM} --time=03:59:00 -p htc \
        --cpus-per-task=16 \
        --wrap="\
set -e; \
if [ ! -d '${RAW_DIR}' ] || [ -z \"\$(ls -A ${RAW_DIR}/*_metadata.csv 2>/dev/null)\" ]; then \
    echo 'ERROR: no pseudobulk metadata found in ${RAW_DIR}' >&2; \
    exit 1; \
fi; \
${SCRIPTS}/run_dglm_for_celltype.sh '${RAW_DIR}' '${FILTERED_DIR}' '${CUTOFF}' '${STAGE_BASE}' '${CHECKPOINTS}'; \
${RSCRIPT} ${SCRIPTS}/dglm_mashr.R --mode combined --checkpoints '${CHECKPOINTS}' --shat_mode sqrt --fast --random_subset_n 3000 --gridmult 2.5; \
${RSCRIPT} ${SCRIPTS}/dglm_fdr_combined.R --checkpoints '${CHECKPOINTS}'; \
echo 'DONE: ${CT}'"
done

echo "Submitted ${#CELL_TYPES[@]} jobs. Check with: squeue -u \$USER"
