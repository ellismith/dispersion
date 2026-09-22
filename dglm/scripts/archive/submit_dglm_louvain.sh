#!/bin/bash
SCRIPT=/scratch/easmit31/dispersion/dglm/scripts/dglm_model.R
BASEDIR=/scratch/easmit31/dispersion/dglm/disp_age__mean_full_louvain
LOGDIR=/scratch/easmit31/dispersion/dglm/logs
RSCRIPT=/home/easmit31/.conda/envs/mashr_env/bin/Rscript
mkdir -p $LOGDIR

declare -A MEM_MAP
MEM_MAP[glutamatergic_neurons]=256G
MEM_MAP[GABAergic_neurons]=128G
MEM_MAP[cerebellar_neurons]=128G

CELL_TYPES=(astrocytes basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons glutamatergic_neurons medium_spiny_neurons midbrain_neurons opc oligodendrocytes vascular_cells)

for CT in "${CELL_TYPES[@]}"; do
    MEM=${MEM_MAP[$CT]:-64G}
    CHECKPOINTS=${BASEDIR}/${CT}
    sbatch \
        --job-name=dglm_${CT}_louvain \
        --output=${LOGDIR}/dglm_${CT}_louvain_%j.out \
        --error=${LOGDIR}/dglm_${CT}_louvain_%j.err \
        --mem=${MEM} --time=03:59:00 -p htc \
        --cpus-per-task=16 \
        --wrap="cd /scratch/easmit31/dispersion/dglm && ${RSCRIPT} ${SCRIPT} --cell_type ${CT} --checkpoints ${CHECKPOINTS} --covariates full --dispersion age"
done
echo "all dglm jobs submitted"
