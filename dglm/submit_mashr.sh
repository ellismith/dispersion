#!/bin/bash

SCRIPT=/scratch/easmit31/variability/dglm/scripts/dglm_mashr.R
CHECKPOINTS=/scratch/easmit31/variability/dglm/checkpoints
LOGDIR=/scratch/easmit31/variability/dglm/logs
RSCRIPT=/home/easmit31/.conda/envs/mashr_env/bin/Rscript
mkdir -p $LOGDIR

CELL_TYPES=(astrocytes basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons glutamatergic_neurons medium_spiny_neurons microglia midbrain_neurons opc oligodendrocytes vascular_cells)

for CT in "${CELL_TYPES[@]}"; do
    sbatch \
        --job-name=mashr_${CT} \
        --output=${LOGDIR}/mashr_${CT}_%j.out \
        --error=${LOGDIR}/mashr_${CT}_%j.err \
        --mem=64G --time=03:59:00 -p htc \
        --wrap="cd /scratch/easmit31/variability/dglm && ${RSCRIPT} ${SCRIPT} --cell_type ${CT} --checkpoints ${CHECKPOINTS}"
done

echo "all jobs submitted"
