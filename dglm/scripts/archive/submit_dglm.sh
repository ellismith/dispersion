#!/bin/bash

SCRIPT=/scratch/easmit31/variability/dglm/scripts/dglm_model.R
CHECKPOINTS=/scratch/easmit31/variability/dglm/checkpoints
LOGDIR=/scratch/easmit31/variability/dglm/logs
RSCRIPT=/home/easmit31/.conda/envs/mashr_env/bin/Rscript
mkdir -p $LOGDIR

declare -A MEM_MAP
MEM_MAP[glutamatergic_neurons]=256G
MEM_MAP[GABAergic_neurons]=128G
MEM_MAP[cerebellar_neurons]=128G

CELL_TYPES=(astrocytes basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons glutamatergic_neurons medium_spiny_neurons microglia midbrain_neurons opc oligodendrocytes vascular_cells)

for CT in "${CELL_TYPES[@]}"; do
    MEM=${MEM_MAP[$CT]:-64G}
    sbatch \
        --job-name=dglm_${CT} \
        --output=${LOGDIR}/dglm_${CT}_%j.out \
        --error=${LOGDIR}/dglm_${CT}_%j.err \
        --mem=${MEM} --time=03:59:00 -p htc \
        --cpus-per-task=16 \
        --wrap="cd /scratch/easmit31/variability/dglm && ${RSCRIPT} ${SCRIPT} --cell_type ${CT} --checkpoints ${CHECKPOINTS}"
done

echo "all jobs submitted"
