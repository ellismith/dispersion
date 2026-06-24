#!/bin/bash

SCRIPT_DIR=/scratch/easmit31/variability/dglm/scripts
OUTDIR=/scratch/easmit31/variability/dglm/checkpoints
LOGDIR=/scratch/easmit31/variability/dglm/logs
mkdir -p $OUTDIR $LOGDIR

CELL_TYPES=(astrocytes basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons glutamatergic_neurons medium_spiny_neurons microglia midbrain_neurons opc oligodendrocytes vascular_cells)

declare -A MEM_MAP
MEM_MAP[glutamatergic_neurons]=256G
MEM_MAP[GABAergic_neurons]=256G
MEM_MAP[cerebellar_neurons]=256G
MEM_MAP[opc]=128G
MEM_MAP[oligodendrocytes]=128G

for CT in "${CELL_TYPES[@]}"; do
    MEM=${MEM_MAP[$CT]:-64G}
    sbatch \
        --job-name=dglm_${CT} \
        --output=${LOGDIR}/dglm_${CT}_%j.out \
        --error=${LOGDIR}/dglm_${CT}_%j.err \
        --mem=${MEM} --time=03:59:00 -p htc \
        --cpus-per-task=16 \
        --wrap="cd /scratch/easmit31/variability/dglm && Rscript ${SCRIPT_DIR}/dglm_model.R --cell_type ${CT} --outdir ${OUTDIR}"
done

echo "all jobs submitted"
