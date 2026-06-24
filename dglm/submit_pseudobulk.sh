#!/bin/bash

SCRIPT=/scratch/easmit31/variability/dglm/scripts/pseudobulk.py
OUTDIR=/scratch/easmit31/variability/dglm/checkpoints
LOGDIR=/scratch/easmit31/variability/dglm/logs
PYTHON=~/.conda/envs/mixed_models/bin/python
mkdir -p $OUTDIR $LOGDIR

declare -A MEM_MAP
MEM_MAP[glutamatergic_neurons]=256G
MEM_MAP[GABAergic_neurons]=128G
MEM_MAP[cerebellar_neurons]=128G

CELL_TYPES=(astrocytes basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons glutamatergic_neurons medium_spiny_neurons microglia midbrain_neurons opc oligodendrocytes vascular_cells)

for CT in "${CELL_TYPES[@]}"; do
    MEM=${MEM_MAP[$CT]:-64G}
    sbatch \
        --job-name=pb_${CT} \
        --output=${LOGDIR}/pb_${CT}_%j.out \
        --error=${LOGDIR}/pb_${CT}_%j.err \
        --mem=${MEM} --time=02:00:00 -p htc \
        --wrap="${PYTHON} ${SCRIPT} --cell_type ${CT} --outdir ${OUTDIR}"
done

echo "all jobs submitted"
