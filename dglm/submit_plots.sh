#!/bin/bash

SCRIPT_DIR=/scratch/easmit31/variability/dglm/scripts
CHECKPOINTS=/scratch/easmit31/variability/dglm/checkpoints
FIGDIR=/scratch/easmit31/variability/dglm/figures
LOGDIR=/scratch/easmit31/variability/dglm/logs
RSCRIPT=/home/easmit31/.conda/envs/mashr_env/bin/Rscript
mkdir -p $LOGDIR

CELL_TYPES=(astrocytes basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons glutamatergic_neurons medium_spiny_neurons microglia midbrain_neurons opc oligodendrocytes vascular_cells)

for CT in "${CELL_TYPES[@]}"; do
    sbatch \
        --job-name=plt_${CT} \
        --output=${LOGDIR}/plt_${CT}_%j.out \
        --error=${LOGDIR}/plt_${CT}_%j.err \
        --mem=32G --time=01:00:00 -p htc \
        --wrap="${RSCRIPT} ${SCRIPT_DIR}/dglm_plot_pval.R \
            --cell_type ${CT} \
            --checkpoints ${CHECKPOINTS} \
            --figdir ${FIGDIR}"
done

echo "all jobs submitted"

for CT in "${CELL_TYPES[@]}"; do
    sbatch \
        --job-name=pltb_${CT} \
        --output=${LOGDIR}/pltb_${CT}_%j.out \
        --error=${LOGDIR}/pltb_${CT}_%j.err \
        --mem=32G --time=01:00:00 -p htc \
        --wrap="${RSCRIPT} ${SCRIPT_DIR}/dglm_plot_beta.R \
            --cell_type ${CT} \
            --checkpoints ${CHECKPOINTS} \
            --figdir ${FIGDIR}"
done
