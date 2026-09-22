#!/bin/bash
SCRIPT=/scratch/easmit31/dispersion/dglm/scripts/dglm_model.R
OUTDIR=/scratch/easmit31/dispersion/dglm/disp_age__mean_full
LOGDIR=/scratch/easmit31/dispersion/dglm/logs
RSCRIPT=/home/easmit31/.conda/envs/mashr_env/bin/Rscript
mkdir -p $LOGDIR
declare -A MEM_MAP
MEM_MAP[glutamatergic_neurons]=256G
MEM_MAP[GABAergic_neurons]=128G
MEM_MAP[cerebellar_neurons]=128G
CELL_TYPES=(astrocytes basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons medium_spiny_neurons microglia midbrain_neurons opc oligodendrocytes vascular_cells)
for CT in "${CELL_TYPES[@]}"; do
    MEM=${MEM_MAP[$CT]:-64G}
    sbatch --job-name=dglm_${CT}_v2 --output=${LOGDIR}/dglm_${CT}_v2_%j.out --error=${LOGDIR}/dglm_${CT}_v2_%j.err \
        --mem=${MEM} --time=03:59:00 -p htc --cpus-per-task=16 \
        --wrap="cd /scratch/easmit31/dispersion/dglm && ${RSCRIPT} ${SCRIPT} --cell_type ${CT} --outdir ${OUTDIR} --covariates full --dispersion age"
done
echo "v2 dglm jobs submitted (12 total)"
