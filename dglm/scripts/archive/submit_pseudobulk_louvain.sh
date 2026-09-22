#!/bin/bash
SCRIPT=/scratch/easmit31/dispersion/dglm/scripts/pseudobulk.py
BASEDIR=/scratch/easmit31/dispersion/dglm/disp_age__mean_full_louvain
LOGDIR=/scratch/easmit31/dispersion/dglm/logs
PYTHON=~/.conda/envs/mixed_models/bin/python
mkdir -p $LOGDIR

declare -A MEM_MAP
MEM_MAP[glutamatergic_neurons]=256G
MEM_MAP[GABAergic_neurons]=128G
MEM_MAP[cerebellar_neurons]=128G

# microglia already done -- excluded from this list
CELL_TYPES=(astrocytes basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons glutamatergic_neurons medium_spiny_neurons midbrain_neurons opc oligodendrocytes vascular_cells)

for CT in "${CELL_TYPES[@]}"; do
    MEM=${MEM_MAP[$CT]:-64G}
    OUTDIR=${BASEDIR}/${CT}
    mkdir -p $OUTDIR
    sbatch \
        --job-name=pb_${CT}_louvain \
        --output=${LOGDIR}/pb_${CT}_louvain_%j.out \
        --error=${LOGDIR}/pb_${CT}_louvain_%j.err \
        --mem=${MEM} --time=02:00:00 -p htc \
        --wrap="${PYTHON} ${SCRIPT} --cell_type ${CT} --outdir ${OUTDIR} --subcluster_col ct_louvain --covariates full"
done
echo "all pseudobulk jobs submitted"
