#!/bin/bash
SCRIPT=/scratch/easmit31/dispersion/dglm/scripts/pseudobulk.py
OUTDIR=/scratch/easmit31/dispersion/dglm/disp_age__mean_full_louvain
LOGDIR=/scratch/easmit31/dispersion/dglm/logs
PYTHON=~/.conda/envs/mixed_models/bin/python
mkdir -p $OUTDIR $LOGDIR
declare -A MEM_MAP
MEM_MAP[glutamatergic_neurons]=256G
MEM_MAP[GABAergic_neurons]=128G
MEM_MAP[cerebellar_neurons]=128G
CELL_TYPES=(astrocytes basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons glutamatergic_neurons medium_spiny_neurons microglia midbrain_neurons opc oligodendrocytes vascular_cells)
for CT in "${CELL_TYPES[@]}"; do
    MEM=${MEM_MAP[$CT]:-128G}
    PART=htc
    TIME=03:59:00
    if [ "$CT" == "glutamatergic_neurons" ]; then PART=highmem; TIME=05:00:00; fi
    sbatch \
        --job-name=pb_${CT}_louvain \
        --output=${LOGDIR}/pb_${CT}_louvain_%j.out \
        --error=${LOGDIR}/pb_${CT}_louvain_%j.err \
        --mem=${MEM} --time=${TIME} -p ${PART} \
        --wrap="${PYTHON} ${SCRIPT} --cell_type ${CT} --outdir ${OUTDIR} --covariates full --min_cells 100 --subcluster_col ct_louvain"
done
echo "louvain pseudobulk jobs submitted (12 total)"
