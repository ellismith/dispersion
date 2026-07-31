#!/bin/bash

H5AD_DIR=/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update
SCRIPT=/scratch/easmit31/variability/gene_variance/run_gene_variance.py
OUTDIR=/scratch/easmit31/variability/gene_variance/results_log
LOGDIR=/scratch/easmit31/variability/gene_variance/logs
PYTHON=~/.conda/envs/mixed_models/bin/python
REGIONS=(ACC CN dlPFC EC HIP IPP lCb M1 MB mdTN NAc)

mkdir -p $OUTDIR $LOGDIR

declare -A H5AD_MAP
H5AD_MAP[astrocytes]=Res1_astrocytes_update.h5ad
H5AD_MAP[basket_cells]=Res1_basket-cells_update.h5ad
H5AD_MAP[cerebellar_neurons]=Res1_cerebellar-neurons_subset.h5ad
H5AD_MAP[ependymal_cells]=Res1_ependymal-cells_new.h5ad
H5AD_MAP[GABAergic_neurons]=Res1_GABAergic-neurons_subset.h5ad
H5AD_MAP[glutamatergic_neurons]=Res1_glutamatergic-neurons_update.h5ad
H5AD_MAP[medium_spiny_neurons]=Res1_medium-spiny-neurons_subset.h5ad
H5AD_MAP[microglia]=Res1_microglia_new.h5ad
H5AD_MAP[midbrain_neurons]=Res1_midbrain-neurons_update.h5ad
H5AD_MAP[opc]=Res1_opc-olig_subset.h5ad
H5AD_MAP[oligodendrocytes]=Res1_opc-olig_subset.h5ad
H5AD_MAP[vascular_cells]=Res1_vascular-cells_subset.h5ad

HTC_TYPES=(astrocytes basket_cells cerebellar_neurons ependymal_cells medium_spiny_neurons microglia midbrain_neurons opc oligodendrocytes vascular_cells)
HIGHMEM_TYPES=(GABAergic_neurons glutamatergic_neurons)

for CT in "${HTC_TYPES[@]}"; do
    for REGION in "${REGIONS[@]}"; do
        sbatch \
            --job-name=gvl_${CT}_${REGION} \
            --output=${LOGDIR}/gvl_${CT}_${REGION}_%j.out \
            --error=${LOGDIR}/gvl_${CT}_${REGION}_%j.err \
            --mem=128G --time=03:59:00 -p htc \
            --wrap="${PYTHON} ${SCRIPT} \
                --h5ad ${H5AD_DIR}/${H5AD_MAP[$CT]} \
                --cell_type ${CT} --region ${REGION} \
                --outdir ${OUTDIR} --min_age 1.0 --min_animals 10"
    done
done

for CT in "${HIGHMEM_TYPES[@]}"; do
    for REGION in "${REGIONS[@]}"; do
        sbatch \
            --job-name=gvl_${CT}_${REGION} \
            --output=${LOGDIR}/gvl_${CT}_${REGION}_%j.out \
            --error=${LOGDIR}/gvl_${CT}_${REGION}_%j.err \
            --mem=256G --time=03:59:00 -p htc \
            --wrap="${PYTHON} ${SCRIPT} \
                --h5ad ${H5AD_DIR}/${H5AD_MAP[$CT]} \
                --cell_type ${CT} --region ${REGION} \
                --outdir ${OUTDIR} --min_age 1.0 --min_animals 10"
    done
done

echo "all jobs submitted"
