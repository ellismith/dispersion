#!/bin/bash

H5AD_DIR=/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update
SCRIPT=/scratch/easmit31/variability/gene_variance/plot_gene_variance.py
INDIR=/scratch/easmit31/variability/gene_variance/results_log
OUTDIR=/scratch/easmit31/variability/gene_variance/plots_log
PYTHON=~/.conda/envs/mixed_models/bin/python
LOGDIR=/scratch/easmit31/variability/gene_variance/logs
REGIONS=(ACC CN dlPFC EC HIP IPP lCb M1 MB mdTN NAc)

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

declare -A MEM_MAP
MEM_MAP[glutamatergic_neurons]=256G
MEM_MAP[GABAergic_neurons]=256G
MEM_MAP[cerebellar_neurons]=256G

for CT in "${!H5AD_MAP[@]}"; do
    H5AD=${H5AD_DIR}/${H5AD_MAP[$CT]}
    MEM=${MEM_MAP[$CT]:-128G}

    # one job per cell type: runs all regions + all_regions plot
    sbatch \
        --job-name=plot_${CT} \
        --output=${LOGDIR}/plot_${CT}_%j.out \
        --error=${LOGDIR}/plot_${CT}_%j.err \
        --mem=${MEM} --time=03:59:00 -p htc \
        --wrap="${PYTHON} ${SCRIPT} \
            --h5ad ${H5AD} --indir ${INDIR} --outdir ${OUTDIR} \
            --cell_type ${CT} --all_regions"

    # per region jobs
    for REGION in "${REGIONS[@]}"; do
        sbatch \
            --job-name=plot_${CT}_${REGION} \
            --output=${LOGDIR}/plot_${CT}_${REGION}_%j.out \
            --error=${LOGDIR}/plot_${CT}_${REGION}_%j.err \
            --mem=${MEM} --time=02:00:00 -p htc \
            --wrap="${PYTHON} ${SCRIPT} \
                --h5ad ${H5AD} --indir ${INDIR} --outdir ${OUTDIR} \
                --cell_type ${CT} --region ${REGION}"
    done
done

echo "all jobs submitted"
