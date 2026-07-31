#!/bin/bash
#SBATCH --job-name=gene_var_mglia_hip
#SBATCH --output=/scratch/easmit31/variability/gene_variance/logs/mglia_hip_%j.out
#SBATCH --error=/scratch/easmit31/variability/gene_variance/logs/mglia_hip_%j.err
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH -p htc

mkdir -p /scratch/easmit31/variability/gene_variance/logs

conda activate mixed_models

python /scratch/easmit31/variability/gene_variance/run_gene_variance.py \
    --h5ad /scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update/Res1_microglia_new.h5ad \
    --cell_type microglia \
    --region HIP \
    --outdir /scratch/easmit31/variability/gene_variance \
    --min_age 1.0 \
    --min_animals 10
