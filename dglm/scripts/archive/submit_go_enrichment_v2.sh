#!/bin/bash
#SBATCH --job-name=go_v2
#SBATCH --output=/scratch/easmit31/dispersion/dglm/logs/go_v2_%j.out
#SBATCH --error=/scratch/easmit31/dispersion/dglm/logs/go_v2_%j.err
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH -p htc

/home/easmit31/.conda/envs/go_env/bin/Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_go_enrichment.R \
    --checkpoints /scratch/easmit31/dispersion/dglm/disp_age__mean_full
