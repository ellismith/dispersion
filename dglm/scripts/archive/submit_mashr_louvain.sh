#!/bin/bash
#SBATCH --job-name=mashr_louvain
#SBATCH --output=/scratch/easmit31/dispersion/dglm/logs/mashr_louvain_%j.out
#SBATCH --error=/scratch/easmit31/dispersion/dglm/logs/mashr_louvain_%j.err
#SBATCH --mem=256G
#SBATCH --time=6:00:00
#SBATCH -p highmem

/home/easmit31/.conda/envs/mashr_env/bin/Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_mashr.R \
    --mode combined \
    --checkpoints /scratch/easmit31/dispersion/dglm/disp_age__mean_full_louvain \
    --fast

/home/easmit31/.conda/envs/mashr_env/bin/Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_fdr_combined.R \
    --checkpoints /scratch/easmit31/dispersion/dglm/disp_age__mean_full_louvain
