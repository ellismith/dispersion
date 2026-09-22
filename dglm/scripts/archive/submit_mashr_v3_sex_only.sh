#!/bin/bash
#SBATCH --job-name=mashr_v3_sex
#SBATCH --output=/scratch/easmit31/dispersion/dglm/logs/mashr_v3_sex_%j.out
#SBATCH --error=/scratch/easmit31/dispersion/dglm/logs/mashr_v3_sex_%j.err
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH -p htc

/home/easmit31/.conda/envs/mashr_env/bin/Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_mashr.R \
    --mode combined \
    --checkpoints /scratch/easmit31/dispersion/dglm/disp_age_sex__mean_full \
    --only_term sex

/home/easmit31/.conda/envs/mashr_env/bin/Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_fdr_combined.R \
    --checkpoints /scratch/easmit31/dispersion/dglm/disp_age_sex__mean_full
