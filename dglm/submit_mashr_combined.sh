#!/bin/bash
#SBATCH --job-name=mashr_combined
#SBATCH --output=/scratch/easmit31/variability/dglm/logs/mashr_combined_%j.out
#SBATCH --error=/scratch/easmit31/variability/dglm/logs/mashr_combined_%j.err
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH -p htc

/home/easmit31/.conda/envs/mashr_env/bin/Rscript /scratch/easmit31/variability/dglm/scripts/dglm_mashr.R \
    --mode combined \
    --checkpoints /scratch/easmit31/variability/dglm/checkpoints

/home/easmit31/.conda/envs/mashr_env/bin/Rscript /scratch/easmit31/variability/dglm/scripts/dglm_fdr_combined.R \
    --checkpoints /scratch/easmit31/variability/dglm/checkpoints
