#!/bin/bash
#SBATCH --job-name=mashr_combined
#SBATCH --output=/scratch/easmit31/dispersion/dglm/logs/mashr_combined_%j.out
#SBATCH --error=/scratch/easmit31/dispersion/dglm/logs/mashr_combined_%j.err
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH -p htc

# Usage: sbatch [--job-name=mashr_v2] submit_mashr_combined.sh <checkpoints_dir>
# e.g.:  sbatch --job-name=mashr_v2 submit_mashr_combined.sh /scratch/easmit31/dispersion/dglm/disp_age__mean_full
#        sbatch --job-name=mashr_v3 submit_mashr_combined.sh /scratch/easmit31/dispersion/dglm/disp_age_sex__mean_full
#
# Fixed from the previous version of this script, which still pointed at the
# stale /scratch/easmit31/variability/dglm/... paths (confirmed abandoned --
# checkpoints/ there is empty) and hardcoded a single --checkpoints dir,
# so it could only ever target one model version at a time. Parameterizing
# it means this one script now covers all model versions, matching how
# dglm_mashr.R itself was unified earlier this session.

if [ -z "$1" ]; then
    echo "Usage: sbatch submit_mashr_combined.sh <checkpoints_dir>"
    exit 1
fi
CHECKPOINTS=$1

/home/easmit31/.conda/envs/mashr_env/bin/Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_mashr.R \
    --mode combined \
    --checkpoints "$CHECKPOINTS"

/home/easmit31/.conda/envs/mashr_env/bin/Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_fdr_combined.R \
    --checkpoints "$CHECKPOINTS"
