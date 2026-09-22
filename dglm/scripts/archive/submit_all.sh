#!/bin/bash
# submit_all.sh
# Submits pseudobulk + DGLM for all three model versions.
# v1: age + sex, per-cell-type mashr
# v2: age + sex, global mashr  
# v3: age + sex + mean_n_umi + n_cells, global mashr

PYSCRIPT=~/.conda/envs/mixed_models/bin/python
RSCRIPT=/home/easmit31/.conda/envs/mashr_env/bin/Rscript
SCRIPTS=/scratch/easmit31/dispersion/dglm/scripts
LOGS=/scratch/easmit31/dispersion/dglm/logs

AGE_SEX=/scratch/easmit31/dispersion/dglm/age_sex
FULL=/scratch/easmit31/dispersion/dglm/age_sex_n_umi_ncells

mkdir -p $AGE_SEX $FULL $LOGS

CTS=(astrocytes basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons
     glutamatergic_neurons medium_spiny_neurons microglia midbrain_neurons
     opc oligodendrocytes vascular_cells)

# ── pseudobulk ────────────────────────────────────────────────────────────
PB_AGE_SEX_JOBS=""
PB_FULL_JOBS=""

echo "=== Submitting pseudobulk ==="
for CT in "${CTS[@]}"; do
    JID=$(sbatch --job-name=pb_as_${CT} \
        --output=${LOGS}/pb_as_${CT}_%j.out \
        --error=${LOGS}/pb_as_${CT}_%j.err \
        --mem=64G --time=02:00:00 -p htc \
        --wrap="$PYSCRIPT $SCRIPTS/pseudobulk.py --cell_type ${CT} --outdir ${AGE_SEX} --covariates age_sex" \
        | awk '{print $4}')
    PB_AGE_SEX_JOBS="${PB_AGE_SEX_JOBS}:${JID}"
    echo "  pb age_sex $CT: $JID"

    JID=$(sbatch --job-name=pb_full_${CT} \
        --output=${LOGS}/pb_full_${CT}_%j.out \
        --error=${LOGS}/pb_full_${CT}_%j.err \
        --mem=64G --time=02:00:00 -p htc \
        --wrap="$PYSCRIPT $SCRIPTS/pseudobulk.py --cell_type ${CT} --outdir ${FULL} --covariates full" \
        | awk '{print $4}')
    PB_FULL_JOBS="${PB_FULL_JOBS}:${JID}"
    echo "  pb full $CT: $JID"
done
PB_AGE_SEX_JOBS=${PB_AGE_SEX_JOBS#:}
PB_FULL_JOBS=${PB_FULL_JOBS#:}

# ── DGLM ─────────────────────────────────────────────────────────────────
DGLM_AGE_SEX_JOBS=""
DGLM_FULL_JOBS=""

echo "=== Submitting DGLM ==="
for CT in "${CTS[@]}"; do
    JID=$(sbatch --job-name=dglm_as_${CT} \
        --output=${LOGS}/dglm_as_${CT}_%j.out \
        --error=${LOGS}/dglm_as_${CT}_%j.err \
        --dependency=afterok:${PB_AGE_SEX_JOBS} \
        --mem=64G --time=3:00:00 -p htc \
        --wrap="$RSCRIPT $SCRIPTS/dglm_model.R --cell_type ${CT} --outdir ${AGE_SEX} --covariates age_sex" \
        | awk '{print $4}')
    DGLM_AGE_SEX_JOBS="${DGLM_AGE_SEX_JOBS}:${JID}"
    echo "  dglm age_sex $CT: $JID"

    JID=$(sbatch --job-name=dglm_full_${CT} \
        --output=${LOGS}/dglm_full_${CT}_%j.out \
        --error=${LOGS}/dglm_full_${CT}_%j.err \
        --dependency=afterok:${PB_FULL_JOBS} \
        --mem=64G --time=3:00:00 -p htc \
        --wrap="$RSCRIPT $SCRIPTS/dglm_model.R --cell_type ${CT} --outdir ${FULL} --covariates full" \
        | awk '{print $4}')
    DGLM_FULL_JOBS="${DGLM_FULL_JOBS}:${JID}"
    echo "  dglm full $CT: $JID"
done
DGLM_AGE_SEX_JOBS=${DGLM_AGE_SEX_JOBS#:}
DGLM_FULL_JOBS=${DGLM_FULL_JOBS#:}

# ── mashr + FDR ───────────────────────────────────────────────────────────
echo "=== Submitting mashr + FDR ==="

# v1: per-cell-type mashr for age_sex
for CT in "${CTS[@]}"; do
    sbatch --job-name=mashr_as_${CT} \
        --output=${LOGS}/mashr_as_${CT}_%j.out \
        --error=${LOGS}/mashr_as_${CT}_%j.err \
        --dependency=afterok:${DGLM_AGE_SEX_JOBS} \
        --mem=64G --time=2:00:00 -p htc \
        --wrap="$RSCRIPT $SCRIPTS/dglm_mashr.R --mode per_ct --cell_type ${CT} --checkpoints ${AGE_SEX}"
    echo "  mashr per_ct age_sex $CT"
done

# v2: global mashr for age_sex
sbatch --job-name=mashr_as_global \
    --output=${LOGS}/mashr_as_global_%j.out \
    --error=${LOGS}/mashr_as_global_%j.err \
    --dependency=afterok:${DGLM_AGE_SEX_JOBS} \
    --mem=128G --time=4:00:00 -p htc \
    --wrap="$RSCRIPT $SCRIPTS/dglm_mashr.R --mode combined --checkpoints ${AGE_SEX} && \
            $RSCRIPT $SCRIPTS/dglm_fdr_combined.R --checkpoints ${AGE_SEX}"
echo "  mashr global age_sex"

# v3: global mashr for full covs
sbatch --job-name=mashr_full_global \
    --output=${LOGS}/mashr_full_global_%j.out \
    --error=${LOGS}/mashr_full_global_%j.err \
    --dependency=afterok:${DGLM_FULL_JOBS} \
    --mem=128G --time=4:00:00 -p htc \
    --wrap="$RSCRIPT $SCRIPTS/dglm_mashr.R --mode combined --checkpoints ${FULL} && \
            $RSCRIPT $SCRIPTS/dglm_fdr_combined.R --checkpoints ${FULL}"
echo "  mashr global full"

echo "=== All jobs submitted ==="
echo "Results will be in:"
echo "  v1 (age+sex, per-ct mashr):    $AGE_SEX"
echo "  v2 (age+sex, global mashr):    $AGE_SEX"
echo "  v3 (full covs, global mashr):  $FULL"
