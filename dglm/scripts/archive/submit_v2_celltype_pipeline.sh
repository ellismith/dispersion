#!/bin/bash
# Whole-cell-type resolution (no louvain subcluster split), both shat
# modes (raw + sqrt). Uses the SAME current pseudobulk.py/filter pipeline
# as the louvain runs. Writes to a directory entirely separate from the
# louvain disp_age__raw_counts/ tree, so nothing overwrites anything.
# Two full checkpoints copies (one per shat mode) since dglm_mashr.R's
# output filename doesn't encode shat_mode.

SCRIPTS=/scratch/easmit31/dispersion/dglm/scripts
BASE=/scratch/easmit31/dispersion/dglm/disp_age__raw_counts_celltype
LOGDIR=/scratch/easmit31/dispersion/dglm/logs
PYTHON=/home/easmit31/.conda/envs/latent_analysis/bin/python
RSCRIPT=/home/easmit31/.conda/envs/mashr_env/bin/Rscript
CUTOFF=0.5

mkdir -p "$BASE" "$BASE/filtered"

CELL_TYPES=(astrocytes basket_cells cerebellar_neurons ependymal_cells GABAergic_neurons glutamatergic_neurons medium_spiny_neurons microglia midbrain_neurons oligodendrocytes opc vascular_cells)

declare -A MEM_MAP
MEM_MAP[glutamatergic_neurons]=256G
MEM_MAP[GABAergic_neurons]=128G
MEM_MAP[cerebellar_neurons]=128G

PB_JOB_IDS=()
for CT in "${CELL_TYPES[@]}"; do
    MEM=${MEM_MAP[$CT]:-64G}
    jid=$(sbatch --parsable \
        --job-name=v2_pb_${CT} \
        --output=${LOGDIR}/v2_pseudobulk_${CT}_%j.out \
        --error=${LOGDIR}/v2_pseudobulk_${CT}_%j.err \
        --mem=${MEM} --time=03:59:00 -p htc \
        --cpus-per-task=16 \
        --wrap="\
${PYTHON} ${SCRIPTS}/pseudobulk.py --cell_type ${CT} --outdir ${BASE}; \
${SCRIPTS}/filter_all_pseudobulk.sh ${BASE} ${BASE}/filtered")
    PB_JOB_IDS+=("$jid")
    echo "Submitted pseudobulk: $CT -> job $jid"
done
PB_DEP=$(IFS=:; echo "${PB_JOB_IDS[*]}")

STAGE_BASE="$BASE/dglm_run_cutoff${CUTOFF}"
CHECKPOINTS="$BASE/dglm_checkpoints_cutoff${CUTOFF}"

DGLM_JID=$(sbatch --parsable \
    --dependency=afterok:${PB_DEP} \
    --job-name=v2_dglm \
    --output=${LOGDIR}/v2_dglm_%j.out \
    --error=${LOGDIR}/v2_dglm_%j.err \
    --mem=64G --time=03:59:00 -p htc \
    --cpus-per-task=16 \
    --wrap="${SCRIPTS}/run_dglm_for_celltype.sh '${BASE}' '${BASE}/filtered' '${CUTOFF}' '${STAGE_BASE}' '${CHECKPOINTS}'")
echo "Submitted dglm (all 12 cell types): job $DGLM_JID (depends on all pseudobulk jobs)"

CHECKPOINTS_RAW="${CHECKPOINTS}_raw"
CHECKPOINTS_SQRT="${CHECKPOINTS}_sqrt"

MASHR_JID=$(sbatch --parsable \
    --dependency=afterok:${DGLM_JID} \
    --job-name=v2_mashr \
    --output=${LOGDIR}/v2_mashr_%j.out \
    --error=${LOGDIR}/v2_mashr_%j.err \
    --mem=64G --time=03:59:00 -p htc \
    --cpus-per-task=16 \
    --wrap="\
mkdir -p '${CHECKPOINTS_RAW}' '${CHECKPOINTS_SQRT}'; \
cp ${CHECKPOINTS}/*_dglm_results.rds '${CHECKPOINTS_RAW}/'; \
cp ${CHECKPOINTS}/*_dglm_results.rds '${CHECKPOINTS_SQRT}/'; \
${RSCRIPT} ${SCRIPTS}/dglm_mashr.R --mode combined --checkpoints '${CHECKPOINTS_RAW}' --shat_mode raw --fast --random_subset_n 3000 --gridmult 2.5; \
${RSCRIPT} ${SCRIPTS}/dglm_fdr_combined.R --checkpoints '${CHECKPOINTS_RAW}'; \
${RSCRIPT} ${SCRIPTS}/dglm_mashr.R --mode combined --checkpoints '${CHECKPOINTS_SQRT}' --shat_mode sqrt --fast --random_subset_n 3000 --gridmult 2.5; \
${RSCRIPT} ${SCRIPTS}/dglm_fdr_combined.R --checkpoints '${CHECKPOINTS_SQRT}'")
echo "Submitted mashr (raw + sqrt, separate dirs): job $MASHR_JID (depends on dglm)"

echo ""
echo "Chain: ${#CELL_TYPES[@]} pseudobulk jobs -> dglm ($DGLM_JID) -> mashr raw+sqrt ($MASHR_JID)"
