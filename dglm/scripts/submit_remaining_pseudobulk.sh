#!/bin/bash
# Submits one Slurm job per remaining (not-yet-pseudobulked) cell type:
# pseudobulk.py (all louvain subclusters x all regions) then the cutoff
# 0.5/1/2 CPM filter, chained in one job so each cell type is fully ready
# for the dglm step once its job finishes.

SCRIPTS=/scratch/easmit31/dispersion/dglm/scripts
BASE=/scratch/easmit31/dispersion/dglm/disp_age__raw_counts
LOGDIR=/scratch/easmit31/dispersion/dglm/logs
PYTHON=/home/easmit31/.conda/envs/latent_analysis/bin/python

declare -A MEM_MAP
MEM_MAP[glutamatergic_neurons]=256G
MEM_MAP[midbrain_neurons]=128G

CELL_TYPES=(glutamatergic_neurons medium_spiny_neurons midbrain_neurons oligodendrocytes opc vascular_cells)

for CT in "${CELL_TYPES[@]}"; do
    MEM=${MEM_MAP[$CT]:-64G}
    mkdir -p "$BASE/$CT" "$BASE/$CT/filtered"

    sbatch \
        --job-name=pseudobulk_${CT} \
        --output=${LOGDIR}/pseudobulk_${CT}_%j.out \
        --error=${LOGDIR}/pseudobulk_${CT}_%j.err \
        --mem=${MEM} --time=04:00:00 -p htc \
        --cpus-per-task=16 \
        --wrap="\
${PYTHON} ${SCRIPTS}/pseudobulk.py --cell_type ${CT} --outdir ${BASE}/${CT} --subcluster_col ct_louvain; \
${SCRIPTS}/filter_all_pseudobulk.sh ${BASE}/${CT} ${BASE}/${CT}/filtered"
done

echo "Submitted ${#CELL_TYPES[@]} jobs."
