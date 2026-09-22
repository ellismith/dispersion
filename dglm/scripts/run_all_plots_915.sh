#!/bin/bash
# Generate all plots for a completed pipeline, organized by figure type.
# Usage: run_all_plots_for_pipeline.sh <pipeline_base_dir> <raw_counts_base_dir> <figures_suffix>
# Example: run_all_plots_for_pipeline.sh disp_age__outlierfilter_300c_QCfiltered disp_age__raw_counts no_outliers

set -e
PIPELINE_BASE=$1
RAW_COUNTS_BASE=$2
FIG_SUFFIX=${3:-pipeline}

SCRIPTS=/scratch/easmit31/dispersion/dglm/scripts
RSCRIPT=/home/easmit31/.conda/envs/mashr_env/bin/Rscript
DGLM_BASE=/scratch/easmit31/dispersion/dglm
LCB="opc microglia astrocytes vascular_cells oligodendrocytes"
CELL_TYPES="GABAergic_neurons microglia astrocytes glutamatergic_neurons cerebellar_neurons medium_spiny_neurons midbrain_neurons basket_cells vascular_cells ependymal_cells opc oligodendrocytes"

# figure output dirs
FIG_ROOT=${DGLM_BASE}/figures_${FIG_SUFFIX}
mkdir -p \
  ${FIG_ROOT}/heatmaps \
  ${FIG_ROOT}/fig3F \
  ${FIG_ROOT}/fig4a \
  ${FIG_ROOT}/bar_charts \
  ${FIG_ROOT}/upset \
  ${FIG_ROOT}/upset_chiou \
  ${FIG_ROOT}/upset_flexible \
  ${FIG_ROOT}/condition_corr \
  ${FIG_ROOT}/top60_scatter \
  ${FIG_ROOT}/prop_increasing

echo "=== Generating plots for pipeline: $PIPELINE_BASE ==="
echo "=== Output root: $FIG_ROOT ==="

for CT in $CELL_TYPES; do
  echo "$LCB" | grep -qw "$CT" && \
    BASE_DIR=${DGLM_BASE}/${PIPELINE_BASE/_QCfiltered/_noLcb_QCfiltered} || \
    BASE_DIR=${DGLM_BASE}/${PIPELINE_BASE}
  MASTER=${BASE_DIR}/${CT}/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv
  RAW_DIR=${DGLM_BASE}/${RAW_COUNTS_BASE}/${CT}
  [ ! -f "$MASTER" ] && echo "SKIP $CT: no master table" && continue
  echo "--- $CT ---"

  $RSCRIPT $SCRIPTS/plot_gene_count_heatmaps.R \
    --base_dir $BASE_DIR --target_celltype $CT --metric genes_tested \
    --figdir ${FIG_ROOT}/heatmaps --per_celltype
  $RSCRIPT $SCRIPTS/plot_gene_count_heatmaps.R \
    --base_dir $BASE_DIR --target_celltype $CT --metric mash_beta_sig \
    --figdir ${FIG_ROOT}/heatmaps --per_celltype

  $RSCRIPT $SCRIPTS/plot_fig4a_style.R \
    --base_dir $BASE_DIR --cell_type $CT --count_mode genes --display split \
    --figdir ${FIG_ROOT}/fig4a
  $RSCRIPT $SCRIPTS/plot_fig4a_style.R \
    --base_dir $BASE_DIR --cell_type $CT --count_mode genes --display net \
    --figdir ${FIG_ROOT}/fig4a

  $RSCRIPT $SCRIPTS/plot_bar_sig_by_region.R \
    --base_dir $BASE_DIR --target_celltype $CT --qthresh 0.05
  for TH in 0.2 0.05; do
    mv ${DGLM_BASE}/figures/bar_sig_by_region_celltype_age_mash_lfsr${TH}.png \
       ${FIG_ROOT}/bar_charts/${CT}_celltype.png 2>/dev/null || true
    mv ${DGLM_BASE}/figures/bar_sig_by_region_subcluster_age_mash_lfsr${TH}.png \
       ${FIG_ROOT}/bar_charts/${CT}_subcluster.png 2>/dev/null || true
  done

  $RSCRIPT $SCRIPTS/fig3F_style_hist.R --master $MASTER \
    --lfsr_thresh 0.05 --title "${CT}" \
    --out ${FIG_ROOT}/fig3F/${CT}.png
  $RSCRIPT $SCRIPTS/fig3F_style_hist.R --master $MASTER \
    --lfsr_thresh 0.05 --split_direction --title "${CT} net direction" \
    --out ${FIG_ROOT}/fig3F/${CT}_netdir.png

  $RSCRIPT $SCRIPTS/upset_plot.R --master $MASTER --group_by region \
    --title "${CT}: UpSet by region" \
    --out ${FIG_ROOT}/upset/${CT}_region.png
  $RSCRIPT $SCRIPTS/upset_plot.R --master $MASTER --group_by subcluster \
    --title "${CT}: UpSet by subcluster" \
    --out ${FIG_ROOT}/upset/${CT}_subcluster.png

  $RSCRIPT $SCRIPTS/plot_upset_chiou_style.R \
    --ct $CT --groupby region --base_dir $BASE_DIR \
    --out ${FIG_ROOT}/upset_chiou/${CT}_region.png
  $RSCRIPT $SCRIPTS/plot_upset_chiou_style.R \
    --ct $CT --groupby subcluster --base_dir $BASE_DIR \
    --out ${FIG_ROOT}/upset_chiou/${CT}_subcluster.png

  $RSCRIPT $SCRIPTS/plot_condition_correlation.R \
    --master $MASTER --title "${CT}: condition correlation (top 30)" \
    --subset_by top_n --top_n 30 \
    --out ${FIG_ROOT}/condition_corr/${CT}_top30.png
  $RSCRIPT $SCRIPTS/plot_condition_correlation.R \
    --master $MASTER --title "${CT}: condition correlation (increasing)" \
    --subset_by direction --direction inc \
    --out ${FIG_ROOT}/condition_corr/${CT}_increasing.png
  $RSCRIPT $SCRIPTS/plot_condition_correlation.R \
    --master $MASTER --title "${CT}: condition correlation (decreasing)" \
    --subset_by direction --direction dec \
    --out ${FIG_ROOT}/condition_corr/${CT}_decreasing.png

  $RSCRIPT $SCRIPTS/plot_top_genes_expr_vs_age_915.R \
    --master $MASTER --raw_counts_dir $RAW_DIR \
    --n_top 30 --lfsr_thresh 0.05 \
    --title "${CT}: top 30 increasing/decreasing genes" \
    --out ${FIG_ROOT}/top60_scatter/${CT}.png

  echo "$CT: done"
done

echo ""
echo "=== Summary figures ==="
$RSCRIPT $SCRIPTS/fig3F_style_hist.R --all --base_dir ${DGLM_BASE}/${PIPELINE_BASE} --lfsr_thresh 0.05 \
  --title "All cell types: subclusters by regions" \
  --out ${FIG_ROOT}/fig3F/ALL_cell_types.png
$RSCRIPT $SCRIPTS/plot_upset_chiou_style.R --ct all --groupby region \
  --base_dir ${DGLM_BASE}/${PIPELINE_BASE} \
  --out ${FIG_ROOT}/upset_chiou/ALL_region.png
$RSCRIPT $SCRIPTS/plot_upset_chiou_style.R --ct all --groupby subcluster \
  --base_dir ${DGLM_BASE}/${PIPELINE_BASE} \
  --out ${FIG_ROOT}/upset_chiou/ALL_subcluster.png
echo "DONE: $FIG_ROOT"
