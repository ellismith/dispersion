#!/bin/bash
CT=$1
MC=$2
SCRIPTS=/scratch/easmit31/dispersion/dglm/scripts
RSCRIPT=/home/easmit31/.conda/envs/mashr_env/bin/Rscript
LCB_TYPES="opc microglia astrocytes vascular_cells oligodendrocytes"

if echo "$LCB_TYPES" | grep -qw "$CT"; then
  BASE_DIR=/scratch/easmit31/dispersion/dglm/disp_age__mincells${MC}_300c_noLcb_QCfiltered
else
  BASE_DIR=/scratch/easmit31/dispersion/dglm/disp_age__mincells${MC}_300c_QCfiltered
fi

MASTER=${BASE_DIR}/${CT}/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv
RAW_DIR=/scratch/easmit31/dispersion/dglm/disp_age__mincells${MC}/${CT}
FIGDIR=${BASE_DIR}/${CT}/figures
mkdir -p $FIGDIR

if [ ! -f "$MASTER" ]; then echo "MISSING master: $MASTER"; exit 1; fi

echo "=== heatmaps ==="
$RSCRIPT $SCRIPTS/plot_gene_count_heatmaps.R --base_dir $BASE_DIR --target_celltype $CT --metric genes_tested --figdir $FIGDIR --per_celltype
$RSCRIPT $SCRIPTS/plot_gene_count_heatmaps.R --base_dir $BASE_DIR --target_celltype $CT --metric mash_beta_sig --figdir $FIGDIR --per_celltype

echo "=== bar charts ==="
$RSCRIPT $SCRIPTS/plot_bar_sig_by_region.R --base_dir $BASE_DIR --target_celltype $CT
mv /scratch/easmit31/dispersion/dglm/figures/bar_sig_by_region_celltype_age_mash_lfsr0.2.png $FIGDIR/bar_sig_by_region_celltype_${CT}_mc${MC}.png 2>/dev/null
mv /scratch/easmit31/dispersion/dglm/figures/bar_sig_by_region_subcluster_age_mash_lfsr0.2.png $FIGDIR/bar_sig_by_region_subcluster_${CT}_mc${MC}.png 2>/dev/null

echo "=== fig4a ==="
$RSCRIPT $SCRIPTS/plot_fig4a_style.R --base_dir $BASE_DIR --cell_type $CT --count_mode genes --display split --figdir $FIGDIR
$RSCRIPT $SCRIPTS/plot_fig4a_style.R --base_dir $BASE_DIR --cell_type $CT --count_mode genes --display net --figdir $FIGDIR

echo "=== fig3F ==="
$RSCRIPT $SCRIPTS/fig3F_style_hist.R --master $MASTER --title "${CT} (mc${MC})" --out $FIGDIR/fig3F_${CT}_mc${MC}.png
$RSCRIPT $SCRIPTS/fig3F_style_hist.R --master $MASTER --split_direction --title "${CT} (mc${MC}) net direction" --out $FIGDIR/fig3F_${CT}_mc${MC}_netdir.png

echo "=== UpSet ==="
$RSCRIPT $SCRIPTS/upset_plot.R --master $MASTER --group_by region --title "${CT} (mc${MC}): UpSet by region" --out $FIGDIR/upset_${CT}_mc${MC}_region.png
$RSCRIPT $SCRIPTS/upset_plot.R --master $MASTER --group_by subcluster --title "${CT} (mc${MC}): UpSet by subcluster" --out $FIGDIR/upset_${CT}_mc${MC}_subcluster.png

echo "=== top60 scatter ==="
$RSCRIPT $SCRIPTS/plot_top_genes_expr_vs_age.R \
  --master $MASTER --raw_counts_dir $RAW_DIR --n_top 30 --lfsr_thresh 0.05 \
  --title "${CT} (mc${MC}): top 30 increasing/decreasing dispersion-age genes" \
  --out $FIGDIR/top60_expr_vs_age_${CT}_mc${MC}.png

echo "=== condition correlation heatmap ==="
$RSCRIPT $SCRIPTS/plot_condition_correlation.R \
  --master $MASTER --title "${CT} (mc${MC}): condition correlation" \
  --out $FIGDIR/condition_corr_${CT}_mc${MC}.png

echo "DONE: ${CT} mc${MC} -- $(ls $FIGDIR | wc -l) files"
