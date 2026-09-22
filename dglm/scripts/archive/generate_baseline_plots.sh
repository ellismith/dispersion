#!/bin/bash
CT=$1
SCRIPTS=/scratch/easmit31/dispersion/dglm/scripts

if [[ "$CT" == "opc" || "$CT" == "microglia" || "$CT" == "astrocytes" || "$CT" == "vascular_cells" || "$CT" == "oligodendrocytes" ]]; then
  BASE_DIR=/scratch/easmit31/dispersion/dglm/disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered
else
  BASE_DIR=/scratch/easmit31/dispersion/dglm/disp_age__FINAL_autosomeX_300c_QCfiltered
fi

MASTER=${BASE_DIR}/${CT}/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv
if [ ! -f "$MASTER" ]; then echo "MISSING: $MASTER"; exit 1; fi

FIGDIR=${BASE_DIR}/${CT}/figures
mkdir -p $FIGDIR

Rscript $SCRIPTS/plot_gene_count_heatmaps.R --base_dir $BASE_DIR --target_celltype $CT --metric genes_tested --figdir $FIGDIR --per_celltype
Rscript $SCRIPTS/plot_gene_count_heatmaps.R --base_dir $BASE_DIR --target_celltype $CT --metric mash_beta_sig --figdir $FIGDIR --per_celltype

Rscript $SCRIPTS/plot_bar_sig_by_region.R --base_dir $BASE_DIR --target_celltype $CT
mv /scratch/easmit31/dispersion/dglm/figures/bar_sig_by_region_celltype_age_mash_lfsr0.2.png $FIGDIR/bar_sig_by_region_celltype_${CT}_baseline.png 2>/dev/null
mv /scratch/easmit31/dispersion/dglm/figures/bar_sig_by_region_subcluster_age_mash_lfsr0.2.png $FIGDIR/bar_sig_by_region_subcluster_${CT}_baseline.png 2>/dev/null

Rscript $SCRIPTS/plot_fig4a_style.R --base_dir $BASE_DIR --cell_type $CT --count_mode genes --display split --figdir $FIGDIR
Rscript $SCRIPTS/plot_fig4a_style.R --base_dir $BASE_DIR --cell_type $CT --count_mode genes --display net --figdir $FIGDIR

Rscript $SCRIPTS/fig3F_style_hist.R --master $MASTER \
  --title "${CT} (baseline): subclusters by regions showing significant effects" \
  --out $FIGDIR/fig3F_${CT}_baseline.png
Rscript $SCRIPTS/fig3F_style_hist.R --master $MASTER --split_direction \
  --title "${CT} (baseline): subclusters by regions, split by direction" \
  --out $FIGDIR/fig3F_${CT}_baseline_netdir.png

Rscript $SCRIPTS/upset_plot.R --master $MASTER --group_by region \
  --title "${CT} (baseline): UpSet by region" \
  --out $FIGDIR/upset_${CT}_baseline_region.png
Rscript $SCRIPTS/upset_plot.R --master $MASTER --group_by subcluster \
  --title "${CT} (baseline): UpSet by subcluster" \
  --out $FIGDIR/upset_${CT}_baseline_subcluster.png

RAW_DIR=/scratch/easmit31/dispersion/dglm/disp_age__raw_counts/${CT}
Rscript $SCRIPTS/plot_top_genes_expr_vs_age.R \
  --master $MASTER --raw_counts_dir $RAW_DIR \
  --title "${CT} (baseline): top increasing/decreasing dispersion-age genes" \
  --out $FIGDIR/top12_expr_vs_age_${CT}_baseline.png

echo "${CT} baseline: $(ls $FIGDIR | wc -l) plot files"
