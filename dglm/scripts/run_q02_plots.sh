#!/bin/bash
# =============================================================================
# run_q02_plots.sh — regenerate figures for the q<0.2 strong-subset run.
# Every call passes an explicit --base_dir / --out into figures_915q02/,
# so nothing in figures_915/ (the q<0.05 run) can be overwritten.
# Failures are echoed ("FAILED: ...") and the script moves on.
# =============================================================================
cd /scratch/easmit31/dispersion/dglm
B=$PWD/disp_age__mask1_915_300c_q0.2          # q<0.2 masters: $B/<CT>/dglm_checkpoints_cutoff0.5/
F=figures_915q02
R=/home/easmit31/.conda/envs/mashr_env/bin/Rscript
CTS="GABAergic_neurons astrocytes glutamatergic_neurons oligodendrocytes ependymal_cells medium_spiny_neurons midbrain_neurons microglia opc vascular_cells"
mkdir -p $F/fig3F $F/upset_chiou $F/bar_charts
fail() { echo "FAILED: $*"; }

# 1. old vs new comparison bars (now including every finished type)
$R scripts/plot_compare_strong_subset.R || fail compare

# 2. top-gene scatters without greyed-out masked points
for CT in $CTS; do
  M=$B/$CT/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv; [ -f $M ] || { echo "no master: $CT"; continue; }
  $R scripts/plot_top_genes_expr_vs_age_915.R --master $M --raw_counts_dir disp_age__raw_counts/$CT \
    --n_top 6 --lfsr_thresh 0.05 --hide_masked --title "$CT (q<0.2 strong subset): top genes" \
    --out $F/top12_${CT}_q02_clean.png || fail "scatter_clean $CT"
done

# 3. fig3F histogram: simple and split-direction (60/40 net_dir rule), min-gene thresholds 1/20/50
for MIN in 1 20 50; do
  $R scripts/fig3F_style_hist.R --all --base_dir $B --lfsr_thresh 0.05 --min_genes $MIN \
    --out $F/fig3F/ALL_cell_types_min${MIN}.png || fail "fig3F min$MIN"
  $R scripts/fig3F_style_hist.R --all --split_direction --base_dir $B --lfsr_thresh 0.05 --min_genes $MIN \
    --out $F/fig3F/ALL_cell_types_netdir60_min${MIN}.png || fail "fig3F netdir min$MIN"
done

# 4. Chiou-style upsets per cell type (region + subcluster), plus >=5%-of-tested strict variants
for CT in $CTS; do
  for G in region subcluster; do
    $R scripts/plot_upset_chiou_style.R --ct $CT --groupby $G --base_dir $B --lfsr 0.05 \
      --out $F/upset_chiou/${CT}_${G}.png || fail "upset $CT $G"
    $R scripts/plot_upset_chiou_style.R --ct $CT --groupby $G --base_dir $B --lfsr 0.05 --min_pct 5 \
      --out $F/upset_chiou/${CT}_${G}_strict5pct.png || fail "upset strict $CT $G"
  done
done

# 5. all-cell-types upset (one row per cell type)
$R scripts/plot_upset_allct.R --base_dir $B --lfsr 0.05 --out $F/upset_chiou/ALL_celltypes.png || fail allct
$R scripts/plot_upset_allct.R --base_dir $B --lfsr 0.05 --min_pct 5 --out $F/upset_chiou/ALL_celltypes_strict5pct.png || fail "allct strict"

# 6. bar charts of sig genes by region (fixed filenames -> one subfolder per cell type)
for CT in $CTS; do
  $R scripts/plot_bar_sig_by_region.R --base_dir $B --target_celltype $CT --sig_col mash_lfsr --qthresh 0.05 \
    --figdir $F/bar_charts/$CT || fail "bars $CT"
done

# 7. louvain upset + summary bars: only if they accept --base_dir and --out (otherwise skip, never guess paths)
for S in plot_upset_louvain plot_summary_bars; do
  if grep -q "base_dir" scripts/$S.R && grep -q -- "--out" scripts/$S.R; then
    $R scripts/$S.R --base_dir $B --out $F/${S}.png || fail $S
  else
    echo "SKIPPED: $S (no --base_dir/--out option; needs a manual call)"
  fi
done

echo "ALL PLOTS DONE"; find $F -name "*.png" | wc -l
