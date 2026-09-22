# Transcriptional Dispersion Aging Analysis — Final Summary

## What this analysis does
Tests whether gene expression variability across animals changes with age,
for each brain cell subtype x brain region combination. Increasing dispersion
= more animal-to-animal variability in old animals. Decreasing = less.

## Final output files (use these)

### Per-cell-type master tables (DGLM + mashr combined)
One row per gene x subcluster x region. Contains mash_beta (effect size),
mash_lfsr (local false sign rate), ensembl_id, symbol, cell_type, region.

**Significance threshold: mash_lfsr < 0.05 (~5% false sign rate)**

| Cell type | Path | n_unique_sig_genes | pct_sig |
|---|---|---|---|
| GABAergic_neurons | disp_age__FINAL_autosomeX_300c_QCfiltered/GABAergic_neurons/.../master_dglm_combined.tsv | 8400 | 9.7% |
| microglia | disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered/microglia/.../master_dglm_combined.tsv | 5474 | 8.2% |
| astrocytes | disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered/astrocytes/.../master_dglm_combined.tsv | 9259 | 10.2% |
| glutamatergic_neurons | disp_age__FINAL_autosomeX_300c_QCfiltered/glutamatergic_neurons/.../master_dglm_combined.tsv | 12755 | 16.6% |
| cerebellar_neurons | disp_age__FINAL_autosomeX_300c_QCfiltered/cerebellar_neurons/.../master_dglm_combined.tsv | 12054 | 36.1%* |
| medium_spiny_neurons | disp_age__FINAL_autosomeX_300c_QCfiltered/medium_spiny_neurons/.../master_dglm_combined.tsv | 7094 | 6.5% |
| midbrain_neurons | disp_age__FINAL_autosomeX_300c_QCfiltered/midbrain_neurons/.../master_dglm_combined.tsv | 10782 | 22.1% |
| basket_cells | disp_age__FINAL_autosomeX_300c_QCfiltered/basket_cells/.../master_dglm_combined.tsv | 10839 | 36.3%* |
| vascular_cells | disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered/vascular_cells/.../master_dglm_combined.tsv | 4634 | 5.2% |
| ependymal_cells | disp_age__FINAL_autosomeX_300c_QCfiltered/ependymal_cells/.../master_dglm_combined.tsv | 2632 | 7.6% |
| opc | disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered/opc/.../master_dglm_combined.tsv | 3699 | 5.3% |
| oligodendrocytes | disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered/oligodendrocytes/.../master_dglm_combined.tsv | 10942 | 27.0% |

*Cerebellar and basket_cells: lCb is their only surviving region — calibration caveat applies.

All paths under: /scratch/easmit31/dispersion/dglm/

### Analysis output tables
- Ranked conditions (top subtype x region by n sig genes):
  /scratch/easmit31/dispersion/dglm/ranked_conditions_baseline_lfsr05.csv

- Dispersion x CCC overlap table:
  /scratch/easmit31/dispersion/dglm/dispersion_ccc_overlap.csv

### Key figures
Per-cell-type figures in: <pipeline_dir>/<CT>/figures/
- fig3F_*.png — subclusters by number of regions showing effects
- upset_*_region.png — UpSet by region combination
- top60_expr_vs_age_*.png — top 30 increasing + 30 decreasing genes scatter
- condition_corr_*.png — condition x condition correlation heatmap

Summary figures in: /scratch/easmit31/dispersion/dglm/figures/
- fig3F_ALL_baseline_lfsr05_CURRENT.png — all 12 cell types combined
- fig3F_ALL_baseline_mingenes_comparison.png — min 1/10/100 genes comparison
- upset_*_baseline_twodots.png — UpSet with separate increasing/decreasing dots
- dispersion_ccc_scatter.png — dispersion vs CCC overlap scatter

## Pipeline (brief)
1. Pseudobulk: sum raw counts per animal x subcluster x region
2. Gene filter: CPM >= 0.5 in >= 50% animals AND filterByExpr()
3. Normalize: TMM + log2(CPM + 0.5)
4. DGLM: mean ~ age + sex + mean_n_umi + n_cells; dispersion ~ age
5. Condition filter: >= 300 cells AND >= 30 animals per combo
6. mashr: shrinkage across conditions (qval_single approach)
7. QC: |beta| > 10 guardrail

## Caveats
- lCb excluded for: opc, microglia, astrocytes, vascular_cells, oligodendrocytes
- Cerebellar and basket_cells: lCb is 100% of their data (calibration caveats)
- Ependymal cells: 100% of significant results from sparse conditions
- Scripts: /scratch/easmit31/dispersion/dglm/scripts/
