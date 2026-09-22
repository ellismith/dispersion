# dispersion

Analysis pipelines for age-associated changes in transcriptional dispersion
(variability) in the rhesus macaque brain aging snRNA-seq dataset (U01).

- **`dglm/`** — DGLM dispersion pipeline adapted from Chiou et al. 2022
  (*Nat Neurosci*), with mashr for effect sharing across cell types and
  regions. Current, final pipeline described below.
- **`gene_variance/`** — OLS-based between-/within-individual variance
  regressions. See [`gene_variance/AGENTS.md`](gene_variance/AGENTS.md).
- **`plotting_scripts/`** — shared style (`plot_style.py`: cell-type
  abbreviations, region order, font sizes) and cross-method comparison plots.

Results and figures live on ASU Sol under `/scratch/easmit31/dispersion/dglm/`
and are not tracked in git.

---

## `dglm/` — current pipeline ("mask1", Sept 2026)

Guiding rule: **quality control drops points, never animals or whole
conditions** (the only condition-level filter is the 300-cell / 30-animal rule).

| Step | Script(s) | What it does |
|---|---|---|
| 1. Pseudobulk | `pseudobulk.py` | **Summed** raw counts per animal × cell type/subcluster × region. No cell or gene filters here. Stale outputs for a label are deleted before writing. |
| 2. Gene filter | `filter_percent_animals.R` (loop: `filter_all_pseudobulk.sh`) | Keep genes with CPM ≥ 0.5 in ≥ 50% of animals; writes an audit trail. |
| 3. DGLM | `run_dglm_for_celltype.sh` → `dglm_model.R` | Per subcluster × region, per gene — see below. |
| 4. Condition filter | `assess_condition_sparsity.py`, `filter_conditions_for_mashr.R` | Keep conditions with ≥ 300 cells and ≥ 30 animals. lCb excluded for astrocytes, microglia, oligodendrocytes, opc, vascular_cells. |
| 5. mashr | `dglm_mashr.R` | Combined mode (all cell type × region conditions pooled). `--shat_mode raw` (Shat = bvar). Strong subset: `qval_single` rule with `--strong_q 0.2` (DGLM q < 0.2 in ≥ 1 condition; one threshold for all cell types). Default `--strong_q 0.05` reproduces the earlier run. |
| 6. FDR / master | `dglm_fdr_combined.R` | Writes `master_dglm_combined.tsv`. Significance = `mash_lfsr < 0.05`. |
| 7. Plots / audit | `run_q02_plots.sh`, `run_all_plots_915.sh`, `funnel_audit.sh [cell_type]` | Figure suite; read-only per-step funnel audit. |

### Step 3 detail (`dglm_model.R`)

Within each subcluster × region, per gene:

1. **Low-point mask** (set to NA before fitting): raw count ≤ 1, or value
   below median − max(3 × MAD, 2). Removes count-0/1 "shelves" that come from
   animals with too few cells, without dropping the animal.
2. **Testability:** gene must be truly detected (count ≥ 2) in ≥ 50% of that
   region's animals, and ≥ 10 animals must remain. The CPM ≥ 0.5 filter alone
   is nearly vacuous at these library sizes (a single read already exceeds
   0.5 CPM), so this real-detection rule is what enforces expression.
3. TMM normalization + log-CPM (`edgeR`, `prior.count = 0.5`).
4. `dglm(family = gaussian, dlink = "log")`; mean ~ covariates, dispersion ~ age.
5. Dispersion p-values from `summary(res$dispersion.fit)`, **not** dglm's own
   `dispersion.summary` (which fixes the Gamma scale at 2 and is
   anti-conservative). BH q-values per condition.

The output field `bvar` is a **standard error**, not a variance (name
inherited from Chiou's original script).

### Result directories (Sol, `dglm/`)

| Directory | Contents |
|---|---|
| `disp_age__raw_counts/` | Pseudobulk + CPM-filtered inputs |
| `disp_age__mask1_915/` | DGLM fits, final model |
| `disp_age__mask1_915_300c/` | Condition-filtered; mashr + FDR with strong q < 0.05 |
| `disp_age__mask1_915_300c_q0.2/` | **Current:** mashr + FDR with strong q < 0.2 |
| `figures_915/`, `figures_915q02/` | Figures for the two trees above |
| `disp_age__FINAL_autosomeX*` | Pre-mask baseline, kept for comparison (copy in `BACKUP_baseline_20260914/`) |
| `archive_pre_normalization_fix/`, `archive_pre_mask1/` | Superseded runs, kept for reference |

Coverage: all 12 parent cell types have DGLM fits; mashr/FDR has been run for
10 (basket_cells and cerebellar_neurons pending).

### Validation

- **Calibration** (after the extraction fix; GABAergic_0 × ACC, 10,170 genes):
  5.92% of genes p < 0.05 for real age vs 6.07% for permuted age.
- **Mask:** checked by hand refits. For example, NAP1L1 in GABAergic_16 × dlPFC
  went from β −0.42 to +0.05 once two count-0 animals were masked, and a saved
  ADAMTS6 β (−0.724) was reproduced by refit (−0.719).
- **Dispersion decreases are not a mask artifact.** In astrocytes_6 × ACC,
  masked and kept animals have the same mean age (10.60 vs 10.70), and the
  strongest decreasing gene keeps its sign with no mask applied.

### Open items

- Sex effects (dispersion ~ age + sex, and the mean-model sex coefficient)
  not yet run on the final model.
- mashr for basket_cells and cerebellar_neurons.
- Small-N regions (e.g. CN, ~12 animals) looked miscalibrated in an earlier
  pre-mask check; not re-checked under the final model.
- Oligodendrocyte files keep an `opc-olig_` filename prefix (harmless).

### Superseded approaches (for context)

Each of these was tried and replaced; the runs are in the archive directories.

- Mean-aggregated pseudobulk was replaced by summed counts.
- Fitting raw counts was replaced by fitting TMM log-CPM.
- The `dispersion.summary` p-value extraction was replaced by `dispersion.fit`
  (see step 3).
- Rejected fixes for low-cell animals: per-animal minimum-cell filters (10 and
  100 cells), a 3SD + <20-cell outlier rule, and a condition filter on median
  cells per animal. All of these dropped animals or whole conditions.
- Rejected ways of adjusting standard errors: the sqrt(bvar) vs raw-bvar Shat
  comparison (settled on raw) and a genomic-control rescaling.
- The `qval_tenth` recurrence rule for the strong subset was replaced by
  `qval_single`.
