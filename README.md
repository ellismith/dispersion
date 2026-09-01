# dispersion

Analysis pipelines for age-associated changes in single-cell transcriptional
variability/dispersion in the macaque brain aging dataset (U01), plus
cross-pipeline comparison plots. Two independent methods are implemented and
compared:

- **`dglm/`** — DGLM-based approach (dispersion GLM: `mean ~ covariates`,
  `dispersion ~ age[+sex]`), adapted from Chiou et al. 2022 (*Nat Neurosci*),
  with mashr for effect-size sharing/shrinkage across regions and cell
  types. See [`dglm/AGENTS.md`](dglm/AGENTS.md) for full pipeline order.
- **`gene_variance/`** — Simpler OLS-based approach: between-individual
  (`|residuals| ~ age + sex` on pseudobulk means) and within-individual
  (per-animal cell-level variance `~ age + sex`) regressions. See
  [`gene_variance/AGENTS.md`](gene_variance/AGENTS.md).

## `dglm/` — pseudobulk and DGLM fitting

`pseudobulk.py` sums raw counts per animal x louvain-subcluster x region
(no gene-level filtering at this stage; self-cleans stale output from
prior runs on every rerun, so reruns into an existing directory are
safe). Gene filtering happens downstream via `filter_percent_animals.R`:
a gene is kept if its CPM is at least a chosen cutoff (0.5 for the current
run) in at least 50% of that combo's animals (raw, non-TMM library
sizes); each run writes a per-gene audit-trail CSV
(`*_filter_stats_cutoff*.csv`) alongside the filtered matrix.

`dglm_model.R` fits DGLM per gene, per region, per subcluster:
`expression ~ age + sex + mean_n_umi + n_cells` (mean submodel),
`dispersion ~ age` (dispersion submodel, `dlink='log'`). Requires a
minimum of 100 cells per animal per subcluster x region (`--min_cells`)
to avoid unstable near-single-cell pseudobulk values. Output columns:
`beta`, `bvar`, `pval`, `qval` (FDR-corrected per condition). `bvar` is
the dispersion submodel's Std. Error directly — confirmed against the
`dglm` package's own printed output — not a variance, despite the name
(inherited from the original script this pipeline is adapted from).

Scripts: `run_all_celltypes.sh` (pseudobulk all subclusters x regions per
cell type), `filter_all_pseudobulk.sh` (CPM filter), `stage_filtered_for_dglm.sh`
+ `run_dglm_for_celltype.sh` (stage filtered output, run DGLM per
subcluster, collect into one shared checkpoints dir per cell type),
`submit_dglm_louvain_pipeline.sh` (submits the full per-cell-type chain as
Slurm jobs).

## `dglm_mashr.R` — mashr fitting

Current approach: for each gene, fit a one-condition-at-a-time adaptive
shrinkage model (`mash_1by1`) across all 127 subcluster x region
conditions; a gene is selected into the "strong subset" if its resulting
LFSR is significant (< 0.05) in at least one condition
(`get_significant_results(m.1by1, thresh=0.05)`). That subset is used to
learn data-driven covariance patterns — PCA (`cov_pca`, top 5 components)
followed by extreme deconvolution (`cov_ed`, which accounts for
measurement noise when estimating the true covariance structure). The
resulting patterns are then applied directly to every gene in one
`mash()` fit — no random-subset step, no null-correlation estimation;
single-stage.

**`--shat_mode raw`** (Shat = `bvar` as-is, no `sqrt()`) is the parameter
used for the current `GABAergic_neurons` run. `bvar` being confirmed as a
genuine Std. Error is the direct justification — mashr's documentation
specifies Std. Error as its expected input, and raw is that value
untransformed. Two things worth flagging about this choice: (1) a
separate calibration check (fraction of raw DGLM z-scores exceeding 1.96,
expected ~5% under the null) found raw gives ~25% — too liberal — so this
parameter is not treated as fully resolved; (2) the resulting strong
subset for `GABAergic_neurons` is large under this setting (14,435 of
15,568 genes, ~93%), consistent with that same calibration concern.
`--shat_mode sqrt` remains the pipeline-wide default for the other 11
cell types, unaffected by this choice.

`dglm_mashr.R --approach {qval_third,qval_recurrence,lfsr_single}
--shat_mode {sqrt,raw}` — `lfsr_single` is the approach described above
and the one in current use; `qval_third` and `qval_recurrence` are
alternative strong-subset rules based on DGLM's own q-values rather than
`mash_1by1`, kept in the script for comparison but not currently in use.
Output filename includes both settings
(`combined_dglm_mashr_results_<approach>_<shat_mode>.rds`) so results
from different runs never collide in the same checkpoints directory.

**Status:** `GABAergic_neurons`, 300-cell/30-animal filtered data
(`disp_age__min300cells_min30animals/`), `--approach lfsr_single
--shat_mode raw` — run in progress. All 12 cell types' main pipeline
results (unfiltered, both dispersion models) still use the prior default
(`--approach qval_recurrence --shat_mode sqrt`), unaffected by this run.

## `plotting_scripts/`

Cross-pipeline comparison plots between `dglm/` and `gene_variance/`:
`plot_concordance_heatmap.py` / `plot_overlap.R` (significant-gene
overlap), `plot_effect_scatter.py` / `plot_rank_correlation.py`
(effect-size/rank correlation), `plot_heatmap_summary.py` (median
standardized effect size; note: significance asterisks are currently
disabled, hardcoded off), `plot_volcano.py` / `plot_mean_vs_variability.py`.
For a DGLM-only volcano at cell-type or subcluster resolution, see
`dglm/scripts/dglm_plot_volcano.R`. Most read from each pipeline's master
FDR TSV — run `dglm_fdr.R`/`dglm_fdr_combined.R` and `fdr_correct.py` first.

## Environments

- `mixed_models` (Python, `gene_variance/`); `latent_analysis` (Python,
  `dglm/` — pseudobulking, plotting, the interactive explorer).
- `mashr_env` (R — DGLM, mashr, FDR pooling, most plots). Missing a
  C/Fortran toolchain (blocks `topGO`/`GO.db`/`biomaRt`) and missing
  `estimate_null_correlation` (only affects `qval_third`/`qval_recurrence`,
  which use `estimate_null_correlation_simple` instead). GO enrichment
  uses a separate `go_env` conda env.

## Cluster conventions

ASU Sol (SLURM). Working directory `/scratch/easmit31/dispersion/`.
**Long-running jobs must use `sbatch`**, not an interactive
`salloc`/`srun` session — those die silently on disconnect even though
`squeue` makes them look like real background jobs (check `sacct`'s
`JobName`; `interacti+` means it wasn't detached). Partition guide: `htc`
(4hr cap, best default availability) for most jobs; `public`/`highmem`
(7-day cap) for jobs needing more time or memory, though `public` can
have longer queue waits if its nodes are reserved/drained. Never
overwrite an mashr/FDR result worth keeping — use a distinct directory
name for any variant run.

## Data

h5ads: `/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update`
(mapped per cell type in `pseudobulk.py`'s `H5AD_MAP`). Human-macaque
ortholog map:
`/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv`.
Autosome+X gene keep-list:
`/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv`.

opc and oligodendrocytes share a raw subcluster label prefix from their
common source clustering but are genuinely separate populations, split
correctly at the pseudobulk stage (opc = subclusters 12,13;
oligodendrocytes = 0-11) — plotting scripts and the explorer relabel these
for display.

## Interactive explorer

`scripts/build_dglm_explorer.py` reads master TSVs directly (pandas, no
R/mashr/`.rds`) and writes `dglm_explorer.html`. Unit selector has three
tiers (all subclusters summed, one cell type, one subcluster); a region
checkbox filter applies across every tab; heatmap tab has a numeric-value
toggle. **Pending:** a "Method" dropdown to show the current
`GABAergic_neurons` run alongside the pipeline-default result once it
finishes — scoped to that cell type only, since no other cell type has
more than one method's results yet.
