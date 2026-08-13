# dispersion

Analysis pipelines for age-associated changes in single-cell transcriptional
variability/dispersion in the macaque brain aging dataset (U01), plus
cross-pipeline comparison plots. Two independent methods are implemented and
compared:

- **`dglm/`** — DGLM-based approach (dispersion GLM: `mean ~ covariates`,
  `dispersion ~ age[+sex]`), adapted from Chiou et al. 2022 (*Nat Neurosci*),
  with mashr for effect-size sharing/shrinkage across regions and cell
  types. **Two real bugs were found and fixed** after initial runs had
  already completed: (1) `dglm_mashr.R` incorrectly applied `sqrt()` to
  `bvar`, a field that is actually a standard error, not a variance,
  distorting every mashr result run before the fix; (2) `pseudobulk.py`
  had no minimum per-animal cell count, letting 1-3-cell "averages"
  destabilize fits. **V2 (`disp_age__mean_full_min100`) is being fully
  regenerated with both fixes; V3 and louvain-subtype resolution are
  paused until V2 is confirmed correct** — do not treat any V3 or
  louvain output as final. See [`dglm/AGENTS.md`](dglm/AGENTS.md) for
  full pipeline order, methodology, and current status.
- **`gene_variance/`** — Simpler OLS-based approach: between-individual
  (`|residuals| ~ age + sex` on pseudobulk means) and within-individual
  (per-animal cell-level variance `~ age + sex`) regressions. See
  [`gene_variance/AGENTS.md`](gene_variance/AGENTS.md) for pipeline order,
  filtering steps, and known gaps.

## `plotting_scripts/`

Cross-pipeline comparison plots between `dglm/` and `gene_variance/`:
- `plot_concordance_heatmap.py`, `plot_overlap.R` — significant-gene overlap
  (Venn diagrams, concordance counts). `plot_overlap.R` accepts either a DGLM
  master TSV (`--dglm_master`, preferred) or the legacy per-cell-type RDS
  directory (`--dglm_dir`, with mashr→raw fallback).
- `plot_effect_scatter.py`, `plot_rank_correlation.py` — effect-size/rank
  correlation between the two methods.
- `plot_heatmap_summary.py` — median standardized effect size per cell type x
  region; color comes from `--dglm_master` (raw betas) or `--gv_master`,
  significance asterisks from a separate `--mashr_tsv` (mash_lfsr) when
  given. Note: asterisk drawing is currently disabled in this script
  (hardcoded off) — check before relying on it.
- `plot_volcano.py`, `plot_mean_vs_variability.py` — per-gene volcano and
  mean-vs-variability-effect plots (cross-pipeline comparison; for a DGLM-only
  volcano at cell-type or louvain resolution, see `dglm/scripts/dglm_plot_volcano.R`).

Most of these read from each pipeline's master FDR TSV, so run
`dglm_fdr.R`/`dglm_fdr_combined.R` and `fdr_correct.py` first.

## Environments

- Python steps (pseudobulking, OLS regressions, outlier detection, most
  plots): conda env `mixed_models`
- R steps (DGLM model fitting, mashr, FDR pooling, some plots): conda env
  `mashr_env` — missing a C/Fortran compiler toolchain, which blocks
  installing a few Bioconductor packages (`topGO`, `GO.db`, `biomaRt`)
  needed for DGLM's GO enrichment step. Being resolved via a separate
  `go_env` conda env built from precompiled conda-forge/bioconda binaries
  (no compilation needed) — `mashr_env` itself is not being modified.

## Cluster conventions

All jobs run on ASU Sol (SLURM, `htc` partition unless noted; large
mashr/louvain runs use `highmem`). Working directory is
`/scratch/easmit31/dispersion/`. Cell types with the largest cell counts
(`glutamatergic_neurons`, `GABAergic_neurons`, `cerebellar_neurons`, and
OPC/oligodendrocyte splits) need bumped memory — see the `MEM_MAP` in each
submit script. Never overwrite existing checkpoints/figures; write new
versioned output directories instead. The DGLM pipeline runs two finalized
model versions (`disp_age__mean_full`, `disp_age_sex__mean_full`; a third,
`disp_age__mean_age_sex`, is deprioritized) under `dglm/<model_version>/`,
plus louvain-subtype resolution under `dglm/<model_version>_louvain/`
(microglia complete for `disp_age__mean_full`; remaining cell types and
the second model version in progress) — see `dglm/AGENTS.md`.

**Long-running jobs must be submitted with `sbatch`**, not run in an
interactive `salloc`/`srun` session or live in a terminal — those die when
the connection drops, even though `squeue` makes them look like real
background jobs. Check `sacct`'s `JobName` column if unsure (`interacti+`
means it wasn't actually detached).

## Data

h5ads per cell type live at
`/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update`, mapped by
cell type in `dglm/scripts/pseudobulk.py`'s `H5AD_MAP` and the equivalent in
`gene_variance` scripts. Human-macaque ortholog mapping:
`/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv`.
