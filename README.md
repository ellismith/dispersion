# dispersion

Analysis pipelines for age-associated changes in single-cell transcriptional
variability/dispersion in the macaque brain aging dataset (U01), plus
cross-pipeline comparison plots. Two independent methods are implemented and
compared:

- **`dglm/`** — DGLM-based approach (dispersion GLM: `mean ~ covariates`,
  `dispersion ~ age`), adapted from Chiou et al. 2022 (Nat Neurosci), with
  mashr for effect-size sharing/shrinkage across regions and cell types.
  See [`dglm/AGENTS.md`](dglm/AGENTS.md) for the full pipeline order,
  environments, and known gotchas.
- **`gene_variance/`** — Simpler OLS-based approach: between-individual
  (`|residuals| ~ age + sex` on pseudobulk means) and within-individual
  (per-animal cell-level variance `~ age + sex`) regressions.

## `plotting_scripts/`

Cross-pipeline comparison plots between `dglm/` and `gene_variance/`
(concordance heatmaps, effect-size scatter, rank correlation, volcano plots,
mean-vs-variability relationship). Most read from each pipeline's master FDR
TSV, so run `dglm_fdr.R`/`dglm_fdr_combined.R` and `fdr_correct.py` first.

## Environments

- Python steps (pseudobulking, OLS regressions, outlier detection, most
  plots): conda env `mixed_models`
- R steps (DGLM model fitting, mashr, FDR pooling, some plots): conda env
  `mashr_env`

## Cluster conventions

All jobs run on ASU Sol (SLURM, `htc` partition unless noted). Working
directory is `/scratch/easmit31/dispersion/` (this repo's scripts default to
an older `/scratch/easmit31/variability/` path in a few places — always pass
`--checkpoints`/`--outdir` explicitly rather than relying on the defaults;
see `dglm/AGENTS.md` for specifics). Cell types with the largest cell counts
(`glutamatergic_neurons`, `GABAergic_neurons`, `cerebellar_neurons`, and
OPC/oligodendrocyte splits) need bumped memory — see the `MEM_MAP` in each
submit script. Never overwrite existing checkpoints/figures; write new
versioned output directories instead. The DGLM pipeline additionally runs
three parallel model versions under `dglm/<model_version>/` — see
`dglm/AGENTS.md`.

## Data

h5ads per cell type live at
`/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update`, mapped by
cell type in `dglm/scripts/_include_options.R` and `gene_variance` scripts.
Human-macaque ortholog mapping:
`/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv`.
