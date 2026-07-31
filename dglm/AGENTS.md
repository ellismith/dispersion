# AGENTS.md — DGLM (transcriptional dispersion) pipeline

Guidance for AI coding agents (and future-me) working in this directory. This
pipeline estimates age-associated changes in transcriptional **dispersion**
(variance beyond what's explained by the mean) per gene, per cell type, per
brain region, in the macaque brain aging dataset. It's a direct adaptation of
the DGLM approach in Chiou et al. 2022 (Nat Neurosci), applied to pseudobulked
snRNA-seq instead of bulk RNA-seq.

## Pipeline order

Each stage reads the previous stage's checkpoint output. Run in this order:

1. **`submit_pseudobulk.sh`** → `scripts/pseudobulk.py`
   Per cell type, computes pseudobulk mean expression per animal per region
   from the h5ad (backed mode, sparse ops). Writes, per cell type:
   - `{cell_type}_{region}_pseudobulk.csv` (genes x animals)
   - `{cell_type}_metadata.csv` (animals x covariates)
   - `{cell_type}_gene_names.csv` (ensembl_id → external_gene_name)

2. **`submit_dglm.sh`** → `scripts/dglm_model.R`
   Fits one DGLM per gene per region: `mean ~ covariates`, `dispersion ~ age`
   (gaussian family, log dispersion link). Parallelized with `parApply`.
   Output: `{cell_type}_dglm_results.rds` — a genes x {beta,bvar,pval,qval} x
   regions array, with per-region BH-FDR, plus a `human_symbols` vector.

3. **`submit_mashr.sh`** (per cell type) and/or **`submit_mashr_combined.sh`**
   (all cell types x regions pooled) → `scripts/dglm_mashr.R`
   Runs mashr for effect-size shrinkage/sharing across conditions. `per_ct`
   mode shares across regions within a cell type; `combined` mode shares
   across cell type x region jointly. Single-region cell types skip mashr
   in `per_ct` mode (saved with `mash=NULL`).

4. **FDR pooling**: `scripts/dglm_fdr.R` (pools per-cell-type `per_ct` mashr
   results → `master_dglm_globalfdr.tsv`) or `scripts/dglm_fdr_combined.R`
   (extracts from `combined_dglm_mashr_results.rds` → `master_dglm_combined.tsv`,
   includes `mash_beta`/`mash_lfsr` columns). `submit_mashr_combined.sh` runs
   both the combined mashr and its FDR step back to back.

5. **`submit_plots.sh`** → `scripts/dglm_plot_pval.R`, `scripts/dglm_plot_beta.R`
   Per-cell-type p-value histograms/QQ plots and beta density/sig-gene-count
   plots. `scripts/dglm_plot_heatmap.R` and `scripts/dglm_plot_fig4.R` are run
   manually (not in a submit script) and prefer the master TSV over
   per-cell-type RDS files when available.

Utility scripts (run manually, not part of the linear chain):
- `scripts/identify_outliers.py` — flags animals by low cell count, expression
  z-score, or a hardcoded known-outlier list (currently `8H2:ACC`).
- `scripts/patch_gene_names.R` — back-fills `human_symbols` on existing
  `*_dglm_mashr_results.rds` files without rerunning the model (ortholog >
  external_gene_name > ensembl_id priority).

## Shared config

All R scripts `source()` `scripts/_include_options.R` first — this is the
single source of truth for cell type list, region list, h5ad file map,
model covariates, mashr thresholds (`fsr.cutoff`, `fraction.shared.cutoff`,
`fraction.unique.cutoff`), colors, and the random seed (42). Change config in
one place, not per-script.

## Paths and environments

- **Actual working directory on Sol: `/scratch/easmit31/dispersion/dglm/`.**
  The scripts in this repo default their `--checkpoints`/paths to an older
  location, `/scratch/easmit31/variability/dglm/` — that's a stale default
  from before the repo was reorganized (variability/ used to hold both this
  work and the separate diversity analysis). In practice every job is
  submitted with an explicit `--checkpoints .../dispersion/dglm/<model_version>`
  overriding the default, so don't trust the hardcoded default paths at the
  top of `submit_*.sh` at face value — check what's actually passed, and
  prefer updating the defaults if you're touching these scripts anyway.
- h5ad source: `/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update`
- Ortholog map: `/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv`
- Python (pseudobulk, outliers): conda env `mixed_models`
- R (dglm, mashr, FDR, plots): conda env `mashr_env`
- SLURM: `htc` partition for everything here; `glutamatergic_neurons`,
  `GABAergic_neurons`, `cerebellar_neurons` need bumped memory (128–256G) —
  see the `MEM_MAP` associative array at the top of each submit script (values
  aren't fully consistent across submit scripts — check the one you're using).

## Model versions

Three parallel model specs are run, each into its own subdirectory of
`/scratch/easmit31/dispersion/dglm/` (set via `_include_options.R` +
`--checkpoints` per run, not a script flag):

| Directory | Mean model | Dispersion model | Role |
|---|---|---|---|
| `disp_age__mean_full` | age + sex + mean_n_umi + n_cells | ~ age | **Main model** |
| `disp_age__mean_age_sex` | age + sex | ~ age | Simpler mean model, no technical covariates |
| `disp_age_sex__mean_full` | age + sex + mean_n_umi + n_cells | ~ age + sex | Tests whether sex also affects dispersion |

After `dglm_fdr.R`/`dglm_fdr_combined.R` produce their generically-named
`master_dglm_globalfdr.tsv` / `master_dglm_combined.tsv` inside a model-version
directory, they get renamed to disambiguate across the three versions, e.g.
`master_disp_age_mean_full_globalfdr.tsv`, `master_disp_age_mean_full_mashr.tsv`
(from `master_dglm_combined.tsv`). `disp_age_sex__mean_full` additionally
splits into age-effect and sex-effect master tables
(`master_disp_age_disp_age_sex_mean_full_globalfdr.tsv`,
`master_disp_sex_mean_full_globalfdr.tsv`) since its dispersion submodel has
two terms. When asked "what's in file X," check the model-version directory
name first — the three versions' outputs otherwise look identical column-wise.

## Known gotcha: repo's `pseudobulk.py` is missing covariate columns the model expects

`_include_options.R`'s `model.covariates = c('age', 'sex', 'mean_n_umi',
'n_cells')` matches the intended "full" mean model (used by
`disp_age__mean_full` and `disp_age_sex__mean_full`, see Model versions
above) — this is a real, intentional model spec, not a typo. But the
`pseudobulk.py` currently checked into this repo only writes `age`, `sex`,
`sequencing_run_id`, `sample_reads`, `n_umi` to the metadata CSV — there's no
`mean_n_umi` or `n_cells` column. `dglm_model.R`'s covariate filter
(`length(unique(m[[cv]])) > 1`) silently drops any covariate missing from the
metadata (`m[[cv]]` is `NULL`, `unique(NULL)` has length 0), so run against
this exact checked-in script, the "full" models would silently collapse to
`age + sex` only. This strongly suggests the committed `pseudobulk.py` is
behind whatever version actually generated the existing
`disp_age__mean_full`/`disp_age_sex__mean_full` results on Sol — before
rerunning either of those two model versions from scratch, find/restore the
version of `pseudobulk.py` that actually emits `mean_n_umi` and `n_cells`
(likely `n_cells` = cells contributing to each animal's pseudobulk row,
`mean_n_umi` = mean UMI count across those cells) and reconcile it with what's
committed here. Either way, check the logged
`message('  covariates: ', ...)` line in the job output to confirm what
actually got used before trusting any rerun's results.

## Conventions (apply repo-wide, not just here)

- Never overwrite existing results — write to a new versioned checkpoint/figure
  directory rather than clobbering one that already has output.
- Before submitting a rerun, verify covariates/config via the logged
  `message()` output (or `run_info`-equivalent) rather than assuming a script
  change took effect as intended (see gotcha above).
- Scripts are delivered/edited as whole-file heredocs or in-place fixes, not
  patched-up "v2" copies — once a bug-fixed version is validated, fix the
  original and delete any intermediate patched versions.
- Don't invent or guess file paths; verify with `ls`/`cat` before assuming.
