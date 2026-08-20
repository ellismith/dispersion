# dispersion

Analysis pipelines for age-associated changes in single-cell transcriptional
variability/dispersion in the macaque brain aging dataset (U01), plus
cross-pipeline comparison plots. Two independent methods are implemented and
compared:

- **`dglm/`** — DGLM-based approach (dispersion GLM: `mean ~ covariates`,
  `dispersion ~ age[+sex]`), adapted from Chiou et al. 2022 (*Nat Neurosci*),
  with mashr for effect-size sharing/shrinkage across regions and cell
  types. See below for current status — **the pipeline itself has three
  confirmed, fixed bugs, plus one open, unresolved statistical question**
  that should be understood before treating any significant-gene count as
  final. See [`dglm/AGENTS.md`](dglm/AGENTS.md) for full pipeline order and
  methodology.
- **`gene_variance/`** — Simpler OLS-based approach: between-individual
  (`|residuals| ~ age + sex` on pseudobulk means) and within-individual
  (per-animal cell-level variance `~ age + sex`) regressions. See
  [`gene_variance/AGENTS.md`](gene_variance/AGENTS.md) for pipeline order,
  filtering steps, and known gaps.

## `dglm/` — three fixed bugs, one open question

**Fixed:**
1. **Shat/sqrt mislabeling.** `dglm_model.R`'s `bvar` output field is the
   dispersion submodel's Std. Error directly (confirmed identical in
   Chiou et al.'s own original script — this naming is inherited, not
   introduced here), not a variance. `dglm_mashr.R` previously applied an
   unconditional `sqrt()` to it; this is now controlled by `--shat_mode`
   (`raw` = Shat=bvar as-is, `sqrt` = Shat=sqrt(bvar), matching Chiou's
   exact original method) so both can be run and compared directly.
2. **No per-animal minimum cell count.** `pseudobulk.py` had no minimum
   cells a single animal had to contribute to a given cell_type/subcluster
   x region pseudobulk value — an animal contributing 1-3 cells produced
   an unstable near-single-cell "average." Fixed with `--min_cells`
   (default 100).
3. **Stale-file bug in `pseudobulk.py`.** The script only ever added
   output files, never cleared old ones from a prior run — on a rerun, a
   region that no longer survived filtering could leave a stale file from
   an earlier run sitting alongside fresh metadata, causing silent
   mismatches downstream. Fixed: `pseudobulk_one()` now clears any
   pre-existing output for a given label before writing fresh files,
   every run.

**Open, unresolved:** neither `--shat_mode raw` nor `--shat_mode sqrt` is
well-calibrated on this dataset. A null-calibration check (fraction of raw
DGLM z-scores exceeding 1.96, which should be ~5% under the null) gives
~25% for raw (too liberal) and ~0.05% for sqrt (Chiou's exact method —
far too conservative, the opposite direction). Diagnosis: `bvar` is
nearly constant across all genes within a given cell type x region
(driven almost entirely by dglm's fixed theoretical dispersion constant
and the shared design, not per-gene fit quality) and shows no correlation
with animal sample size (Spearman ρ=-0.002 across 92 conditions) — ruling
out "just needs more animals" as an explanation. This is likely a
structural property of `dglm`'s dispersion-submodel SE calculation, not a
one-line fixable bug. **Do not treat either `raw` or `sqrt` mashr output
as a final, settled significance threshold.**

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
  plots): conda env `mixed_models` for `gene_variance/`; `latent_analysis`
  for `dglm/`'s Python steps (pseudobulking, plotting, the explorer).
- R steps (DGLM model fitting, mashr, FDR pooling, most plots): conda env
  `mashr_env` — missing a C/Fortran compiler toolchain, blocking
  `topGO`/`GO.db`/`biomaRt`. GO enrichment instead uses a separate `go_env`
  conda env built from precompiled conda-forge/bioconda binaries.

## Cluster conventions

All jobs run on ASU Sol (SLURM). Working directory is
`/scratch/easmit31/dispersion/`. **Long-running jobs must be submitted with
`sbatch`**, not run in an interactive `salloc`/`srun` session or live in a
terminal — those die when the connection drops, even though `squeue` makes
them look like real background jobs. Check `sacct`'s `JobName` column if
unsure (`interacti+` means it wasn't actually detached).

Partition guide: `htc` (4hr cap, ~200+ nodes, best default availability) for
most jobs; `highmem` (7-day cap, 11 nodes) or `public` (7-day cap, ~160
nodes, ~500GB/node) for jobs needing more memory or more than 4 hours —
`public` tends to clear faster than `highmem` for large-memory jobs since it
has far more nodes, as long as the memory request fits under a single
node's ~500GB ceiling. `pseudobulk.py` self-cleans stale output on every
run, so reruns into an existing directory are safe. Never overwrite an
mashr/FDR result you want to keep — use a clearly distinct directory name
for any variant run.

## Data

h5ads per cell type live at
`/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update`, mapped by
cell type in `dglm/scripts/pseudobulk.py`'s `H5AD_MAP` and the equivalent in
`gene_variance` scripts. Human-macaque ortholog mapping:
`/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv`.

## `dglm/` — pseudobulk, filtering, and per-cell-type mashr

`pseudobulk.py` sums raw counts per animal x louvain-subcluster x region
(no gene-level filtering at this stage). Gene filtering happens downstream
via `filter_percent_animals.R`: a gene is kept if its CPM is at least a
chosen cutoff (0.5 used for the current run) in at least 50% of that
combo's animals (raw, non-TMM library sizes); each run also writes a
per-gene audit-trail CSV (`*_filter_stats_cutoff*.csv`) alongside the
filtered matrix. `filter_all_pseudobulk.sh` wraps this over every
pseudobulk file in a directory.

mashr is run per cell type, not pooled across all 12 at once — each cell
type's louvain subclusters x regions form one independent mash fit
(`dglm_mashr.R --mode combined --checkpoints <cell-type-scoped dir>`),
using `--shat_mode sqrt` (matches Chiou et al.'s exact method). Two
dispersion models are run per cell type: `dispersion ~ age` and
`dispersion ~ age + sex` — pipeline scripts and directories use `age` /
`age_sex` (or `agesex`) to distinguish them. Each cell type's results live
under `dglm/disp_age__raw_counts/<cell_type>/`:
`dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv` for the age-only
model, `dglm_checkpoints_cutoff0.5_agesex/master_dglm_age_disp_combined.tsv`
(age term) and `master_dglm_sex_disp_combined.tsv` (sex term) for the
age+sex model. Master TSV columns: `ensembl_id, symbol, cell_type, region,
beta, bvar, pvalue, mash_beta, mash_lfsr, qvalue`.

Pipeline scripts: `run_all_celltypes.sh` / `submit_remaining_pseudobulk.sh`
(pseudobulk all subclusters x regions per cell type), `stage_filtered_for_dglm.sh`
+ `run_dglm_for_celltype.sh` (stages filtered output, runs dglm per
subcluster, collects into one shared checkpoints dir per cell type),
`submit_dglm_louvain_pipeline.sh` / `submit_v3_agesex_louvain_pipeline.sh`
(submit the full dglm → combined mashr → master TSV chain per cell type as
separate Slurm jobs, for the age-only and age+sex models respectively).
Diagnostic/QC plots: `plot_cpm_valley.py` (per-cell-type and per-region CPM
density with fixed threshold lines), `plot_genes_before_after.py` /
`plot_genes_before_after_per_region.py` (gene counts surviving each cutoff,
with outlier subclusters labeled).

opc and oligodendrocytes share a raw subcluster label prefix from their
common source clustering, but are genuinely separate populations, correctly
split at the pseudobulk stage (opc = subclusters 12,13; oligodendrocytes =
subclusters 0-11) — the plotting scripts and explorer relabel these for
display (`opc_12`, `oligodendrocytes_0`, etc.).

**Status:** all 12 cell types complete for both the age-only and age+sex
dispersion models.

### Plotting scripts (`scripts/`, run as Slurm jobs, not interactively)

- `plot_gene_count_heatmaps.R` — unit x region heatmaps. `--metric` selects
  what's shown: `genes_tested`, `genes_significant`, or mash-based effect
  metrics (`pct_sig`, `n_increase`, `n_decrease`, `n_net`, `mash_beta_sig`,
  `mash_beta_mag_sig`, `raw_beta_all`, `raw_beta_sig`). Produces
  subcluster-level and cell-type-level (aggregated) heatmaps together; add
  `--per_celltype` for one additional PNG per cell type. `--dispersion age`
  or `--dispersion age_sex` selects the model; output filenames include it.
- `plot_bar_sig_by_region.R` — significant-gene counts by region, split by
  direction (increase/decrease), same subcluster + cell-type dual-level
  structure.
- `dglm_plot_volcano.R` — per-subcluster volcano plots (effect size vs.
  significance), faceted by region.

### Interactive explorer

`scripts/build_dglm_explorer.py` reads the master TSVs directly (via
pandas) and writes `dglm_explorer.html` in one step:
```
conda run -n latent_analysis python scripts/build_dglm_explorer.py --out dglm_explorer.html
```
The unit selector has three tiers: all subclusters summed together, one
cell type (bar tab offers grouped small-multiples or summed across its
subclusters; heatmap offers grouped, mean-per-cell-type, or
median-per-cell-type; volcano shows grouped small-multiples), or one
specific subcluster. Subcluster/cell-type ordering is numeric (0,1,2…10,11).
The heatmap tab has a toggle to show numeric values on cells.

### Not part of the active pipeline

`scripts/archive/pool_and_filter_celltype.py` filters CPM at the cell-type
level (pooling all of a cell type's subcluster x region pseudobulk before
filtering) rather than per subcluster x region. Output sits in
`dglm/disp_age__celltype_pooled/` (gitignored). Not connected to
DGLM/mashr — the per-subcluster x region filtering described above is what's
actually used.
