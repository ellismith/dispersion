# dispersion

Analysis pipelines for age-associated changes in single-cell transcriptional
variability/dispersion in the macaque brain aging dataset (U01), plus
cross-pipeline comparison plots. Two independent methods are implemented and
compared:

- **`dglm/`** — DGLM-based approach (dispersion GLM: `mean ~ covariates`,
  `dispersion ~ age[+sex]`), adapted from Chiou et al. 2022 (*Nat Neurosci*),
  with mashr for effect-size sharing/shrinkage across regions and cell
  types. See below for current status. See [`dglm/AGENTS.md`](dglm/AGENTS.md)
  for full pipeline order and methodology.
- **`gene_variance/`** — Simpler OLS-based approach: between-individual
  (`|residuals| ~ age + sex` on pseudobulk means) and within-individual
  (per-animal cell-level variance `~ age + sex`) regressions. See
  [`gene_variance/AGENTS.md`](gene_variance/AGENTS.md) for pipeline order,
  filtering steps, and known gaps.

## `dglm/` — three fixed bugs

**Fixed:**
1. **Shat/sqrt mislabeling.** `dglm_model.R`'s `bvar` output field is the
   dispersion submodel's Std. Error directly — confirmed twice, directly
   from the `dglm` package's own output (on synthetic test data and on a
   real gene from the pipeline; the coefficient table's own printed header
   reads `"Std. Error"`, not a variance). `bvar` is a variable name this
   pipeline chose — it is not a label the package assigns, and does not
   reflect a different underlying quantity. `dglm_mashr.R` previously
   applied an unconditional `sqrt()` to it; this is now controlled by
   `--shat_mode` (`raw` = Shat=bvar as-is, `sqrt` = Shat=sqrt(bvar)) so
   both can be run and compared directly. See "Shat calibration" below —
   neither mode is well-calibrated on this dataset, and this remains
   unresolved.
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

### Shat calibration — open, unresolved

Neither `--shat_mode raw` nor `--shat_mode sqrt` is well-calibrated on this
dataset, and switching between them does not resolve the underlying issue.
A null-calibration check (fraction of raw DGLM z-scores exceeding 1.96,
which should be ~5% under the null) gives ~25% for raw (too liberal) and
~0.05% for sqrt (far too conservative). A separate, independent check on
`GABAergic_neurons` (`mash_1by1`, no `cov_ed`) found the same pattern from
the other direction: under `raw`, 92-95% of all genes cleared a
significance threshold — consistent with "too liberal," not a fix.

Diagnosis: `bvar` is nearly constant across all genes within a given cell
type x region (driven almost entirely by dglm's fixed theoretical
dispersion constant and the shared design, not per-gene fit quality) and
shows no correlation with animal sample size (Spearman ρ=-0.002 across 92
conditions) — ruling out "just needs more animals" as an explanation. This
is likely a structural property of `dglm`'s dispersion-submodel SE
calculation at this sample size (~30-50 animals/condition) — the same
class of problem empirical-Bayes shrinkage methods (mashr, limma,
DESeq2, etc.) exist to work around, not a one-line fixable bug. **Do not
treat either `raw` or `sqrt` mashr output as a final, settled significance
threshold** on its own.

## `dglm_mashr.R` — strong-subset selection, `--approach`

mashr's data-driven covariance step (`cov_pca`/`cov_ed`) needs a "strong
subset" of genes with reliable signal to learn shared effect patterns
from; those patterns are then applied to every gene. Three strong-subset
definitions were tested on `GABAergic_neurons` (300-cell/30-animal
filtered data, 127 conditions, ~15,479 genes after the autosome+X filter):

1. **`--approach qval_third`**: DGLM q<0.05 in ≥1/3 of all conditions.
   The original rule this pipeline is adapted from, calibrated for a
   coarser ~15-region resolution; at this pipeline's 127 subcluster x
   region conditions, ≥1/3 means ≥43 — **7 genes**. Too few to learn a
   reliable covariance pattern from.
2. **`--approach lfsr_single`**: a gene counts if mashr's own
   `mash_1by1` LFSR is significant in at least one condition. 47 genes at
   the standard thresh=0.05, 525 at thresh=0.1. **Directly checked: 100%
   of genes selected this way are significant in exactly one condition;
   0% recur in a second.** This holds at every threshold tested —
   loosening the threshold only adds more single-condition genes; it does
   not find genes with real cross-condition recurrence, because the rule
   never checks for that.
3. **`--approach qval_recurrence`** (current default): DGLM q-value
   checked for recurrence across many conditions directly, instead of
   mashr's LFSR. Real signal confirmed independently: the overall raw
   q-value distribution shows genuine enrichment above the null (11.8% of
   all tests pass q<0.05, vs. 5% expected by chance). At q<0.01 in ≥10 of
   127 conditions: **1,428 genes**, with confirmed real recurrence (unlike
   `lfsr_single`).

`dglm_mashr.R` takes `--approach {qval_third,qval_recurrence,lfsr_single}`
(default `qval_recurrence`):
- `qval_third`/`qval_recurrence`: full two-stage fit — a random 50% gene
  subset is used to estimate the null correlation
  (`estimate_null_correlation_simple`, since this environment's mashr
  build lacks `estimate_null_correlation`) and fit canonical + data-driven
  covariances together; the resulting mixture is then applied to every
  gene (`fixg=TRUE`). Note: if the strong subset is a large fraction of
  all genes, the random subset will overlap it heavily and the
  null-correlation estimate becomes unreliable — at 1,428/15,479 genes
  (~9%) this is not a practical risk for the current default, but should
  be reconsidered before using a much larger strong subset.
- `lfsr_single`: simpler single-stage fit — data-driven covariances only,
  fit directly on all genes at once, no random-subset/null-correlation
  step.

Output filename includes both the approach and shat_mode
(`combined_dglm_mashr_results_<approach>_<shat_mode>.rds`) so results from
different settings never collide in the same checkpoints directory.

**Status, `GABAergic_neurons` only:** `qval_recurrence` + `sqrt` is the
current default and has been run on the 300-cell/30-animal filtered data.
`qval_third` and `lfsr_single`, both with `--shat_mode raw`, are in
progress for direct comparison. Not yet applied to any other cell type or
to the unfiltered (194-condition) data.

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
  `topGO`/`GO.db`/`biomaRt`, and missing `estimate_null_correlation`
  (using `estimate_null_correlation_simple` instead — see `--approach`
  above). GO enrichment instead uses a separate `go_env` conda env built
  from precompiled conda-forge/bioconda binaries.

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
`public` can have longer queue waits if its nodes are reserved/drained, so
`htc` is usually still worth trying first even for larger-memory jobs,
switching to `public` only if `htc` doesn't fit the time budget.
`pseudobulk.py` self-cleans stale output on every run, so reruns into an
existing directory are safe. Never overwrite an mashr/FDR result you want
to keep — use a clearly distinct directory name for any variant run.

## Data

h5ads per cell type live at
`/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update`, mapped by
cell type in `dglm/scripts/pseudobulk.py`'s `H5AD_MAP` and the equivalent in
`gene_variance` scripts. Human-macaque ortholog mapping:
`/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv`.
Autosome+X gene keep-list (used by `dglm_mashr.R` and related scripts):
`/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv`.

## `dglm/` — pseudobulk, filtering, and per-cell-type mashr

`pseudobulk.py` sums raw counts per animal x louvain-subcluster x region
(no gene-level filtering at this stage). Gene filtering happens downstream
via `filter_percent_animals.R`: a gene is kept if its CPM is at least a
chosen cutoff (0.5 used for the current run) in at least 50% of that
combo's animals (raw, non-TMM library sizes); each run also writes a
per-gene audit-trail CSV (`*_filter_stats_cutoff*.csv`) alongside the
filtered matrix.

mashr is run per cell type, not pooled across all 12 at once — each cell
type's louvain subclusters x regions form one independent mash fit
(`dglm_mashr.R --checkpoints <cell-type-scoped dir>`). Two dispersion
models are run per cell type: `dispersion ~ age` and `dispersion ~ age +
sex`. Each cell type's results live under
`dglm/disp_age__raw_counts/<cell_type>/`. Master TSV columns: `ensembl_id,
symbol, cell_type, region, beta, bvar, pvalue, mash_beta, mash_lfsr,
qvalue`.

opc and oligodendrocytes share a raw subcluster label prefix from their
common source clustering, but are genuinely separate populations, correctly
split at the pseudobulk stage (opc = subclusters 12,13; oligodendrocytes =
subclusters 0-11) — the plotting scripts and explorer relabel these for
display.

**Status:** all 12 cell types complete for both the age-only and age+sex
dispersion models (using `--approach qval_recurrence --shat_mode sqrt`,
the pipeline default). The strong-subset investigation above has only been
applied to `GABAergic_neurons` so far, as a targeted deep-dive — not yet
extended to the other 11 cell types.

### Interactive explorer

`scripts/build_dglm_explorer.py` reads the master TSVs directly (via
pandas, no R/mashr/`.rds` involved) and writes `dglm_explorer.html`. The
unit selector has three tiers (all subclusters summed, one cell type, one
subcluster); a region checkbox filter restricts which regions appear
across every tab; the heatmap tab has a toggle to show numeric values on
cells.

**Pending:** a "Method" dropdown to compare the different `--approach` x
`--shat_mode` combinations for `GABAergic_neurons`, once the
`qval_third`/`raw` and `lfsr_single`/`raw` runs finish. Scoped to
`GABAergic_neurons` only for now, since no other cell type has more than
one method's results.

### Not part of the active pipeline

`scripts/archive/pool_and_filter_celltype.py` filters CPM at the cell-type
level rather than per subcluster x region. Not connected to DGLM/mashr.

A condition-level sparsity filter (`scripts/assess_condition_sparsity.py`
+ `scripts/filter_conditions_for_mashr.R`) was trialed on
`GABAergic_neurons`, dropping subcluster x region combinations below a
minimum cell/animal count before mashr. Findings: this genuinely helps
some cell types (e.g. `microglia`, `vascular_cells`) but barely helps
others (`ependymal_cells`, `glutamatergic_neurons` stay badly sparse even
at loose thresholds) — see
`dglm/disp_age__condition_sparsity_diagnostics/threshold_comparison_all_celltypes.csv`
for real numbers across all 12 cell types at two thresholds. Not currently
applied to the main pipeline; superseded for now by the strong-subset
investigation above, which addresses a related but distinct problem
(which genes teach mashr its covariance patterns, vs. which conditions
mashr sees at all).
