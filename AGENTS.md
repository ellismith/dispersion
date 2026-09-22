# AGENTS.md — DGLM (transcriptional dispersion) pipeline

Guidance for AI coding agents (and future-me) working in this directory. This
pipeline estimates age-associated changes in transcriptional **dispersion**
(variance beyond what's explained by the mean) per gene, per cell type (or
louvain subcluster), per brain region, in the macaque brain aging dataset.
It's adapted from the DGLM approach in Chiou et al. 2022 (*Nat Neurosci*,
github.com/CayoBiobankResearchUnit/brain_transcriptome_aging_bulk), applied
to pseudobulked snRNA-seq instead of their bulk RNA-seq across 15 regions
(no cell-type dimension in their data — that's ours to add).

Working directory on Sol: **`/scratch/easmit31/dispersion/dglm/`**.

## Pipeline order

1. **`submit_pseudobulk.sh`** → `scripts/pseudobulk.py --outdir <dir>
   --covariates {age_sex,full} [--subcluster_col ct_louvain]`
   Per cell type, computes pseudobulk mean expression per animal per region
   from the h5ad (backed mode, sparse ops).
   - `--covariates full` (default): metadata gets `age, sex, mean_n_umi,
     n_cells`. `--covariates age_sex`: just `age, sex`.
   - `--subcluster_col <col>` (optional): instead of pseudobulking the whole
     cell type as one unit, splits by every distinct value of that obs
     column (e.g. `ct_louvain` → `microglia_0` .. `microglia_16`) and
     pseudobulks each separately, using the subcluster's own value as the
     output label. `--min_subcluster_cells` (default 100) skips any
     subcluster too small to bother with. Downstream scripts don't need to
     know the difference — a subcluster label just slots in as a
     `--cell_type` value like any other.
   - opc/oligodendrocytes share one h5ad (`Res1_opc-olig_subset.h5ad`),
     split by the `louvain` column (clusters `12`/`13` = opc, else
     oligodendrocytes) — this filter always applies first, before any
     `--subcluster_col` split.

   Writes, per unit (cell type or subcluster):
   `{label}_{region}_pseudobulk.csv`, `{label}_metadata.csv`,
   `{label}_gene_names.csv`.

2. **`submit_dglm.sh`** → `scripts/dglm_model.R --outdir <dir>
   --covariates {age_sex,full} --dispersion {age,age_sex}`
   Fits one DGLM per gene per region: `mean ~ covariates`,
   `dispersion ~ predictor(s)` (gaussian family, log dispersion link).
   The three model versions map onto these two flags:
   - `disp_age__mean_age_sex`: `--covariates age_sex --dispersion age`
   - `disp_age__mean_full`: `--covariates full --dispersion age` (default)
   - `disp_age_sex__mean_full`: `--covariates full --dispersion age_sex`

   `--dispersion age` produces a 4-stat array (`beta`/`bvar`/`pval`/`qval`
   per region). `--dispersion age_sex` produces an 8-stat array
   (`beta_age`/`bvar_age`/`pval_age`/`qval_age` +
   `beta_sex`/`bvar_sex`/`pval_sex`/`qval_sex`). `dglm_mashr.R` auto-detects
   which shape it's looking at — no flag needed there.

   Validates required metadata columns are present (`stop()`s with a clear
   error otherwise) and saves `run_info` (covariates, dispersion mode, date,
   cell_type, outdir) into the results RDS, so you can always check after
   the fact what a saved result was fit with.

   **This one script replaces what used to be two separate scripts**
   (`dglm_model.R` for the standard case, `dglm_model_age_sex_disp.R` for
   the two-term case) — those are obsolete now.

3. **`submit_mashr.sh`** (per cell type) / **`submit_mashr_combined.sh`**
   → `scripts/dglm_mashr.R --mode {per_ct,combined} --checkpoints <dir>
   [--cell_type <ct>] [--fast]`
   `combined` mode is the true structural match to Chiou et al. (they only
   ever ran one mashr across all 15 regions, no cell-type dimension) —
   `per_ct` is our own extension with no analog in the paper.

   **Discovers units dynamically** from whatever `*_dglm_results.rds` files
   actually exist in `--checkpoints` (not a hardcoded 12-cell-type list) —
   this is what makes `combined` mode work unchanged on a louvain
   checkpoints dir full of `microglia_0`..`microglia_16` instead of the 12
   canonical cell types, with zero code changes needed.

   **Auto-detects single-term vs two-term** dglm output from column names
   and runs mashr once (single term) or **twice** (once per term: age, sex)
   accordingly, saving separate output files per term.

   **`--fast`**: skips data-driven covariance estimation (`cov_pca`/
   `cov_ed`), uses canonical covariances only. Large condition counts (100+,
   e.g. a full louvain combined run) make the data-driven step scale badly
   — component count grows much faster than condition count, and a
   176-condition run that should've taken roughly double a working
   84-condition run (44 min) instead ran 2+ hours without finishing. Use
   `--fast` for anything with many conditions; it's a real, documented
   accuracy/speed tradeoff (canonical-only misses correlation patterns
   data-driven covariances would catch), not a hack.

   **This one script replaces what used to be two separate scripts**
   (`dglm_mashr.R` for the standard case, a since-removed
   `dglm_mashr_age_sex_disp.R` for the two-term case) — the two-term case
   is now just auto-detected within the same script.

   Uses **two separate cutoffs** (see Shared config) — matching a detail in
   Chiou et al.'s actual script: `strong.subset.qval.cutoff` (0.05) selects
   which genes are "strong" enough to drive `cov_pca()`/`cov_ed()`
   estimation; `fsr.cutoff` (0.2) is only the *final* LFSR significance
   threshold. These feed different steps and should not be collapsed into
   one value.

   **Matches Chiou et al.'s actual math as closely as this environment
   allows** (confirmed by cloning their real repo, not going from memory):
   - No z-score scaling of `Bhat`/`Shat` (their script doesn't; removing
     ours also fixed a real bug — `mash_beta` had been written to master
     TSVs in scaled units with nothing multiplying it back out)
   - Null correlation: their exact call is `estimate_null_correlation(temp,
     temp.U.c)`, but that function **does not exist** in mashr 0.2.79
     (installed here) — confirmed directly via `ls('package:mashr')`, only
     `estimate_null_correlation_simple()` and `mash_estimate_corr_em()` are
     exported. Uses `_simple()` for this reason, not as an unexamined
     shortcut.

   Output filenames carry both cutoff params:
   `{cell_type}_dglm_mashr_results[_<term>]_strong{strong.subset.qval.cutoff}_lfsr{fsr.cutoff}.rds`
   (per_ct) / `combined_dglm_mashr_results[_<term>]_strong..._lfsr....rds`
   (combined).

4. **FDR pooling**: `scripts/dglm_fdr.R` (per-cell-type → `master_dglm_globalfdr.tsv`)
   or `scripts/dglm_fdr_combined.R` (from the combined mashr RDS →
   `master_dglm_combined.tsv`, includes `mash_beta`/`mash_lfsr`). Both
   construct the params-suffixed mashr filename themselves from
   `_include_options.R`, so a cutoff change just needs that file updated,
   not a new flag. For the two-term (age+sex dispersion) case, use
   `dglm_fdr_age_sex_disp_combined.R` instead — produces two master TSVs
   (`master_dglm_age_disp_combined.tsv`, `master_dglm_sex_disp_combined.tsv`).

   **IMPORTANT: rerun this after any mashr rerun before trusting plots.**
   The master TSV is a separate artifact from the mashr RDS — regenerating
   the RDS does not automatically update the TSV.

5. **Plots**, all under `scripts/`:
   - `dglm_plot_pval.R`, `dglm_plot_beta.R` — per cell type. `dglm_plot_beta.R`
     was redesigned to drop its old per_ct-mashr dependency entirely:
     DGLM-only plots (density/bar count) read straight from
     `{cell_type}_dglm_results.rds` (no mashr needed at all), mashr plots
     filter the *combined* master TSV by `cell_type` — one shared mashr
     source instead of needing a separate per-cell-type mashr rerun.
     `--dglm_qthresh` (default 0.05) is separate from `fsr.cutoff` (mashr
     plots only) — don't reuse one threshold for both.
   - `dglm_plot_heatmap.R` — `--master_tsv` (required), `--sig_col
     {qvalue,mash_lfsr}`, `--qthresh`, `--sex_mode`. Cell-type acronym labels
     from `_include_options.R`'s `cell.type.labels`.
   - `dglm_plot_fig4.R` — Panel a (per-region/per-celltype sig-gene bar
     chart, the paper's Fig 4a/b analog) defaults to `master_dglm_combined.tsv`
     + `mash_lfsr<0.2`, auto-scales y-axis to data magnitude (not a fixed
     range), `--by_celltype` facets use free_y. Panel b (per-animal z-score
     jitter for one example gene) defaults to only showing regions where
     that gene is significant for the given cell type — `--regions`/
     `--all_regions` to override.
   - `dglm_plot_volcano.R` — **mashr-based** volcano (`mash_beta` vs.
     `-log10(mash_lfsr)`), faceted by region, for one `--cell_type` (works
     identically for a cell type or a louvain subcluster label). Top-N
     significant genes by `|beta|` labeled.
   - `dglm_plot_volcano_raw.R` — **early-look only, pre-mashr.** Reads
     `{cell_type}_dglm_results.rds` directly, no mashr required at all — use
     to look at results the moment `dglm_model.R` finishes. **Read the
     header before trusting the y-axis** — see "Known limitation" below.
   - `dglm_plot_summary_scatter.R` — one point per unit (cell type or
     subcluster), x = number of *distinct genes* significant in ≥⅓ of that
     unit's own tested regions (matches Chiou et al.'s `mash.sig` criterion
     exactly — their `≥5` of 15 regions is also ⅓, just adapted to work when
     units have varying numbers of tested regions instead of always 15;
     **not** a raw sum of gene×region hits, which conflates "many genes each
     hitting once" with "few genes hitting everywhere"), y = mean `|beta|`
     among those genes, color = % increasing vs. decreasing. `ggrepel` is
     optional in both this and `dglm_plot_volcano.R`/`dglm_plot_volcano_raw.R`
     — falls back to plain (non-repelled) `geom_text` if not installed, so a
     missing package never blocks the plot.
   - `dglm_go_enrichment.R` — adapted from Chiou et al.'s `dglm_topgo.R`
     (topGO Fisher's exact + KS tests, region-specificity permutation test).
     **Written but never run** — blocked on missing R packages (`topGO`,
     `GO.db`, `biomaRt`) in `mashr_env`, which itself is blocked on a
     missing C/Fortran compiler toolchain in that conda env. Their exact
     `get.ks.pval`/`get.ks.score`/`get.fisher.pval`/`get.fisher.score`
     helper functions (needed to override topGO's internal S4 methods) were
     added to `_include_options.R`, pulled from their real repo rather than
     guessed. Runs the paper's "union across regions" + "per-region"
     analysis once per cell type (their script has no cell-type dimension
     to loop over) — their raw-DGLM-only block is commented out in their
     own script and was not replicated, since it was never actually run by
     them either.

Utility scripts: `identify_outliers.py`, `patch_gene_names.R` (unchanged,
see prior notes).

## Known limitation: raw per-region q-values are not real significance tests

`dglm()`'s dispersion-submodel standard error comes out **identical for
every gene** within a given cell_type × region (confirmed empirically: two
totally unrelated synthetic response vectors, same design matrix, produce
bit-identical SEs to full floating-point precision). This is a property of
the `dglm` package's method itself — SE derived from the design matrix
alone, not each gene's actual data — present regardless of input, and
equally true of Chiou et al.'s own results since they use the identical
package. **Not a bug, does not need fixing.**

Practical consequence: **within one region**, ranking genes by the raw
per-region `qval` (or the pooled `qvalue` column in master TSVs) is
mathematically equivalent to ranking by `|beta|` alone — not real per-gene
statistical discrimination, despite the name. `mash_lfsr` does not have
this problem, since mashr compares each gene's pattern **across** regions,
where SE genuinely does vary (different animal counts per region).

**Rule of thumb: `mash_lfsr` is the only trustworthy significance measure.**
Use raw `qvalue`/`qval` only for an early look before mashr has run, and
read it as "sorted by effect size," not a real significance test.

## Shared config

All R scripts `source()` `scripts/_include_options.R` first — single source
of truth for cell type list, region list, h5ad file map, model covariates,
mashr thresholds (`fsr.cutoff` = 0.2, final LFSR significance;
`strong.subset.qval.cutoff` = 0.05, strong-subset selection — see step 3 for
why these differ), the paper's exact topGO helper functions, plotting
colors, `cell.type.labels` acronyms, and the random seed (42). Mashr output
filenames embed both cutoff values, so changing either produces a new
filename rather than silently overwriting old results.

## Paths and environments

- Working directory: `/scratch/easmit31/dispersion/dglm/`. Scripts default
  to this correctly now (an earlier stale `variability/` default has been
  fixed).
- h5ad source: `/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update`
  — filenames follow `Res1_<cell-type>_<suffix>.h5ad`, mapped in
  `pseudobulk.py`'s `H5AD_MAP` (suffix varies: `_update`, `_new`, `_subset`
  — don't assume a pattern, check the map). Subtype labels come from the
  `ct_louvain` obs column (confirmed present, cleanly self-prefixed e.g.
  `microglia_16` — other candidate columns like `cell_subcluster_assign`/
  `tacco_cell_subcluster` had cross-contamination even within a
  single-cell-type file and are not the right column to use).
- Ortholog map: `/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv`
- Python (pseudobulk, outliers): conda env `mixed_models`
- R (dglm, mashr, FDR, plots): conda env `mashr_env` — missing a C/Fortran
  compiler toolchain (blocks GO enrichment's packages; unresolved, low
  priority)
- SLURM: `htc` partition; `glutamatergic_neurons`, `GABAergic_neurons`,
  `cerebellar_neurons` need bumped memory — check `MEM_MAP` per submit
  script, not consistent across all of them.

## Model versions — current status

| Directory | Mean model | Dispersion model | Status |
|---|---|---|---|
| `disp_age__mean_full` | full | ~ age | **Correct and verified** — rerun with all current fixes (no z-scaling, real null-correlation, correct strong-subset cutoff), confirmed via log output + fresh file timestamps |
| `disp_age__mean_age_sex` | age_sex | ~ age | Not yet rerun with corrected cutoffs |
| `disp_age_sex__mean_full` | full | ~ age + sex | Not yet rerun with corrected cutoffs |

`disp_age_sex__mean_full`'s two-term output now works with the unified
`dglm_mashr.R` (auto-detected) — no longer needs a separate dedicated
mashr script.

## Louvain / subtype-level analysis

In progress. `dglm/disp_age__mean_full_louvain/` — microglia's 17
`ct_louvain` subclusters, pseudobulked and DGLM-fit (all succeeded, 7-11
regions each depending on subcluster size). Combined mashr on all 17 at
once (176 conditions) proved much more expensive than the reference
84-condition run — use `--fast` for this scale (see step 3). Not yet
expanded to other cell types; that's a straightforward repeat of the same
steps once microglia's results look right.

## Conventions (repo-wide)

- Never overwrite existing results — write to a new versioned
  checkpoint/figure directory rather than clobbering one that already has
  output.
- Before submitting a rerun, confirm `--covariates`/`--dispersion` match
  intent, and check `run_info` (or job log) after a run to see what
  actually happened, not what was assumed.
- SLURM jobs that need to survive a dropped connection **must** be
  submitted with `sbatch`, not run interactively (`salloc`/`srun` sessions,
  or live terminal execution) — an interactive allocation dies when the
  connection drops, same as a live terminal command, even though it looks
  like a "real job" in `squeue`. Check `JobName` in `sacct` — `interacti+`
  means it was never actually detached.
- Long-running jobs should log real progress (component counts, timing per
  step), not just a start/end message — silence for 2+ hours with no way to
  tell "slow but working" from "hung" is a real cost, not a minor
  inconvenience.
- Scripts are delivered/edited as whole-file heredocs or in-place fixes,
  not patched-up "v2" copies.
- Don't invent or guess file paths; verify with `ls`/`cat` before assuming.
