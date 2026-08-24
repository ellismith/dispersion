cat /scratch/easmit31/dispersion/README.md
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
as a final, settled significance threshold** — both are being generated
and compared side-by-side (see `generate_dglm_explorer.py`'s Methods tab)
pending further input.

**Current status (V2, `disp_age__mean_full`, dispersion~age):** all four
combinations of (Shat mode: raw/sqrt) x (resolution: cell-type/louvain)
are being generated for direct comparison:
- raw + cell-type: complete (pseudobulk → DGLM → mashr → FDR → all plots)
- sqrt + cell-type: mashr in progress (`disp_age__mean_full_SQRT`)
- raw + louvain: DGLM complete for all 12 cell types; mashr in progress
  (`disp_age__mean_full_louvain`)
- sqrt + louvain: DGLM complete (symlinked); mashr in progress
  (`disp_age__mean_full_louvain_SQRT`)

**V3** (`disp_age_sex__mean_full`, dispersion~age+sex) is deprioritized
for now — its DGLM output has been updated with all three fixes above,
but its mashr results predate that and are stale. Not being worked on
until V2's four columns are settled.

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
  `mashr_env` — missing a C/Fortran compiler toolchain, blocking
  `topGO`/`GO.db`/`biomaRt`. GO enrichment instead uses a separate `go_env`
  conda env built from precompiled conda-forge/bioconda binaries.

## Cluster conventions

All jobs run on ASU Sol (SLURM). **Prefer the `htc` partition even for
large-memory (256G) jobs when the job is expected to finish well under 4
hours** — `htc` has consistently had better availability than `highmem`
this project, and 256G mashr/louvain runs have completed in ~2.5hrs.
Reach for `highmem` only if a genuinely longer time limit is needed.
Working directory is `/scratch/easmit31/dispersion/`. Cell types with the
largest cell counts (`glutamatergic_neurons`, `GABAergic_neurons`,
`cerebellar_neurons`, and OPC/oligodendrocyte splits) need bumped memory —
see the `MEM_MAP` in each submit script. `pseudobulk.py` now self-cleans
stale output on every run (see Bug 3 above), so reruns into an existing
directory are safe. Never overwrite an mashr/FDR result you want to keep —
use a clearly distinct directory name for any variant run (e.g. the
`_SQRT` suffix convention used for the sqrt Shat-mode comparison).

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

## `dglm/` — pseudobulk & filtering rebuild (applies within either V2 or V3)

A further round of fixes on top of the three above, specific to `pseudobulk.py`
and gene-level filtering:

**Fixed:**
1. **`pseudobulk.py` was averaging, not summing, raw counts per animal.**
   Per-animal pseudobulk values were the *mean* expression across a cell
   type/subcluster's cells, not the *sum* — this breaks `edgeR::cpm()`'s
   library-size normalization downstream, since CPM assumes summed raw
   counts. Fixed: now sums raw counts per animal x louvain-subcluster x
   region, computed via a single sparse indicator-matrix multiply per
   region rather than a per-animal loop (same math, faster at scale).
2. **`--min_cells` and the old gene-presence filter removed entirely** from
   `pseudobulk.py` — no gene-level filtering happens at the pseudobulk
   stage anymore. Filtering moved downstream to an explicit, inspectable
   CPM-based rule: a gene is kept if its CPM is **at least** a chosen
   cutoff (0.5 used for the current run) in **at least 50%** of that
   combo's animals — both comparisons inclusive (`>=`). Raw (non-TMM)
   library sizes. See `filter_percent_animals.R`; each run also writes a
   per-gene audit-trail CSV (`*_filter_stats_cutoff*.csv`) alongside the
   filtered matrix.
3. **Mashr run per cell type, not pooled across all 12 at once.** A
   deliberate change from the prior single-pooled-mashr approach — each
   cell type's louvain subclusters x regions form one independent mash
   fit (`dglm_mashr.R --mode combined --checkpoints <cell-type-scoped dir>`).
   `--shat_mode sqrt` only going forward (matches Chiou et al.'s exact
   method) — the raw-vs-sqrt calibration comparison above remains an open
   question, but this round commits to sqrt rather than continuing to
   generate both.

**New scripts:** `run_all_celltypes.sh` / `submit_remaining_pseudobulk.sh`
(pseudobulk all subclusters x regions per cell type), `filter_percent_animals.R`
+ `filter_all_pseudobulk.sh` (the CPM filter, cutoffs 0.5/1/2), `stage_filtered_for_dglm.sh`
+ `run_dglm_for_celltype.sh` (stages filtered output into `dglm_model.R`'s
expected naming, runs dglm per subcluster, collects into one shared
checkpoints dir per cell type), `submit_dglm_louvain_pipeline.sh` (submits
the full dglm → combined mashr → master TSV chain per cell type as
separate Slurm jobs). Diagnostic/QC plots: `plot_cpm_valley.py` (per-cell-type
and per-region CPM density with fixed threshold lines), `plot_genes_before_after.py`
/ `plot_genes_before_after_per_region.py` (gene counts surviving each cutoff,
with outlier subclusters labeled). Explorer rebuild: `build_dglm_explorer_data.R`
+ `build_dglm_explorer_html.py` (regenerates `dglm/dglm_explorer.html` — the
Model dropdown is now a single option, since this pipeline no longer produces
whole-cell-type or raw-shat variants).

**Cluster note:** `htc` caps jobs at 240 minutes (`--time=03:59:00` max, not
4:00:00). Long-running steps should be `sbatch`'d, not run in an interactive
session — confirmed firsthand this round that an interactive-session process
(pseudobulk and a mashr run both) can die silently with no error and no
trace if the session is interrupted, even mid-run.

**Status as of this round:** microglia's per-cell-type combined mashr in
progress; astrocytes pseudobulk+filter done, dglm/mashr not yet started; the
other 10 cell types' pseudobulk+filter+dglm+mashr chain submitted and running
in parallel via `submit_dglm_louvain_pipeline.sh`.
