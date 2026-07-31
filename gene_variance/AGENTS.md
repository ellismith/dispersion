# AGENTS.md — gene_variance (OLS) pipeline

Guidance for AI coding agents (and future-me) working in this directory. This
is a simpler, independent alternative to the `dglm/` pipeline: instead of a
dispersion GLM, it uses plain OLS regressions on per-animal summary
statistics to test whether gene expression variability changes with age, per
gene, per cell type, per brain region. Used as a cross-check against `dglm/`
(see `plotting_scripts/`).

Working directory on Sol: **`/scratch/easmit31/dispersion/gene_variance/`**
(note: `submit_all.sh`, `submit_plots.sh`, and `run_mglia_hip.sh` currently
still default/hardcode paths under the older
`/scratch/easmit31/variability/gene_variance/` location — this is the same
stale-path issue that was fixed in `dglm/`'s scripts but hasn't been fixed
here yet. Update these three scripts' paths before relying on their defaults,
or override with explicit args.).

## Two variance modes, computed together

`run_gene_variance.py` is run once per cell type x region and computes both
of the following in a single pass (not separate scripts):

- **Between-individual**: pseudobulk mean expression per animal (log-
  normalized), then `|residuals(mean ~ age + sex)| ~ age + sex` — i.e. do the
  *animal-level* means get more spread out with age, after accounting for
  the age/sex trend in the mean itself.
- **Within-individual**: per-animal variance across that animal's own cells
  (log-normalized), `log1p(variance) ~ age + sex` — i.e. does cell-to-cell
  heterogeneity *within* an animal increase with age.

Both regressions share the same OLS covariates (`age`, plus sex as a
`pd.get_dummies` column) and the same per-gene summary stats computed from
**raw counts** (`.raw` if present, else `.X`): `mean_expr`, `var_between`,
`var_within` — these ride along in the output purely for plotting/filtering,
they aren't part of either regression.

## Pipeline order

1. **`submit_all.sh`** → `run_gene_variance.py --h5ad <path> --cell_type <ct>
   --region <region> --outdir results_log/ --min_age 1.0 --min_animals 10`
   One SLURM job per (cell type x region) pair. Writes, per cell type x
   region:
   - `{cell_type}_{region}_between.tsv`
   - `{cell_type}_{region}_within.tsv`
   - `{cell_type}_{region}_args.json` (records the exact args that job ran with)

   Filtering before regression: region filter → age filter (`>= min_age`) →
   opc/oligodendrocyte louvain split (clusters `12`,`13` = opc, everything
   else = oligodendrocytes; both come from the same `Res1_opc-olig_subset.h5ad`)
   → optional HVG filter (`--n_hvgs`, off by default) → per-gene
   `min_animals` expressing filter. Genes not expressed in enough animals
   after all upstream filters are dropped before regression, not before.

   `HTC_TYPES` get 128G, `HIGHMEM_TYPES` (`GABAergic_neurons`,
   `glutamatergic_neurons`) get 256G — see the `MEM_MAP`-equivalent
   arrays at the top of `submit_all.sh`. Note this doesn't include
   `cerebellar_neurons` at high-mem here, unlike `dglm/`'s submit scripts —
   check `submit_all.sh` directly rather than assuming parity with `dglm/`.

2. **`fdr_correct.py --indir results_log --outdir <dir>`**
   Pools every `*_between.tsv` and every `*_within.tsv` across all cell
   types/regions, drops NaN p-values, runs a single BH-FDR correction
   *within each mode* (between and within are corrected separately from each
   other). Outputs `master_between.tsv`, `master_within.tsv` — these are
   what `plotting_scripts/` and `dglm/`'s comparison plots read from
   (`--gv_master`).

3. **Plots** (`submit_plots.sh`, per cell type):
   - `plot_gene_variance.py` — per-gene scatter, x = log10(mean expression),
     y = age slope, colored by significance; `--region` for one region or
     `--all_regions` to overlay all 11.
   - `plot_heatmap.py` — count of significant genes (q<0.05) per cell type x
     region, 4 heatmaps (between/within x increase/decrease). `--louvain`
     mode (subtype-level rows within one cell type) is **not implemented
     yet** — it raises `NotImplementedError` if you pass it; the louvain
     column isn't in the master TSVs to support it.

`run_mglia_hip.sh` is a one-off manual rerun script for a single cell
type/region pair (microglia x HIP) — not part of the batch pipeline, useful
as a template if you need to rerun just one pair.

## Paths and environments

- h5ad source: `/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update`
  (cell-type → filename map is duplicated in `submit_all.sh` and
  `submit_plots.sh` — if you add/rename a cell type's h5ad, update both)
- Ortholog map (for `human_symbol` column): `/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv`
- Python env: `mixed_models` (same as `dglm/`'s pseudobulking step)
- SLURM: `htc` partition, per-job memory per the `HTC_TYPES`/`HIGHMEM_TYPES`
  split in `submit_all.sh` (see note above on high-mem cell types differing
  from `dglm/`)

## Known gaps / things to check before trusting output

- Regression failures are silently set to NaN per-gene (`slope, pval,
  intercept = np.nan, np.nan, np.nan`) with only an aggregate warning count
  printed — if you need to know *which* genes failed, that's not currently
  saved anywhere; you'd have to rerun with added logging.
- `--n_hvgs` (HVG filtering) exists as a flag but isn't used by
  `submit_all.sh`'s default invocation — full gene set is fit, not just HVGs,
  unless someone runs the script manually with that flag.
- `plot_heatmap.py --louvain` is a stub — don't rely on it.

## Conventions (same as repo-wide)

- Never overwrite existing results — write to a new versioned output
  directory rather than clobbering one that already has output.
- Don't invent or guess file paths; verify with `ls`/`cat` before assuming,
  especially given the stale-path caveat above.
