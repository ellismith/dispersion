"""
plot_raw_dispersion_heatmap.py

Cell type x region summary heatmaps of RAW transcriptional variability --
computed directly from pseudobulk expression, independent of any DGLM/age
model fit. This is deliberately NOT derived from dglm_model.R's output:
the DGLM only ever fits dispersion AS A FUNCTION of age
(log(dispersion) = intercept + beta*age), and the saved arrays only
contain that age coefficient, not a standalone baseline dispersion value
-- so "the dispersion values themselves" (as opposed to the age EFFECT on
them, which plot_heatmap_ct_region.py already covers) has to come from
the raw data, not the model.

Since this doesn't depend on any model fit, it only needs to be computed
ONCE -- not separately per model version (V2/V3 share the same underlying
pseudobulk data).

Three metrics, one heatmap each:
  cv        : coefficient of variation (SD/mean) per gene, median across
              genes -- fully normalizes for expression level, so cell
              types/regions aren't just ranked by how highly expressed
              their genes happen to be
  fano      : Fano factor (variance/mean) per gene, median across genes --
              standard for count-like data, partial normalization
  log_var   : log10(variance) per gene, median across genes -- no
              normalization, reflects absolute magnitude (will mostly
              track expression level)

Genes are filtered to the same autosome+X keep-list used throughout the
rest of the pipeline this session (autosome_x_genes.csv), for consistency
-- otherwise this would be the one place Y/MT/scaffold genes sneak back
into a summary.

Discovers which {cell_type}_{region}_pseudobulk.csv files actually exist
in --checkpoints rather than assuming a full cell-type x region grid (same
pattern dglm_mashr.R uses for combined mode) -- not every cell type has
usable data in every region.

Usage:
  python plot_raw_dispersion_heatmap.py \
      --checkpoints /scratch/easmit31/dispersion/dglm/disp_age__mean_full \
      --outdir /scratch/easmit31/dispersion/dglm/disp_age__mean_full/figures
"""
import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import CELL_TYPES, CT_LABELS, REGIONS, TICK_FS, LABEL_FS, CBAR_FS, TITLE_FS

parser = argparse.ArgumentParser()
parser.add_argument('--checkpoints', default='/scratch/easmit31/dispersion/dglm/disp_age__mean_full',
                    help='dir containing {cell_type}_{region}_pseudobulk.csv files')
parser.add_argument('--outdir', default='/scratch/easmit31/dispersion/dglm/disp_age__mean_full/figures')
parser.add_argument('--gene_keep_list', default='/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv',
                    help='autosome+X gene keep-list, same filter used throughout the rest of the pipeline. Pass "" to skip.')
parser.add_argument('--min_mean', type=float, default=1e-6,
                    help='genes with mean expression below this are excluded from cv/fano for that cell_type x region (avoids divide-by-near-zero blowing up the median)')
parser.add_argument('--scale_factor', type=float, default=1e6,
                    help='rescale pseudobulk expression (raw proportions, ~1e-3 to 1e-2) by this factor before computing stats -- default 1e6 treats them as CPM (counts-per-million), the standard RNA-seq convention. CV is scale-invariant (unaffected). Fano and log_var scale with this -- 1e6 is what shifts log10(variance) from negative (log10 of a tiny raw proportion) into positive, human-readable territory, not an arbitrary correction.')
parser.add_argument('--outfmt', default='png')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

keep_genes = None
if args.gene_keep_list:
    keep_genes = set(pd.read_csv(args.gene_keep_list)['ensembl_gene_id'])
    print(f'Loaded gene keep-list: {len(keep_genes)} autosome+X genes')

# discover which cell_type x region pseudobulk files actually exist
pb_files = glob.glob(os.path.join(args.checkpoints, '*_pseudobulk.csv'))
print(f'Found {len(pb_files)} pseudobulk files in {args.checkpoints}')

# {cell_type}_{region}_pseudobulk.csv -- cell_type can itself contain
# underscores (e.g. GABAergic_neurons), so match against the known
# CELL_TYPES/REGIONS lists rather than a naive split
combos = []
for f in pb_files:
    base = os.path.basename(f)[:-len('_pseudobulk.csv')]
    matched = False
    for ct in CELL_TYPES:
        if base.startswith(ct + '_'):
            region = base[len(ct) + 1:]
            if region in REGIONS:
                combos.append((ct, region, f))
                matched = True
                break
    if not matched:
        print(f'  WARNING: could not parse cell_type/region from {base} -- skipping')
print(f'Parsed {len(combos)} valid cell_type x region combos')

METRICS = {
    'cv': dict(
        label='Median coefficient of variation\n(SD / mean, per gene)',
        fname='heatmap_raw_cv',
        cmap='viridis',
        fn=lambda mean, var: np.sqrt(var) / mean,
    ),
    'fano': dict(
        label='Median Fano factor\n(variance / mean, per gene)',
        fname='heatmap_raw_fano',
        cmap='viridis',
        fn=lambda mean, var: var / mean,
    ),
    'log_var': dict(
        label='Median log10(variance)\n(per gene, no normalization)',
        fname='heatmap_raw_log_var',
        cmap='viridis',
        fn=lambda mean, var: np.log10(var + 1e-300),
    ),
}

mats = {m: pd.DataFrame(np.nan, index=CELL_TYPES, columns=REGIONS) for m in METRICS}

for ct, region, f in combos:
    pb = pd.read_csv(f, index_col=0)
    if keep_genes is not None:
        pb = pb[pb.index.isin(keep_genes)]
    if pb.shape[0] == 0 or pb.shape[1] < 2:
        continue
    pb = pb * args.scale_factor
    mean = pb.mean(axis=1)
    var  = pb.var(axis=1, ddof=1)
    valid = mean > args.min_mean
    for metric, m in METRICS.items():
        vals = m['fn'](mean[valid], var[valid])
        vals = vals.replace([np.inf, -np.inf], np.nan).dropna()
        if len(vals) > 0:
            mats[metric].loc[ct, region] = vals.median()

def make_heatmap(mat, title, fname, cmap):
    row_labels = [CT_LABELS.get(ct, ct) for ct in CELL_TYPES]
    n_rows, n_cols = mat.shape
    fig, ax = plt.subplots(figsize=(n_cols * 0.85 + 1.5, n_rows * 0.7 + 1.5))

    vals = mat.values.astype(float)
    vmax = np.nanpercentile(vals, 95)
    vmin = np.nanpercentile(vals, 5)
    im = ax.imshow(vals, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(REGIONS, rotation=45, ha='right', fontsize=TICK_FS)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, rotation=0, fontsize=TICK_FS)
    ax.set_xlabel('Region', fontsize=LABEL_FS)
    ax.set_ylabel('Cell type', fontsize=LABEL_FS)
    ax.set_title(title, fontsize=TITLE_FS, pad=10)

    cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.ax.tick_params(labelsize=CBAR_FS)
    cbar.set_label(title.split('\n')[0], fontsize=CBAR_FS)

    plt.tight_layout()
    out = os.path.join(args.outdir, f'{fname}.{args.outfmt}')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')

for metric, m in METRICS.items():
    print(f'Plotting {metric}...')
    make_heatmap(mats[metric], m['label'], m['fname'], m['cmap'])

print('done.')
