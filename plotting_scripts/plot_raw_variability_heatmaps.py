#!/usr/bin/env python3
"""
plot_raw_variability_heatmaps.py

Raw between-animal variability level heatmaps:
  1. OLS: median log10(var_between) per CT x region
  2. Population centroid: mean distance to population centroid per CT x region

Usage:
  python plot_raw_variability_heatmaps.py
"""
import os, sys, argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from plot_style import (CELL_TYPES, CT_LABELS,
                        CENTROID_CT_ORDER, CENTROID_CT_LABELS,
                        REGIONS, TICK_FS, LABEL_FS, TITLE_FS, CBAR_FS)

parser = argparse.ArgumentParser()
parser.add_argument('--ols_master',       default='/scratch/easmit31/dispersion/gene_variance/results_log/master_between.tsv')
parser.add_argument('--centroid_summary', default='/scratch/easmit31/factor_analysis/sex_effects/population_centroid_mean_dist_summary.csv')
parser.add_argument('--outdir',           default='/scratch/easmit31/dispersion/plotting_scripts/figures/raw_variability')
parser.add_argument('--outfmt',           default='png')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

def make_heatmap(mat, row_labels, col_labels, title, fname, cbar_label):
    n_rows, n_cols = mat.shape
    fig, ax = plt.subplots(figsize=(n_cols * 0.85 + 1.5, n_rows * 0.7 + 1.5))
    vals = mat.values.astype(float)
    vmin = np.nanpercentile(vals, 5)
    vmax = np.nanpercentile(vals, 95)
    im = ax.imshow(vals, cmap='Reds', vmin=vmin, vmax=vmax, aspect='auto')
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_labels, rotation=45, ha='right', fontsize=TICK_FS)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, rotation=0, fontsize=TICK_FS)
    ax.set_xlabel('Region', fontsize=LABEL_FS)
    ax.set_ylabel('Cell type', fontsize=LABEL_FS)
    ax.set_title(title, fontsize=TITLE_FS, pad=10)
    cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.ax.tick_params(labelsize=CBAR_FS)
    cbar.set_label(cbar_label, fontsize=CBAR_FS)
    plt.tight_layout()
    out = os.path.join(args.outdir, f'{fname}_final.{args.outfmt}')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')

# ── 1. OLS: median log10(var_between) ────────────────────────────────────────
print('Loading OLS...')
ols = pd.read_csv(args.ols_master, sep='\t')
ols = ols[ols['var_between'].notna() & (ols['var_between'] > 0)]
ols['log10_var'] = np.log10(ols['var_between'])

mat = pd.DataFrame(np.nan, index=CELL_TYPES, columns=REGIONS)
for ct in CELL_TYPES:
    for region in REGIONS:
        sub = ols[(ols['cell_type'] == ct) & (ols['region'] == region)]
        if len(sub) == 0:
            continue
        mat.loc[ct, region] = sub['log10_var'].median()

row_labels = [CT_LABELS.get(ct, ct) for ct in CELL_TYPES]
make_heatmap(mat, row_labels, REGIONS,
             title='Between-individual gene expression variance\nMedian log10(var_between) per CT x region',
             fname='ols_raw_var_between',
             cbar_label='Median log10(var_between)')

# ── 2. Population centroid: mean distance ─────────────────────────────────────
print('Loading population centroid...')
cent = pd.read_csv(args.centroid_summary)

present_cts = [c for c in CENTROID_CT_ORDER if c in cent['cell_type'].unique()]
row_labels  = [CENTROID_CT_LABELS.get(ct, ct) for ct in present_cts]

mat = cent.pivot(index='cell_type', columns='region', values='mean_dist') \
          .reindex(index=present_cts, columns=REGIONS)

make_heatmap(mat, row_labels, REGIONS,
             title='Between-animal variability — population centroid distance\nMean across animals and louvains per CT x region',
             fname='population_centroid_raw_mean_dist',
             cbar_label='Mean distance to population centroid (PCA units)')

print('done.')
