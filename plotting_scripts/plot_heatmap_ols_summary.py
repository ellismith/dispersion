"""
plot_heatmap_ols_summary.py

Cell type x region summary heatmap for gene_variance OLS results.
Color = median age slope per cell type x region.
Asterisk = >10% of genes with raw p<0.05 (no FDR, OLS global FDR is too conservative).

Usage:
  python plot_heatmap_ols_summary.py
  python plot_heatmap_ols_summary.py --pthresh 0.05 --frac_thresh 0.10
"""
import argparse
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_style import CELL_TYPES, CT_LABELS, REGIONS, TICK_FS, LABEL_FS, CBAR_FS, TITLE_FS

parser = argparse.ArgumentParser()
parser.add_argument('--master',      default='/scratch/easmit31/dispersion/gene_variance/results_log/master_between.tsv')
parser.add_argument('--outdir',      default='/scratch/easmit31/dispersion/gene_variance/figures')
parser.add_argument('--pthresh',     type=float, default=0.05,
                    help='Raw p-value threshold for significance (no FDR)')
parser.add_argument('--frac_thresh', type=float, default=0.10,
                    help='Fraction of genes at pthresh required for asterisk')
parser.add_argument('--no_sig', action='store_true',
                    help='Suppress significance asterisks entirely')
parser.add_argument('--outfmt',      default='png')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

print('Loading OLS results...')
df = pd.read_csv(args.master, sep='\t')
df = df[df['age_slope'].notna() & df['pvalue'].notna()]
df['sig'] = df['pvalue'] < args.pthresh
print(f'Loaded {len(df)} tests')

mat     = pd.DataFrame(np.nan,  index=CELL_TYPES, columns=REGIONS)
sig_mat = pd.DataFrame(False,   index=CELL_TYPES, columns=REGIONS)

for ct in CELL_TYPES:
    for region in REGIONS:
        sub = df[(df['cell_type'] == ct) & (df['region'] == region)]
        if len(sub) == 0:
            continue
        mat.loc[ct, region]     = sub['age_slope'].median()
        sig_mat.loc[ct, region] = sub['sig'].mean() > args.frac_thresh

row_labels = [CT_LABELS.get(ct, ct) for ct in CELL_TYPES]
n_rows, n_cols = mat.shape

vals = mat.values.astype(float)
vmax = np.nanpercentile(np.abs(vals), 95)
vmax = max(vmax, 1e-4)

fig, ax = plt.subplots(figsize=(n_cols * 0.85 + 1.5, n_rows * 0.7 + 1.5))
im = ax.imshow(vals, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')

if not args.no_sig:
    for i in range(n_rows):
        for j in range(n_cols):
            if sig_mat.iloc[i, j]:
                ax.text(j, i, '*', ha='center', va='center',
                        fontsize=14, color='black', fontweight='bold')

ax.set_xticks(range(n_cols))
ax.set_xticklabels(REGIONS, rotation=45, ha='right', fontsize=TICK_FS)
ax.set_yticks(range(n_rows))
ax.set_yticklabels(row_labels, rotation=0, fontsize=TICK_FS)
ax.set_xlabel('Region', fontsize=LABEL_FS)
ax.set_ylabel('Cell type', fontsize=LABEL_FS)
sig_note = '' if args.no_sig else f'\n* = >{int(args.frac_thresh*100)}% genes p<{args.pthresh}'
ax.set_title(f'Gene variance OLS\nMedian age slope (between-individual){sig_note}',
             fontsize=TITLE_FS, pad=10)

cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
cbar.ax.tick_params(labelsize=CBAR_FS)
cbar.set_label('Median age slope', fontsize=CBAR_FS)

plt.tight_layout()
out = os.path.join(args.outdir, f'heatmap_ols_median_slope_final.{args.outfmt}')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {out}')
