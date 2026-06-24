"""
plot_mean_vs_variability.py

Scatter of mean expression age effect (x) vs dispersion age effect (y).
Each point = one gene x cell type x region.
Shows relationship between mean and variability changes with age.

Requires gene_variance master TSV (has both mean slope and variability slope)
or DGLM master TSV paired with a mean effect TSV.

For gene_variance: mean effect = OLS slope of mean ~ age, variability = OLS slope of |residuals| ~ age
For DGLM: mean effect from mean submodel, dispersion effect from dispersion submodel

Usage:
  python plot_mean_vs_variability.py --pipeline gv
  python plot_mean_vs_variability.py --pipeline gv --cell_type microglia
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import os

parser = argparse.ArgumentParser()
parser.add_argument('--pipeline',    choices=['gv','dglm'], required=True)
parser.add_argument('--gv_master',   default='/scratch/easmit31/variability/gene_variance/results_log/master_between.tsv')
parser.add_argument('--dglm_master', default='/scratch/easmit31/variability/dglm/checkpoints/master_dglm_globalfdr.tsv')
parser.add_argument('--outdir',      default='/scratch/easmit31/variability/plotting_scripts/figures')
parser.add_argument('--qthresh',     type=float, default=0.05)
parser.add_argument('--cell_type',   default=None)
parser.add_argument('--region',      default=None)
parser.add_argument('--outfmt',      default='png')
parser.add_argument('--top_n_label', type=int, default=8)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

if args.pipeline == 'gv':
    df = pd.read_csv(args.gv_master, sep='\t')
    df = df[df['age_slope'].notna() & df['pvalue'].notna()]
    # gene_variance master has mean_expr but not a separate mean age slope
    # use age_slope as variability effect, mean_expr as x proxy
    # check columns
    print('Columns:', df.columns.tolist())
    if 'mean_age_slope' in df.columns:
        df['x_effect'] = df['mean_age_slope']
        x_label = 'Mean expression age slope'
    else:
        # fallback: use mean_expr (average expression level)
        df['x_effect'] = df['mean_expr']
        x_label = 'Mean expression level'
    df['y_effect']  = df['age_slope']
    y_label = 'Variability age slope'
    title_prefix = 'gene_variance'
else:
    df = pd.read_csv(args.dglm_master, sep='\t')
    df = df[df['beta'].notna() & df['pvalue'].notna() & (df['beta'].abs() <= 100)]
    if 'mean_beta' in df.columns:
        df['x_effect'] = df['mean_beta']
        x_label = 'Mean expression age effect (beta)'
    else:
        df['x_effect'] = df['beta']
        x_label = 'Dispersion age effect (beta)'
    df['y_effect']  = df['beta'] / np.sqrt(df['bvar'].clip(lower=1e-10))
    y_label = 'Standardized dispersion age effect'
    title_prefix = 'DGLM'

df['sig'] = df['qvalue'] < args.qthresh

if args.cell_type:
    df = df[df['cell_type'] == args.cell_type]
if args.region:
    df = df[df['region'] == args.region]

if len(df) == 0:
    print('No data after filtering')
    exit(0)

print(f'n genes: {len(df)}, sig: {df["sig"].sum()}')

# color: red=sig increasing, blue=sig decreasing, gray=ns
colors = np.where(~df['sig'], '#AAAAAA',
         np.where(df['y_effect'] > 0, '#d73027', '#4575b4'))

fig, ax = plt.subplots(figsize=(7, 6))

# nonsig first
mask_ns = ~df['sig']
ax.scatter(df.loc[mask_ns, 'x_effect'], df.loc[mask_ns, 'y_effect'],
           c='#AAAAAA', s=3, alpha=0.3, linewidths=0, zorder=1)

# sig on top
mask_s = df['sig']
ax.scatter(df.loc[mask_s, 'x_effect'], df.loc[mask_s, 'y_effect'],
           c=colors[mask_s], s=6, alpha=0.7, linewidths=0, zorder=2)

ax.axhline(0, color='black', linewidth=0.5, alpha=0.4)
ax.axvline(0, color='black', linewidth=0.5, alpha=0.4)

# label top genes
if args.top_n_label > 0 and df['sig'].any():
    top = df[df['sig']].nlargest(args.top_n_label, 'y_effect')
    bot = df[df['sig']].nsmallest(args.top_n_label, 'y_effect')
    for _, row in pd.concat([top, bot]).iterrows():
        sym = str(row.get('symbol', ''))
        if not sym or sym == 'nan': continue
        ax.annotate(sym, (row['x_effect'], row['y_effect']),
                    fontsize=6, ha='left', va='bottom',
                    xytext=(3,3), textcoords='offset points',
                    color='#222222',
                    path_effects=[pe.withStroke(linewidth=2, foreground='white')])

ax.set_xlabel(x_label, fontsize=11)
ax.set_ylabel(y_label, fontsize=11)

parts = [title_prefix]
if args.cell_type: parts.append(args.cell_type)
if args.region:    parts.append(args.region)
parts.append(f'q<{args.qthresh}')
ax.set_title(' — '.join(parts), fontsize=10)

from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#d73027', markersize=6, label='Sig increasing variability'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#4575b4', markersize=6, label='Sig decreasing variability'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#AAAAAA', markersize=6, label='Not significant'),
], fontsize=8, framealpha=0.8)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

suffix = ''
if args.cell_type: suffix += f'_{args.cell_type}'
if args.region:    suffix += f'_{args.region}'
out = os.path.join(args.outdir, f'mean_vs_variability_{args.pipeline}{suffix}.{args.outfmt}')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {out}')
print('done.')
