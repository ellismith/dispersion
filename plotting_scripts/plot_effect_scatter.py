"""
plot_effect_scatter.py

Effect size correlation scatter between gene_variance and DGLM.
x = gene_variance standardized age slope
y = DGLM standardized beta (beta/sqrt(bvar))
One point per gene x cell type x region (matched on symbol|ct|region key).

Color coding:
  - red:    sig in both, increasing
  - blue:   sig in both, decreasing
  - orange: sig in gene_variance only
  - purple: sig in DGLM only
  - gray:   not sig in either

Usage:
  python plot_effect_scatter.py
  python plot_effect_scatter.py --cell_type microglia --region HIP
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import stats
import os

parser = argparse.ArgumentParser()
parser.add_argument('--gv_master',   default='/scratch/easmit31/variability/gene_variance/results_log/master_between.tsv')
parser.add_argument('--dglm_master', default='/scratch/easmit31/variability/dglm/checkpoints/master_dglm_globalfdr.tsv')
parser.add_argument('--outdir',      default='/scratch/easmit31/variability/plotting_scripts/figures')
parser.add_argument('--q_gv',        type=float, default=0.05)
parser.add_argument('--q_dglm',      type=float, default=0.05)
parser.add_argument('--cell_type',   default=None)
parser.add_argument('--region',      default=None)
parser.add_argument('--outfmt',      default='png')
parser.add_argument('--top_n_label', type=int, default=8)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# ── load ──────────────────────────────────────────────────────────────────
print('Loading gene_variance results')
gv = pd.read_csv(args.gv_master, sep='\t')
gv = gv[gv['age_slope'].notna() & gv['pvalue'].notna()]
gv['std_effect'] = gv['age_slope'] / gv['age_slope'].std()
gv['sig']        = gv['qvalue'] < args.q_gv
gv['symbol']     = gv['human_symbol'].fillna(gv['ensembl_id'])
gv['key']        = gv['symbol'] + '|' + gv['cell_type'] + '|' + gv['region']

print('Loading DGLM results')
dglm = pd.read_csv(args.dglm_master, sep='\t')
dglm = dglm[dglm['beta'].notna() & dglm['pvalue'].notna() & (dglm['beta'].abs() <= 100)]
dglm['std_effect'] = dglm['beta'] / np.sqrt(dglm['bvar'].clip(lower=1e-10))
dglm['sig']        = dglm['qvalue'] < args.q_dglm
dglm['key']        = dglm['symbol'] + '|' + dglm['cell_type'] + '|' + dglm['region']

# ── filter ────────────────────────────────────────────────────────────────
if args.cell_type:
    gv   = gv[gv['cell_type'] == args.cell_type]
    dglm = dglm[dglm['cell_type'] == args.cell_type]
if args.region:
    gv   = gv[gv['region'] == args.region]
    dglm = dglm[dglm['region'] == args.region]

# ── merge on key ──────────────────────────────────────────────────────────
merged = pd.merge(
    gv[['key','symbol','cell_type','region','std_effect','sig']],
    dglm[['key','std_effect','sig']],
    on='key',
    suffixes=('_gv','_dglm')
)
print(f'Matched gene x ct x region pairs: {len(merged)}')

if len(merged) == 0:
    print('No overlapping keys — check symbol matching')
    exit(0)

# ── color by sig status ───────────────────────────────────────────────────
def get_color(row):
    if row['sig_gv'] and row['sig_dglm']:
        return '#d73027' if row['std_effect_gv'] > 0 else '#4575b4'
    elif row['sig_gv']:
        return '#fc8d59'
    elif row['sig_dglm']:
        return '#9970ab'
    else:
        return '#AAAAAA'

merged['color'] = merged.apply(get_color, axis=1)
merged['alpha'] = merged['color'].apply(lambda c: 0.3 if c == '#AAAAAA' else 0.7)
merged['size']  = merged['color'].apply(lambda c: 4 if c == '#AAAAAA' else 8)

# ── spearman correlation ──────────────────────────────────────────────────
rho, pval = stats.spearmanr(merged['std_effect_gv'], merged['std_effect_dglm'])
print(f'Spearman r = {rho:.3f}, p = {pval:.2e}')

# ── plot ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 7))

for color, zorder in [('#AAAAAA',1),('#fc8d59',2),('#9970ab',2),('#4575b4',3),('#d73027',3)]:
    mask = merged['color'] == color
    if mask.sum() == 0:
        continue
    ax.scatter(
        merged.loc[mask, 'std_effect_gv'],
        merged.loc[mask, 'std_effect_dglm'],
        c=color,
        s=merged.loc[mask, 'size'],
        alpha=merged.loc[mask, 'alpha'].iloc[0],
        linewidths=0,
        zorder=zorder
    )

ax.axhline(0, color='black', linewidth=0.5, alpha=0.3)
ax.axvline(0, color='black', linewidth=0.5, alpha=0.3)

lim = max(abs(merged['std_effect_gv'].quantile(0.99)),
          abs(merged['std_effect_dglm'].quantile(0.99)))
ax.plot([-lim, lim], [-lim, lim], color='gray', linewidth=0.8,
        linestyle='--', alpha=0.5, zorder=0)

# label top genes sig in both
both_sig = merged[merged['sig_gv'] & merged['sig_dglm']]
if len(both_sig) > 0 and args.top_n_label > 0:
    top = both_sig.nlargest(args.top_n_label // 2, 'std_effect_gv')
    bot = both_sig.nsmallest(args.top_n_label // 2, 'std_effect_gv')
    for _, row in pd.concat([top, bot]).iterrows():
        sym = str(row.get('symbol', ''))
        if not sym or sym == 'nan':
            continue
        ax.annotate(sym, (row['std_effect_gv'], row['std_effect_dglm']),
                    fontsize=6, ha='left', va='bottom',
                    xytext=(3, 3), textcoords='offset points',
                    color='#222222',
                    path_effects=[pe.withStroke(linewidth=2, foreground='white')])

ax.set_xlabel('gene_variance standardized age slope', fontsize=11)
ax.set_ylabel('DGLM standardized beta', fontsize=11)

parts = ['Effect size correlation']
if args.cell_type:
    parts.append(args.cell_type)
if args.region:
    parts.append(args.region)
ax.set_title(' — '.join(parts) + f'\nSpearman r = {rho:.3f}', fontsize=10)

from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#d73027',
           markersize=6, label='Sig both (increasing)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#4575b4',
           markersize=6, label='Sig both (decreasing)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#fc8d59',
           markersize=6, label='Sig gene_variance only'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#9970ab',
           markersize=6, label='Sig DGLM only'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#AAAAAA',
           markersize=6, label='Not significant'),
], fontsize=8, framealpha=0.8)

# clip axes to 99th percentile
xlim = np.percentile(np.abs(merged['std_effect_gv']), 99) * 1.1
ylim = np.percentile(np.abs(merged['std_effect_dglm']), 99) * 1.1
ax.set_xlim(-xlim, xlim)
ax.set_ylim(-ylim, ylim)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

suffix = ''
if args.cell_type:
    suffix += f'_{args.cell_type}'
if args.region:
    suffix += f'_{args.region}'
out = os.path.join(args.outdir, f'effect_scatter{suffix}.{args.outfmt}')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {out}')
print('done.')
