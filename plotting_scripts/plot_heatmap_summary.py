"""
plot_heatmap_summary.py

Heatmap of median standardized effect size per cell type x region.
Matches lochNESS/centroid style: RdBu_r colormap, asterisks for significance.

Supports:
  - gene_variance between-individual pipeline (--pipeline gv)
  - DGLM pipeline (--pipeline dglm)
  - cell type level (default) or subtype level (--by_subtype, requires --cell_type)

Usage:
  python plot_heatmap_summary.py --pipeline gv
  python plot_heatmap_summary.py --pipeline dglm
  python plot_heatmap_summary.py --pipeline gv --by_subtype --cell_type microglia
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

parser = argparse.ArgumentParser()
parser.add_argument('--pipeline',    choices=['gv','dglm'], required=True)
parser.add_argument('--gv_master',   default='/scratch/easmit31/variability/gene_variance/results_log/master_between.tsv')
parser.add_argument('--dglm_master', default='/scratch/easmit31/variability/dglm/checkpoints/master_dglm_globalfdr.tsv')
parser.add_argument('--outdir',      default='/scratch/easmit31/variability/plotting_scripts/figures')
parser.add_argument('--qthresh',     type=float, default=0.05)
parser.add_argument('--by_subtype',  action='store_true')
parser.add_argument('--cell_type',   default=None, help='required for --by_subtype')
parser.add_argument('--outfmt',      default='png')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

CELL_TYPES = ['astrocytes','basket_cells','cerebellar_neurons','ependymal_cells',
              'GABAergic_neurons','glutamatergic_neurons','medium_spiny_neurons',
              'microglia','midbrain_neurons','opc','oligodendrocytes','vascular_cells']

REGIONS = ['ACC','CN','dlPFC','EC','HIP','IPP','lCb','M1','MB','mdTN','NAc']

CT_LABELS = {
    'astrocytes': 'AST', 'basket_cells': 'BC', 'cerebellar_neurons': 'CER',
    'ependymal_cells': 'EPEN', 'GABAergic_neurons': 'INH',
    'glutamatergic_neurons': 'EXC', 'medium_spiny_neurons': 'MSN',
    'microglia': 'MGL', 'midbrain_neurons': 'MBN', 'opc': 'OPC',
    'oligodendrocytes': 'OLIG', 'vascular_cells': 'VASC'
}

# ── load data ─────────────────────────────────────────────────────────────
if args.pipeline == 'gv':
    print('Loading gene_variance results')
    df = pd.read_csv(args.gv_master, sep='\t')
    df = df[df['age_slope'].notna() & df['qvalue'].notna()]
    df['effect']    = df['age_slope']
    df['std_effect'] = df['age_slope'] / df['age_slope'].std()
    label = 'gene_variance (OLS)'
    effect_label = 'Median age slope'
else:
    print('Loading DGLM results')
    df = pd.read_csv(args.dglm_master, sep='\t')
    df = df[df['beta'].notna() & df['qvalue'].notna() & (df['beta'].abs() <= 100)]
    df['std_effect'] = df['beta'] / df['bvar'].apply(lambda x: np.sqrt(x) if x > 0 else np.nan)
    df['effect']     = df['beta']
    label = 'DGLM (dispersion)'
    effect_label = 'Median standardized beta'

df['sig'] = df['qvalue'] < args.qthresh

# ── filter by cell type if by_subtype ────────────────────────────────────
if args.by_subtype:
    if args.cell_type is None:
        raise ValueError('--cell_type required with --by_subtype')
    df = df[df['cell_type'] == args.cell_type]
    row_var  = 'louvain' if 'louvain' in df.columns else 'cell_type'
    row_label = f'{args.cell_type} subtypes'
else:
    row_var   = 'cell_type'
    row_label = 'Cell type'

# ── compute median std_effect per row x region ────────────────────────────
rows    = CELL_TYPES if not args.by_subtype else sorted(df[row_var].unique())
row_labels = [CT_LABELS.get(r, r) for r in rows] if not args.by_subtype else rows

mat_effect = pd.DataFrame(np.nan, index=rows, columns=REGIONS)
mat_sig    = pd.DataFrame(False,  index=rows, columns=REGIONS)

for row in rows:
    for region in REGIONS:
        sub = df[(df[row_var] == row) & (df['region'] == region)]
        if len(sub) == 0:
            continue
        mat_effect.loc[row, region] = sub['std_effect'].median()
        # sig: majority of genes significant
        mat_sig.loc[row, region] = sub['sig'].mean() > 0.1

# ── plot ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(len(REGIONS)*0.85 + 1.5, len(rows)*0.7 + 1.5))

vmax = np.nanpercentile(np.abs(mat_effect.values), 95)
vmax = max(vmax, 0.1)

im = ax.imshow(mat_effect.values.astype(float),
               cmap='RdBu_r', vmin=-vmax, vmax=vmax,
               aspect='auto')

# asterisks for significant cells
for i, row in enumerate(rows):
    for j, region in enumerate(REGIONS):
        if mat_sig.loc[row, region]:
            ax.text(j, i, '*', ha='center', va='center',
                    fontsize=12, color='black', fontweight='bold')

ax.set_xticks(range(len(REGIONS)))
ax.set_xticklabels(REGIONS, fontsize=10)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels(row_labels, fontsize=10)

cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
cbar.set_label(effect_label, fontsize=9)

suffix = f'_{args.cell_type}' if args.by_subtype else ''
title  = f'{label}\nMedian standardized effect size per {row_label.lower()} × region (q<{args.qthresh})'
ax.set_title(title, fontsize=10, pad=10)

plt.tight_layout()
out = os.path.join(args.outdir, f'heatmap_summary_{args.pipeline}{suffix}.{args.outfmt}')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {out}')
print('done.')
