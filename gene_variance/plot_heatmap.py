"""
plot_heatmap.py

Summary heatmaps of gene variance results across cell types x regions.
Produces 4 heatmaps: between positive, between negative, within positive, within negative.
Values = number of significant genes at q<0.05.

Optional louvain mode: rows = louvain clusters within a cell type instead of cell types.
Requires running fdr_correct.py first to produce master_between.tsv and master_within.tsv.

Usage:
    # default: all cell types x regions
    python plot_heatmap.py --indir PATH --outdir PATH

    # specific cell type only
    python plot_heatmap.py --indir PATH --outdir PATH --cell_type microglia

    # louvain resolution (must specify cell type)
    python plot_heatmap.py --indir PATH --outdir PATH --cell_type microglia --louvain
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ── args ──────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument('--indir',     required=True,  help='directory containing master TSVs')
parser.add_argument('--outdir',    required=True,  help='output directory for plots')
parser.add_argument('--cell_type', default=None,   help='optional: restrict to one cell type')
parser.add_argument('--louvain',   action='store_true', help='louvain resolution (requires --cell_type)')
parser.add_argument('--qthresh',   type=float, default=0.05)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

if args.louvain and args.cell_type is None:
    raise ValueError("--louvain requires --cell_type")

REGIONS = ['ACC', 'CN', 'dlPFC', 'EC', 'HIP', 'IPP', 'lCb', 'M1', 'MB', 'mdTN', 'NAc']

def make_heatmap(pivot, title, outpath, cmap, vmax):
    # only keep regions that exist as columns
    cols = [r for r in REGIONS if r in pivot.columns]
    pivot = pivot[cols]

    fig, ax = plt.subplots(figsize=(len(cols) * 0.8 + 2, len(pivot) * 0.5 + 2))
    sns.heatmap(
        pivot,
        ax=ax,
        cmap=cmap,
        vmin=0,
        vmax=vmax,
        annot=True,
        fmt='.0f',
        linewidths=0.5,
        linecolor='lightgray',
        cbar_kws={'label': 'n significant genes'}
    )
    ax.set_title(title, fontsize=13, pad=10)
    ax.set_xlabel('Region', fontsize=11)
    ax.set_ylabel('', fontsize=11)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  saved: {outpath}")

for mode in ['between', 'within']:
    print(f"\nprocessing {mode}")
    master = pd.read_csv(os.path.join(args.indir, f'master_{mode}.tsv'), sep='\t')

    # filter to cell type if specified
    if args.cell_type is not None:
        master = master[master['cell_type'] == args.cell_type]
        print(f"  filtered to {args.cell_type}: {len(master)} rows")

    # louvain mode: use louvain column as row grouper
    # need to add louvain to master — it's in the per-gene TSVs but not in master yet
    # for now, louvain mode reads directly from per-cell-type TSVs
    if args.louvain:
        raise NotImplementedError("louvain mode requires louvain column — add --louvain support after confirming louvain is in TSVs")

    sig = master[master['qvalue'] < args.qthresh]

    row_col = 'cell_type'

    # positive
    pos = sig[sig['age_slope'] > 0].groupby([row_col, 'region']).size().reset_index(name='n')
    pivot_pos = pos.pivot(index=row_col, columns='region', values='n').fillna(0)

    # negative
    neg = sig[sig['age_slope'] < 0].groupby([row_col, 'region']).size().reset_index(name='n')
    pivot_neg = neg.pivot(index=row_col, columns='region', values='n').fillna(0)

    vmax = max(pivot_pos.values.max(), pivot_neg.values.max())
    label = args.cell_type if args.cell_type else 'all_celltypes'

    make_heatmap(
        pivot=pivot_pos,
        title=f'{mode} variance increases with age (q<{args.qthresh})',
        outpath=os.path.join(args.outdir, f'heatmap_{mode}_positive_{label}.png'),
        cmap='Reds',
        vmax=vmax
    )

    make_heatmap(
        pivot=pivot_neg,
        title=f'{mode} variance decreases with age (q<{args.qthresh})',
        outpath=os.path.join(args.outdir, f'heatmap_{mode}_negative_{label}.png'),
        cmap='Blues',
        vmax=vmax
    )

print("\ndone.")
