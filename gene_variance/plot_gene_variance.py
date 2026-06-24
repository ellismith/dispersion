"""
plot_gene_variance.py

Produces scatter plots for gene variance results.
x = log10(mean expression across animals)
y = age slope from regression (positive = variance increases with age)
color = red (sig positive), blue (sig negative), gray (ns)

Two modes:
  1. Single region (--region)
  2. All regions overlaid (--all_regions)

Usage:
    python plot_gene_variance.py --h5ad PATH --indir PATH --outdir PATH
                                 --cell_type STR [--region STR] [--all_regions]
"""

import argparse
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser()
parser.add_argument('--h5ad',                 required=True)
parser.add_argument('--indir',                required=True)
parser.add_argument('--outdir',               required=True)
parser.add_argument('--cell_type',            required=True)
parser.add_argument('--region',               default=None)
parser.add_argument('--all_regions',          action='store_true')
parser.add_argument('--min_age',              type=float, default=1.0)
parser.add_argument('--min_animals',          type=int,   default=10)
parser.add_argument('--qthresh',              type=float, default=0.05)
args = parser.parse_args()

if not args.region and not args.all_regions:
    raise ValueError("must specify --region or --all_regions")

os.makedirs(args.outdir, exist_ok=True)

REGIONS = ['ACC', 'CN', 'dlPFC', 'EC', 'HIP', 'IPP', 'lCb', 'M1', 'MB', 'mdTN', 'NAc']

# ── load master TSVs for qvalues ──────────────────────────────────────────
master_between = pd.read_csv(os.path.join(args.indir, 'master_between.tsv'), sep='\t')
master_within  = pd.read_csv(os.path.join(args.indir, 'master_within.tsv'),  sep='\t')

def get_reg(master, cell_type, region):
    sub = master[(master['cell_type'] == cell_type) & (master['region'] == region)]
    return sub.set_index('ensembl_id')[['age_slope', 'qvalue', 'mean_expr']]

# ── scatter helper ────────────────────────────────────────────────────────
def scatter(ax, reg_df, qthresh):
    if reg_df.empty:
        return 0, 0, 0

    x     = np.log10(reg_df['mean_expr'].values + 1e-10)
    slope = reg_df['age_slope'].values
    qval  = reg_df['qvalue'].values

    sig_pos = (qval < qthresh) & (slope > 0)
    sig_neg = (qval < qthresh) & (slope < 0)
    ns      = ~(sig_pos | sig_neg)

    ax.scatter(x[ns],      slope[ns],      c='lightgray', s=2, alpha=0.4, linewidths=0)
    ax.scatter(x[sig_neg], slope[sig_neg], c='blue',      s=4, alpha=0.8, linewidths=0)
    ax.scatter(x[sig_pos], slope[sig_pos], c='red',       s=4, alpha=0.8, linewidths=0)

    return int(sig_pos.sum()), int(sig_neg.sum()), int(ns.sum())

# ── finalize plot ─────────────────────────────────────────────────────────
def finalize_plot(ax, n_pos, n_neg, n_ns, title, outpath, qthresh):
    ax.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.scatter([], [], c='red',       s=15, label=f'increases with age (n={n_pos})')
    ax.scatter([], [], c='blue',      s=15, label=f'decreases with age (n={n_neg})')
    ax.scatter([], [], c='lightgray', s=15, label=f'ns (n={n_ns})')
    ax.set_xlabel('log10(mean expression)', fontsize=11)
    ax.set_ylabel('age slope (variance ~ age)', fontsize=11)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=8, markerscale=2)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  saved: {outpath}")

# ── determine regions ─────────────────────────────────────────────────────
regions = REGIONS if args.all_regions else [args.region]

# ── main loop ─────────────────────────────────────────────────────────────
for mode in ['between', 'within']:
    master = master_between if mode == 'between' else master_within

    if args.all_regions:
        fig, ax = plt.subplots(figsize=(6, 5))
        total_pos, total_neg, total_ns = 0, 0, 0

    for region in regions:
        tsv = os.path.join(args.indir, f'{args.cell_type}_{region}_{mode}.tsv')
        if not os.path.exists(tsv):
            continue
        try:
            df = pd.read_csv(tsv, sep='\t')
            if df.empty:
                continue
        except Exception:
            continue

        reg_df = get_reg(master, args.cell_type, region)
        if reg_df.empty:
            continue

        if not args.all_regions:
            fig, ax = plt.subplots(figsize=(6, 5))
            n_pos, n_neg, n_ns = scatter(ax, reg_df, args.qthresh)
            finalize_plot(
                ax, n_pos, n_neg, n_ns,
                title   = f'{args.cell_type} {region} {mode}-individual (q<{args.qthresh})',
                outpath = os.path.join(args.outdir, f'{args.cell_type}_{region}_{mode}.png'),
                qthresh = args.qthresh
            )
        else:
            n_pos, n_neg, n_ns = scatter(ax, reg_df, args.qthresh)
            total_pos += n_pos
            total_neg += n_neg
            total_ns  += n_ns

    if args.all_regions:
        finalize_plot(
            ax, total_pos, total_neg, total_ns,
            title   = f'{args.cell_type} all regions {mode}-individual (q<{args.qthresh})',
            outpath = os.path.join(args.outdir, f'{args.cell_type}_allregions_{mode}.png'),
            qthresh = args.qthresh
        )

print("done.")
