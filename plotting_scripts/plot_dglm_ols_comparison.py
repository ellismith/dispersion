"""
plot_dglm_ols_comparison.py

Comparison plots between DGLM and OLS gene variance results.

Produces per cell type x region:
  1. Overlap heatmap — fraction of sig genes shared between DGLM and OLS
  2. Direction consistency heatmap — % of shared sig genes agreeing on direction
  3. Rank correlation heatmap — Spearman r of effect sizes across all genes

Usage:
  python plot_dglm_ols_comparison.py
  python plot_dglm_ols_comparison.py --dglm_q 0.05 --ols_p 0.05 --frac_thresh 0.10
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
from scipy.stats import spearmanr
from plot_style import CELL_TYPES, CT_LABELS, REGIONS, TICK_FS, LABEL_FS, TITLE_FS, CBAR_FS

parser = argparse.ArgumentParser()
parser.add_argument('--dglm_master', default='/scratch/easmit31/dispersion/dglm/disp_age__mean_full/master_disp_age_mean_full_globalfdr.tsv')
parser.add_argument('--ols_master',  default='/scratch/easmit31/dispersion/gene_variance/results_log/master_between.tsv')
parser.add_argument('--outdir',      default='/scratch/easmit31/dispersion/plotting_scripts/figures/dglm_ols_comparison')
parser.add_argument('--dglm_q',      type=float, default=0.05,  help='DGLM significance threshold (qvalue)')
parser.add_argument('--ols_p',       type=float, default=0.05,  help='OLS significance threshold (raw pvalue)')
parser.add_argument('--frac_thresh', type=float, default=0.10,  help='Min fraction sig genes for asterisk')
parser.add_argument('--outfmt',      default='png')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

print('Loading DGLM...')
dglm = pd.read_csv(args.dglm_master, sep='\t')
dglm = dglm[dglm['beta'].notna() & dglm['pvalue'].notna() & (dglm['beta'].abs() <= 100)]
dglm['std_beta'] = dglm['beta'] / dglm['bvar'].apply(lambda x: np.sqrt(max(x, 1e-10)))
dglm['sig'] = dglm['qvalue'] < args.dglm_q
print(f'  {len(dglm)} tests')

print('Loading OLS...')
ols = pd.read_csv(args.ols_master, sep='\t')
ols = ols[ols['age_slope'].notna() & ols['pvalue'].notna()]
ols['sig'] = ols['pvalue'] < args.ols_p
ols = ols.rename(columns={'human_symbol': 'symbol'})
print(f'  {len(ols)} tests')

# ── Merge on ensembl_id + cell_type + region ──────────────────────────────────
merged = pd.merge(
    dglm[['ensembl_id', 'symbol', 'cell_type', 'region', 'std_beta', 'beta', 'sig']],
    ols[['ensembl_id', 'cell_type', 'region', 'age_slope', 'sig']],
    on=['ensembl_id', 'cell_type', 'region'],
    suffixes=('_dglm', '_ols')
)
print(f'Merged: {len(merged)} genes in common')

# ── Build matrices ────────────────────────────────────────────────────────────
overlap_mat   = pd.DataFrame(np.nan, index=CELL_TYPES, columns=REGIONS)
direction_mat = pd.DataFrame(np.nan, index=CELL_TYPES, columns=REGIONS)
spearman_mat  = pd.DataFrame(np.nan, index=CELL_TYPES, columns=REGIONS)
overlap_sig   = pd.DataFrame(False,  index=CELL_TYPES, columns=REGIONS)
direction_sig  = pd.DataFrame(False,  index=CELL_TYPES, columns=REGIONS)
spearman_sig  = pd.DataFrame(False,  index=CELL_TYPES, columns=REGIONS)

for ct in CELL_TYPES:
    for region in REGIONS:
        sub = merged[(merged['cell_type'] == ct) & (merged['region'] == region)]
        if len(sub) < 10:
            continue

        # 1. overlap: jaccard of sig gene sets
        sig_dglm = set(sub.loc[sub['sig_dglm'], 'ensembl_id'])
        sig_ols  = set(sub.loc[sub['sig_ols'],  'ensembl_id'])
        union    = sig_dglm | sig_ols
        if len(union) > 0:
            overlap_mat.loc[ct, region] = len(sig_dglm & sig_ols) / len(union)
            # sig if both methods have >frac_thresh sig genes
            overlap_sig.loc[ct, region] = (
                sub['sig_dglm'].mean() > args.frac_thresh and
                sub['sig_ols'].mean() > args.frac_thresh
            )

        # 2. direction consistency among genes sig in both
        both_sig = sub[sub['sig_dglm'] & sub['sig_ols']]
        if len(both_sig) >= 5:
            agree = ((both_sig['beta'] > 0) == (both_sig['age_slope'] > 0)).mean()
            direction_mat.loc[ct, region] = agree
            direction_sig.loc[ct, region] = agree > 0.8

        # 3. spearman r of effect sizes across all genes
        if len(sub) >= 10:
            r, p = spearmanr(sub['std_beta'], sub['age_slope'])
            spearman_mat.loc[ct, region] = r
            spearman_sig.loc[ct, region] = p < 0.05

# ── Plot function ─────────────────────────────────────────────────────────────
def make_heatmap(mat, sig, title, fname, cmap, vmin, vmax, cbar_label):
    row_labels = [CT_LABELS.get(ct, ct) for ct in CELL_TYPES]
    n_rows, n_cols = mat.shape
    fig, ax = plt.subplots(figsize=(n_cols * 0.85 + 1.5, n_rows * 0.7 + 1.5))

    vals = mat.values.astype(float)
    im = ax.imshow(vals, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')

    for i in range(n_rows):
        for j in range(n_cols):
            if sig.iloc[i, j]:
                ax.text(j, i, '*', ha='center', va='center',
                        fontsize=14, color='black', fontweight='bold')

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(REGIONS, rotation=45, ha='right', fontsize=TICK_FS)
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

# ── Plot 1: overlap (Jaccard) ─────────────────────────────────────────────────
make_heatmap(
    overlap_mat, overlap_sig,
    title=f'DGLM vs OLS — sig gene overlap (Jaccard)\n(DGLM q<{args.dglm_q}, OLS p<{args.ols_p})\n* = both methods >{int(args.frac_thresh*100)}% sig',
    fname='overlap_jaccard',
    cmap='Reds', vmin=0, vmax=0.3,
    cbar_label='Jaccard index'
)

# ── Plot 2: direction consistency ─────────────────────────────────────────────
make_heatmap(
    direction_mat, direction_sig,
    title=f'DGLM vs OLS — direction consistency\namong genes sig in both\n* = >80% agreement',
    fname='direction_consistency',
    cmap='RdYlGn', vmin=0.5, vmax=1.0,
    cbar_label='Fraction agreeing on direction'
)

# ── Plot 3: rank correlation ──────────────────────────────────────────────────
make_heatmap(
    spearman_mat, spearman_sig,
    title=f'DGLM vs OLS — effect size rank correlation\n(Spearman r, all genes)\n* = p<0.05',
    fname='rank_correlation',
    cmap='RdBu_r', vmin=-0.5, vmax=0.5,
    cbar_label='Spearman r'
)

print('done.')
