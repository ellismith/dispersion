"""
plot_heatmap_ct_region.py

Flexible cell type x region summary heatmap for DGLM dispersion results.
--metric options:
  median_std_beta   : median standardized beta (beta/sqrt(bvar))
  median_beta       : median raw beta
  neg_log10_pval    : -log10(median pvalue)
  frac_sig          : fraction of genes significant at --qthresh
  net_direction     : (n_increasing - n_decreasing) / n_total sig genes

Usage:
  python plot_heatmap_ct_region.py --metric median_std_beta
  python plot_heatmap_ct_region.py --metric median_std_beta --covariate "sex (M vs F ref)"
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
parser.add_argument('--dglm_master', default='/scratch/easmit31/dispersion/dglm/disp_age__mean_full/master_disp_age_mean_full_globalfdr.tsv')
parser.add_argument('--outdir',      default='/scratch/easmit31/dispersion/dglm/figures/disp_age__mean_full')
parser.add_argument('--metric',      default='all',
                    choices=['median_std_beta','median_beta','neg_log10_pval',
                             'frac_sig','net_direction','all'])
parser.add_argument('--qthresh',     type=float, default=0.05,
                    help='FDR threshold for significance asterisks')
parser.add_argument('--frac_thresh', type=float, default=0.10,
                    help='Min fraction of sig genes required for asterisk')
parser.add_argument('--no_sig', action='store_true',
                    help='Suppress significance asterisks entirely')
parser.add_argument('--covariate',   default='age',
                    help='Label for covariate being plotted, used in title/cbar (e.g. age, sex)')
parser.add_argument('--outfmt',      default='png')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

print('Loading results...')
df = pd.read_csv(args.dglm_master, sep='\t')
df = df[df['beta'].notna() & df['pvalue'].notna() & (df['beta'].abs() <= 100)]
df['std_beta'] = df['beta'] / df['bvar'].apply(lambda x: np.sqrt(max(x, 1e-10)))
df['sig']      = df['qvalue'] < args.qthresh
print(f'Loaded {len(df)} tests')

METRICS = {
    'median_std_beta': dict(
        label=f'Median standardized beta [{args.covariate}]\n(beta/sqrt(bvar))',
        fname='heatmap_median_std_beta',
        cmap='RdBu_r', diverging=True,
        fn=lambda sub: sub['std_beta'].median()
    ),
    'median_beta': dict(
        label=f'Median raw beta [{args.covariate}]\n(DGLM dispersion coefficient)',
        fname='heatmap_median_beta',
        cmap='RdBu_r', diverging=True,
        fn=lambda sub: sub['beta'].median()
    ),
    'neg_log10_pval': dict(
        label=f'-log10(median p-value) [{args.covariate}]',
        fname='heatmap_neg_log10_pval',
        cmap='Reds', diverging=False,
        fn=lambda sub: -np.log10(sub['pvalue'].median() + 1e-300)
    ),
    'frac_sig': dict(
        label=f'Fraction of genes significant [{args.covariate}]\n(q<{args.qthresh})',
        fname=f'heatmap_frac_sig_q{args.qthresh}',
        cmap='Reds', diverging=False,
        fn=lambda sub: sub['sig'].mean()
    ),
    'net_direction': dict(
        label=f'Net direction among sig genes [{args.covariate}]\n(n_inc - n_dec) / n_sig (q<{args.qthresh})',
        fname=f'heatmap_net_direction_q{args.qthresh}',
        cmap='RdBu_r', diverging=True,
        fn=lambda sub: (
            (sub.loc[sub['sig'], 'beta'] > 0).sum() -
            (sub.loc[sub['sig'], 'beta'] < 0).sum()
        ) / max(sub['sig'].sum(), 1)
    ),
}

metrics_to_run = list(METRICS.keys()) if args.metric == 'all' else [args.metric]

def make_heatmap(mat, sig_mat, title, fname, cmap, diverging):
    row_labels = [CT_LABELS.get(ct, ct) for ct in CELL_TYPES]
    n_rows, n_cols = mat.shape
    fig, ax = plt.subplots(figsize=(n_cols * 0.85 + 1.5, n_rows * 0.7 + 1.5))

    vals = mat.values.astype(float)
    if diverging:
        vmax = np.nanpercentile(np.abs(vals), 95)
        vmax = max(vmax, 0.01)
        im = ax.imshow(vals, cmap=cmap, vmin=-vmax, vmax=vmax, aspect='auto')
    else:
        vmax = np.nanpercentile(vals, 95)
        im = ax.imshow(vals, cmap=cmap, vmin=0, vmax=vmax, aspect='auto')

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
    sig_note = '' if args.no_sig else f'\n* = >{int(args.frac_thresh*100)}% genes q<{args.qthresh}'
    ax.set_title(f'{title}{sig_note}',
                 fontsize=TITLE_FS, pad=10)

    cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.ax.tick_params(labelsize=CBAR_FS)
    cbar.set_label(title.split('\n')[0], fontsize=CBAR_FS)

    plt.tight_layout()
    out = os.path.join(args.outdir, f'{fname}_final.{args.outfmt}')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out}')

for metric in metrics_to_run:
    m = METRICS[metric]
    print(f'Computing {metric}...')
    mat     = pd.DataFrame(np.nan,  index=CELL_TYPES, columns=REGIONS)
    sig_mat = pd.DataFrame(False,   index=CELL_TYPES, columns=REGIONS)
    for ct in CELL_TYPES:
        for region in REGIONS:
            sub = df[(df['cell_type'] == ct) & (df['region'] == region)]
            if len(sub) == 0:
                continue
            mat.loc[ct, region]     = m['fn'](sub)
            sig_mat.loc[ct, region] = sub['sig'].mean() > args.frac_thresh
    make_heatmap(mat, sig_mat, m['label'], m['fname'], m['cmap'], m['diverging'])

print('done.')
