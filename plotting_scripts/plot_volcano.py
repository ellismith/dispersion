"""
plot_volcano.py

Volcano plot of age- or sex-associated dispersion changes.
x = standardized effect size, y = -log10(pvalue)
Significant genes colored red (increasing) or blue (decreasing), gray otherwise.

Best used with --cell_type and --region to get one point per gene.
Without filtering, each gene appears multiple times (once per ct x region).

Modes:
  --pipeline gv    : gene_variance between-individual
  --pipeline dglm  : DGLM dispersion
  --pipeline both  : only genes significant in both

Predictors:
  --predictor age  : age effect (default)
  --predictor sex  : sex effect (requires sex columns in master TSV)

Usage:
  python plot_volcano.py --pipeline dglm --cell_type microglia --region HIP
  python plot_volcano.py --pipeline gv --cell_type astrocytes --region ACC
  python plot_volcano.py --pipeline both --cell_type microglia --region HIP
  python plot_volcano.py --pipeline dglm --predictor sex --cell_type microglia --region HIP
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import os

parser = argparse.ArgumentParser()
parser.add_argument('--pipeline',    choices=['gv','dglm','both'], required=True)
parser.add_argument('--predictor',   choices=['age','sex'], default='age')
parser.add_argument('--gv_master',   default='/scratch/easmit31/variability/gene_variance/results_log/master_between.tsv')
parser.add_argument('--dglm_master', default='/scratch/easmit31/variability/dglm/checkpoints/master_dglm_globalfdr.tsv')
parser.add_argument('--outdir',      default='/scratch/easmit31/variability/plotting_scripts/figures')
parser.add_argument('--qthresh',     type=float, default=0.05)
parser.add_argument('--cell_type',   default=None)
parser.add_argument('--region',      default=None)
parser.add_argument('--outfmt',      default='png')
parser.add_argument('--top_n_label', type=int, default=10)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# ── load helpers ──────────────────────────────────────────────────────────
def load_gv(path, predictor='age'):
    df = pd.read_csv(path, sep='\t')
    df = df[df['pvalue'].notna() & df['age_slope'].notna()]
    if predictor == 'sex':
        slope_col = 'sex_slope'   if 'sex_slope'   in df.columns else 'age_slope'
        pval_col  = 'sex_pvalue'  if 'sex_pvalue'  in df.columns else 'pvalue'
        qval_col  = 'sex_qvalue'  if 'sex_qvalue'  in df.columns else 'qvalue'
    else:
        slope_col = 'age_slope'
        pval_col  = 'pvalue'
        qval_col  = 'qvalue'
    df['std_effect'] = df[slope_col] / df[slope_col].std()
    df['neg_log10p'] = -np.log10(df[pval_col].clip(lower=1e-300))
    df['qvalue']     = df[qval_col]
    df['pipeline']   = 'gv'
    return df

def load_dglm(path, predictor='age'):
    df = pd.read_csv(path, sep='\t')
    df = df[df['pvalue'].notna() & df['beta'].notna() & (df['beta'].abs() <= 100)]
    df['std_effect'] = df['beta'] / np.sqrt(df['bvar'].clip(lower=1e-10))
    df['neg_log10p'] = -np.log10(df['pvalue'].clip(lower=1e-300))
    df['pipeline']   = 'dglm'
    return df

# ── load ──────────────────────────────────────────────────────────────────
if args.pipeline == 'gv':
    df = load_gv(args.gv_master, args.predictor)
    title_prefix = 'gene_variance (OLS)'
elif args.pipeline == 'dglm':
    df = load_dglm(args.dglm_master, args.predictor)
    title_prefix = 'DGLM (dispersion)'
else:
    gv   = load_gv(args.gv_master, args.predictor)
    dglm = load_dglm(args.dglm_master, args.predictor)
    gv['key']   = gv['symbol']   + '|' + gv['cell_type']   + '|' + gv['region']
    dglm['key'] = dglm['symbol'] + '|' + dglm['cell_type'] + '|' + dglm['region']
    gv_sig_keys   = set(gv.loc[gv['qvalue']   < args.qthresh, 'key'])
    dglm_sig_keys = set(dglm.loc[dglm['qvalue'] < args.qthresh, 'key'])
    overlap_keys  = gv_sig_keys & dglm_sig_keys
    df = dglm[dglm['key'].isin(overlap_keys)].copy()
    title_prefix = f'Overlap (sig in both, n={len(overlap_keys)})'

# ── filter to cell type x region ──────────────────────────────────────────
if args.cell_type:
    df = df[df['cell_type'] == args.cell_type]
if args.region:
    df = df[df['region'] == args.region]

if len(df) == 0:
    print('No data after filtering — exiting')
    exit(0)

df['sig'] = df['qvalue'] < args.qthresh

print(f'Total genes: {len(df)}')
print(f'Sig: {df["sig"].sum()}')
print(f'Sig increasing: {(df["sig"] & (df["std_effect"] > 0)).sum()}')
print(f'Sig decreasing: {(df["sig"] & (df["std_effect"] < 0)).sum()}')

# ── color ─────────────────────────────────────────────────────────────────
colors = np.where(~df['sig'], '#AAAAAA',
         np.where(df['std_effect'] > 0, '#d73027', '#4575b4'))

# ── plot ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))

# nonsig first
mask_ns = ~df['sig']
ax.scatter(df.loc[mask_ns, 'std_effect'], df.loc[mask_ns, 'neg_log10p'],
           c='#AAAAAA', s=4, alpha=0.3, linewidths=0, zorder=1)

# sig on top
mask_s = df['sig']
ax.scatter(df.loc[mask_s, 'std_effect'], df.loc[mask_s, 'neg_log10p'],
           c=colors[mask_s], s=8, alpha=0.8, linewidths=0, zorder=2)

# significance threshold line
if df['sig'].any():
    thresh_p = df.loc[df['sig'], 'pvalue'].max()
    ax.axhline(-np.log10(thresh_p), color='black', linewidth=0.8,
               linestyle='--', alpha=0.5)

ax.axvline(0, color='black', linewidth=0.5, alpha=0.3)

# label top genes
if args.top_n_label > 0 and df['sig'].any():
    top = df[df['sig']].nlargest(args.top_n_label, 'std_effect')
    bot = df[df['sig']].nsmallest(args.top_n_label, 'std_effect')
    for _, row in pd.concat([top, bot]).iterrows():
        sym = str(row.get('symbol', ''))
        if not sym or sym == 'nan': continue
        ax.annotate(sym, (row['std_effect'], row['neg_log10p']),
                    fontsize=7, ha='left', va='bottom',
                    xytext=(3, 3), textcoords='offset points',
                    color='#222222',
                    path_effects=[pe.withStroke(linewidth=2, foreground='white')])

ax.set_xlabel(f'Standardized effect size ({args.predictor})', fontsize=11)
ax.set_ylabel('-log$_{10}$(p-value)', fontsize=11)

parts = [title_prefix]
if args.cell_type: parts.append(args.cell_type)
if args.region:    parts.append(args.region)
parts.append(f'{args.predictor} effect, q<{args.qthresh}')
ax.set_title(' — '.join(parts), fontsize=10)

from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#d73027',
           markersize=6, label='Increasing dispersion'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#4575b4',
           markersize=6, label='Decreasing dispersion'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#AAAAAA',
           markersize=6, label='Not significant'),
], fontsize=8, framealpha=0.8)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

suffix = ''
if args.predictor != 'age': suffix += f'_{args.predictor}'
if args.cell_type: suffix += f'_{args.cell_type}'
if args.region:    suffix += f'_{args.region}'
out = os.path.join(args.outdir, f'volcano_{args.pipeline}{suffix}.{args.outfmt}')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {out}')
print('done.')
