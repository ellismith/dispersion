"""
plot_heatmap_summary.py

Color  = median standardized DGLM beta per cell type x region (raw, not mashr)
Asterisks = mashr LFSR < lfsr_thresh if --mashr_tsv given, else global FDR q < qthresh

Usage:
  python plot_heatmap_summary.py --pipeline dglm \
      --dglm_master master_globalfdr.tsv \
      --mashr_tsv master_mashr.tsv
  python plot_heatmap_summary.py --pipeline gv --gv_master master_between.tsv
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser()
parser.add_argument('--pipeline',    choices=['gv','dglm'], required=True)
parser.add_argument('--gv_master',   default='/scratch/easmit31/dispersion/gene_variance/results_log/master_between.tsv')
parser.add_argument('--dglm_master', default=None, help='for color: raw DGLM betas')
parser.add_argument('--mashr_tsv',   default=None, help='for asterisks: mashr LFSR')
parser.add_argument('--outdir',      default='/scratch/easmit31/dispersion/plotting_scripts/figures')
parser.add_argument('--qthresh',     type=float, default=0.05)
parser.add_argument('--lfsr_thresh', type=float, default=0.2)
parser.add_argument('--outfmt',      default='png')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

CELL_TYPES = ['astrocytes','basket_cells','cerebellar_neurons','ependymal_cells',
              'GABAergic_neurons','glutamatergic_neurons','medium_spiny_neurons',
              'microglia','midbrain_neurons','opc','oligodendrocytes','vascular_cells']
CT_LABELS  = {'astrocytes':'AST','basket_cells':'BC','cerebellar_neurons':'CER',
               'ependymal_cells':'EPEN','GABAergic_neurons':'INH',
               'glutamatergic_neurons':'EXC','medium_spiny_neurons':'MSN',
               'microglia':'MGL','midbrain_neurons':'MBN','opc':'OPC',
               'oligodendrocytes':'OLIG','vascular_cells':'VASC'}
REGIONS = ['ACC','CN','dlPFC','EC','HIP','IPP','lCb','M1','MB','mdTN','NAc']

# ── color: load raw effect sizes ──────────────────────────────────────────
if args.pipeline == 'gv':
    df = pd.read_csv(args.gv_master, sep='\t')
    df = df[df['age_slope'].notna()]
    df['std_effect'] = df['age_slope'] / df['age_slope'].std()
    label        = 'gene_variance (OLS)'
    effect_label = 'Median standardized age slope'
else:
    df = pd.read_csv(args.dglm_master, sep='\t')
    df = df[df['beta'].notna() & (df['beta'].abs() <= 100)]
    df['std_effect'] = df['beta'] / df['bvar'].apply(lambda x: np.sqrt(max(x,1e-10)))
    label        = 'DGLM (dispersion)'
    effect_label = 'Median standardized beta'

# ── asterisks: mashr LFSR or global FDR ──────────────────────────────────
if args.mashr_tsv is not None:
    mdf = pd.read_csv(args.mashr_tsv, sep='\t')
    mdf = mdf[mdf['mash_lfsr'].notna()]
    mdf['sig'] = mdf['mash_lfsr'] < args.lfsr_thresh
    sig_source = mdf
    sig_col    = 'sig'
    sig_label  = f'mashr lfsr<{args.lfsr_thresh}'
else:
    df['sig']  = df['qvalue'] < args.qthresh
    sig_source = df
    sig_col    = 'sig'
    sig_label  = f'q<{args.qthresh}'

# ── compute per cell type x region ───────────────────────────────────────
mat_effect = pd.DataFrame(np.nan,  index=CELL_TYPES, columns=REGIONS)
mat_sig    = pd.DataFrame(False,   index=CELL_TYPES, columns=REGIONS)

for ct in CELL_TYPES:
    for region in REGIONS:
        sub = df[(df['cell_type']==ct) & (df['region']==region)]
        if len(sub) > 0:
            mat_effect.loc[ct, region] = sub['std_effect'].median()
        sig_sub = sig_source[(sig_source['cell_type']==ct) & (sig_source['region']==region)]
        if len(sig_sub) > 0:
            mat_sig.loc[ct, region] = sig_sub[sig_col].any()

# ── plot ──────────────────────────────────────────────────────────────────
row_labels = [CT_LABELS.get(ct, ct) for ct in CELL_TYPES]
fig, ax    = plt.subplots(figsize=(len(REGIONS)*0.85+1.5, len(CELL_TYPES)*0.7+1.5))

vmax = np.nanpercentile(np.abs(mat_effect.values.astype(float)), 95)
vmax = max(vmax, 0.01)

im = ax.imshow(mat_effect.values.astype(float),
               cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')

for i, ct in enumerate(CELL_TYPES):
    for j, region in enumerate(REGIONS):
        if False:
            ax.text(j, i, '*', ha='center', va='center',
                    fontsize=12, color='black', fontweight='bold')

ax.set_xticks(range(len(REGIONS)))
ax.set_xticklabels(REGIONS, fontsize=10)
ax.set_yticks(range(len(CELL_TYPES)))
ax.set_yticklabels(row_labels, fontsize=10)
plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02, label=effect_label)
ax.set_title(f'{label}\n{effect_label} per cell type x region (* = {sig_label})', fontsize=10)
plt.tight_layout()

out = os.path.join(args.outdir, f'heatmap_summary_{args.pipeline}.{args.outfmt}')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {out}')
