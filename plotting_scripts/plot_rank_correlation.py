"""
plot_rank_correlation.py

Spearman rank correlation of -log10(pvalue) between gene_variance and DGLM
per cell type x region. Summarized as heatmap (cell types x regions).

Usage:
  python plot_rank_correlation.py
  python plot_rank_correlation.py --qthresh 0.05
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os

parser = argparse.ArgumentParser()
parser.add_argument('--gv_master',   default='/scratch/easmit31/variability/gene_variance/results_log/master_between.tsv')
parser.add_argument('--dglm_master', default='/scratch/easmit31/variability/dglm/checkpoints/master_dglm_globalfdr.tsv')
parser.add_argument('--outdir',      default='/scratch/easmit31/variability/plotting_scripts/figures')
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

# ── load ──────────────────────────────────────────────────────────────────
print('Loading gene_variance')
gv = pd.read_csv(args.gv_master, sep='\t')
gv = gv[gv['pvalue'].notna()]
gv['neg_log10p'] = -np.log10(gv['pvalue'].clip(lower=1e-300))
gv['symbol']     = gv['human_symbol'].fillna(gv['ensembl_id'])
gv['key']        = gv['symbol'] + '|' + gv['cell_type'] + '|' + gv['region']

print('Loading DGLM')
dglm = pd.read_csv(args.dglm_master, sep='\t')
dglm = dglm[dglm['pvalue'].notna() & (dglm['beta'].abs() <= 100)]
dglm['neg_log10p'] = -np.log10(dglm['pvalue'].clip(lower=1e-300))
dglm['key']        = dglm['symbol'] + '|' + dglm['cell_type'] + '|' + dglm['region']

# ── compute spearman r per cell type x region ─────────────────────────────
mat_r   = pd.DataFrame(np.nan, index=CELL_TYPES, columns=REGIONS)
mat_p   = pd.DataFrame(np.nan, index=CELL_TYPES, columns=REGIONS)
mat_n   = pd.DataFrame(0,      index=CELL_TYPES, columns=REGIONS)

for ct in CELL_TYPES:
    for region in REGIONS:
        gv_sub   = gv[(gv['cell_type']==ct) & (gv['region']==region)]
        dglm_sub = dglm[(dglm['cell_type']==ct) & (dglm['region']==region)]
        merged   = pd.merge(gv_sub[['key','neg_log10p']],
                            dglm_sub[['key','neg_log10p']],
                            on='key', suffixes=('_gv','_dglm'))
        if len(merged) < 10:
            continue
        rho, pval = stats.spearmanr(merged['neg_log10p_gv'],
                                    merged['neg_log10p_dglm'])
        mat_r.loc[ct, region] = rho
        mat_p.loc[ct, region] = pval
        mat_n.loc[ct, region] = len(merged)

print('Spearman r range:', mat_r.stack().min(), mat_r.stack().max())

# ── plot ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(len(REGIONS)*0.85+2, len(CELL_TYPES)*0.7+2))

vals  = mat_r.values.astype(float)
im    = ax.imshow(vals, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')

for i, ct in enumerate(CELL_TYPES):
    for j, region in enumerate(REGIONS):
        r = mat_r.loc[ct, region]
        n = mat_n.loc[ct, region]
        p = mat_p.loc[ct, region]
        if np.isnan(r):
            continue
        # asterisk if p < 0.05
        star = '*' if not np.isnan(p) and p < 0.05 else ''
        ax.text(j, i, f'{r:.2f}{star}', ha='center', va='center',
                fontsize=7, color='black' if abs(r) < 0.7 else 'white')

ax.set_xticks(range(len(REGIONS)))
ax.set_xticklabels(REGIONS, fontsize=10)
ax.set_yticks(range(len(CELL_TYPES)))
ax.set_yticklabels([CT_LABELS.get(ct, ct) for ct in CELL_TYPES], fontsize=10)

cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
cbar.set_label('Spearman r (-log10 p)', fontsize=9)

ax.set_title('Rank correlation: gene_variance vs DGLM\nper cell type × region', fontsize=11)
plt.tight_layout()

out = os.path.join(args.outdir, f'rank_correlation_heatmap.{args.outfmt}')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {out}')
print('done.')
