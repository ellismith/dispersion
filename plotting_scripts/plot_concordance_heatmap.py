"""
plot_concordance_heatmap.py

Heatmap showing n genes significant in:
  - both pipelines
  - gene_variance only
  - DGLM only
per cell type x region. Three separate heatmaps.

Usage:
  python plot_concordance_heatmap.py
  python plot_concordance_heatmap.py --q_gv 0.05 --q_dglm 0.05
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser()
parser.add_argument('--gv_master',   default='/scratch/easmit31/variability/gene_variance/results_log/master_between.tsv')
parser.add_argument('--dglm_master', default='/scratch/easmit31/variability/dglm/checkpoints/master_dglm_globalfdr.tsv')
parser.add_argument('--outdir',      default='/scratch/easmit31/variability/plotting_scripts/figures')
parser.add_argument('--q_gv',        type=float, default=0.05)
parser.add_argument('--q_dglm',      type=float, default=0.05)
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
gv['symbol'] = gv['human_symbol'].fillna(gv['ensembl_id'])
gv['key']    = gv['symbol'] + '|' + gv['cell_type'] + '|' + gv['region']
gv['sig']    = gv['qvalue'] < args.q_gv

print('Loading DGLM')
dglm = pd.read_csv(args.dglm_master, sep='\t')
dglm = dglm[dglm['pvalue'].notna() & (dglm['beta'].abs() <= 100)]
dglm['key'] = dglm['symbol'] + '|' + dglm['cell_type'] + '|' + dglm['region']
dglm['sig'] = dglm['qvalue'] < args.q_dglm

gv_sig_keys   = set(gv.loc[gv['sig'], 'key'])
dglm_sig_keys = set(dglm.loc[dglm['sig'], 'key'])

# ── compute counts per cell type x region ────────────────────────────────
mat_both  = pd.DataFrame(np.nan, index=CELL_TYPES, columns=REGIONS)
mat_gv    = pd.DataFrame(np.nan, index=CELL_TYPES, columns=REGIONS)
mat_dglm  = pd.DataFrame(np.nan, index=CELL_TYPES, columns=REGIONS)

for ct in CELL_TYPES:
    for region in REGIONS:
        gv_sub   = gv[(gv['cell_type']==ct) & (gv['region']==region)]
        dglm_sub = dglm[(dglm['cell_type']==ct) & (dglm['region']==region)]
        if len(gv_sub) == 0 or len(dglm_sub) == 0:
            continue
        all_keys  = set(gv_sub['key']) | set(dglm_sub['key'])
        gv_keys   = set(gv_sub.loc[gv_sub['sig'], 'key'])
        dglm_keys = set(dglm_sub.loc[dglm_sub['sig'], 'key'])
        both_keys = gv_keys & dglm_keys
        mat_both.loc[ct, region]  = len(both_keys)
        mat_gv.loc[ct, region]    = len(gv_keys - dglm_keys)
        mat_dglm.loc[ct, region]  = len(dglm_keys - gv_keys)

# ── plot three heatmaps side by side ─────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(len(REGIONS)*2.2+2, len(CELL_TYPES)*0.65+2))

panels = [
    (mat_both, 'Both pipelines', '#8B2FC9', 'Purples'),
    (mat_gv,   'gene_variance only', '#fc8d59', 'Oranges'),
    (mat_dglm, 'DGLM only', '#9970ab', 'RdPu'),
]

for ax, (mat, title, color, cmap) in zip(axes, panels):
    vals = mat.values.astype(float)
    vmax = np.nanpercentile(vals, 95)
    if vmax == 0: vmax = 1

    im = ax.imshow(vals, cmap=cmap, vmin=0, vmax=vmax, aspect='auto')

    for i, ct in enumerate(CELL_TYPES):
        for j, region in enumerate(REGIONS):
            v = mat.loc[ct, region]
            if np.isnan(v):
                continue
            ax.text(j, i, str(int(v)), ha='center', va='center',
                    fontsize=6.5, color='black' if v < vmax*0.7 else 'white')

    ax.set_xticks(range(len(REGIONS)))
    ax.set_xticklabels(REGIONS, fontsize=9, rotation=45, ha='right')
    ax.set_yticks(range(len(CELL_TYPES)))
    ax.set_yticklabels([CT_LABELS.get(ct, ct) for ct in CELL_TYPES], fontsize=9)
    ax.set_title(title, fontsize=10, pad=8)

    cbar = plt.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label('n genes', fontsize=8)

fig.suptitle(f'Gene overlap per cell type × region\n(gene_variance q<{args.q_gv}, DGLM q<{args.q_dglm})',
             fontsize=11, y=1.01)
plt.tight_layout()

out = os.path.join(args.outdir, f'concordance_heatmap.{args.outfmt}')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {out}')
print('done.')
