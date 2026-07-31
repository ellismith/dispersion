"""
plot_pca_dglm_correlation.py

Correlates PCA population centroid age effect (mean Pearson r across louvains)
with DGLM dispersion age effect (median standardized beta) per cell type x region.
One point per ct x region, Spearman r across all points.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import stats
import glob, os

POP_DIR   = '/scratch/easmit31/factor_analysis/population_centroid_outputs'
DGLM_FILE = '/scratch/easmit31/dispersion/dglm/disp_age__mean_full/master_disp_age_mean_full_globalfdr.tsv'
OUT_DIR   = '/scratch/easmit31/dispersion/plotting_scripts/figures'
os.makedirs(OUT_DIR, exist_ok=True)

# cell type name mapping PCA -> DGLM
CT_MAP = {
    'astrocytes':          'astrocytes',
    'GABAergic-neurons':   'GABAergic_neurons',
    'GABAergic':           'GABAergic_neurons',
    'glutamatergic-neurons': 'glutamatergic_neurons',
    'microglia':           'microglia',
    'oligodendrocytes':    'oligodendrocytes',
    'ependymal-cells':     'ependymal_cells',
    'basket-cells':        'basket_cells',
    'cerebellar-neurons':  'cerebellar_neurons',
    'medium-spiny-neurons':'medium_spiny_neurons',
    'midbrain-neurons':    'midbrain_neurons',
    'opc':                 'opc',
    'opc-olig':            None,  # skip
    'vascular-cells':      'vascular_cells',
}

CT_COLORS = {
    'astrocytes':'#1b9e77', 'GABAergic_neurons':'#d95f02',
    'glutamatergic_neurons':'#7570b3', 'microglia':'#e7298a',
    'oligodendrocytes':'#66a61e', 'ependymal_cells':'#e6ab02',
    'basket_cells':'#a6761d', 'cerebellar_neurons':'#666666',
    'medium_spiny_neurons':'#ed1c24', 'midbrain_neurons':'#00aeef',
    'opc':'#86328c', 'vascular_cells':'#4575b4',
}

# ── load PCA centroid data ────────────────────────────────────────────────
pca_rows = []
for f in glob.glob(f'{POP_DIR}/*_population_centroid_summary.csv'):
    bn     = os.path.basename(f).replace('_population_centroid_summary.csv','')
    parts  = bn.split('_')
    ct_raw = parts[0]
    region = parts[1]
    ct_dglm = CT_MAP.get(ct_raw)
    if ct_dglm is None:
        continue
    df = pd.read_csv(f)
    if 'r_mean_dist' not in df.columns:
        continue
    pca_rows.append(dict(cell_type=ct_dglm, region=region,
                         pca_r=df['r_mean_dist'].mean()))

pca = pd.DataFrame(pca_rows)
print(f'PCA: {len(pca)} ct x region pairs')

# ── load DGLM — compute median std beta per ct x region ──────────────────
print('Loading DGLM...')
dglm = pd.read_csv(DGLM_FILE, sep='\t')
dglm = dglm[dglm['beta'].notna() & (dglm['beta'].abs() <= 100)]
dglm['std_effect'] = dglm['beta'] / dglm['bvar'].apply(lambda x: np.sqrt(max(x,1e-10)))
dglm_sum = dglm.groupby(['cell_type','region'])['std_effect'].median().reset_index()
dglm_sum.columns = ['cell_type','region','dglm_beta']
print(f'DGLM: {len(dglm_sum)} ct x region pairs')

# ── merge ─────────────────────────────────────────────────────────────────
merged = pd.merge(pca, dglm_sum, on=['cell_type','region'])
print(f'Matched: {len(merged)} ct x region pairs')

rho, pval = stats.spearmanr(merged['pca_r'], merged['dglm_beta'])
print(f'Spearman r = {rho:.3f}, p = {pval:.3e}')

# ── plot ──────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7,6))

for ct in merged['cell_type'].unique():
    sub = merged[merged['cell_type']==ct]
    ax.scatter(sub['pca_r'], sub['dglm_beta'],
               c=CT_COLORS.get(ct,'gray'), s=60, alpha=0.8, label=ct, zorder=2)
    for _, row in sub.iterrows():
        ax.annotate(row['region'], (row['pca_r'], row['dglm_beta']),
                    fontsize=6, ha='left', xytext=(3,3), textcoords='offset points',
                    path_effects=[pe.withStroke(linewidth=2, foreground='white')])

ax.axhline(0, color='black', linewidth=0.5, alpha=0.3)
ax.axvline(0, color='black', linewidth=0.5, alpha=0.3)

ax.set_xlabel('PCA population centroid age effect\n(mean Pearson r across louvains)', fontsize=11)
ax.set_ylabel('DGLM dispersion age effect\n(median standardized beta)', fontsize=11)
ax.set_title(f'PCA centroid vs DGLM dispersion\nSpearman r = {rho:.3f}, p = {pval:.3e}', fontsize=11)

ax.legend(fontsize=7, bbox_to_anchor=(1.01,1), loc='upper left', framealpha=0.8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

out = os.path.join(OUT_DIR, 'pca_centroid_vs_dglm_dispersion.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {out}')
print('done.')
