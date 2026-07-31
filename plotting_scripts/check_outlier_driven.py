"""
check_outlier_driven.py

For top DGLM dispersion genes, checks whether variance is driven by outlier animals
or is broadly distributed. Computes:
  - CV (coefficient of variation) per gene per ct x region
  - Kurtosis (high kurtosis = outlier driven)
  - Plots distribution of z-scores for top genes
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os, glob

DGLM_FILE  = '/scratch/easmit31/dispersion/dglm/disp_age__mean_full/master_disp_age_mean_full_globalfdr.tsv'
PB_DIR     = '/scratch/easmit31/dispersion/dglm/disp_age__mean_full'
META_FILE  = '/scratch/easmit31/u01_metadata_age_sex.csv'
OUT_DIR    = '/scratch/easmit31/dispersion/plotting_scripts/figures'
os.makedirs(OUT_DIR, exist_ok=True)

# top genes by n sig ct x region
print('Loading DGLM results...')
dglm = pd.read_csv(DGLM_FILE, sep='\t')
dglm = dglm[dglm['beta'].notna() & (dglm['beta'].abs() <= 100) & (dglm['qvalue'] < 0.05)]

top_genes = dglm.groupby('symbol').size().sort_values(ascending=False).head(10).index.tolist()
print('Top genes:', top_genes)

meta = pd.read_csv(META_FILE)
age_map = meta.set_index('animal_id')['age']

# for each top gene, look at microglia HIP as example
ct     = 'microglia'
region = 'HIP'
pb_file = os.path.join(PB_DIR, f'{ct}_{region}_pseudobulk.csv')
pb = pd.read_csv(pb_file, index_col=0)

# get ensembl ids for top genes
gene_map = dglm[dglm['symbol'].isin(top_genes)][['ensembl_id','symbol']].drop_duplicates()
gene_map = gene_map[gene_map['ensembl_id'].isin(pb.index)]

print(f'\nGenes found in {ct} {region} pseudobulk: {len(gene_map)}')

results = []
for _, row in gene_map.iterrows():
    vals = pb.loc[row['ensembl_id']].values.astype(float)
    ages = [age_map.get(a, np.nan) for a in pb.columns]
    # exclude infants
    mask = np.array(ages) >= 1.0
    vals = vals[mask]
    ages_filt = np.array(ages)[mask]
    z = (vals - vals.mean()) / (vals.std() + 1e-10)
    results.append(dict(
        symbol   = row['symbol'],
        cv       = vals.std() / (vals.mean() + 1e-10),
        kurtosis = stats.kurtosis(z),
        max_abs_z = np.max(np.abs(z)),
        n_animals = len(vals)
    ))

res_df = pd.DataFrame(results)
print('\nOutlier diagnostics:')
print(res_df.to_string(index=False))

# plot z-score distributions for top 6 genes
genes_to_plot = gene_map.head(6)
fig, axes = plt.subplots(2, 3, figsize=(12, 7))
axes = axes.flatten()

for i, (_, row) in enumerate(genes_to_plot.iterrows()):
    vals = pb.loc[row['ensembl_id']].values.astype(float)
    ages = np.array([age_map.get(a, np.nan) for a in pb.columns])
    mask = ages >= 1.0
    vals = vals[mask]
    ages_filt = ages[mask]
    z = (vals - vals.mean()) / (vals.std() + 1e-10)

    ax = axes[i]
    ax.scatter(ages_filt, z, c=ages_filt, cmap='RdYlBu_r', s=30, alpha=0.8)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.axhline(2, color='red', linewidth=0.5, linestyle='--', alpha=0.5)
    ax.axhline(-2, color='red', linewidth=0.5, linestyle='--', alpha=0.5)
    ax.set_title(f'{row["symbol"]}\nkurt={stats.kurtosis(z):.2f}, max|z|={np.max(np.abs(z)):.2f}',
                 fontsize=9)
    ax.set_xlabel('Age')
    ax.set_ylabel('Z-score')

plt.suptitle(f'Top dispersion genes — {ct} {region}', fontsize=11)
plt.tight_layout()
out = os.path.join(OUT_DIR, f'outlier_check_{ct}_{region}.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'\nSaved: {out}')
print('done.')
