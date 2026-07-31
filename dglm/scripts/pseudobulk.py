"""
pseudobulk.py

Computes pseudobulk expression per animal x region for a given cell type.
Saves:
  - {outdir}/{cell_type}_{region}_pseudobulk.csv  (genes x animals)
  - {outdir}/{cell_type}_metadata.csv             (animals x covariates)
  - {outdir}/{cell_type}_gene_names.csv           (ensembl_id -> external_gene_name)

Usage:
  python pseudobulk.py --cell_type microglia --outdir /path/to/outdir
  python pseudobulk.py --cell_type microglia --outdir /path/to/outdir --covariates age_sex
  python pseudobulk.py --cell_type microglia --outdir /path/to/outdir --covariates full
"""

import argparse
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import os

parser = argparse.ArgumentParser()
parser.add_argument('--cell_type',   required=True)
parser.add_argument('--outdir',      required=True)
parser.add_argument('--h5ad_dir',    default='/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update')
parser.add_argument('--min_age',     type=float, default=1.0)
parser.add_argument('--min_animals', type=int,   default=10)
parser.add_argument('--covariates',  type=str,   default='full', choices=['age_sex','full'],
                    help='age_sex: age+sex only; full: age+sex+mean_n_umi+n_cells')
args = parser.parse_args()

H5AD_MAP = {
    'astrocytes':            'Res1_astrocytes_update.h5ad',
    'basket_cells':          'Res1_basket-cells_update.h5ad',
    'cerebellar_neurons':    'Res1_cerebellar-neurons_subset.h5ad',
    'ependymal_cells':       'Res1_ependymal-cells_new.h5ad',
    'GABAergic_neurons':     'Res1_GABAergic-neurons_subset.h5ad',
    'glutamatergic_neurons': 'Res1_glutamatergic-neurons_update.h5ad',
    'medium_spiny_neurons':  'Res1_medium-spiny-neurons_subset.h5ad',
    'microglia':             'Res1_microglia_new.h5ad',
    'midbrain_neurons':      'Res1_midbrain-neurons_update.h5ad',
    'opc':                   'Res1_opc-olig_subset.h5ad',
    'oligodendrocytes':      'Res1_opc-olig_subset.h5ad',
    'vascular_cells':        'Res1_vascular-cells_subset.h5ad',
}

OPC_LOUVAIN = {'12', '13'}
REGIONS     = ['ACC','CN','dlPFC','EC','HIP','IPP','lCb','M1','MB','mdTN','NAc']

os.makedirs(args.outdir, exist_ok=True)

print(f"loading {args.cell_type} | covariates: {args.covariates}")
h5ad_path = os.path.join(args.h5ad_dir, H5AD_MAP[args.cell_type])
adata     = ad.read_h5ad(h5ad_path, backed='r')

gene_names            = adata.var[['external_gene_name']].copy()
gene_names.index.name = 'ensembl_id'
gene_names.to_csv(os.path.join(args.outdir, f"{args.cell_type}_gene_names.csv"))

obs        = adata.obs.copy()
obs['_idx'] = np.arange(len(obs))
obs        = obs[obs['age'].astype(float) >= args.min_age]

if args.cell_type == 'opc':
    obs = obs[obs['louvain'].astype(str).isin(OPC_LOUVAIN)]
elif args.cell_type == 'oligodendrocytes':
    obs = obs[~obs['louvain'].astype(str).isin(OPC_LOUVAIN)]

genes        = adata.var_names.tolist()
all_metadata = []

for region in REGIONS:
    obs_r   = obs[obs['region'] == region]
    animals = obs_r['animal_id'].unique()
    if len(animals) < 5:
        continue

    pb_rows = []
    pb_animals = []

    for animal in animals:
        cell_idx = obs_r[obs_r['animal_id'] == animal]['_idx'].values
        X_sub    = adata.X[cell_idx, :]
        if not sparse.issparse(X_sub):
            X_sub = sparse.csr_matrix(X_sub)
        else:
            X_sub = X_sub.tocsr()
        mean_vec = np.asarray(X_sub.mean(axis=0)).flatten()
        pb_rows.append(mean_vec)
        pb_animals.append(animal)

        row      = obs_r[obs_r['animal_id'] == animal]
        meta_row = {
            'animal_id': animal,
            'age':       float(row['age'].iloc[0]),
            'sex':       row['sex'].iloc[0],
            'region':    region,
        }
        if args.covariates == 'full':
            meta_row['mean_n_umi'] = float(row['n_umi'].mean())
            meta_row['n_cells']    = len(row)
        all_metadata.append(meta_row)

    pb_mat     = np.vstack(pb_rows)
    expressed  = (pb_mat > 0).sum(axis=0)
    keep       = expressed >= args.min_animals
    pb_keep    = pb_mat[:, keep]
    genes_keep = [g for g, k in zip(genes, keep) if k]

    if len(genes_keep) == 0:
        continue

    df = pd.DataFrame(pb_keep.T, index=genes_keep, columns=pb_animals)
    df.to_csv(os.path.join(args.outdir, f"{args.cell_type}_{region}_pseudobulk.csv"))
    print(f"  {region}: {len(genes_keep)} genes x {len(pb_animals)} animals")

meta_df = pd.DataFrame(all_metadata).drop_duplicates(subset=['animal_id','region'])
meta_df.to_csv(os.path.join(args.outdir, f"{args.cell_type}_metadata.csv"), index=False)
print(f"metadata columns: {list(meta_df.columns)}")
print("done.")
