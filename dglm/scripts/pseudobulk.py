"""
pseudobulk.py

Reads h5ad for a given cell type, computes pseudobulk mean per animal per region
using sparse matrix operations throughout. Saves:
  - {outdir}/{cell_type}_{region}_pseudobulk.csv  (genes x animals)
  - {outdir}/{cell_type}_metadata.csv              (animals x covariates)
  - {outdir}/{cell_type}_gene_names.csv            (ensembl_id -> external_gene_name)

Usage:
    python pseudobulk.py --cell_type microglia --outdir /path/to/output
"""

import argparse
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import os

parser = argparse.ArgumentParser()
parser.add_argument('--cell_type',    required=True)
parser.add_argument('--outdir',       required=True)
parser.add_argument('--h5ad_dir',     default='/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update')
parser.add_argument('--min_age',      type=float, default=1.0)
parser.add_argument('--min_animals',  type=int,   default=10)
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
REGIONS     = ['ACC', 'CN', 'dlPFC', 'EC', 'HIP', 'IPP', 'lCb', 'M1', 'MB', 'mdTN', 'NAc']

os.makedirs(args.outdir, exist_ok=True)

# ── load in backed mode ───────────────────────────────────────────────────
h5ad_path = os.path.join(args.h5ad_dir, H5AD_MAP[args.cell_type])
print(f"loading {h5ad_path} (backed mode)")
adata = ad.read_h5ad(h5ad_path, backed='r')
print(f"  shape: {adata.shape}")

# ── save gene name mapping immediately ────────────────────────────────────
gene_names = adata.var[['external_gene_name']].copy()
gene_names.index.name = 'ensembl_id'
gene_names_path = os.path.join(args.outdir, f"{args.cell_type}_gene_names.csv")
gene_names.to_csv(gene_names_path)
print(f"  gene names saved: {gene_names_path}")

obs = adata.obs.copy()
obs['_idx'] = np.arange(len(obs))

# age filter
obs = obs[obs['age'].astype(float) >= args.min_age]
print(f"  after age filter: {len(obs)} cells")

# opc/oligo split
if args.cell_type == 'opc':
    obs = obs[obs['louvain'].astype(str).isin(OPC_LOUVAIN)]
    print(f"  after opc louvain filter: {len(obs)} cells")
elif args.cell_type == 'oligodendrocytes':
    obs = obs[~obs['louvain'].astype(str).isin(OPC_LOUVAIN)]
    print(f"  after oligo louvain filter: {len(obs)} cells")

genes       = adata.var_names.tolist()
all_metadata = []

for region in REGIONS:
    print(f"\nprocessing region {region}")
    obs_r   = obs[obs['region'] == region]
    animals = obs_r['animal_id'].unique()
    print(f"  n animals: {len(animals)}")

    if len(animals) < 5:
        print(f"  too few animals, skipping")
        continue

    pb_rows    = []
    pb_animals = []

    for animal in animals:
        cell_idx = obs_r[obs_r['animal_id'] == animal]['_idx'].values

        X_sub = adata.X[cell_idx, :]
        if not sparse.issparse(X_sub):
            X_sub = sparse.csr_matrix(X_sub)
        else:
            X_sub = X_sub.tocsr()

        mean_vec = np.asarray(X_sub.mean(axis=0)).flatten()
        pb_rows.append(mean_vec)
        pb_animals.append(animal)

        row = obs_r[obs_r['animal_id'] == animal]
        all_metadata.append({
            'animal_id':         animal,
            'age':               float(row['age'].iloc[0]),
            'sex':               row['sex'].iloc[0],
            'sequencing_run_id': row['sequencing_run_id'].iloc[0],
            'sample_reads':      float(row['sample_reads'].mean()),
            'n_umi':             float(row['n_umi'].mean()),
            'region':            region,
        })

    pb_mat = np.vstack(pb_rows)

    expressed  = (pb_mat > 0).sum(axis=0)
    keep       = expressed >= args.min_animals
    pb_keep    = pb_mat[:, keep]
    genes_keep = [g for g, k in zip(genes, keep) if k]
    print(f"  genes after min_animals filter: {len(genes_keep)}")

    if len(genes_keep) == 0:
        print(f"  no genes passed filter, skipping")
        continue

    df = pd.DataFrame(pb_keep.T, index=genes_keep, columns=pb_animals)
    out_path = os.path.join(args.outdir, f"{args.cell_type}_{region}_pseudobulk.csv")
    df.to_csv(out_path)
    print(f"  saved: {out_path}")

meta_df   = pd.DataFrame(all_metadata).drop_duplicates(subset=['animal_id', 'region'])
meta_path = os.path.join(args.outdir, f"{args.cell_type}_metadata.csv")
meta_df.to_csv(meta_path, index=False)
print(f"\nmetadata saved: {meta_path}")
print("done.")
