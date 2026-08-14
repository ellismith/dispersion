"""
pseudobulk.py

Computes RAW SUMMED-COUNT pseudobulk expression per animal x region for a
given cell type, via a single sparse indicator-matrix multiply per region
(animals x cells indicator @ cells x genes) rather than a per-animal loop --
same summed-count math, much faster at scale. No per-animal min-cell filter,
no gene-presence filter -- all gene-level filtering happens downstream via
CPM-based group-median filtering (filter_group_median_cpm.R), not here.

Saves:
  - {outdir}/{label}_{region}_pseudobulk.csv  (genes x animals, RAW COUNTS)
  - {outdir}/{label}_metadata.csv             (animals x covariates)
  - {outdir}/{label}_gene_names.csv           (ensembl_id -> external_gene_name)

{label} is normally --cell_type. With --subcluster_col set, this instead
loops over every distinct value of that obs column WITHIN the given cell
type (e.g. ct_louvain values like 'astrocytes_8') and pseudobulks each one
separately, with {label} = that subcluster's value.

Usage:
  python pseudobulk.py --cell_type microglia --outdir /path/to/outdir
  python pseudobulk.py --cell_type microglia --outdir /path/to/outdir --subcluster_col ct_louvain
  # fast single-combo test run:
  python pseudobulk.py --cell_type astrocytes --outdir /path/to/test_outdir \
      --subcluster_col ct_louvain --target_subcluster astrocytes_8 --target_region ACC
"""

import argparse
import re
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import os

parser = argparse.ArgumentParser()
parser.add_argument('--cell_type',      required=True)
parser.add_argument('--outdir',         required=True)
parser.add_argument('--h5ad_dir',       default='/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update')
parser.add_argument('--min_age',        type=float, default=1.0)
parser.add_argument('--min_animals',    type=int,   default=5,
                    help='minimum animals a region needs to be included at all (region-level cutoff only -- '
                         'NOT a gene filter; gene-level filtering happens downstream)')
parser.add_argument('--covariates',     type=str,   default='full', choices=['age_sex','full'],
                    help='age_sex: age+sex only; full: age+sex+mean_n_umi+n_cells')
parser.add_argument('--subcluster_col', type=str,   default=None,
                    help='obs column to split into subtype-level pseudobulks (e.g. ct_louvain). '
                         'If unset, pseudobulks the whole cell type as one unit (original behavior).')
parser.add_argument('--min_subcluster_cells', type=int, default=100,
                    help='skip a subcluster entirely if it has fewer than this many cells total '
                         '(across all animals/regions) -- not worth pseudobulking at all below this')
parser.add_argument('--target_region', type=str, default=None,
                    help='restrict to a single region (skip all others) -- for fast iteration/testing')
parser.add_argument('--target_subcluster', type=str, default=None,
                    help='restrict to a single --subcluster_col value (skip all others) -- for fast iteration/testing')
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

def sanitize(label):
    """Make a subcluster value safe to use in filenames."""
    return re.sub(r'[^A-Za-z0-9_.-]', '_', str(label))

def pseudobulk_one(obs, adata, label, gene_names_df):
    import glob as _glob
    for _f in _glob.glob(os.path.join(args.outdir, f"{label}_*_pseudobulk.csv")) + \
               [os.path.join(args.outdir, f"{label}_metadata.csv"),
                os.path.join(args.outdir, f"{label}_gene_names.csv")]:
        if os.path.exists(_f):
            os.remove(_f)

    genes        = adata.var_names.tolist()
    all_metadata = []
    wrote_any    = False

    regions_to_run = [args.target_region] if args.target_region else REGIONS

    for region in regions_to_run:
        obs_r = obs[obs['region'] == region]
        if obs_r.empty:
            continue

        unique_animals = obs_r['animal_id'].unique()
        n_animals = len(unique_animals)
        if n_animals < args.min_animals:
            print(f"  [{label}] {region}: only {n_animals} animals (< --min_animals={args.min_animals}), skipping")
            continue

        # --- vectorized pseudobulk: one sparse matmul per region instead of a per-animal loop ---
        animal_pos    = {a: i for i, a in enumerate(unique_animals)}
        animal_codes  = obs_r['animal_id'].map(animal_pos).values  # per-cell animal index
        n_cells_region = len(obs_r)

        # indicator: n_animals x n_cells_region, indicator[a, c] = 1 if cell c belongs to animal a
        indicator = sparse.csr_matrix(
            (np.ones(n_cells_region), (animal_codes, np.arange(n_cells_region))),
            shape=(n_animals, n_cells_region)
        )

        cell_idx  = obs_r['_idx'].values
        X_region  = X_FULL[cell_idx, :]
        if not sparse.issparse(X_region):
            X_region = sparse.csr_matrix(X_region)
        else:
            X_region = X_region.tocsr()

        # single sparse matmul: (n_animals x n_cells) @ (n_cells x n_genes) -> n_animals x n_genes
        # equivalent to summing each animal's cell rows, computed for all animals at once
        pb_mat = indicator.dot(X_region)
        if sparse.issparse(pb_mat):
            pb_mat = np.asarray(pb_mat.todense())

        pb_animals = list(unique_animals)
        # -----------------------------------------------------------------------------------

        meta_g = obs_r.groupby('animal_id', sort=False).agg(
            age=('age', 'first'), sex=('sex', 'first'), n_cells=('animal_id', 'size')
        )
        if args.covariates == 'full':
            meta_g['mean_n_umi'] = obs_r.groupby('animal_id', sort=False)['n_umi'].mean()
        meta_g = meta_g.reindex(pb_animals)

        for animal in pb_animals:
            meta_row = {
                'animal_id': animal,
                'age':       float(meta_g.loc[animal, 'age']),
                'sex':       meta_g.loc[animal, 'sex'],
                'region':    region,
                'n_cells':   int(meta_g.loc[animal, 'n_cells']),
            }
            if args.covariates == 'full':
                meta_row['mean_n_umi'] = float(meta_g.loc[animal, 'mean_n_umi'])
            all_metadata.append(meta_row)

        df = pd.DataFrame(pb_mat.T, index=genes, columns=pb_animals)
        df.to_csv(os.path.join(args.outdir, f"{label}_{region}_pseudobulk.csv"))
        print(f"  [{label}] {region}: {len(genes)} genes x {len(pb_animals)} animals (raw counts)")
        wrote_any = True

    if wrote_any:
        meta_df = pd.DataFrame(all_metadata).drop_duplicates(subset=['animal_id','region'])
        meta_df.to_csv(os.path.join(args.outdir, f"{label}_metadata.csv"), index=False)
        gene_names_df.to_csv(os.path.join(args.outdir, f"{label}_gene_names.csv"))
    return wrote_any

print(f"loading {args.cell_type} | covariates: {args.covariates}" +
      (f" | subcluster_col: {args.subcluster_col}" if args.subcluster_col else ""))
h5ad_path = os.path.join(args.h5ad_dir, H5AD_MAP[args.cell_type])
adata     = ad.read_h5ad(h5ad_path, backed='r')
print(f"Loading full expression matrix for {args.cell_type} into memory "
      f"(one-time read, avoids repeated backed random-access per animal)...")
X_FULL = adata.X[:]
print("Done loading into memory.")
gene_names            = adata.var[['external_gene_name']].copy()
gene_names.index.name = 'ensembl_id'
obs         = adata.obs.copy()
obs['_idx'] = np.arange(len(obs))
obs         = obs[obs['age'].astype(float) >= args.min_age]
if args.cell_type == 'opc':
    obs = obs[obs['louvain'].astype(str).isin(OPC_LOUVAIN)]
elif args.cell_type == 'oligodendrocytes':
    obs = obs[~obs['louvain'].astype(str).isin(OPC_LOUVAIN)]

if args.subcluster_col is None:
    pseudobulk_one(obs, adata, args.cell_type, gene_names)
else:
    if args.subcluster_col not in obs.columns:
        raise SystemExit(f"--subcluster_col '{args.subcluster_col}' not found in obs columns: "
                          f"{list(obs.columns)}")

    subcluster_counts = obs[args.subcluster_col].value_counts()
    print(f"Found {len(subcluster_counts)} distinct '{args.subcluster_col}' values in {args.cell_type}")

    n_written = 0
    for subcluster_value, n_cells in subcluster_counts.items():
        if args.target_subcluster and str(subcluster_value) != args.target_subcluster:
            continue
        if n_cells < args.min_subcluster_cells:
            print(f"  skipping {subcluster_value}: only {n_cells} cells total (< --min_subcluster_cells={args.min_subcluster_cells})")
            continue

        label     = sanitize(subcluster_value)
        obs_sub   = obs[obs[args.subcluster_col] == subcluster_value]
        wrote_any = pseudobulk_one(obs_sub, adata, label, gene_names)
        if wrote_any:
            n_written += 1
        else:
            print(f"  {label}: no region had enough animals -- nothing written")

    print(f"done: {n_written}/{len(subcluster_counts)} subclusters produced usable output")
print("done.")
