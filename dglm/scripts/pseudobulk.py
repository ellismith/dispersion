"""
pseudobulk.py

Computes pseudobulk expression per animal x region for a given cell type.
Saves:
  - {outdir}/{label}_{region}_pseudobulk.csv  (genes x animals)
  - {outdir}/{label}_metadata.csv             (animals x covariates)
  - {outdir}/{label}_gene_names.csv           (ensembl_id -> external_gene_name)

{label} is normally --cell_type. With --subcluster_col set, this instead
loops over every distinct value of that obs column WITHIN the given cell
type (e.g. ct_louvain values like 'microglia_16') and pseudobulks each one
separately, with {label} = that subcluster's value. Downstream scripts
(dglm_model.R, dglm_mashr.R) don't need to know the difference -- they just
take a --cell_type string and an --outdir, so a subcluster label slots in
exactly like a normal cell type would.

Usage:
  python pseudobulk.py --cell_type microglia --outdir /path/to/outdir
  python pseudobulk.py --cell_type microglia --outdir /path/to/outdir --covariates age_sex
  python pseudobulk.py --cell_type microglia --outdir /path/to/outdir --covariates full
  python pseudobulk.py --cell_type microglia --outdir /path/to/outdir --subcluster_col ct_louvain
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
parser.add_argument('--min_animals',    type=int,   default=10)
parser.add_argument('--covariates',     type=str,   default='full', choices=['age_sex','full'],
                    help='age_sex: age+sex only; full: age+sex+mean_n_umi+n_cells')
parser.add_argument('--subcluster_col', type=str,   default=None,
                    help='obs column to split into subtype-level pseudobulks (e.g. ct_louvain). '
                         'If unset, pseudobulks the whole cell type as one unit (original behavior).')
parser.add_argument('--min_cells',           type=int, default=100,
                    help='minimum cells a single animal must contribute to a given cell_type/subcluster x region '
                         'to be included -- an animal-level pseudobulk mean built from very few cells (e.g. 1-3) '
                         'is essentially that cell\'s raw value, not a real average, and can destabilize the '
                         'dispersion submodel fit for that gene/condition')
parser.add_argument('--min_subcluster_cells', type=int, default=100,
                    help='skip a subcluster entirely if it has fewer than this many cells total '
                         '(across all animals/regions) -- not worth pseudobulking at all below this')
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
    # Always clear any pre-existing output for this label before writing --
    # this script previously only ever ADDED files, never removed stale
    # ones from a prior run. If a region no longer survives filtering on a
    # rerun, its OLD file from a previous run stayed behind and silently
    # mismatched against the fresh metadata.csv (duplicate/NA row-names
    # crash in dglm_model.R downstream). Confirmed root cause across THREE
    # separate cell types this session (ependymal_cells, glutamatergic_neurons,
    # opc) -- fixing here once, systemically, instead of per-cell-type.
    import glob as _glob
    for _f in _glob.glob(os.path.join(args.outdir, f"{label}_*_pseudobulk.csv")) + \
               [os.path.join(args.outdir, f"{label}_metadata.csv"),
                os.path.join(args.outdir, f"{label}_gene_names.csv")]:
        if os.path.exists(_f):
            os.remove(_f)
    """Runs the full per-region pseudobulk + metadata pipeline for one unit
    (either the whole cell type, or one subcluster), writing files under
    {label}_*. Returns True if anything was written, False if this unit had
    no viable region with enough animals/genes."""
    genes        = adata.var_names.tolist()
    all_metadata = []
    wrote_any    = False

    for region in REGIONS:
        obs_r   = obs[obs['region'] == region]
        animals = obs_r['animal_id'].unique()
        if len(animals) < 5:
            continue

        pb_rows    = []
        pb_animals = []

        n_skipped_low_cells = 0
        for animal in animals:
            cell_idx = obs_r[obs_r['animal_id'] == animal]['_idx'].values
            n_cells_animal = len(cell_idx)
            if n_cells_animal < args.min_cells:
                # too few cells for this animal's pseudobulk mean to be
                # trustworthy -- confirmed this matters: a 1-cell and a
                # 3-cell animal alone (out of 47) were enough to destabilize
                # a dispersion-submodel sex coefficient to beta=43 with an
                # otherwise perfectly normal-looking Shat, invisible to any
                # downstream sanity check
                n_skipped_low_cells += 1
                continue
            X_sub    = X_FULL[cell_idx, :]
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
                'n_cells':   n_cells_animal,
            }
            if args.covariates == 'full':
                meta_row['mean_n_umi'] = float(row['n_umi'].mean())
            all_metadata.append(meta_row)

        if n_skipped_low_cells > 0:
            print(f"  [{label}] {region}: skipped {n_skipped_low_cells} animal(s) with "
                  f"< --min_cells={args.min_cells} cells for this condition")

        if len(pb_rows) == 0:
            print(f"  [{label}] {region}: no animals left after --min_cells filter, skipping")
            continue

        pb_mat     = np.vstack(pb_rows)
        expressed  = (pb_mat > 0).sum(axis=0)
        keep       = expressed >= args.min_animals
        pb_keep    = pb_mat[:, keep]
        genes_keep = [g for g, k in zip(genes, keep) if k]

        if len(genes_keep) == 0:
            continue

        df = pd.DataFrame(pb_keep.T, index=genes_keep, columns=pb_animals)
        df.to_csv(os.path.join(args.outdir, f"{label}_{region}_pseudobulk.csv"))
        print(f"  [{label}] {region}: {len(genes_keep)} genes x {len(pb_animals)} animals")
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

# Load the full expression matrix into memory ONCE, up front, into a
# SEPARATE variable -- not adata.X = adata.X[:], which on a backed
# AnnData tries to WRITE the new value back into the (read-only) h5ad
# file and crashes with "no write intent on file". X_FULL is used for
# all per-animal slicing instead of adata.X, avoiding both the repeated
# backed random-access reads AND the accidental write-back.
print(f"Loading full expression matrix for {args.cell_type} into memory "
      f"(one-time read, avoids repeated backed random-access per animal)...")
X_FULL = adata.X[:]
print("Done loading into memory.")

gene_names            = adata.var[['external_gene_name']].copy()
gene_names.index.name = 'ensembl_id'

obs         = adata.obs.copy()
obs['_idx'] = np.arange(len(obs))
obs         = obs[obs['age'].astype(float) >= args.min_age]

# opc/oligodendrocytes share one h5ad -- this filter identifies which parent
# cell type's cells we're even looking at. Applies regardless of
# --subcluster_col, and always BEFORE any subcluster split.
if args.cell_type == 'opc':
    obs = obs[obs['louvain'].astype(str).isin(OPC_LOUVAIN)]
elif args.cell_type == 'oligodendrocytes':
    obs = obs[~obs['louvain'].astype(str).isin(OPC_LOUVAIN)]

if args.subcluster_col is None:
    # original behavior: whole cell type as one unit
    pseudobulk_one(obs, adata, args.cell_type, gene_names)
else:
    if args.subcluster_col not in obs.columns:
        raise SystemExit(f"--subcluster_col '{args.subcluster_col}' not found in obs columns: "
                          f"{list(obs.columns)}")

    subcluster_counts = obs[args.subcluster_col].value_counts()
    print(f"Found {len(subcluster_counts)} distinct '{args.subcluster_col}' values in {args.cell_type}")

    n_written = 0
    for subcluster_value, n_cells in subcluster_counts.items():
        if n_cells < args.min_subcluster_cells:
            print(f"  skipping {subcluster_value}: only {n_cells} cells total (< --min_subcluster_cells={args.min_subcluster_cells})")
            continue

        label     = sanitize(subcluster_value)
        obs_sub   = obs[obs[args.subcluster_col] == subcluster_value]
        wrote_any = pseudobulk_one(obs_sub, adata, label, gene_names)
        if wrote_any:
            n_written += 1
        else:
            print(f"  {label}: no region had enough animals/genes -- nothing written")

    print(f"done: {n_written}/{len(subcluster_counts)} subclusters produced usable output")

print("done.")
