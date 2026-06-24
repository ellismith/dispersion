"""
run_gene_variance.py

For a given cell type x region h5ad, computes:
  1. Between-individual: pseudobulk mean per animal, then |residuals| ~ age + sex OLS
  2. Within-individual: per-animal variance across cells, regressed ~ age + sex OLS

Fixes vs v1:
  - sparse matrix ops throughout (no toarray unless necessary)
  - raw gene alignment via dict lookup with missing gene guard
  - log1p transform of variance before regression
  - robust sex encoding via pd.get_dummies
  - minimum animals check
  - regression failure counting
  - args metadata saved alongside outputs

Usage:
    python run_gene_variance.py --h5ad PATH --cell_type STR --region STR --outdir PATH
                                [--min_age FLOAT] [--min_animals INT] [--n_hvgs INT]
"""

import argparse
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import statsmodels.api as sm
import os
import json
import sys

# ── args ──────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument('--h5ad',        required=True)
parser.add_argument('--cell_type',   required=True)
parser.add_argument('--region',      required=True)
parser.add_argument('--outdir',      required=True)
parser.add_argument('--min_age',     type=float, default=1.0)
parser.add_argument('--min_animals', type=int,   default=10)
parser.add_argument('--n_hvgs',      type=int,   default=None)
parser.add_argument('--orthologs',   default='/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# ── save args metadata ────────────────────────────────────────────────────
meta = vars(args)
meta_path = os.path.join(args.outdir, f"{args.cell_type}_{args.region}_args.json")
with open(meta_path, 'w') as f:
    json.dump(meta, f, indent=2)
print(f"args saved: {meta_path}")

# ── load orthologs ────────────────────────────────────────────────────────
print("loading ortholog map")
orth = pd.read_csv(args.orthologs)
orth = orth[orth['Human homology type'] == 'ortholog_one2one'][['Gene stable ID', 'Human gene name']].drop_duplicates()
orth = orth.set_index('Gene stable ID')['Human gene name'].to_dict()

# ── load ──────────────────────────────────────────────────────────────────
print(f"loading {args.h5ad}")
adata = ad.read_h5ad(args.h5ad)
print(f"  loaded: {adata.shape}")

# ── filter: region ────────────────────────────────────────────────────────
adata = adata[adata.obs['region'] == args.region].copy()
print(f"  after region filter ({args.region}): {adata.shape}")

# ── filter: age ───────────────────────────────────────────────────────────
adata = adata[adata.obs['age'].astype(float) >= args.min_age].copy()
print(f"  after age filter (>={args.min_age}): {adata.shape}")

# ── opc/oligo split ───────────────────────────────────────────────────────
if args.cell_type == 'opc':
    adata = adata[adata.obs['louvain'].astype(str).isin(['12', '13'])].copy()
    print(f"  after opc louvain filter (12,13): {adata.shape}")
elif args.cell_type == 'oligodendrocytes':
    adata = adata[~adata.obs['louvain'].astype(str).isin(['12', '13'])].copy()
    print(f"  after oligo louvain filter (not 12,13): {adata.shape}")

# ── optional: HVG filter ──────────────────────────────────────────────────
if args.n_hvgs is not None:
    import scanpy as sc
    print(f"  computing top {args.n_hvgs} HVGs")
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvgs, flavor='seurat_v3')
    adata = adata[:, adata.var['highly_variable']].copy()
    print(f"  after HVG filter: {adata.shape}")

# ── get log-normalized matrix (keep sparse) ───────────────────────────────
X_log = adata.X  # keep sparse

genes      = adata.var_names.tolist()
animal_ids = adata.obs['animal_id'].values
ages       = adata.obs['age'].astype(float).values

# ── sex encoding ──────────────────────────────────────────────────────────
sex_dummies = pd.get_dummies(adata.obs['sex'], drop_first=True).values.astype(float)
print(f"  sex categories: {pd.get_dummies(adata.obs['sex'], drop_first=True).columns.tolist()}")

# ── raw counts for visualization ──────────────────────────────────────────
if adata.raw is not None:
    print("  using .raw for visualization columns")
    raw_genes   = list(adata.raw.var_names)
    raw_idx_map = {g: i for i, g in enumerate(raw_genes)}
    # guard: only keep genes present in raw
    genes_in_raw = [g for g in genes if g in raw_idx_map]
    missing      = len(genes) - len(genes_in_raw)
    if missing > 0:
        print(f"  warning: {missing} genes not found in .raw, dropping")
    raw_gene_idx = np.array([raw_idx_map[g] for g in genes_in_raw])
    # subset adata to genes present in raw
    adata        = adata[:, genes_in_raw].copy()
    X_log        = adata.X
    genes        = genes_in_raw
    use_raw      = True
else:
    print("  .raw not found, using .X for visualization columns")
    use_raw      = False

animals = np.unique(animal_ids)
print(f"  n animals: {len(animals)}")

# ── minimum animals check ─────────────────────────────────────────────────
if len(animals) < 5:
    print(f"  ERROR: only {len(animals)} animals after filtering, skipping")
    sys.exit(0)

# ── pseudobulk: mean and variance per animal (sparse-aware) ───────────────
print("computing pseudobulk mean and variance per animal")
n_animals = len(animals)
n_genes   = len(genes)

mean_mat_log = np.zeros((n_animals, n_genes))
var_mat_log  = np.zeros((n_animals, n_genes))
mean_mat_raw = np.zeros((n_animals, n_genes))
var_mat_raw  = np.zeros((n_animals, n_genes))
age_vec      = np.zeros(n_animals)
sex_mat      = np.zeros((n_animals, sex_dummies.shape[1]))

for i, animal in enumerate(animals):
    mask = animal_ids == animal

    # log-normalized (sparse-aware)
    X_sub = X_log[mask]
    if sparse.issparse(X_sub):
        mean_mat_log[i] = np.asarray(X_sub.mean(axis=0)).flatten()
        # variance: E[x^2] - E[x]^2
        X_sub_sq        = X_sub.copy()
        X_sub_sq.data **= 2
        mean_sq         = np.asarray(X_sub_sq.mean(axis=0)).flatten()
        var_mat_log[i]  = mean_sq - mean_mat_log[i] ** 2
    else:
        mean_mat_log[i] = X_sub.mean(axis=0)
        var_mat_log[i]  = X_sub.var(axis=0)

    # raw counts for visualization
    if use_raw:
        raw_sub         = adata.raw.X[mask][:, raw_gene_idx]
        if sparse.issparse(raw_sub):
            mean_mat_raw[i] = np.asarray(raw_sub.mean(axis=0)).flatten()
            raw_sub_sq      = raw_sub.copy()
            raw_sub_sq.data **= 2
            mean_sq_raw     = np.asarray(raw_sub_sq.mean(axis=0)).flatten()
            var_mat_raw[i]  = mean_sq_raw - mean_mat_raw[i] ** 2
        else:
            mean_mat_raw[i] = raw_sub.mean(axis=0)
            var_mat_raw[i]  = raw_sub.var(axis=0)
    else:
        mean_mat_raw[i] = mean_mat_log[i]
        var_mat_raw[i]  = var_mat_log[i]

    age_vec[i]  = ages[mask][0]
    sex_mat[i]  = sex_dummies[mask][0]

# ── filter: min animals expressing gene ───────────────────────────────────
expressed  = (mean_mat_log > 0).sum(axis=0)
keep       = expressed >= args.min_animals
mean_mat_log = mean_mat_log[:, keep]
var_mat_log  = var_mat_log[:, keep]
mean_mat_raw = mean_mat_raw[:, keep]
var_mat_raw  = var_mat_raw[:, keep]
genes_kept   = [g for g, k in zip(genes, keep) if k]
print(f"  genes after min_animals filter: {len(genes_kept)} / {len(genes)}")

# ── precompute per-gene summary stats from raw (for plotting) ─────────────
gene_mean_expr   = mean_mat_raw.mean(axis=0)
gene_var_between = mean_mat_raw.var(axis=0)
gene_var_within  = var_mat_raw.mean(axis=0)

# ── regression helper ─────────────────────────────────────────────────────
def run_variance_regression(matrix, age_vec, sex_mat, genes,
                             gene_mean_expr, gene_var_between, gene_var_within, mode):
    results      = []
    n_failures   = 0
    covariates   = np.column_stack([age_vec, sex_mat])
    X_cov        = sm.add_constant(covariates)

    for j, gene in enumerate(genes):
        y = matrix[:, j]
        try:
            if mode == 'between':
                mean_fit = sm.OLS(y, X_cov).fit()
                y2       = np.log1p(np.abs(mean_fit.resid))
            else:
                y2       = np.log1p(y)  # log1p of per-animal variance

            var_fit   = sm.OLS(y2, X_cov).fit()
            slope     = var_fit.params[1]
            pval      = var_fit.pvalues[1]
            intercept = var_fit.params[0]
        except Exception:
            slope, pval, intercept = np.nan, np.nan, np.nan
            n_failures += 1

        results.append({
            'ensembl_id':   gene,
            'human_symbol': orth.get(gene, np.nan),
            'mean_expr':    gene_mean_expr[j],
            'var_between':  gene_var_between[j],
            'var_within':   gene_var_within[j],
            'age_slope':    slope,
            'pvalue':       pval,
            'intercept':    intercept
        })

    if n_failures > 0:
        print(f"  warning: {n_failures} regression failures (set to nan)")

    return pd.DataFrame(results)

# ── run: between-individual ───────────────────────────────────────────────
print("running between-individual regression")
between_df = run_variance_regression(
    mean_mat_log, age_vec, sex_mat, genes_kept,
    gene_mean_expr, gene_var_between, gene_var_within, mode='between'
)
out_between = os.path.join(args.outdir, f"{args.cell_type}_{args.region}_between.tsv")
between_df.to_csv(out_between, sep='\t', index=False)
print(f"  saved: {out_between}")

# ── run: within-individual ────────────────────────────────────────────────
print("running within-individual regression")
within_df = run_variance_regression(
    var_mat_log, age_vec, sex_mat, genes_kept,
    gene_mean_expr, gene_var_between, gene_var_within, mode='within'
)
out_within = os.path.join(args.outdir, f"{args.cell_type}_{args.region}_within.tsv")
within_df.to_csv(out_within, sep='\t', index=False)
print(f"  saved: {out_within}")

print("done.")
