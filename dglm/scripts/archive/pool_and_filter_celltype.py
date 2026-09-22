#!/usr/bin/env python3
"""
Pools ALL of a cell type's raw pseudobulk files (every subcluster x region)
into one combined raw-count matrix (genes x animals, summed), computes CPM
using that pooled matrix's own library sizes, and applies the same filter
rule (CPM>=0.5 in >=50% of animals) ONCE per cell type -- instead of once
per narrow subcluster x region slice. Saves into a completely separate
directory; does not touch the existing per-slice filtered results.

Usage:
    python pool_and_filter_celltype.py
"""
import glob
import os
import pandas as pd

RAW_BASE = '/scratch/easmit31/dispersion/dglm/disp_age__raw_counts'
OUT_BASE = '/scratch/easmit31/dispersion/dglm/disp_age__celltype_pooled'
CPM_CUTOFF = 0.5
ANIMAL_FRAC = 0.5

cell_type_dirs = sorted(d for d in glob.glob(os.path.join(RAW_BASE, '*')) if os.path.isdir(d))

results = []
for d in cell_type_dirs:
    ct = os.path.basename(d)
    files = sorted(glob.glob(os.path.join(d, f'{ct}_*_pseudobulk.csv')))
    if not files:
        continue
    pooled = None
    for f in files:
        df = pd.read_csv(f, index_col=0)
        pooled = df if pooled is None else pooled.add(df, fill_value=0)
    pooled = pooled.fillna(0)

    lib_size = pooled.sum(axis=0)
    cpm = pooled.div(lib_size, axis=1) * 1e6
    n_animals = pooled.shape[1]
    pass_frac = (cpm >= CPM_CUTOFF).sum(axis=1)
    keep = pass_frac >= (ANIMAL_FRAC * n_animals)
    n_kept = int(keep.sum())

    pooled.to_csv(os.path.join(OUT_BASE, f'{ct}_pooled_pseudobulk.csv'))
    kept_genes = pooled.index[keep]
    with open(os.path.join(OUT_BASE, f'{ct}_pooled_filtered_genes.txt'), 'w') as fh:
        fh.write('\n'.join(kept_genes))

    results.append((ct, len(files), n_animals, pooled.shape[0], n_kept))
    print(f'{ct}: {len(files)} files pooled | {n_animals} animals | {pooled.shape[0]} raw genes | {n_kept} kept at CPM>={CPM_CUTOFF} in >={int(ANIMAL_FRAC*100)}% animals')

print()
print('=== summary ===')
print(f'{"cell_type":25s} {"n_files":8s} {"n_animals":10s} {"raw_genes":10s} {"kept_genes":10s}')
for ct, nf, na, rg, kg in results:
    print(f'{ct:25s} {nf:<8d} {na:<10d} {rg:<10d} {kg:<10d}')
