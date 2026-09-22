#!/usr/bin/env python3
"""
assess_condition_sparsity.py

Diagnostic step, run BEFORE any actual filtering: for a given cell type,
computes n_cells (summed) and n_animals (row count) per subcluster x region
directly from each subcluster's *_metadata.csv (one row per animal x region,
n_cells already given per animal per region -- no need to touch raw h5ad or
pseudobulk counts). Flags which subcluster x region combinations would fail
a >=300 cells AND >=30 animals threshold, and reports how many of the
CURRENT mashr conditions (from the existing master TSV) that corresponds to
-- i.e. how sparse the post-filter matrix would actually be.

Does NOT modify or filter anything -- read-only assessment. Output goes to
a new, separate directory (disp_age__min300cells_min30animals/), not mixed
with existing results.

Usage:
    python assess_condition_sparsity.py --cell_type GABAergic_neurons
"""
import argparse
import glob
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--cell_type', required=True)
parser.add_argument('--base_dir', default='/scratch/easmit31/dispersion/dglm/disp_age__raw_counts')
parser.add_argument('--out_dir', default='/scratch/easmit31/dispersion/dglm/disp_age__min300cells_min30animals')
parser.add_argument('--min_cells', type=int, default=300)
parser.add_argument('--min_animals', type=int, default=30)
args = parser.parse_args()

os.makedirs(args.out_dir, exist_ok=True)
ct_dir = os.path.join(args.base_dir, args.cell_type)

meta_files = sorted(glob.glob(os.path.join(ct_dir, '*_metadata.csv')))
print(f'Found {len(meta_files)} subcluster metadata file(s) for {args.cell_type}')

rows = []
for f in meta_files:
    base = os.path.basename(f)
    subcluster = base[:-len('_metadata.csv')]
    if subcluster.startswith('opc-olig_'):
        subcluster = f'{args.cell_type}_' + subcluster[len('opc-olig_'):]
    df = pd.read_csv(f)
    grp = df.groupby('region').agg(n_cells=('n_cells', 'sum'), n_animals=('animal_id', 'nunique')).reset_index()
    grp['subcluster'] = subcluster
    rows.append(grp)

if not rows:
    raise SystemExit(f'No metadata files found for {args.cell_type} under {ct_dir}')

result = pd.concat(rows, ignore_index=True)
result['pass_cells'] = result['n_cells'] >= args.min_cells
result['pass_animals'] = result['n_animals'] >= args.min_animals
result['pass_both'] = result['pass_cells'] & result['pass_animals']
result = result[['subcluster', 'region', 'n_cells', 'n_animals', 'pass_cells', 'pass_animals', 'pass_both']]
result = result.sort_values(['subcluster', 'region']).reset_index(drop=True)

out_csv = os.path.join(args.out_dir, f'{args.cell_type}_condition_sparsity.csv')
result.to_csv(out_csv, index=False)
print(f'Saved full per-condition table: {out_csv}')

n_total = len(result)
n_pass = int(result['pass_both'].sum())
n_fail = n_total - n_pass
print()
print(f'=== {args.cell_type}: subcluster x region combinations from metadata ===')
print(f'Total: {n_total} | would PASS (>={args.min_cells} cells AND >={args.min_animals} animals): {n_pass} | would be REMOVED: {n_fail}')

# cross-check against the conditions actually used in the current mashr run
master_tsv = os.path.join(ct_dir, 'dglm_checkpoints_cutoff0.5', 'master_dglm_combined.tsv')
if os.path.exists(master_tsv):
    current = pd.read_csv(master_tsv, sep='\t', usecols=['cell_type', 'region'])
    current_conditions = current.drop_duplicates()[['cell_type', 'region']]
    current_conditions.columns = ['subcluster', 'region']
    merged = current_conditions.merge(result[['subcluster', 'region', 'pass_both']], on=['subcluster', 'region'], how='left')
    n_current = len(merged)
    n_current_pass = int(merged['pass_both'].sum())
    print()
    print(f'=== cross-check against CURRENT mashr conditions (master_dglm_combined.tsv) ===')
    print(f'Current conditions actually used in mashr: {n_current}')
    print(f'Of those, would still pass the new filter: {n_current_pass}')
    print(f'Of those, would be REMOVED by the new filter: {n_current - n_current_pass}')
else:
    print(f'(no existing master TSV found at {master_tsv} to cross-check against)')

print()
print('Removed combinations (first 20):')
print(result[~result['pass_both']].head(20).to_string(index=False))
