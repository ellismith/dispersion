#!/usr/bin/env python3
"""
build_autosome_x_gene_list.py

One-time extraction of the autosome+X gene keep-list, matching Chiou et
al.'s filter (chromosome_name %in% c(1:20,'X')) without needing biomaRt
(biomaRt is blocked in mashr_env by the missing compiler toolchain -- same
blocker as GO enrichment).

Reads var['chr'] directly from each h5ad in H5AD_MAP (extracted from
pseudobulk.py's source text, NOT imported -- pseudobulk.py runs its own
argparse at module level, unguarded by __main__, so a normal import drags
in its whole CLI and would error or worse, try to actually execute the
pipeline). Extracting the dict literal directly avoids executing any of
pseudobulk.py's code, while still reading the real H5AD_MAP so it can't
drift out of sync with a retyped copy. Checks that
chr assignment is IDENTICAL for any gene seen in more than one h5ad file
(same reference genome expected across all cell types) -- asserts rather
than assumes, since this filter is only trustworthy if that holds.

Confirmed from a live check (astrocytes h5ad):
  var index / 'ensembl_gene_id' column: ENSMMUG... IDs (matches master TSV's
  ensembl_id column, so the keep-list will actually line up downstream)
  chr values: plain strings, autosomes '1'-'20', plus 'X', 'Y', 'MT', and
  ~200 unplaced scaffold contigs (e.g. 'QNVO02000432.1') -- no 'chr' prefix.

Usage:
    conda activate mixed_models
    python3 build_autosome_x_gene_list.py
Writes:
    /scratch/easmit31/dispersion/dglm/autosome_x_genes.csv
    (one column: ensembl_gene_id -- genes to KEEP)
"""
import ast
import os
import re

import anndata as ad
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Extract H5AD_MAP from pseudobulk.py's SOURCE TEXT without importing/running
# the file -- pseudobulk.py calls parse_args() at module level (not guarded
# by if __name__=='__main__'), so a normal import fails on missing required
# args, or worse, could go on to execute real pipeline logic with garbage
# args. This regex expects a flat dict (string keys/values, no nesting) --
# if pseudobulk.py's H5AD_MAP structure ever changes, this will raise
# rather than silently grab the wrong text.
pb_path = os.path.join(SCRIPT_DIR, 'pseudobulk.py')
with open(pb_path) as f:
    pb_src = f.read()
m = re.search(r'H5AD_MAP\s*=\s*(\{[^}]*\})', pb_src)
if not m:
    raise RuntimeError(
        f'Could not find a flat H5AD_MAP = {{...}} dict in {pb_path} -- '
        'its structure may have changed. Check by hand before proceeding.'
    )
H5AD_MAP = ast.literal_eval(m.group(1))
print(f'Extracted H5AD_MAP from source: {len(H5AD_MAP)} cell type(s)')

H5AD_DIR = '/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update'
OUT_PATH = '/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv'

AUTOSOMES = [str(i) for i in range(1, 21)]
KEEP_CHRS = set(AUTOSOMES + ['X'])

seen = {}  # ensembl_gene_id -> chr, to check consistency across files
unique_files = sorted(set(H5AD_MAP.values()))
print(f'Checking {len(unique_files)} unique h5ad file(s) for chr consistency...')

for fname in unique_files:
    path = os.path.join(H5AD_DIR, fname)
    if not os.path.exists(path):
        print(f'  MISSING (skipping): {path}')
        continue
    print(f'  reading {fname} (backed mode)...')
    a = ad.read_h5ad(path, backed='r')
    if 'ensembl_gene_id' not in a.var.columns or 'chr' not in a.var.columns:
        raise ValueError(f'{fname}: expected var columns ensembl_gene_id/chr not found -- got {a.var.columns.tolist()}')
    for gene_id, chrom in zip(a.var['ensembl_gene_id'], a.var['chr']):
        if gene_id in seen and seen[gene_id] != chrom:
            raise ValueError(
                f'Inconsistent chr for {gene_id}: {seen[gene_id]!r} in an earlier '
                f'file vs {chrom!r} in {fname} -- reference genome may differ '
                f'across h5ads. STOP and investigate before trusting this filter.'
            )
        seen[gene_id] = chrom

print(f'Total unique genes seen across all files: {len(seen)}')

chr_counts = pd.Series(list(seen.values())).value_counts()
print('Top chr value counts:')
print(chr_counts.head(25))

keep = sorted([g for g, c in seen.items() if c in KEEP_CHRS])
dropped_y = sum(1 for c in seen.values() if c == 'Y')
dropped_mt = sum(1 for c in seen.values() if c == 'MT')
dropped_other = len(seen) - len(keep) - dropped_y - dropped_mt

print(f'Keeping {len(keep)} genes (autosomes 1-20 + X)')
print(f'Dropping {len(seen) - len(keep)} genes total: '
      f'{dropped_y} on Y, {dropped_mt} on MT, {dropped_other} on unplaced/other scaffolds')

out = pd.DataFrame({'ensembl_gene_id': keep})
out.to_csv(OUT_PATH, index=False)
print(f'Saved: {OUT_PATH}')
