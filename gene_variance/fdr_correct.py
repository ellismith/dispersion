"""
fdr_correct.py

Pools all between and within TSVs, applies single BH-FDR correction, saves master tables.

Usage:
    python fdr_correct.py --indir PATH --outdir PATH
"""

import argparse
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import glob
import os

parser = argparse.ArgumentParser()
parser.add_argument('--indir',  required=True)
parser.add_argument('--outdir', required=True)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

for mode in ['between', 'within']:
    print(f"\nprocessing {mode}")
    files = glob.glob(os.path.join(args.indir, f'*_{mode}.tsv'))
    print(f"  found {len(files)} TSVs")

    dfs = []
    skipped = 0
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t')
            if df.empty:
                skipped += 1
                continue
        except Exception as e:
            print(f"  skipping {os.path.basename(f)}: {e}")
            skipped += 1
            continue

        basename  = os.path.basename(f).replace(f'_{mode}.tsv', '')
        parts     = basename.split('_')
        region    = parts[-1]
        cell_type = '_'.join(parts[:-1])

        df['cell_type'] = cell_type
        df['region']    = region
        dfs.append(df)

    print(f"  skipped {skipped} files, loaded {len(dfs)}")

    master = pd.concat(dfs, ignore_index=True)
    print(f"  total tests: {len(master)}")

    n_before = len(master)
    master   = master.dropna(subset=['pvalue'])
    print(f"  dropped {n_before - len(master)} nan pvalues")

    _, qvalues, _, _ = multipletests(master['pvalue'].values, method='fdr_bh')
    master['qvalue'] = qvalues

    out = os.path.join(args.outdir, f'master_{mode}.tsv')
    master.to_csv(out, sep='\t', index=False)
    print(f"  saved: {out}")
    print(f"  sig at q<0.05: {(master['qvalue'] < 0.05).sum()}")
    print(f"  sig positive:  {((master['qvalue'] < 0.05) & (master['age_slope'] > 0)).sum()}")
    print(f"  sig negative:  {((master['qvalue'] < 0.05) & (master['age_slope'] < 0)).sum()}")

print("\ndone.")
