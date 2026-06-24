"""
identify_outliers.py

Identifies outlier animals per cell type x region based on:
  1. Low cell count (< --min_cells threshold)
  2. Expression outliers (pseudobulk mean z-score > --zscore_thresh SD from group)
  3. Known outliers from --known_outliers list

Produces:
  - outlier_report.tsv: all flagged animal x ct x region combinations
  - outlier_summary.png: heatmap of n outlier flags per animal x region

Usage:
  python identify_outliers.py
  python identify_outliers.py --min_cells 20 --zscore_thresh 3
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import os
from scipy import sparse, stats

parser = argparse.ArgumentParser()
parser.add_argument('--h5ad_dir',      default='/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update')
parser.add_argument('--checkpoints',   default='/scratch/easmit31/variability/dglm/checkpoints')
parser.add_argument('--outdir',        default='/scratch/easmit31/variability/dglm/checkpoints')
parser.add_argument('--min_cells',     type=int,   default=10, help='min cells per animal x ct x region')
parser.add_argument('--zscore_thresh', type=float, default=3.0, help='SD threshold for expression outlier')
parser.add_argument('--known_outliers',type=str,   default='8H2:ACC',
                    help='comma-separated animal:region pairs to always flag')
parser.add_argument('--outfmt',        default='png')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

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
CELL_TYPES  = list(H5AD_MAP.keys())
REGIONS     = ['ACC','CN','dlPFC','EC','HIP','IPP','lCb','M1','MB','mdTN','NAc']

# parse known outliers
known = {}
for pair in args.known_outliers.split(','):
    if ':' in pair:
        animal, region = pair.strip().split(':')
        known.setdefault(region, set()).add(animal)

all_flags = []

for ct in CELL_TYPES:
    h5ad_path = os.path.join(args.h5ad_dir, H5AD_MAP[ct])
    print(f'\nLoading {ct}')
    adata = ad.read_h5ad(h5ad_path, backed='r')
    obs   = adata.obs.copy()
    obs['_idx'] = np.arange(len(obs))

    # age filter
    obs = obs[obs['age'].astype(float) >= 1.0]

    # opc/oligo split
    if ct == 'opc':
        obs = obs[obs['louvain'].astype(str).isin(OPC_LOUVAIN)]
    elif ct == 'oligodendrocytes':
        obs = obs[~obs['louvain'].astype(str).isin(OPC_LOUVAIN)]

    for region in REGIONS:
        obs_r   = obs[obs['region'] == region]
        animals = obs_r['animal_id'].unique()
        if len(animals) < 5:
            continue

        # load pseudobulk CSV for expression outlier detection
        pb_file = os.path.join(args.checkpoints, f'{ct}_{region}_pseudobulk.csv')
        if os.path.exists(pb_file):
            pb = pd.read_csv(pb_file, index_col=0)
        else:
            pb = None

        for animal in animals:
            obs_a  = obs_r[obs_r['animal_id'] == animal]
            n_cells = len(obs_a)
            flags  = []

            # flag 1: low cell count
            if n_cells < args.min_cells:
                flags.append(f'low_cells({n_cells})')

            # flag 2: known outlier
            if region in known and animal in known[region]:
                flags.append('known_outlier')

            # flag 3: expression outlier
            if pb is not None and animal in pb.columns:
                animal_expr = pb[animal].values
                group_mean  = pb.drop(columns=[animal]).mean(axis=1).values
                group_std   = pb.drop(columns=[animal]).std(axis=1).values
                z_scores    = (animal_expr - group_mean) / (group_std + 1e-10)
                mean_abs_z  = np.nanmean(np.abs(z_scores))
                if mean_abs_z > args.zscore_thresh:
                    flags.append(f'expression_outlier(mean_abs_z={mean_abs_z:.2f})')

            if flags:
                all_flags.append({
                    'animal_id': animal,
                    'cell_type': ct,
                    'region':    region,
                    'n_cells':   n_cells,
                    'flags':     ';'.join(flags)
                })

    adata.file.close()

# save report
report = pd.DataFrame(all_flags)
out_tsv = os.path.join(args.outdir, 'outlier_report.tsv')
report.to_csv(out_tsv, sep='\t', index=False)
print(f'\nTotal flagged: {len(report)}')
print(f'Saved: {out_tsv}')

# summary heatmap: n flags per animal x region
if len(report) > 0:
    animals_all = sorted(report['animal_id'].unique())
    summary = pd.DataFrame(0, index=animals_all, columns=REGIONS)
    for _, row in report.iterrows():
        summary.loc[row['animal_id'], row['region']] += 1

    fig, ax = plt.subplots(figsize=(len(REGIONS)*0.8+2, len(animals_all)*0.3+2))
    im = ax.imshow(summary.values, cmap='Reds', aspect='auto',
                   vmin=0, vmax=summary.values.max())
    for i, animal in enumerate(animals_all):
        for j, region in enumerate(REGIONS):
            v = summary.loc[animal, region]
            if v > 0:
                ax.text(j, i, str(v), ha='center', va='center', fontsize=6)
    ax.set_xticks(range(len(REGIONS)))
    ax.set_xticklabels(REGIONS, fontsize=9, rotation=45, ha='right')
    ax.set_yticks(range(len(animals_all)))
    ax.set_yticklabels(animals_all, fontsize=7)
    plt.colorbar(im, ax=ax, shrink=0.5, label='n flags')
    ax.set_title(f'Outlier flags per animal x region\n(min_cells={args.min_cells}, zscore_thresh={args.zscore_thresh})',
                 fontsize=10)
    plt.tight_layout()
    out_fig = os.path.join(args.outdir, f'outlier_summary.{args.outfmt}')
    plt.savefig(out_fig, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_fig}')

print('done.')
EOFcat > /scratch/easmit31/variability/dglm/scripts/identify_outliers.py << 'EOF'
"""
identify_outliers.py

Identifies outlier animals per cell type x region based on:
  1. Low cell count (< --min_cells threshold)
  2. Expression outliers (pseudobulk mean z-score > --zscore_thresh SD from group)
  3. Known outliers from --known_outliers list

Produces:
  - outlier_report.tsv: all flagged animal x ct x region combinations
  - outlier_summary.png: heatmap of n outlier flags per animal x region

Usage:
  python identify_outliers.py
  python identify_outliers.py --min_cells 20 --zscore_thresh 3
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import os
from scipy import sparse, stats

parser = argparse.ArgumentParser()
parser.add_argument('--h5ad_dir',      default='/scratch/nsnyderm/u01/intermediate_files/cell-class_h5ad_update')
parser.add_argument('--checkpoints',   default='/scratch/easmit31/variability/dglm/checkpoints')
parser.add_argument('--outdir',        default='/scratch/easmit31/variability/dglm/checkpoints')
parser.add_argument('--min_cells',     type=int,   default=10, help='min cells per animal x ct x region')
parser.add_argument('--zscore_thresh', type=float, default=3.0, help='SD threshold for expression outlier')
parser.add_argument('--known_outliers',type=str,   default='8H2:ACC',
                    help='comma-separated animal:region pairs to always flag')
parser.add_argument('--outfmt',        default='png')
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

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
CELL_TYPES  = list(H5AD_MAP.keys())
REGIONS     = ['ACC','CN','dlPFC','EC','HIP','IPP','lCb','M1','MB','mdTN','NAc']

# parse known outliers
known = {}
for pair in args.known_outliers.split(','):
    if ':' in pair:
        animal, region = pair.strip().split(':')
        known.setdefault(region, set()).add(animal)

all_flags = []

for ct in CELL_TYPES:
    h5ad_path = os.path.join(args.h5ad_dir, H5AD_MAP[ct])
    print(f'\nLoading {ct}')
    adata = ad.read_h5ad(h5ad_path, backed='r')
    obs   = adata.obs.copy()
    obs['_idx'] = np.arange(len(obs))

    # age filter
    obs = obs[obs['age'].astype(float) >= 1.0]

    # opc/oligo split
    if ct == 'opc':
        obs = obs[obs['louvain'].astype(str).isin(OPC_LOUVAIN)]
    elif ct == 'oligodendrocytes':
        obs = obs[~obs['louvain'].astype(str).isin(OPC_LOUVAIN)]

    for region in REGIONS:
        obs_r   = obs[obs['region'] == region]
        animals = obs_r['animal_id'].unique()
        if len(animals) < 5:
            continue

        # load pseudobulk CSV for expression outlier detection
        pb_file = os.path.join(args.checkpoints, f'{ct}_{region}_pseudobulk.csv')
        if os.path.exists(pb_file):
            pb = pd.read_csv(pb_file, index_col=0)
        else:
            pb = None

        for animal in animals:
            obs_a  = obs_r[obs_r['animal_id'] == animal]
            n_cells = len(obs_a)
            flags  = []

            # flag 1: low cell count
            if n_cells < args.min_cells:
                flags.append(f'low_cells({n_cells})')

            # flag 2: known outlier
            if region in known and animal in known[region]:
                flags.append('known_outlier')

            # flag 3: expression outlier
            if pb is not None and animal in pb.columns:
                animal_expr = pb[animal].values
                group_mean  = pb.drop(columns=[animal]).mean(axis=1).values
                group_std   = pb.drop(columns=[animal]).std(axis=1).values
                z_scores    = (animal_expr - group_mean) / (group_std + 1e-10)
                mean_abs_z  = np.nanmean(np.abs(z_scores))
                if mean_abs_z > args.zscore_thresh:
                    flags.append(f'expression_outlier(mean_abs_z={mean_abs_z:.2f})')

            if flags:
                all_flags.append({
                    'animal_id': animal,
                    'cell_type': ct,
                    'region':    region,
                    'n_cells':   n_cells,
                    'flags':     ';'.join(flags)
                })

    adata.file.close()

# save report
report = pd.DataFrame(all_flags)
out_tsv = os.path.join(args.outdir, 'outlier_report.tsv')
report.to_csv(out_tsv, sep='\t', index=False)
print(f'\nTotal flagged: {len(report)}')
print(f'Saved: {out_tsv}')

# summary heatmap: n flags per animal x region
if len(report) > 0:
    animals_all = sorted(report['animal_id'].unique())
    summary = pd.DataFrame(0, index=animals_all, columns=REGIONS)
    for _, row in report.iterrows():
        summary.loc[row['animal_id'], row['region']] += 1

    fig, ax = plt.subplots(figsize=(len(REGIONS)*0.8+2, len(animals_all)*0.3+2))
    im = ax.imshow(summary.values, cmap='Reds', aspect='auto',
                   vmin=0, vmax=summary.values.max())
    for i, animal in enumerate(animals_all):
        for j, region in enumerate(REGIONS):
            v = summary.loc[animal, region]
            if v > 0:
                ax.text(j, i, str(v), ha='center', va='center', fontsize=6)
    ax.set_xticks(range(len(REGIONS)))
    ax.set_xticklabels(REGIONS, fontsize=9, rotation=45, ha='right')
    ax.set_yticks(range(len(animals_all)))
    ax.set_yticklabels(animals_all, fontsize=7)
    plt.colorbar(im, ax=ax, shrink=0.5, label='n flags')
    ax.set_title(f'Outlier flags per animal x region\n(min_cells={args.min_cells}, zscore_thresh={args.zscore_thresh})',
                 fontsize=10)
    plt.tight_layout()
    out_fig = os.path.join(args.outdir, f'outlier_summary.{args.outfmt}')
    plt.savefig(out_fig, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_fig}')

print('done.')
