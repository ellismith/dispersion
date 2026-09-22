#!/usr/bin/env python3
"""
Plot log10(median CPM) distributions from pseudobulk CSVs (genes x animals).
Fixed CPM thresholds (0.5, 1, 2) are marked as vertical lines.

Modes:
  --csv PATH                            single cell_type x region file
  --cell_type NAME --aggregate_regions  pool all region files for one cell type
  --all_cell_types                      grid of subplots, one per cell type,
                                         each aggregated across all its regions
  --all_cell_types --per_region         one PNG per cell type, each a grid of
                                         subplots (one per region for that type)

--input_dir defaults to ../disp_age__mean_full relative to this script's own
location, so it works regardless of which directory you run it from.

Usage:
  python plot_cpm_valley.py --csv <path> [--out OUT.png] [--label LABEL]
  python plot_cpm_valley.py --cell_type GABAergic_neurons --aggregate_regions [--out OUT.png]
  python plot_cpm_valley.py --all_cell_types [--out OUT.png]
  python plot_cpm_valley.py --all_cell_types --per_region [--out_dir DIR]
"""
import argparse
import glob
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_INPUT_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "disp_age__mean_full"))

CELL_TYPES = [
    "astrocytes", "ependymal_cells", "GABAergic_neurons", "microglia",
    "oligodendrocytes", "opc", "vascular_cells", "basket_cells",
    "cerebellar_neurons", "medium_spiny_neurons", "midbrain_neurons",
    "glutamatergic_neurons",
]

FIXED_THRESHOLDS = [0.5, 1, 2]
HIST_COLOR = "#4d4d4d"  # darker gray


def cpm_normalize(df):
    return df.div(df.sum(axis=0), axis=1) * 1e6


def median_log10_cpm_from_csv(csv_path):
    df = pd.read_csv(csv_path, index_col=0)
    cpm = cpm_normalize(df)
    median_cpm = cpm.median(axis=1)
    median_cpm = median_cpm[median_cpm > 0]
    return np.log10(median_cpm)


def find_region_files(input_dir, cell_type):
    pattern = os.path.join(input_dir, f"{cell_type}_*_pseudobulk.csv")
    return sorted(glob.glob(pattern))


def region_from_filename(path, cell_type):
    base = os.path.basename(path)
    base = base[len(cell_type) + 1:]  # strip "<cell_type>_"
    return base.replace("_pseudobulk.csv", "")


def discover_cell_types(input_dir):
    return [ct for ct in CELL_TYPES if find_region_files(input_dir, ct)]


def median_log10_cpm_aggregated(input_dir, cell_type):
    files = find_region_files(input_dir, cell_type)
    if not files:
        pattern = os.path.join(input_dir, f"{cell_type}_*_pseudobulk.csv")
        raise FileNotFoundError(f"No files matched pattern: {pattern}")
    cpm_frames = [cpm_normalize(pd.read_csv(f, index_col=0)) for f in files]
    combined = pd.concat(cpm_frames, axis=1)
    median_cpm = combined.median(axis=1)
    median_cpm = median_cpm[median_cpm > 0]
    return np.log10(median_cpm), len(files)


def draw_distribution(ax, log10_vals, title):
    ax.hist(log10_vals, bins=80, density=True, color=HIST_COLOR, alpha=0.85)
    for thresh, color in zip(FIXED_THRESHOLDS, ["tab:blue", "tab:green", "tab:orange"]):
        ax.axvline(np.log10(thresh), color=color, linestyle="--", lw=1.2,
                   label=f"{thresh} CPM")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("log10(median CPM)")
    ax.set_ylabel("density")


def grid_shape(n, ncols=4):
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    return nrows, ncols


def plot_cell_type_per_region(input_dir, cell_type, out_dir):
    files = find_region_files(input_dir, cell_type)
    if not files:
        raise FileNotFoundError(f"No region files found for '{cell_type}' in {input_dir}")

    nrows, ncols = grid_shape(len(files))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.2 * nrows))
    axes = np.array(axes).reshape(-1)

    for i, f in enumerate(files):
        log10_vals = median_log10_cpm_from_csv(f)
        region = region_from_filename(f, cell_type)
        draw_distribution(axes[i], log10_vals.values,
                           f"{region}\n(genes={len(log10_vals)})")

    for j in range(len(files), len(axes)):
        axes[j].axis("off")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9)
    fig.suptitle(f"{cell_type}: median-CPM distribution by region")
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])

    out_path = os.path.join(out_dir, f"{cell_type}_per_region_cpm.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", default=None)
    parser.add_argument("--input_dir", default=DEFAULT_INPUT_DIR,
                         help=f"Directory of pseudobulk CSVs (default: {DEFAULT_INPUT_DIR})")
    parser.add_argument("--cell_type", default=None)
    parser.add_argument("--aggregate_regions", action="store_true")
    parser.add_argument("--all_cell_types", action="store_true")
    parser.add_argument("--per_region", action="store_true",
                         help="With --all_cell_types: one PNG per cell type, subplots per region")
    parser.add_argument("--out", default=None)
    parser.add_argument("--out_dir", default=None,
                         help="Output directory for --all_cell_types --per_region (default: cwd)")
    parser.add_argument("--label", default=None)
    args = parser.parse_args()

    args.input_dir = os.path.abspath(args.input_dir)
    if args.csv:
        args.csv = os.path.abspath(args.csv)
    if args.out:
        args.out = os.path.abspath(args.out)
    out_dir = os.path.abspath(args.out_dir) if args.out_dir else os.getcwd()

    if args.all_cell_types:
        if not os.path.isdir(args.input_dir):
            parser.error(f"--input_dir does not exist: {args.input_dir}")
        cell_types = discover_cell_types(args.input_dir)
        if not cell_types:
            parser.error(f"No pseudobulk files found in {args.input_dir}")

        if args.per_region:
            os.makedirs(out_dir, exist_ok=True)
            for ct in cell_types:
                out_path = plot_cell_type_per_region(args.input_dir, ct, out_dir)
                print(f"Saved: {out_path}")
            return

        n = len(cell_types)
        nrows, ncols = grid_shape(n)
        fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.2 * nrows))
        axes = np.array(axes).reshape(-1)

        for i, ct in enumerate(cell_types):
            log10_vals, n_regions = median_log10_cpm_aggregated(args.input_dir, ct)
            draw_distribution(axes[i], log10_vals.values,
                               f"{ct}\n(regions={n_regions}, genes={len(log10_vals)})")

        for j in range(len(cell_types), len(axes)):
            axes[j].axis("off")

        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9)
        fig.suptitle("Median-CPM distributions by cell type (aggregated across regions)")
        fig.tight_layout(rect=[0, 0.04, 1, 0.96])

        out_path = args.out or os.path.abspath("all_celltypes_cpm_aggregated.png")
        fig.savefig(out_path, dpi=150)
        print(f"Saved: {out_path}")
        return

    if args.cell_type:
        if not args.aggregate_regions:
            parser.error("--cell_type requires --aggregate_regions (use --csv for a single cell_type x region file)")
        if not os.path.isdir(args.input_dir):
            parser.error(f"--input_dir does not exist: {args.input_dir}")
        log10_vals, n_regions = median_log10_cpm_aggregated(args.input_dir, args.cell_type)
        label = args.label or args.cell_type
        out_path = args.out or os.path.abspath(f"{label}_aggregated_cpm.png")
        fig, ax = plt.subplots(figsize=(8, 5))
        draw_distribution(ax, log10_vals.values,
                           f"{label}: median-CPM across all regions (regions={n_regions}, genes={len(log10_vals)})")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        print(f"Saved: {out_path}")
        return

    if args.csv:
        if not os.path.isfile(args.csv):
            parser.error(f"--csv file does not exist: {args.csv}")
        log10_vals = median_log10_cpm_from_csv(args.csv)
        label = args.label or os.path.basename(args.csv).replace("_pseudobulk.csv", "")
        out_path = args.out or os.path.abspath(f"{label}_cpm.png")
        fig, ax = plt.subplots(figsize=(8, 5))
        draw_distribution(ax, log10_vals.values,
                           f"{label}: median-CPM distribution across genes (n={len(log10_vals)})")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out_path, dpi=150)
        print(f"Saved: {out_path}")
        return

    parser.error("Must specify one of: --csv, --cell_type with --aggregate_regions, or --all_cell_types")


if __name__ == "__main__":
    main()
