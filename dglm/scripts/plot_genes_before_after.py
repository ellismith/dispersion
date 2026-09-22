#!/usr/bin/env python3
"""
For every subcluster x region combo in a raw pseudobulk directory, counts
genes before filtering and genes remaining after filter_percent_animals.R
at cutoffs 0.5, 1, 2. Plots the distribution (one box per cutoff, across
all combos) as a single PNG.

Usage:
  python plot_genes_before_after.py --raw_dir <dir> --filtered_dir <dir> --out <out.png> [--label "astrocytes"]
"""
import argparse
import glob
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def count_genes(csv_path):
    """Fast line count minus header, without loading the full file."""
    with open(csv_path) as f:
        n = sum(1 for _ in f)
    return n - 1  # subtract header row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw_dir", required=True, help="Directory with *_pseudobulk.csv (unfiltered)")
    parser.add_argument("--filtered_dir", required=True, help="Directory with *_filtered_cutoff{0.5,1,2}.csv")
    parser.add_argument("--out", required=True)
    parser.add_argument("--label", default=None, help="Cell type label for plot title")
    args = parser.parse_args()

    raw_files = sorted(glob.glob(os.path.join(args.raw_dir, "*_pseudobulk.csv")))
    if not raw_files:
        raise SystemExit(f"No *_pseudobulk.csv files found in {args.raw_dir}")

    cutoffs = [0.5, 1, 2]
    before_counts = []
    after_counts = {c: [] for c in cutoffs}
    missing = []

    for raw_path in raw_files:
        label = os.path.basename(raw_path).replace("_pseudobulk.csv", "")
        n_before = count_genes(raw_path)
        before_counts.append(n_before)

        for c in cutoffs:
            filt_path = os.path.join(args.filtered_dir, f"{label}_filtered_cutoff{c}.csv")
            if os.path.isfile(filt_path):
                after_counts[c].append(count_genes(filt_path))
            else:
                missing.append(f"{label} (cutoff {c})")

    if missing:
        print(f"Warning: {len(missing)} combo/cutoff filtered files not found, e.g.: {missing[:5]}")

    data = [before_counts] + [after_counts[c] for c in cutoffs]
    tick_labels = ["before\nfilter"] + [f"cutoff\n{c}" for c in cutoffs]

    fig, ax = plt.subplots(figsize=(7, 5))
    bp = ax.boxplot(data, labels=tick_labels, patch_artist=True, showmeans=True)
    for patch in bp['boxes']:
        patch.set_facecolor("#888888")
        patch.set_alpha(0.6)

    n_combos = len(raw_files)
    title_label = args.label or os.path.basename(args.raw_dir.rstrip("/"))
    ax.set_title(f"{title_label}: genes before/after filtering (n={n_combos} louvain x region combos)")
    ax.set_ylabel("number of genes")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved: {args.out}")
    print(f"medians -- before: {sorted(before_counts)[len(before_counts)//2]}, "
          + ", ".join(f"cutoff {c}: {sorted(after_counts[c])[len(after_counts[c])//2]}" for c in cutoffs))


if __name__ == "__main__":
    main()
