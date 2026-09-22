#!/usr/bin/env python3
"""
Same before/after-filtering gene counts as plot_genes_before_after.py, but
broken into one subplot per region -- each subplot's boxes are across that
region's louvain subclusters only. Outliers (matplotlib's standard 1.5*IQR
rule, computed within each region separately) are drawn as colored points
labeled with the SUBCLUSTER they came from (region is already fixed by the
subplot itself, so labeling with region would be redundant).

Usage:
  python plot_genes_before_after_per_region.py --raw_dir <dir> --filtered_dir <dir> --out <out.png> [--label "astrocytes"]
"""
import argparse
import glob
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REGIONS = ['ACC','CN','dlPFC','EC','HIP','IPP','lCb','M1','MB','mdTN','NAc']
CUTOFFS = [0.5, 1, 2]
OUTLIER_COLOR = "#D62728"


def count_genes(csv_path):
    with open(csv_path) as f:
        n = sum(1 for _ in f)
    return n - 1


def region_from_filename(path):
    base = os.path.basename(path).replace("_pseudobulk.csv", "")
    for region in REGIONS:
        if base.endswith(f"_{region}"):
            return region, base[:-(len(region) + 1)]  # (region, subcluster label)
    return None, None


def iqr_outliers(values):
    values = np.asarray(values)
    q1, q3 = np.percentile(values, [25, 75])
    iqr = q3 - q1
    lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    return (values < lo) | (values > hi)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw_dir", required=True)
    parser.add_argument("--filtered_dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--label", default=None)
    args = parser.parse_args()

    raw_files = sorted(glob.glob(os.path.join(args.raw_dir, "*_pseudobulk.csv")))
    if not raw_files:
        raise SystemExit(f"No *_pseudobulk.csv files found in {args.raw_dir}")

    by_region = {r: {"before": [], **{c: [] for c in CUTOFFS}} for r in REGIONS}

    for raw_path in raw_files:
        region, subcluster = region_from_filename(raw_path)
        if region is None:
            continue
        label = os.path.basename(raw_path).replace("_pseudobulk.csv", "")
        by_region[region]["before"].append((subcluster, count_genes(raw_path)))
        for c in CUTOFFS:
            filt_path = os.path.join(args.filtered_dir, f"{label}_filtered_cutoff{c}.csv")
            if os.path.isfile(filt_path):
                by_region[region][c].append((subcluster, count_genes(filt_path)))

    used_regions = [r for r in REGIONS if len(by_region[r]["before"]) > 0]
    n = len(used_regions)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.6 * nrows))
    axes = np.array(axes).reshape(-1)

    tick_labels = ["before"] + [f"{c}" for c in CUTOFFS]

    for i, region in enumerate(used_regions):
        stages = [by_region[region]["before"]] + [by_region[region][c] for c in CUTOFFS]
        data = [[v for _, v in stage] for stage in stages]
        n_combos = len(stages[0])

        bp = axes[i].boxplot(data, tick_labels=tick_labels, patch_artist=True,
                              showmeans=True, showfliers=False)
        for patch in bp['boxes']:
            patch.set_facecolor("#888888")
            patch.set_alpha(0.6)

        for x_pos, stage in enumerate(stages, start=1):
            if len(stage) < 4:
                continue
            stage_subclusters, stage_values = zip(*stage)
            stage_values = np.asarray(stage_values)
            outlier_mask = iqr_outliers(stage_values)
            outlier_subclusters = np.array(stage_subclusters)[outlier_mask]
            outlier_values = stage_values[outlier_mask]
            for subcluster_name, val in zip(outlier_subclusters, outlier_values):
                axes[i].scatter(x_pos, val, color=OUTLIER_COLOR, zorder=5, s=25)
                axes[i].annotate(subcluster_name.split('_')[-1], (x_pos, val),
                                  textcoords="offset points", xytext=(5, 0),
                                  fontsize=6, color=OUTLIER_COLOR)

        axes[i].set_title(f"{region} (n={n_combos} louvains)", fontsize=9)
        axes[i].set_ylabel("n genes")
        axes[i].tick_params(axis='x', labelsize=8)

    for j in range(len(used_regions), len(axes)):
        axes[j].axis("off")

    title_label = args.label or os.path.basename(args.raw_dir.rstrip("/"))
    fig.suptitle(f"{title_label}: genes before/after filtering, by region\n"
                 f"(box=median, green triangle=mean, red=outlier (1.5x IQR within region), "
                 f"labeled by louvain subcluster; x-axis: before filter, cutoff 0.5/1/2)")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(args.out, dpi=150)
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
