#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument(
    "--input",
    default="/scratch/easmit31/dispersion/dglm/scripts/mean_expression_age_effects.csv"
)
parser.add_argument(
    "--outdir",
    default="/scratch/easmit31/dispersion/dglm/scripts/mean_expression_volcanoes"
)
parser.add_argument(
    "--facet",
    action="store_true",
    help="Make one plot per cell type with regions as separate panels"
)
parser.add_argument(
    "--top",
    type=int,
    default=5,
    help="Number of top increasing/decreasing genes to label"
)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

df = pd.read_csv(args.input)
df["neglog10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
df["significant"] = df["qvalue"] < 0.05

def clean_name(x):
    return str(x).replace("/", "_").replace(" ", "_")

def plot_one(d, title, outfile):
    fig, ax = plt.subplots(figsize=(7, 6))

    sns.scatterplot(
        data=d,
        x="beta",
        y="neglog10p",
        hue="significant",
        palette={False: "lightgray", True: "firebrick"},
        s=12,
        linewidth=0,
        legend=False,
        ax=ax
    )

    ax.axvline(0, color="black", linewidth=0.8)
    ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8)

    # Label strongest positive and negative effects among significant genes
    sig = d[d["significant"]].copy()
    up = sig[sig["beta"] > 0].nlargest(args.top, "beta")
    down = sig[sig["beta"] < 0].nsmallest(args.top, "beta")

    for _, r in pd.concat([up, down]).iterrows():
        label = r["symbol"] if pd.notna(r["symbol"]) else r["ensembl_id"]
        ax.annotate(
            label,
            (r["beta"], r["neglog10p"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7
        )

    ax.set_xlabel("Age effect on mean expression (β)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(title)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()

if args.facet:
    # One figure per cell type, regions as panels
    for cell_type, dct in df.groupby("cell_type"):
        regions = sorted(dct["region"].dropna().unique())
        n = len(regions)
        ncols = min(3, n)
        nrows = int(np.ceil(n / ncols))

        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(6 * ncols, 5 * nrows),
            squeeze=False
        )

        for ax, region in zip(axes.flat, regions):
            d = dct[dct["region"] == region]

            sns.scatterplot(
                data=d,
                x="beta",
                y="neglog10p",
                hue="significant",
                palette={False: "lightgray", True: "firebrick"},
                s=10,
                linewidth=0,
                legend=False,
                ax=ax
            )

            ax.axvline(0, color="black", linewidth=0.8)
            ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8)
            ax.set_title(region)
            ax.set_xlabel("Age effect (β)")
            ax.set_ylabel("-log10(p-value)")

        for ax in axes.flat[n:]:
            ax.axis("off")

        fig.suptitle(cell_type, fontsize=14)
        plt.tight_layout()

        outfile = os.path.join(
            args.outdir,
            f"{clean_name(cell_type)}_mean_expression_volcano_faceted.png"
        )
        plt.savefig(outfile, dpi=300, bbox_inches="tight")
        plt.close()

else:
    # Default: one plot per cell type x region
    for (cell_type, region), d in df.groupby(["cell_type", "region"]):
        title = f"{cell_type} × {region}"

        outfile = os.path.join(
            args.outdir,
            f"{clean_name(cell_type)}__{clean_name(region)}_mean_expression_volcano.png"
        )

        plot_one(d, title, outfile)

print(f"Saved plots to: {args.outdir}")
