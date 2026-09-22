#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

GO = "/scratch/easmit31/dispersion/dglm/scripts/go_enrichment_results/mean_expression/dglm_mean_expression_GO_results.csv"
OUT = "/scratch/easmit31/dispersion/dglm/scripts/go_enrichment_results/mean_expression/plots"

os.makedirs(OUT, exist_ok=True)

d = pd.read_csv(GO)

d = d[d["qvalue"] < 0.05].copy()

if len(d) == 0:
    raise SystemExit("No GO terms with q < 0.05")

d["neglog10q"] = -np.log10(d["qvalue"].clip(lower=1e-300))

for ct in sorted(d["cell_type"].unique()):

    x = d[d["cell_type"] == ct].copy()

    # Top 5 terms per subtype x region x direction
    x = (
        x.sort_values("qvalue")
         .groupby(["region", "direction"], group_keys=False)
         .head(5)
    )

    x["label"] = (
        x["region"] + " | " +
        x["direction"].str.replace("_with_age", "", regex=False) +
        " | " +
        x["go_name"]
    )

    x = x.sort_values("neglog10q")

    fig, ax = plt.subplots(figsize=(11, max(6, 0.32 * len(x))))

    ax.barh(
        np.arange(len(x)),
        x["neglog10q"]
    )

    ax.set_yticks(np.arange(len(x)))
    ax.set_yticklabels(x["label"], fontsize=8)

    ax.set_xlabel("-log10(FDR)")
    ax.set_title(f"{ct}: mean-expression age effects — GO Biological Process")

    ax.axvline(-np.log10(0.05), linestyle="--", linewidth=1)

    plt.tight_layout()

    outfile = os.path.join(
        OUT,
        f"mean_expression_GO_{ct}.png"
    )

    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()

    print(outfile)
