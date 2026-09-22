#!/usr/bin/env python3

import os
import pandas as pd

INPUT = "/scratch/easmit31/dispersion/dglm/disp_age__raw_counts"
OUTPUT = "/scratch/easmit31/dispersion/dglm/scripts/mean_expression_age_effects.csv"

files = []

for root, _, filenames in os.walk(INPUT):
    for f in filenames:
        if f == "master_dglm_combined.tsv":
            files.append(os.path.join(root, f))

if not files:
    raise FileNotFoundError(f"No master_dglm_combined.tsv files found under {INPUT}")

dfs = []

for f in files:
    df = pd.read_csv(f, sep="\t")

    required = ["ensembl_id", "symbol", "cell_type", "region", "beta", "pvalue", "qvalue"]
    missing = [x for x in required if x not in df.columns]

    if missing:
        print(f"Skipping {f}: missing {missing}")
        continue

    x = df[required].copy()
    x["original_subtype"] = x["cell_type"].astype(str)
    dfs.append(x)

results = pd.concat(dfs, ignore_index=True)

results = results.drop_duplicates(
    subset=["original_subtype", "region", "ensembl_id"]
)

def relabel_subtype(x):
    if x.startswith("opc-olig_"):
        n = x.replace("opc-olig_", "")
        if n in {"12", "13"}:
            return f"opc_{n}"
        return f"olig_{n}"
    return x

results["cell_type"] = results["original_subtype"].map(relabel_subtype)

results["parent_cell_type"] = results["cell_type"].map(
    lambda x:
        "opc" if x.startswith("opc_") else
        "olig" if x.startswith("olig_") else
        x
)

results["direction"] = results["beta"].map(
    lambda x:
        "increases_with_age" if x > 0 else
        "decreases_with_age" if x < 0 else
        "no_change"
)

results["rank"] = pd.NA

up = results["direction"] == "increases_with_age"
down = results["direction"] == "decreases_with_age"

results.loc[up, "rank"] = (
    results.loc[up]
    .groupby(["cell_type", "region"])["beta"]
    .rank(method="first", ascending=False)
)

results.loc[down, "rank"] = (
    results.loc[down]
    .groupby(["cell_type", "region"])["beta"]
    .rank(method="first", ascending=True)
)

results = results.sort_values(
    ["cell_type", "region", "direction", "rank"]
)

results.to_csv(OUTPUT, index=False)

print(f"Wrote {len(results):,} rows to:")
print(OUTPUT)
print()
print("Subtype labels:")
print(
    results[["original_subtype", "cell_type", "parent_cell_type"]]
    .drop_duplicates()
    .sort_values(["parent_cell_type", "cell_type"])
    .to_string(index=False)
)
