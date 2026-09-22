#!/usr/bin/env Rscript
# Filters a raw-count pseudobulk matrix: keep gene if CPM >= cutoff in at
# least n_min animals (not just any single animal -- that degenerate case
# is why the earlier group=animal_id version couldn't distinguish 0.5 vs 1).
#
# Usage:
#   Rscript filter_manual_cutoff.R <pseudobulk_csv> <outdir> [n_min] [label]

suppressMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: filter_manual_cutoff.R <pseudobulk_csv> <outdir> [n_min] [label]")
csv_path <- args[1]
outdir   <- args[2]
n_min    <- if (length(args) >= 3) as.numeric(args[3]) else 5
label    <- if (length(args) >= 4) args[4] else tools::file_path_sans_ext(basename(csv_path))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

counts <- read.csv(csv_path, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
cat(sprintf("Loaded %s: %d genes x %d animals\n", csv_path, nrow(counts), ncol(counts)))
cat(sprintf("n_min = %d (gene must clear cutoff CPM in at least this many animals)\n", n_min))

y <- DGEList(counts = counts)
cpms <- edgeR::cpm(y)

for (cutoff in c(0.5, 1)) {
  keep <- rowSums(cpms >= cutoff) >= n_min
  y_filt <- y[keep, , keep.lib.sizes = FALSE]
  cat(sprintf("cutoff=%s: kept %d / %d genes\n", cutoff, sum(keep), length(keep)))

  out_df <- as.data.frame(y_filt$counts)
  out_path <- file.path(outdir, sprintf("%s_filtered_cutoff%s.csv", label, cutoff))
  write.csv(out_df, out_path, row.names = TRUE)
  cat(sprintf("Saved: %s\n", out_path))
}
