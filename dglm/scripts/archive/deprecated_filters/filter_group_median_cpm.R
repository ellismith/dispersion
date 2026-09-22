#!/usr/bin/env Rscript
# Keep gene if its median CPM across all animals exceeds the cutoff
# (equivalent to: >= 50% of animals express it above cutoff CPM).
# One group containing every animal -- filter_group_median_cpm exactly as
# specified, with group = single pooled group.
#
# Usage:
#   Rscript filter_group_median_cpm.R <pseudobulk_csv> <outdir> [label]

suppressMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: filter_group_median_cpm.R <pseudobulk_csv> <outdir> [label]")
csv_path <- args[1]
outdir   <- args[2]
label    <- if (length(args) >= 3) args[3] else tools::file_path_sans_ext(basename(csv_path))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

counts <- read.csv(csv_path, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
cat(sprintf("Loaded %s: %d genes x %d animals\n", csv_path, nrow(counts), ncol(counts)))

y <- DGEList(counts = counts)

filter_group_median_cpm <- function(y, group, cutoff = 1) {
  cpms <- edgeR::cpm(y)
  group_medians <- sapply(
    unique(group),
    function(g) apply(cpms[, group == g, drop = FALSE], 1, median)
  )
  apply(group_medians > cutoff, 1, any)
}

# single group containing every animal -- median CPM across all animals
group <- rep("all", ncol(counts))

for (cutoff in c(0.5, 1)) {
  keep <- filter_group_median_cpm(y, group = group, cutoff = cutoff)
  y_filt <- y[keep, , keep.lib.sizes = FALSE]
  cat(sprintf("cutoff=%s: kept %d / %d genes\n", cutoff, sum(keep), length(keep)))

  out_df <- as.data.frame(y_filt$counts)
  out_path <- file.path(outdir, sprintf("%s_filtered_cutoff%s.csv", label, cutoff))
  write.csv(out_df, out_path, row.names = TRUE)
  cat(sprintf("Saved: %s\n", out_path))
}
