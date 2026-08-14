#!/usr/bin/env Rscript
#
# filter_percent_animals.R
#
# Gene-level low-expression filter for a raw-count pseudobulk matrix
# (genes x animals). A gene is kept if its CPM exceeds the cutoff in at
# least 50% of animals in this file.
#
# LOGIC:
#   n_above <- rowSums(cpms > cutoff)      # how many animals clear the cutoff, per gene
#   keep    <- n_above >= ceiling(0.5 * n_animals)
#
# This is mathematically equivalent, for an ODD number of animals, to
# "median CPM across all animals > cutoff" -- the median of n (odd) values
# is literally the middle-ranked value, so if it's above the cutoff, at
# least half the values must also be above it, and vice versa. Verified
# directly on real data (astrocytes/louvain-8/ACC test case): computing
# both the proportion rule and the median rule independently produced
# ZERO discordant genes. For an even number of animals the two can
# technically diverge at the exact boundary, but that's not the case we're
# running here. Kept as the proportion-rule form since that's what was
# explicitly requested, not because it differs from the median rule.
#
# CPM is computed from RAW (not TMM-normalized) library sizes
# (cpm(y, normalized.lib.sizes = FALSE)) -- calcNormFactors() has not been
# run anywhere in this pipeline yet. That's a deliberate, separate decision
# for later (TMM normalization happens downstream of this filtering step,
# not as part of it).
#
# OUTPUT per cutoff:
#   {label}_filtered_cutoff{cutoff}.csv       -- filtered raw-count matrix (genes x animals)
#   {label}_filter_stats_cutoff{cutoff}.csv   -- audit trail: one row per gene, with
#                                                 n_animals_above, median_cpm, and keep,
#                                                 so any future "why did this gene get
#                                                 dropped/kept" question can be answered
#                                                 by reading this file directly instead of
#                                                 re-deriving it from the raw counts.
#
# Usage:
#   Rscript filter_percent_animals.R <pseudobulk_csv> <outdir> [label]

suppressMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: filter_percent_animals.R <pseudobulk_csv> <outdir> [label]")
csv_path <- args[1]
outdir   <- args[2]
label    <- if (length(args) >= 3) args[3] else tools::file_path_sans_ext(basename(csv_path))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- load raw counts (genes x animals) ---
counts <- read.csv(csv_path, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
n_animals <- ncol(counts)
min_animals <- ceiling(0.5 * n_animals)  # ceiling so "50%" is unambiguous for odd n_animals

cat(sprintf("Loaded %s: %d genes x %d animals\n", csv_path, nrow(counts), n_animals))
cat(sprintf("50%% expression criterion requires >= %d animals\n", min_animals))

# --- CPM, from raw library sizes (no TMM at this stage) ---
y <- DGEList(counts = counts)
cpms <- edgeR::cpm(y, normalized.lib.sizes = FALSE)

# also compute the median-CPM version for the audit trail, since it's the
# same underlying quantity most people reading the stats file will expect
# to see, even though the KEEP decision itself is made via the proportion rule
median_cpm <- apply(cpms, 1, median)

cutoffs <- c(0.5, 1, 2)

for (cutoff in cutoffs) {
  n_above <- rowSums(cpms > cutoff)
  keep <- n_above >= min_animals

  cat(sprintf("cutoff=%s: kept %d / %d genes\n", cutoff, sum(keep), length(keep)))

  # --- filtered raw-count matrix ---
  y_filt <- y[keep, , keep.lib.sizes = FALSE]
  out_df <- as.data.frame(y_filt$counts)
  out_path <- file.path(outdir, sprintf("%s_filtered_cutoff%s.csv", label, cutoff))
  write.csv(out_df, out_path, row.names = TRUE)
  cat(sprintf("Saved: %s\n", out_path))

  # --- audit trail: per-gene stats behind this cutoff's keep decision ---
  filter_stats <- data.frame(
    gene = rownames(counts),
    n_animals_above_cutoff = n_above,
    median_cpm = median_cpm,
    keep = keep
  )
  stats_path <- file.path(outdir, sprintf("%s_filter_stats_cutoff%s.csv", label, cutoff))
  write.csv(filter_stats, stats_path, row.names = FALSE)
  cat(sprintf("Saved: %s\n", stats_path))
}
