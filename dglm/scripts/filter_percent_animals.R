#!/usr/bin/env Rscript
#
# filter_percent_animals.R
#
# Gene-level low-expression filter for a raw-count pseudobulk matrix
# (genes x animals). A gene is kept if its CPM is AT LEAST the cutoff
# (>=, not >) in at least 50% of animals in this file.
#
# LOGIC:
#   n_above <- rowSums(cpms >= cutoff)     # how many animals meet/clear the cutoff, per gene
#   keep    <- n_above >= ceiling(0.5 * n_animals)
#
# Both comparisons are inclusive (>=): a gene with CPM exactly equal to the
# cutoff in a given animal counts as meeting it, and a gene meeting the
# cutoff in exactly 50% of animals (for even n_animals) counts as passing.
#
# This is mathematically equivalent, for an ODD number of animals, to
# "median CPM across all animals >= cutoff" -- the median of n (odd) values
# is literally the middle-ranked value, so if it's at or above the cutoff,
# at least half the values must also be at or above it, and vice versa.
# Verified directly on real data (astrocytes/louvain-8/ACC test case):
# computing both the proportion rule and the median rule independently
# produced ZERO discordant genes. Kept as the proportion-rule form since
# that's what was explicitly requested, not because it differs from the
# median rule.
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
#                                                 n_animals_above_cutoff (>=), median_cpm,
#                                                 and keep, so any future "why did this gene
#                                                 get dropped/kept" question can be answered
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
min_animals <- ceiling(0.5 * n_animals)  # ceiling so "at least 50%" is unambiguous for odd n_animals

cat(sprintf("Loaded %s: %d genes x %d animals\n", csv_path, nrow(counts), n_animals))
cat(sprintf("at-least-50%% expression criterion requires >= %d animals\n", min_animals))

# --- CPM, from raw library sizes (no TMM at this stage) ---
y <- DGEList(counts = counts)
cpms <- edgeR::cpm(y, normalized.lib.sizes = FALSE)

# median-CPM version for the audit trail (same underlying quantity, useful
# to have logged even though the KEEP decision is made via the proportion rule)
median_cpm <- apply(cpms, 1, median)

# Standard edgeR low-expression filter, computed once and used as an
# ADDITIONAL required condition alongside the CPM-proportion rule (not a
# replacement for it). Catches genes that clear the CPM-proportion rule
# only because pseudobulk library sizes are small (e.g. a raw count of 1
# can exceed CPM=0.5 in a low-depth animal), which filterByExpr correctly
# rejects using its own count- and library-size-aware logic. Verified
# directly: on two real genes found via visual inspection to be almost
# entirely raw-count 0/1/2 (SNORD81, MADCAM1), both passed the CPM-proportion
# rule alone but were correctly rejected by filterByExpr on its own.
keep_fbe <- suppressWarnings(edgeR::filterByExpr(y))

cutoffs <- c(0.5, 1, 2)

for (cutoff in cutoffs) {
  n_above <- rowSums(cpms >= cutoff)   # >= : "at least" the cutoff
  keep_cpm <- n_above >= min_animals   # >= : "at least" 50% of animals -- UNCHANGED, original rule
  keep <- keep_cpm & keep_fbe          # additionally require filterByExpr

  cat(sprintf("cutoff=%s: kept %d / %d genes (CPM rule alone would keep %d; filterByExpr additionally removed %d)\n",
              cutoff, sum(keep), length(keep), sum(keep_cpm), sum(keep_cpm & !keep_fbe)))

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
    keep_cpm_rule = keep_cpm,
    keep_filterByExpr = keep_fbe,
    keep = keep
  )
  stats_path <- file.path(outdir, sprintf("%s_filter_stats_cutoff%s.csv", label, cutoff))
  write.csv(filter_stats, stats_path, row.names = FALSE)
  cat(sprintf("Saved: %s\n", stats_path))
}
