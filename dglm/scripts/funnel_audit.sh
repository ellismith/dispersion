#!/bin/bash
# =============================================================================
# funnel_audit.sh — READ-ONLY audit of the final DGLM dispersion pipeline
# (mask1_915). Prints data dimensions at every stage and the key parameters
# as written in the code on disk. Writes nothing.
#
# Usage: bash scripts/funnel_audit.sh [cell_type]   (default GABAergic_neurons)
#
# Where the data lives (confirmed for GABAergic):
#   DGLM inputs  : <run_info$outdir>/filtered/<subcluster>_<region>_filtered_cutoff0.5.csv
#                  (raw, pre-CPM: <outdir>/<subcluster>_<region>_pseudobulk.csv)
#   DGLM results : disp_age__mask1_915/<CT>/dglm_checkpoints_cutoff0.5/*_dglm_results.rds
#   mashr input  : disp_age__mask1_915_300c/<CT>/dglm_checkpoints_cutoff0.5/*_dglm_results.rds
#   mashr result : same dir, combined_dglm_mashr_results_qval_single_raw_fill1000.rds
#   final table  : same dir, master_dglm_combined.tsv
# =============================================================================

CT=${1:-GABAergic_neurons}
cd /scratch/easmit31/dispersion/dglm
S=disp_age__mask1_915/${CT}/dglm_checkpoints_cutoff0.5
F=disp_age__mask1_915_300c/${CT}/dglm_checkpoints_cutoff0.5
echo "CELL TYPE: $CT"; echo "DGLM results: $S"; echo "mashr/final: $F"

# -----------------------------------------------------------------------------
# Parameters as written in the code (line numbers cited in the doc)
#   dglm_model.R  : TMM + logCPM (prior 0.5); count<=1 mask; MAD low-point cut;
#                   >=50% detection per region; <10 animals -> skip; gaussian
#                   mean model with log dispersion link; p from dispersion.fit;
#                   qval = BH per subcluster x region
#   dglm_mashr.R  : qval_single strong subset (line 118)
#   dglm_fdr_combined.R : master qvalue = BH across the whole master (line 140)
# -----------------------------------------------------------------------------
echo; echo "=== PARAMS: dglm_model.R ==="
grep -n "calcNormFactors\|prior.count\|counts <= 1\|madv\|lowpt\|detected\|nrow(d) < 10\|family=\|dlink\|dispersion.fit)\|p.adjust" scripts/dglm_model.R | head -20
echo; echo "=== PARAMS: dglm_mashr.R (strong subset) ==="
sed -n '117,118p' scripts/dglm_mashr.R
echo; echo "=== PARAMS: dglm_fdr_combined.R (master qvalue) ==="
sed -n '140p' scripts/dglm_fdr_combined.R

Rscript - "$S" "$F" "$CT" << 'REOF'
args <- commandArgs(trailingOnly=TRUE); S <- args[1]; F <- args[2]; CT <- args[3]
q3 <- function(v) sprintf("min %s | median %s | max %s", min(v), median(v), max(v))
nlines <- function(f) length(readLines(f)) - 1                          # rows minus header = genes
ncols  <- function(f) length(strsplit(readLines(f, n=1), ",")[[1]]) - 1  # columns minus gene id = animals

# STAGES 1-4 — per subcluster x region: raw genes -> CPM-filtered genes -> genes tested.
# The input folder is read from run_info$outdir in each rds, i.e. what DGLM actually used.
fs <- list.files(S, pattern="_dglm_results.rds$", full.names=TRUE)
rows <- list(); meta <- list()
for (f in fs) {
  x <- readRDS(f); a <- x$array; ri <- x$run_info; od <- ri$outdir
  meta[[ri$cell_type]] <- data.frame(sub=ri$cell_type, outdir=od, date=format(ri$date))
  for (r in dimnames(a)[[3]]) {
    stem <- paste0(ri$cell_type, "_", r)
    filt <- file.path(od, "filtered", paste0(stem, "_filtered_cutoff0.5.csv"))
    if (!file.exists(filt)) filt <- file.path(od, paste0(stem, "_filtered_cutoff0.5.csv"))
    raw  <- file.path(od, paste0(stem, "_pseudobulk.csv"))
    rows[[stem]] <- data.frame(sub=ri$cell_type, region=r,
      raw_genes = if (file.exists(raw))  nlines(raw)  else NA,
      cpm_genes = if (file.exists(filt)) nlines(filt) else NA,
      animals   = if (file.exists(filt)) ncols(filt)  else NA,
      tested    = sum(!is.na(a[,1,r])))
  }
}
meta <- do.call(rbind, meta); d <- do.call(rbind, rows)
cat("\n=== RUN INFO (", nrow(meta), "subclusters ) ===\n")
cat("mean model:", paste(x$run_info$covariates, collapse=" + "), " | dispersion:", x$run_info$disp_covariates, "\n")
cat("response:", x$run_info$response, "\n")
cat("input folders:", paste(unique(dirname(meta$outdir)), collapse="; "), "\n")
cat("fit dates:", min(meta$date), "to", max(meta$date), "\n")

cat("\n=== 1. raw pseudobulk (pre-CPM) ===\n")
cat("files found:", sum(!is.na(d$raw_genes)), "of", nrow(d), " genes:", q3(na.omit(d$raw_genes)), "\n")
cat("\n=== 2. after CPM>=0.5 in >=50% animals ===\n")
cat("files found:", sum(!is.na(d$cpm_genes)), "of", nrow(d), " genes:", q3(na.omit(d$cpm_genes)),
    "\n animals per condition:", q3(na.omit(d$animals)), "\n")
cat("\n=== 3-4. after mask + 50% detection + 10-animal floor -> genes tested in DGLM ===\n")
cat("subcluster x region combos in rds:", nrow(d), " with >=1 gene tested:", sum(d$tested > 0), "\n")
cat("genes tested (combos with results):", q3(d$tested[d$tested > 0]), "\n")

# STAGE 5 — condition filter: >=300 cells AND >=30 animals (lCb also out for 5 glial types)
sc <- file.path("disp_age__mask1_915/sparsity_diagnostics", paste0(CT, "_condition_sparsity.csv"))
s5 <- read.csv(sc); s5 <- s5[!is.na(s5$subcluster) & !is.na(s5$region),]
k <- s5[s5$pass_both %in% c("True", TRUE),]
cat("\n=== 5. condition filter (", sc, ") ===\n")
cat("total:", nrow(s5), " pass cells:", sum(s5$pass_cells %in% c("True",TRUE)),
    " pass animals:", sum(s5$pass_animals %in% c("True",TRUE)), " pass both:", nrow(k), "\n")
cat("kept: n_cells", q3(k$n_cells), "\n      n_animals", q3(k$n_animals), "\n")

# STAGE 6 — strong subset, recomputed from mashr's input rds with dglm_mashr.R line 118:
# gene has per-condition DGLM qval < 0.05 in >= 1 kept condition
ff <- list.files(F, pattern="_dglm_results.rds$", full.names=TRUE)
hits <- c(); nq <- 0; nc <- 0
for (f in ff) { a <- readRDS(f)$array; q <- a[,4,,drop=FALSE]
  nc <- nc + sum(apply(q, 3, function(z) any(!is.na(z))))
  nq <- nq + sum(q < 0.05, na.rm=TRUE)
  hits <- c(hits, rownames(a)[apply(q, 1, function(z) any(z < 0.05, na.rm=TRUE))]) }
cat("\n=== 6. mashr input:", length(ff), "rds,", nc, "conditions with data ===\n")
cat("DGLM per-condition qval<0.05 rows:", nq, "  strong subset (unique genes):", length(unique(hits)), "\n")

# mashr result object: print its structure so the doc reports what is stored (read-only)
mf <- file.path(F, "combined_dglm_mashr_results_qval_single_raw_fill1000.rds")
if (file.exists(mf)) { cat("\n=== 6b. mashr result object ===\n"); str(readRDS(mf), max.level=1) }

# STAGE 7 — FDR + final master. qvalue = global BH (dglm_fdr_combined.R line 140).
# Significance = mash_lfsr < 0.05; direction = sign(mash_beta); genes keyed by ensembl_id.
# One gene can go up in one condition and down in another, so inc + dec > unique sig.
mm <- read.delim(file.path(F, "master_dglm_combined.tsv")); s <- mm[mm$mash_lfsr < 0.05,]
cat("\n=== 7. master table ===\n")
cat("rows:", nrow(mm), " conditions:", length(unique(paste(mm$cell_type, mm$region))),
    " unique genes:", length(unique(mm$ensembl_id)), "\n")
cat("global-BH qvalue<0.05 rows:", sum(mm$qvalue < 0.05, na.rm=TRUE), "\n")
cat("sig rows (lfsr<0.05):", nrow(s), " (", round(100*nrow(s)/nrow(mm), 2), "% )\n")
cat("unique sig genes:", length(unique(s$ensembl_id)),
    " inc:", length(unique(s$ensembl_id[s$mash_beta > 0])),
    " dec:", length(unique(s$ensembl_id[s$mash_beta < 0])), "\n")
cat("conditions with >=1 sig gene:", length(unique(paste(s$cell_type, s$region))), "\n")
REOF
