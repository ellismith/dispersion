#!/usr/bin/env Rscript
# Pre-mashr gate: for every condition that survived DGLM, verify that
# NO animal in the fit has <100 cells, and report how many raw-zero
# values remain among kept animals. Read-only.
args = commandArgs(trailingOnly=TRUE)
staged = args[1]   # dglm_run dir with subcluster folders
worst = data.frame()
for (d in list.dirs(staged, recursive=FALSE)) {
  sb = basename(d)
  meta = read.csv(file.path(d, paste0(sb, "_metadata.csv")))
  for (f in list.files(d, "_pseudobulk.csv$", full.names=TRUE)) {
    reg = sub(paste0("^",sb,"_"), "", sub("_pseudobulk.csv$","",basename(f)))
    m = meta[meta$region==reg,]
    kept = m[!is.na(m$n_cells) & m$n_cells >= 100,]
    if (nrow(kept) < 5) next   # condition dies anyway
    cnt = read.csv(f, row.names=1, check.names=FALSE)
    cnt = as.matrix(cnt[, kept$animal_id, drop=FALSE])
    zero_frac = rowMeans(cnt == 0)
    worst = rbind(worst, data.frame(cond=paste(sb,reg), n_animals=nrow(kept),
      min_cells=min(kept$n_cells), genes_gt20pct_zero=sum(zero_frac>0.2),
      genes_total=nrow(cnt), worst_zero_frac=max(zero_frac)))
  }
}
worst = worst[order(-worst$worst_zero_frac),]
cat("Conditions surviving min100 + >=5 animals:", nrow(worst), "\n")
cat("Conditions where ANY kept animal has <100 cells (must be 0):", sum(worst$min_cells<100), "\n\n")
cat("Worst 15 by residual zero fraction (zeros that remain even after min100):\n")
print(head(worst,15), row.names=FALSE)
