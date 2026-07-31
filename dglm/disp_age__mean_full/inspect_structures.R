# Quick structure dump - one cell type's dglm results + the combined mashr object
# Run: Rscript inspect_structures.R > inspect_output.txt

cat("=== Single cell-type dglm results (astrocytes) ===\n")
dglm_res <- readRDS("astrocytes_dglm_results.rds")
cat("class:", class(dglm_res), "\n")
if (is.list(dglm_res) && !is.data.frame(dglm_res)) {
  cat("names:", paste(names(dglm_res), collapse=", "), "\n")
  str(dglm_res, max.level = 2, list.len = 20)
} else {
  cat("dim:", paste(dim(dglm_res), collapse=" x "), "\n")
  cat("colnames:", paste(colnames(dglm_res), collapse=", "), "\n")
  print(head(dglm_res, 3))
}

cat("\n\n=== Combined mashr object ===\n")
mash_obj <- readRDS("combined_mashr_disp_age_mean_full.rds")
cat("class:", class(mash_obj), "\n")
cat("names:", paste(names(mash_obj), collapse=", "), "\n")
str(mash_obj, max.level = 2, list.len = 20)

if ("result" %in% names(mash_obj)) {
  cat("\nresult names:", paste(names(mash_obj$result), collapse=", "), "\n")
  cat("PosteriorMean dim:", paste(dim(mash_obj$result$PosteriorMean), collapse=" x "), "\n")
  cat("PosteriorMean rownames head:", paste(head(rownames(mash_obj$result$PosteriorMean)), collapse=", "), "\n")
  cat("PosteriorMean colnames:", paste(colnames(mash_obj$result$PosteriorMean), collapse=", "), "\n")
}
