suppressMessages({
  library(dglm)
  library(edgeR)
  library(parallel)
})

RAW = "/scratch/easmit31/dispersion/dglm/disp_age__raw_counts"
autosome_x = read.csv("/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv", stringsAsFactors=FALSE)$ensembl_gene_id

fit_new = function(gene_idx, expr_mat, meta) {
  d = meta
  d$e = as.numeric(expr_mat[gene_idx, ])
  f = try(dglm(e ~ age+sex+mean_n_umi+n_cells, ~age, family=gaussian(), dlink="log", data=d), silent=TRUE)
  if (inherits(f, "try-error")) return(NA)
  tb = try(coef(summary(f$dispersion.fit)), silent=TRUE)
  if (inherits(tb,"try-error") || !("age" %in% rownames(tb))) return(NA)
  tb["age","Pr(>|t|)"]
}

for (ct in c("cerebellar_neurons", "basket_cells")) {
  ct_dir = file.path(RAW, ct)
  pb_files = list.files(ct_dir, pattern="_pseudobulk.csv$", full.names=TRUE)
  chosen = NULL
  for (f in pb_files) {
    hdr = read.csv(f, nrows=0, check.names=FALSE)
    if (ncol(hdr) >= 20) { chosen = f; break }
  }
  if (is.null(chosen)) chosen = pb_files[1]

  base = sub("_pseudobulk.csv$", "", basename(chosen))
  parts = strsplit(base, "_")[[1]]
  region = parts[length(parts)]
  subcluster = sub(paste0("_", region, "$"), "", base)

  meta_file = file.path(ct_dir, paste0(subcluster, "_metadata.csv"))
  counts = read.csv(chosen, row.names=1, check.names=FALSE)
  counts = counts[intersect(rownames(counts), autosome_x), , drop=FALSE]
  counts = as.matrix(counts); storage.mode(counts) = "numeric"

  meta = read.csv(meta_file, stringsAsFactors=FALSE)
  meta = meta[meta$region == region, ]
  meta = meta[match(colnames(counts), meta$animal_id), ]
  meta$sex = as.factor(meta$sex); meta$age = as.numeric(meta$age)
  meta$mean_n_umi = as.numeric(meta$mean_n_umi); meta$n_cells = as.numeric(meta$n_cells)

  keep = rowSums(cpm(counts) >= 0.5) >= (0.5 * ncol(counts))
  counts = counts[keep, , drop=FALSE]

  y = DGEList(counts=counts); y = calcNormFactors(y, method="TMM")
  logCPM = cpm(y, log=TRUE, prior.count=0.5)

  n = nrow(logCPM); n_test = min(n, 3000); idx = sample(1:n, n_test)

  real_p = unlist(mclapply(idx, fit_new, expr_mat=logCPM, meta=meta, mc.cores=16))
  set.seed(1)
  meta_perm = meta; meta_perm$age = sample(meta$age)
  perm_p = unlist(mclapply(idx, fit_new, expr_mat=logCPM, meta=meta_perm, mc.cores=16))

  cat(ct, "(", subcluster, region, "| n_animals=", ncol(counts), "):\n")
  cat("  Real age  % p<0.05:", round(100*mean(real_p<0.05, na.rm=TRUE),1), "\n")
  cat("  Permuted age % p<0.05:", round(100*mean(perm_p<0.05, na.rm=TRUE),1), "\n\n")
}
