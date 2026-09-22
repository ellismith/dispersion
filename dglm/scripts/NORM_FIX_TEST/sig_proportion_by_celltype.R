suppressMessages({
  library(dglm)
  library(edgeR)
  library(parallel)
})

cell_types = c("GABAergic_neurons","microglia","astrocytes","glutamatergic_neurons",
               "cerebellar_neurons","medium_spiny_neurons","midbrain_neurons",
               "basket_cells","vascular_cells","ependymal_cells","opc","oligodendrocytes")

RAW = "/scratch/easmit31/dispersion/dglm/disp_age__raw_counts"
autosome_x = read.csv("/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv", stringsAsFactors=FALSE)$ensembl_gene_id

fit_gene = function(gene_idx, expr_mat, meta, use_new) {
  d = meta
  d$e = as.numeric(expr_mat[gene_idx, ])
  f = try(dglm(e ~ age+sex+mean_n_umi+n_cells, ~age, family=gaussian(), dlink="log", data=d), silent=TRUE)
  if (inherits(f, "try-error")) return(NA)
  if (use_new) {
    tb = try(coef(summary(f$dispersion.fit)), silent=TRUE)
    if (inherits(tb,"try-error") || !("age" %in% rownames(tb))) return(NA)
    return(tb["age","Pr(>|t|)"])
  } else {
    tb = try(coef(summary(f)$dispersion.summary), silent=TRUE)
    if (inherits(tb,"try-error") || !("age" %in% rownames(tb))) return(NA)
    return(tb["age",4])
  }
}

results = list()

for (ct in cell_types) {
  ct_dir = file.path(RAW, ct)
  pb_files = list.files(ct_dir, pattern="_pseudobulk.csv$", full.names=TRUE)
  if (length(pb_files) == 0) next

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
  if (!file.exists(meta_file)) next

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
  if (nrow(counts) < 50 || ncol(counts) < 10) next

  y = DGEList(counts=counts); y = calcNormFactors(y, method="TMM")
  logCPM = cpm(y, log=TRUE, prior.count=0.5)

  n = nrow(logCPM)
  n_test = min(n, 3000)
  idx = sample(1:n, n_test)

  old_real = unlist(mclapply(idx, fit_gene, expr_mat=logCPM, meta=meta, use_new=FALSE, mc.cores=16))
  new_real = unlist(mclapply(idx, fit_gene, expr_mat=logCPM, meta=meta, use_new=TRUE, mc.cores=16))

  results[[ct]] = data.frame(
    cell_type = ct, subcluster = subcluster, region = region, n_genes_tested = n_test, n_animals = ncol(counts),
    old_pct_sig = round(100*mean(old_real<0.05, na.rm=TRUE), 1),
    new_pct_sig = round(100*mean(new_real<0.05, na.rm=TRUE), 1)
  )
  cat(ct, "(", subcluster, region, "):", "old=", results[[ct]]$old_pct_sig, "% new=", results[[ct]]$new_pct_sig, "%\n")
}

final = do.call(rbind, results)
write.csv(final, "/scratch/easmit31/dispersion/dglm/sig_proportion_by_celltype.csv", row.names=FALSE)
cat("\nSaved.\n")
