#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
  library(dglm)
})

e.this = read.csv("/scratch/easmit31/dispersion/dglm/disp_age__raw_counts/GABAergic_neurons/filtered/GABAergic_neurons_0_ACC_filtered_cutoff0.5.csv", row.names=1, check.names=FALSE)
meta = read.csv("/scratch/easmit31/dispersion/dglm/disp_age__raw_counts/GABAergic_neurons/GABAergic_neurons_0_metadata.csv", stringsAsFactors=FALSE)
meta = meta[meta$region=="ACC",]
meta = meta[match(colnames(e.this), meta$animal_id),]
rownames(meta) = meta$animal_id
meta$sex = as.factor(meta$sex)
meta$age = as.numeric(meta$age)
meta$mean_n_umi = as.numeric(meta$mean_n_umi)
meta$n_cells = as.numeric(meta$n_cells)

y = DGEList(counts = as.matrix(e.this))
y = calcNormFactors(y, method="TMM")
logCPM = cpm(y, log=TRUE, prior.count=0.5)

fit_one_gene = function(expr, dat) {
  dat$e = as.numeric(expr)
  out = tryCatch(
    suppressWarnings(
      dglm(e ~ age + sex + mean_n_umi + n_cells, ~ age,
           family=gaussian(), dlink="log", data=dat)
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(c(beta=NA, se=NA, z=NA, converged=FALSE))

  beta = unname(out$dispersion.fit$coefficients["age"])
  se = tryCatch(
    unname(summary(out$dispersion.fit)$coefficients["age","Std. Error"]),
    error = function(e) NA_real_
  )
  converged = isTRUE(out$fit$converged) && isTRUE(out$dispersion.fit$converged)
  z = beta/se
  c(beta=beta, se=se, z=z, converged=converged)
}

run_all_genes = function(age_vector) {
  dat = data.frame(age=as.numeric(age_vector), sex=meta$sex,
                    mean_n_umi=as.numeric(meta$mean_n_umi), n_cells=as.numeric(meta$n_cells))
  ans = t(vapply(seq_len(nrow(logCPM)), function(i) fit_one_gene(logCPM[i,], dat), numeric(4)))
  ans = as.data.frame(ans)
  ans$gene = rownames(logCPM)
  ans$converged = as.logical(ans$converged)
  ans
}

set.seed(1)
real = run_all_genes(meta$age)
perm = run_all_genes(sample(meta$age))

summarize_run = function(x, label) {
  valid = x[x$converged & is.finite(x$beta) & is.finite(x$se) & x$se>0 & is.finite(x$z), , drop=FALSE]
  cat("\n", label, "\n", sep="")
  cat("All genes: ", nrow(x), "\n", sep="")
  cat("Converged valid genes: ", nrow(valid), "\n", sep="")
  cat("Percent |z| > 1.96 among converged valid genes: ",
      round(100*mean(abs(valid$z) > 1.96), 2), "%\n", sep="")
  invisible(valid)
}

summarize_run(real, "REAL AGE")
summarize_run(perm, "PERMUTED AGE")
