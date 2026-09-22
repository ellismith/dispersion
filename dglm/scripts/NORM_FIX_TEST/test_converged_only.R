library(dglm)
library(edgeR)

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

# Fit one gene, counting warnings, returning z and a "clean" flag
fit_one = function(d) {
  n_warn = 0
  z = NA
  withCallingHandlers({
    f = try(dglm(e ~ age+sex+mean_n_umi+n_cells, ~age, family=gaussian(), dlink="log", data=d), silent=TRUE)
    if (!inherits(f,"try-error")) {
      tb = coef(summary(f)$dispersion.summary)
      z = tb["age",1]/tb["age",2]
    }
  }, warning = function(w) { n_warn <<- n_warn + 1; invokeRestart("muffleWarning") })
  c(z=z, n_warn=n_warn)
}

run_combo = function(meta_use) {
  out = t(sapply(1:nrow(logCPM), function(i){
    d = meta_use; d$e = as.numeric(logCPM[i,])
    fit_one(d)
  }))
  out
}

cat("=== REAL age ===\n")
real = run_combo(meta)
clean_real = real[,"n_warn"]==0 & !is.na(real[,"z"])
cat("Total genes:", nrow(real), "| Clean (0 warnings):", sum(clean_real), "\n")
cat("ALL genes pct |z|>1.96:", round(100*mean(abs(real[,"z"])>1.96, na.rm=TRUE),1), "%\n")
cat("CLEAN-only pct |z|>1.96:", round(100*mean(abs(real[clean_real,"z"])>1.96),1), "%\n")

cat("\n=== PERMUTED age (set.seed(1)) ===\n")
meta_perm = meta
set.seed(1)
meta_perm$age = sample(meta_perm$age)
perm = run_combo(meta_perm)
clean_perm = perm[,"n_warn"]==0 & !is.na(perm[,"z"])
cat("Total genes:", nrow(perm), "| Clean (0 warnings):", sum(clean_perm), "\n")
cat("ALL genes pct |z|>1.96:", round(100*mean(abs(perm[,"z"])>1.96, na.rm=TRUE),1), "%\n")
cat("CLEAN-only pct |z|>1.96:", round(100*mean(abs(perm[clean_perm,"z"])>1.96),1), "%\n")

cat("\n=== SUMMARY (clean fits only) ===\n")
cat("Real age:      ", round(100*mean(abs(real[clean_real,"z"])>1.96),1), "%\n")
cat("Permuted age:  ", round(100*mean(abs(perm[clean_perm,"z"])>1.96),1), "%\n")
cat("(target under null: ~5%)\n")
