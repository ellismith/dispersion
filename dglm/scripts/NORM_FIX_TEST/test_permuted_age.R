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

fit_all_z = function(meta_use, expr_mat) {
  sapply(1:nrow(expr_mat), function(i){
    d = meta_use
    d$e = as.numeric(expr_mat[i,])
    f = try(dglm(e ~ age+sex+mean_n_umi+n_cells, ~age, family=gaussian(), dlink="log", data=d), silent=TRUE)
    if(inherits(f,"try-error")) return(NA)
    tb = coef(summary(f)$dispersion.summary)
    tb["age",1]/tb["age",2]
  })
}

cat("=== REAL age, TMM+logCPM ===\n")
z_real = fit_all_z(meta, logCPM)
cat("pct |z|>1.96 =", round(100*mean(abs(z_real)>1.96, na.rm=TRUE),1), "%\n")
w_real = warnings()
cat("\n--- warnings (real) ---\n")
print(w_real)

cat("\n=== PERMUTED age (set.seed(1)), TMM+logCPM ===\n")
meta_perm = meta
set.seed(1)
meta_perm$age = sample(meta_perm$age)
z_perm = fit_all_z(meta_perm, logCPM)
cat("pct |z|>1.96 =", round(100*mean(abs(z_perm)>1.96, na.rm=TRUE),1), "%\n")
w_perm = warnings()
cat("\n--- warnings (permuted) ---\n")
print(w_perm)

cat("\n=== SUMMARY ===\n")
cat("Real age:      ", round(100*mean(abs(z_real)>1.96, na.rm=TRUE),1), "%\n")
cat("Permuted age:  ", round(100*mean(abs(z_perm)>1.96, na.rm=TRUE),1), "%\n")
cat("(target under null: ~5%)\n")
