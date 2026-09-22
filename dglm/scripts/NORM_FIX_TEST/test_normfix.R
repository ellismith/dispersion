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

gene_idx = 50
gene_name = rownames(e.this)[gene_idx]

# --- OLD: raw counts ---
d_old = meta
d_old$e = as.numeric(e.this[gene_idx,])
fit_old = dglm(e ~ age + sex + mean_n_umi + n_cells, ~ age, family=gaussian(), dlink="log", data=d_old)
old_row = coef(summary(fit_old)$dispersion.summary)["age", c(1,2,4)]

# --- NEW: TMM + log-CPM ---
y = DGEList(counts = as.matrix(e.this))
y = calcNormFactors(y, method="TMM")
logCPM = cpm(y, log=TRUE, prior.count=0.5)

d_new = meta
d_new$e = as.numeric(logCPM[gene_idx,])
fit_new = dglm(e ~ age + sex + mean_n_umi + n_cells, ~ age, family=gaussian(), dlink="log", data=d_new)
new_row = coef(summary(fit_new)$dispersion.summary)["age", c(1,2,4)]

cat("=== Gene:", gene_name, "| GABAergic_neurons_0, ACC ===\n\n")
cat("OLD (raw counts):      beta=", old_row[1], " SE=", old_row[2], " p=", old_row[3], "\n")
cat("NEW (TMM+logCPM):      beta=", new_row[1], " SE=", new_row[2], " p=", new_row[3], "\n")

# quick calibration check across ALL genes in this combo, both ways
z_old = sapply(1:nrow(e.this), function(i){
  d=meta; d$e=as.numeric(e.this[i,])
  f=try(dglm(e ~ age+sex+mean_n_umi+n_cells, ~age, family=gaussian(), dlink="log", data=d), silent=TRUE)
  if(inherits(f,"try-error")) return(NA)
  tb=coef(summary(f)$dispersion.summary)
  tb["age",1]/tb["age",2]
})
z_new = sapply(1:nrow(logCPM), function(i){
  d=meta; d$e=as.numeric(logCPM[i,])
  f=try(dglm(e ~ age+sex+mean_n_umi+n_cells, ~age, family=gaussian(), dlink="log", data=d), silent=TRUE)
  if(inherits(f,"try-error")) return(NA)
  tb=coef(summary(f)$dispersion.summary)
  tb["age",1]/tb["age",2]
})

cat("\n=== Calibration across all", nrow(e.this), "genes in this combo (expect ~5% under null) ===\n")
cat("OLD raw counts: pct |z|>1.96 =", round(100*mean(abs(z_old)>1.96, na.rm=TRUE),1), "%\n")
cat("NEW TMM+logCPM: pct |z|>1.96 =", round(100*mean(abs(z_new)>1.96, na.rm=TRUE),1), "%\n")
