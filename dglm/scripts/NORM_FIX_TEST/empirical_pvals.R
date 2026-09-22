library(dglm)
library(edgeR)
library(parallel)

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

n_genes = nrow(logCPM)
n_perm = 100
n_cores = 16

fit_beta_one = function(gene_idx, age_vec) {
  d = meta
  d$age = age_vec
  d$e = as.numeric(logCPM[gene_idx,])
  f = try(dglm(e ~ age+sex+mean_n_umi+n_cells, ~age, family=gaussian(), dlink="log", data=d), silent=TRUE)
  if (inherits(f,"try-error")) return(NA)
  tb = coef(summary(f)$dispersion.summary)
  unname(tb["age",1])
}

cat("=== REAL age: fitting all", n_genes, "genes ===\n")
real_beta = unlist(mclapply(1:n_genes, fit_beta_one, age_vec=meta$age, mc.cores=n_cores))
names(real_beta) = rownames(logCPM)
cat("Done. NAs:", sum(is.na(real_beta)), "\n")

cat("\n=== Running", n_perm, "permutations of age ===\n")
perm_beta = matrix(NA, nrow=n_genes, ncol=n_perm)
for (p in 1:n_perm) {
  set.seed(p)
  age_perm = sample(meta$age)
  perm_beta[,p] = unlist(mclapply(1:n_genes, fit_beta_one, age_vec=age_perm, mc.cores=n_cores))
  cat("  permutation", p, "of", n_perm, "done\n")
}

cat("\n=== Computing empirical p-values ===\n")
emp_p = sapply(1:n_genes, function(i) {
  if (is.na(real_beta[i])) return(NA)
  perm_i = perm_beta[i,]
  perm_i = perm_i[!is.na(perm_i)]
  if (length(perm_i)==0) return(NA)
  (1 + sum(abs(perm_i) >= abs(real_beta[i]))) / (1 + length(perm_i))
})

emp_q = p.adjust(emp_p, method="BH")

results = data.frame(
  gene = rownames(logCPM),
  beta = real_beta,
  empirical_p = emp_p,
  empirical_q = emp_q
)
write.csv(results, "/scratch/easmit31/dispersion/dglm/scripts/NORM_FIX_TEST/empirical_pval_results.csv", row.names=FALSE)

cat("\n=== SUMMARY ===\n")
cat("Genes tested:", n_genes, "\n")
cat("Significant at empirical_q < 0.05:", sum(emp_q < 0.05, na.rm=TRUE),
    "(", round(100*mean(emp_q<0.05, na.rm=TRUE),2), "%)\n")
cat("Significant at empirical_q < 0.2:", sum(emp_q < 0.2, na.rm=TRUE),
    "(", round(100*mean(emp_q<0.2, na.rm=TRUE),2), "%)\n")
cat("\nSaved: empirical_pval_results.csv\n")
