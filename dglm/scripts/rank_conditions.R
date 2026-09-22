#!/usr/bin/env Rscript
# Rank all subtype x region conditions by number of unique significant genes,
# split by direction (increasing vs decreasing dispersion with age).
# Usage: Rscript rank_conditions.R --pipeline baseline|mc10 --lfsr 0.05 --out <csv>

suppressMessages(library(optparse))
option_list = list(
  make_option('--pipeline', type='character', default='baseline'),
  make_option('--lfsr', type='double', default=0.05),
  make_option('--out', type='character', default='ranked_conditions.csv')
)
opt = parse_args(OptionParser(option_list=option_list))

lcb = c("opc","microglia","astrocytes","vascular_cells","oligodendrocytes")
cell_types = c("GABAergic_neurons","microglia","astrocytes","glutamatergic_neurons",
               "cerebellar_neurons","medium_spiny_neurons","midbrain_neurons",
               "basket_cells","vascular_cells","ependymal_cells","opc","oligodendrocytes")

all_rows = list()
for (ct in cell_types) {
  if (opt$pipeline == "baseline") {
    base = if(ct %in% lcb) "disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered" else "disp_age__FINAL_autosomeX_300c_QCfiltered"
  } else {
    base = if(ct %in% lcb) "disp_age__mincells10_300c_noLcb_QCfiltered" else "disp_age__mincells10_300c_QCfiltered"
  }
  f = file.path("/scratch/easmit31/dispersion/dglm", base, ct, "dglm_checkpoints_cutoff0.5", "master_dglm_combined.tsv")
  if (!file.exists(f)) { cat("SKIP:", ct, "\n"); next }
  df = read.csv(f, sep="\t")
  sig = df[df$mash_lfsr < opt$lfsr, ]

  conditions = unique(sig[, c("cell_type","region")])
  for (i in seq_len(nrow(conditions))) {
    subcl = conditions$cell_type[i]; reg = conditions$region[i]
    rows = sig[sig$cell_type==subcl & sig$region==reg, ]
    n_inc = length(unique(rows$ensembl_id[rows$mash_beta > 0]))
    n_dec = length(unique(rows$ensembl_id[rows$mash_beta <= 0]))
    n_total = length(unique(rows$ensembl_id))
    median_beta_inc = if(n_inc>0) round(median(rows$mash_beta[rows$mash_beta>0], na.rm=TRUE), 4) else NA
    median_beta_dec = if(n_dec>0) round(median(rows$mash_beta[rows$mash_beta<=0], na.rm=TRUE), 4) else NA
    all_rows[[paste(subcl,reg)]] = data.frame(
      cell_type_broad=ct, subcluster=subcl, region=reg,
      n_sig_total=n_total, n_increasing=n_inc, n_decreasing=n_dec,
      pct_increasing=round(100*n_inc/n_total,1),
      median_beta_increasing=median_beta_inc,
      median_beta_decreasing=median_beta_dec
    )
  }
}

results = do.call(rbind, all_rows)

cat("=== TOP 20 CONDITIONS: MOST INCREASING GENES ===\n")
top_inc = results[order(-results$n_increasing), ][1:min(20,nrow(results)), ]
print(top_inc[, c("cell_type_broad","subcluster","region","n_increasing","median_beta_increasing")], row.names=FALSE)

cat("\n=== TOP 20 CONDITIONS: MOST DECREASING GENES ===\n")
top_dec = results[order(-results$n_decreasing), ][1:min(20,nrow(results)), ]
print(top_dec[, c("cell_type_broad","subcluster","region","n_decreasing","median_beta_decreasing")], row.names=FALSE)

write.csv(results[order(-results$n_sig_total), ], opt$out, row.names=FALSE)
cat("\nFull ranked table saved to:", opt$out, "\n")
cat("Total conditions:", nrow(results), "\n")
