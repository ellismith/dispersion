#!/usr/bin/env Rscript
y_genes = c("SRY","DDX3Y","UTY","USP9Y","ZFY","KDM5D","EIF1AY")
cell_types = c("GABAergic_neurons","microglia","astrocytes","glutamatergic_neurons",
               "cerebellar_neurons","medium_spiny_neurons","midbrain_neurons",
               "basket_cells","vascular_cells","ependymal_cells","opc","oligodendrocytes")

for (mc in c(10, 100)) {
  base = paste0("/scratch/easmit31/dispersion/dglm/disp_age__mincells", mc, "_300c_QCfiltered")
  cat("\n########## min_cells_per_animal =", mc, "##########\n")
  for (ct in cell_types) {
    f = file.path(base, ct, "dglm_checkpoints_cutoff0.5", "master_dglm_combined.tsv")
    if (!file.exists(f)) { cat(sprintf("%-25s NOT YET DONE\n", ct)); next }
    df = read.csv(f, sep="\t")
    dups = sum(duplicated(df[,c("cell_type","region","ensembl_id")]))
    na_beta = sum(is.na(df$mash_beta))
    lfsr_ok = all(df$mash_lfsr >= 0 & df$mash_lfsr <= 1, na.rm=TRUE)
    beta_sane = max(abs(df$mash_beta), na.rm=TRUE) < 10
    mt = sum(grepl("^MT-", df$symbol, ignore.case=TRUE))
    yg = length(intersect(y_genes, unique(df$symbol)))
    pct_sig = round(100*sum(df$mash_lfsr<0.2,na.rm=TRUE)/nrow(df),1)
    status = if (dups==0 && na_beta==0 && lfsr_ok && beta_sane && mt==0 && yg==0) "PASS" else "FAIL"
    cat(sprintf("%-25s %-5s rows=%-9d sig=%5.1f%% dups=%d na=%d MT=%d Y=%d\n",
                ct, status, nrow(df), pct_sig, dups, na_beta, mt, yg))
  }
}
