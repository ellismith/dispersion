#!/usr/bin/env Rscript
# Compares original (no per-animal cell filter) vs min_cells_per_animal=10
# vs =100, for every cell type that has all three finished.

y_genes = c("SRY","DDX3Y","UTY","USP9Y","ZFY","KDM5D","EIF1AY")

cell_types = c("GABAergic_neurons","microglia","astrocytes","glutamatergic_neurons",
               "cerebellar_neurons","medium_spiny_neurons","midbrain_neurons",
               "basket_cells","vascular_cells","ependymal_cells","opc","oligodendrocytes")
lcb_affected = c("opc","microglia","astrocytes","vascular_cells","oligodendrocytes")

get_orig_dir = function(ct) if (ct %in% lcb_affected) "disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered" else "disp_age__FINAL_autosomeX_300c_QCfiltered"

cat(sprintf("%-25s %12s %12s %12s\n", "cell_type", "orig_%_sig", "mc10_%_sig", "mc100_%_sig"))

for (ct in cell_types) {
  f_orig = file.path("/scratch/easmit31/dispersion/dglm", get_orig_dir(ct), ct, "dglm_checkpoints_cutoff0.5", "master_dglm_combined.tsv")
  f_mc10 = file.path("/scratch/easmit31/dispersion/dglm/disp_age__mincells10_300c_QCfiltered", ct, "dglm_checkpoints_cutoff0.5", "master_dglm_combined.tsv")
  f_mc100 = file.path("/scratch/easmit31/dispersion/dglm/disp_age__mincells100_300c_QCfiltered", ct, "dglm_checkpoints_cutoff0.5", "master_dglm_combined.tsv")

  pct = function(f) {
    if (!file.exists(f)) return(NA)
    df = read.csv(f, sep="\t")
    round(100*sum(df$mash_lfsr<0.2,na.rm=TRUE)/nrow(df), 1)
  }

  cat(sprintf("%-25s %11s%% %11s%% %11s%%\n", ct,
              ifelse(is.na(pct(f_orig)),"pend",pct(f_orig)),
              ifelse(is.na(pct(f_mc10)),"pend",pct(f_mc10)),
              ifelse(is.na(pct(f_mc100)),"pend",pct(f_mc100))))
}
