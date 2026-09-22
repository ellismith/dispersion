#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(patchwork)})

ct_colors = c(
  "GABAergic_neurons"="#E07B39","glutamatergic_neurons"="#E8A857",
  "cerebellar_neurons"="#F2C97E","medium_spiny_neurons"="#C4A862",
  "midbrain_neurons"="#A07840","basket_cells"="#8B5E3C",
  "microglia"="#3B7DD8","astrocytes"="#5BA4CF",
  "opc"="#2E5FA3","oligodendrocytes"="#1A3A6B",
  "vascular_cells"="#7FBEDB","ependymal_cells"="#B8D9E8"
)
region_colors = c(
  "ACC"="#5555b0","CN"="#2caa4d","dlPFC"="#345296","EC"="#532b70",
  "HIP"="#84db55","IPP"="#6c9cc0","lCb"="#e6cd69","M1"="#467aa1",
  "MB"="#fc7462","mdTN"="#117867","NAc"="#148e34"
)

lcb = c("opc","microglia","astrocytes","vascular_cells","oligodendrocytes")
cell_types = c("GABAergic_neurons","glutamatergic_neurons","cerebellar_neurons",
               "medium_spiny_neurons","midbrain_neurons","basket_cells",
               "microglia","astrocytes","opc","oligodendrocytes","vascular_cells","ependymal_cells")

rows = list()
for (ct in cell_types) {
  base = if(ct %in% lcb) "disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered" else "disp_age__FINAL_autosomeX_300c_QCfiltered"
  f = paste0("/scratch/easmit31/dispersion/dglm/",base,"/",ct,"/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv")
  if(!file.exists(f)) next
  df = read.csv(f,sep="\t")
  sig = df[df$mash_lfsr<0.05,]
  for (cond in unique(paste(sig$cell_type, sig$region))) {
    parts = strsplit(cond," ")[[1]]
    subcl = parts[1]; reg = parts[2]
    s = sig[sig$cell_type==subcl & sig$region==reg,]
    n_inc = length(unique(s$ensembl_id[s$mash_beta>0]))
    n_dec = length(unique(s$ensembl_id[s$mash_beta<=0]))
    total = n_inc + n_dec
    if(total == 0) next
    rows[[cond]] = data.frame(
      cell_type=ct, region=reg,
      prop_inc=n_inc/total,
      prop_dec=n_dec/total,
      n_sig=total
    )
  }
}
df_plot = do.call(rbind, rows)

# cell type order: neuronal then non-neuronal
neuronal = c("GABAergic_neurons","glutamatergic_neurons","cerebellar_neurons",
             "medium_spiny_neurons","midbrain_neurons","basket_cells")
ct_order = c(neuronal, setdiff(cell_types, neuronal))
df_plot$cell_type = factor(df_plot$cell_type, levels=ct_order)
df_plot$region = factor(df_plot$region, levels=names(region_colors))

# short CT labels
ct_labels = c(
  "GABAergic_neurons"="GABA","glutamatergic_neurons"="Glut",
  "cerebellar_neurons"="Cereb","medium_spiny_neurons"="MSN",
  "midbrain_neurons"="Midbrain","basket_cells"="Basket",
  "microglia"="Microglia","astrocytes"="Astro",
  "opc"="OPC","oligodendrocytes"="Oligo",
  "vascular_cells"="Vasc","ependymal_cells"="Epend"
)

p_ct = ggplot(df_plot, aes(x=cell_type, y=prop_inc, fill=cell_type)) +
  geom_violin(scale="width", alpha=0.8, color=NA) +
  geom_boxplot(width=0.15, fill="white", outlier.size=0.5, color="grey30") +
  geom_hline(yintercept=0.5, linetype="dashed", color="grey50") +
  scale_fill_manual(values=ct_colors, guide="none") +
  scale_x_discrete(labels=ct_labels) +
  scale_y_continuous(limits=c(0,1), labels=scales::percent) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        panel.grid.major.x=element_blank()) +
  labs(x=NULL, y="Proportion of significant genes\nthat increase dispersion with age",
       title="By cell type")

p_reg = ggplot(df_plot, aes(x=region, y=prop_inc, fill=region)) +
  geom_violin(scale="width", alpha=0.8, color=NA) +
  geom_boxplot(width=0.15, fill="white", outlier.size=0.5, color="grey30") +
  geom_hline(yintercept=0.5, linetype="dashed", color="grey50") +
  scale_fill_manual(values=region_colors, guide="none") +
  scale_y_continuous(limits=c(0,1), labels=scales::percent) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        panel.grid.major.x=element_blank()) +
  labs(x=NULL, y=NULL, title="By region")

combined = p_ct | p_reg
combined = combined + plot_annotation(
  title="Proportion of significant genes showing increasing vs decreasing dispersion with age",
  subtitle="Each point = one subcluster x region condition. Dashed line = 50% (equal inc/dec). Baseline pipeline, lfsr<0.05.")

ggsave("/scratch/easmit31/dispersion/dglm/figures/prop_increasing_violin.png",
       combined, width=16, height=7, dpi=150)
cat("Saved.\n")
