#!/usr/bin/env Rscript
suppressMessages({library(optparse); library(ggplot2); library(patchwork)})

option_list = list(
  make_option('--pipeline', type='character', default='mc10'),
  make_option('--lfsr_thresh', type='double', default=0.05),
  make_option('--out', type='character', default='upset_combined.png')
)
opt = parse_args(OptionParser(option_list=option_list))

region_colors = c("ACC"="#5555b0","CN"="#2caa4d","dlPFC"="#345296","EC"="#532b70",
                  "HIP"="#84db55","IPP"="#6c9cc0","lCb"="#e6cd69","M1"="#467aa1",
                  "MB"="#fc7462","mdTN"="#117867","NAc"="#148e34")
ct_colors = c("GABAergic_neurons"="#E07B39","glutamatergic_neurons"="#E8A857",
              "cerebellar_neurons"="#F2C97E","medium_spiny_neurons"="#C4A862",
              "midbrain_neurons"="#A07840","basket_cells"="#8B5E3C",
              "microglia"="#3B7DD8","astrocytes"="#5BA4CF","opc"="#2E5FA3",
              "oligodendrocytes"="#1A3A6B","vascular_cells"="#7FBEDB",
              "ependymal_cells"="#B8D9E8")

lcb_affected = c("opc","microglia","astrocytes","vascular_cells","oligodendrocytes")
cell_types = c("GABAergic_neurons","glutamatergic_neurons","cerebellar_neurons",
               "medium_spiny_neurons","midbrain_neurons","basket_cells",
               "microglia","astrocytes","opc","oligodendrocytes","vascular_cells","ependymal_cells")

mc = sub("mc","",opt$pipeline)
all_data = list()

for (ct in cell_types) {
  if (opt$pipeline == "baseline") {
    base = if (ct %in% lcb_affected) "disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered" else "disp_age__FINAL_autosomeX_300c_QCfiltered"
  } else {
    base = if (ct %in% lcb_affected) paste0("disp_age__mincells",mc,"_300c_noLcb_QCfiltered") else paste0("disp_age__mincells",mc,"_300c_QCfiltered")
  }
  f = file.path("/scratch/easmit31/dispersion/dglm", base, ct, "dglm_checkpoints_cutoff0.5", "master_dglm_combined.tsv")
  if (!file.exists(f)) { cat("SKIP:", ct, "\n"); next }
  d = read.csv(f, sep="\t")
  d$cell_type_broad = ct
  all_data[[ct]] = d
}
df = do.call(rbind, all_data)
sig = df[df$mash_lfsr < opt$lfsr_thresh, ]

# extract louvain number
sig$louvain_num = as.integer(sub(".*_([0-9]+)$","\\1", sig$cell_type))
sig$row_label = paste0(sig$louvain_num, "_", sig$region)

# build combos sorted: by cell_type order, then louvain_num, then region
combos = unique(sig[, c("cell_type_broad","cell_type","louvain_num","region","row_label")])
combos = combos[order(match(combos$cell_type_broad, cell_types), combos$louvain_num, combos$region), ]
combos$row_label = factor(combos$row_label, levels=rev(unique(combos$row_label)))
combos$ct_color = ct_colors[combos$cell_type_broad]

# gene counts
combos$n_inc = sapply(seq_len(nrow(combos)), function(i) {
  length(unique(sig$ensembl_id[sig$cell_type==combos$cell_type[i] & sig$region==combos$region[i] & sig$mash_beta>0]))
})
combos$n_dec = sapply(seq_len(nrow(combos)), function(i) {
  length(unique(sig$ensembl_id[sig$cell_type==combos$cell_type[i] & sig$region==combos$region[i] & sig$mash_beta<=0]))
})

all_regions = c("ACC","CN","dlPFC","EC","HIP","IPP","lCb","M1","MB","mdTN","NAc")

# increasing bar
p_inc = ggplot(combos, aes(y=row_label, x=n_inc, fill=cell_type_broad)) +
  geom_col() +
  scale_fill_manual(values=ct_colors, guide="none") +
  theme_minimal(base_size=5) +
  theme(axis.text.y=element_text(size=3.5), panel.grid.major.y=element_blank()) +
  scale_x_reverse() +
  labs(x="Increasing genes", y=NULL)

# decreasing bar
p_dec = ggplot(combos, aes(y=row_label, x=n_dec, fill=cell_type_broad)) +
  geom_col() +
  scale_fill_manual(values=ct_colors, guide="none") +
  theme_minimal(base_size=5) +
  theme(axis.text.y=element_blank(), panel.grid.major.y=element_blank()) +
  labs(x="Decreasing genes", y=NULL)

# dot matrix: two sets of dots — increasing (red) and decreasing (blue)
dot_rows = list()
for (i in seq_len(nrow(combos))) {
  for (r in all_regions) {
    has_inc = any(sig$cell_type==combos$cell_type[i] & sig$region==r & sig$mash_beta>0)
    has_dec = any(sig$cell_type==combos$cell_type[i] & sig$region==r & sig$mash_beta<=0)
    if (r == combos$region[i]) {
      dot_rows[[paste(i,r,"inc")]] = data.frame(row_label=combos$row_label[i], region=r, direction="Increasing", present=has_inc)
      dot_rows[[paste(i,r,"dec")]] = data.frame(row_label=combos$row_label[i], region=r, direction="Decreasing", present=has_dec)
    } else {
      dot_rows[[paste(i,r,"inc")]] = data.frame(row_label=combos$row_label[i], region=r, direction="Increasing", present=FALSE)
      dot_rows[[paste(i,r,"dec")]] = data.frame(row_label=combos$row_label[i], region=r, direction="Decreasing", present=FALSE)
    }
  }
}
dot_df = do.call(rbind, dot_rows)
dot_df$row_label = factor(dot_df$row_label, levels=levels(combos$row_label))
dot_df$region = factor(dot_df$region, levels=all_regions)
dot_df$direction = factor(dot_df$direction, levels=c("Increasing","Decreasing"))
dot_df$x_pos = as.numeric(dot_df$region) + ifelse(dot_df$direction=="Increasing", -0.2, 0.2)

p_mid = ggplot(dot_df[dot_df$present,], aes(x=x_pos, y=row_label, color=direction)) +
  geom_point(size=1.0) +
  geom_point(data=dot_df[!dot_df$present,], aes(x=x_pos, y=row_label), color="grey90", size=0.5, alpha=0.3) +
  scale_color_manual(values=c("Increasing"="firebrick","Decreasing"="steelblue"), name=NULL) +
  scale_x_continuous(breaks=seq_along(all_regions), labels=all_regions) +
  theme_minimal(base_size=6) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=6),
        axis.text.y=element_blank(),
        panel.grid.major=element_line(color="grey92"),
        legend.position="top") +
  labs(x="Region", y=NULL)

combined = p_inc + p_mid + p_dec +
  plot_layout(widths=c(2, 1.5, 2)) +
  plot_annotation(title=paste0("All cell types (", opt$pipeline, "): unique significant genes per subcluster x region"),
                   subtitle=paste0("Rows sorted by cell type (neuronal then non-neuronal), louvain#, region. lfsr<", opt$lfsr_thresh, ". Two dots per region: red=increasing, blue=decreasing."))

h = max(12, nrow(combos) * 0.11)
ggsave(opt$out, combined, width=16, height=h, dpi=150)
cat("Saved:", opt$out, "| rows:", nrow(combos), "\n")
