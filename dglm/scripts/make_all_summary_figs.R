suppressMessages({library(ggplot2); library(patchwork)})

lcb = c("opc","microglia","astrocytes","vascular_cells","oligodendrocytes")
cell_types = c("GABAergic_neurons","glutamatergic_neurons","cerebellar_neurons",
               "medium_spiny_neurons","midbrain_neurons","basket_cells",
               "microglia","astrocytes","opc","oligodendrocytes","vascular_cells","ependymal_cells")
neuronal = c("GABAergic_neurons","glutamatergic_neurons","cerebellar_neurons",
             "medium_spiny_neurons","midbrain_neurons","basket_cells")
all_regions = c("ACC","CN","dlPFC","EC","HIP","IPP","lCb","M1","MB","mdTN","NAc")

load_sig = function(ct) {
  base = if(ct %in% lcb) "disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered" else "disp_age__FINAL_autosomeX_300c_QCfiltered"
  f = paste0("/scratch/easmit31/dispersion/dglm/",base,"/",ct,"/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv")
  if(!file.exists(f)) return(NULL)
  df = read.csv(f,sep="\t")
  df[df$mash_lfsr<0.05,]
}

make_upset = function(ct, top_n=30) {
  sig = load_sig(ct)
  if(is.null(sig)) return(NULL)
  subclusters = unique(sig$cell_type)
  combos_inc = sapply(subclusters, function(s) paste(sort(unique(sig$region[sig$cell_type==s & sig$mash_beta>0])), collapse=","))
  combos_dec = sapply(subclusters, function(s) paste(sort(unique(sig$region[sig$cell_type==s & sig$mash_beta<=0])), collapse=","))
  get_genes = function(tab, combos, direction) {
    sapply(names(tab), function(combo) {
      regs = strsplit(combo,",")[[1]]
      subs = subclusters[combos==combo]
      length(unique(sig$ensembl_id[sig$cell_type %in% subs & (if(direction=="inc") sig$mash_beta>0 else sig$mash_beta<=0) & sig$region %in% regs]))
    })
  }
  inc_tab = table(combos_inc[combos_inc!=""])
  dec_tab = table(combos_dec[combos_dec!=""])
  inc_genes = get_genes(inc_tab, combos_inc, "inc")
  dec_genes = get_genes(dec_tab, combos_dec, "dec")
  all_combos = union(names(inc_genes), names(dec_genes))
  total = setNames(rep(0,length(all_combos)), all_combos)
  total[names(inc_genes)] = total[names(inc_genes)] + inc_genes
  total[names(dec_genes)] = total[names(dec_genes)] + dec_genes
  combo_order = names(sort(total, decreasing=TRUE))[1:min(top_n,length(total))]
  all_combos_named = setNames(rep(0L, length(combo_order)), combo_order)
  inc_vals = all_combos_named; inc_vals[names(inc_genes)[names(inc_genes) %in% combo_order]] = inc_genes[names(inc_genes) %in% combo_order]
  dec_vals = all_combos_named; dec_vals[names(dec_genes)[names(dec_genes) %in% combo_order]] = dec_genes[names(dec_genes) %in% combo_order]
  bar_df = data.frame(
    combo=factor(combo_order, levels=combo_order),
    n_inc=as.integer(inc_vals),
    n_dec=as.integer(dec_vals)
  )
  p_top = ggplot(bar_df, aes(x=combo, y=n_inc)) +
    geom_col(fill="firebrick") +
    theme_minimal(base_size=9) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), panel.grid.major.x=element_blank()) +
    labs(y="Unique genes\n(increasing)")
  p_bottom = ggplot(bar_df, aes(x=combo, y=-n_dec)) +
    geom_col(fill="steelblue") +
    theme_minimal(base_size=9) +
    theme(axis.text.x=element_text(angle=90,hjust=1,size=6), panel.grid.major.x=element_blank()) +
    labs(y="Unique genes\n(decreasing)", x="Region combination")
  dot_rows = list()
  for (cn in combo_order) {
    inc_regs = if(!is.na(inc_genes[cn])&&inc_genes[cn]>0) strsplit(cn,",")[[1]] else character(0)
    dec_regs = if(!is.na(dec_genes[cn])&&dec_genes[cn]>0) strsplit(cn,",")[[1]] else character(0)
    for (r in all_regions) {
      dot_rows[[paste(cn,r,"inc")]] = data.frame(combo=cn, region=r, direction="Increasing", present=r %in% inc_regs)
      dot_rows[[paste(cn,r,"dec")]] = data.frame(combo=cn, region=r, direction="Decreasing", present=r %in% dec_regs)
    }
  }
  dot_df = do.call(rbind, dot_rows)
  dot_df$combo = factor(dot_df$combo, levels=combo_order)
  dot_df$region = factor(dot_df$region, levels=all_regions)
  dot_df$direction = factor(dot_df$direction, levels=c("Increasing","Decreasing"))
  p_mid = ggplot(dot_df, aes(x=combo, y=region, color=direction)) +
    geom_point(aes(alpha=present, size=present, shape=direction)) +
    scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.07), guide="none") +
    scale_size_manual(values=c("TRUE"=2.5,"FALSE"=0.8), guide="none") +
    scale_shape_manual(values=c("Increasing"=17,"Decreasing"=16), name=NULL) +
    scale_color_manual(values=c("Increasing"="firebrick","Decreasing"="steelblue"), name=NULL) +
    theme_minimal(base_size=9) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
          panel.grid.major=element_line(color="grey92"), legend.position="top") +
    labs(y="Region")
  p_top / p_mid / p_bottom + plot_layout(heights=c(1.5,2.5,1.5)) +
    plot_annotation(title=paste0(ct," (baseline, lfsr<0.05): UpSet by region combination"),
                    subtitle="Bars=unique sig genes. Triangles=increasing, circles=decreasing dispersion with age.")
}

make_fig3F = function(min_genes) {
  all_rows = list()
  for (ct in cell_types) {
    sig = load_sig(ct)
    if(is.null(sig)) next
    for (subcl in unique(sig$cell_type)) {
      sub_sig = sig[sig$cell_type==subcl,]
      n_regions = sum(sapply(unique(sub_sig$region), function(r)
        length(unique(sub_sig$ensembl_id[sub_sig$region==r])) >= min_genes))
      if (n_regions > 0)
        all_rows[[paste(ct,subcl)]] = data.frame(n_regions=n_regions,
          class=ifelse(ct %in% neuronal,"Neuronal","Non-neuronal"))
    }
  }
  per_sub = do.call(rbind, all_rows)
  per_sub$class = factor(per_sub$class, levels=c("Neuronal","Non-neuronal"))
  ggplot(per_sub, aes(x=n_regions, fill=class)) +
    geom_histogram(binwidth=1, color="white", position="stack") +
    scale_fill_manual(values=c("Neuronal"="#E07B39","Non-neuronal"="#3B7DD8"), name=NULL) +
    scale_x_continuous(breaks=1:max(per_sub$n_regions)) +
    scale_y_continuous(breaks=function(x) unique(floor(pretty(x)))) +
    theme_minimal(base_size=14) +
    theme(legend.position="top", panel.grid.major.x=element_blank()) +
    labs(x="Number of regions", y="Number of subclusters",
         title=paste0("Min ", min_genes, " gene(s) per region"))
}

cat("Building UpSet...\n")
for (ct in cell_types) {
  p = make_upset(ct)
  if(!is.null(p)) {
    ggsave(paste0("/scratch/easmit31/dispersion/dglm/figures/upset_",ct,"_baseline_twodots.png"),
           p, width=16, height=12, dpi=150)
    cat(ct, "UpSet saved.\n")
  }
}

cat("Building fig3F...\n")
p1=make_fig3F(1); p10=make_fig3F(10); p100=make_fig3F(100)
combined = p1/p10/p100 +
  plot_annotation(title="All cell types (baseline, lfsr<0.05): subclusters by regions",
                  subtitle="Min unique significant genes per region")
ggsave("/scratch/easmit31/dispersion/dglm/figures/fig3F_ALL_baseline_mingenes_comparison.png",
       combined, width=12, height=16, dpi=150)
cat("fig3F saved.\n")
