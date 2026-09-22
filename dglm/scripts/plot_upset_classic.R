#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(patchwork); library(dplyr); library(optparse)})

option_list = list(
  make_option('--ct', type='character', default='all'),
  make_option('--pipeline', type='character', default='baseline'),
  make_option('--lfsr', type='double', default=0.05),
  make_option('--top_n', type='integer', default=30),
  make_option('--groupby', type='character', default='region'),
  make_option('--out', type='character', default=NULL)
)
opt = parse_args(OptionParser(option_list=option_list))

region_colors = c("ACC"="#5555b0","CN"="#2caa4d","dlPFC"="#345296","EC"="#532b70",
  "HIP"="#84db55","IPP"="#6c9cc0","lCb"="#e6cd69","M1"="#467aa1",
  "MB"="#fc7462","mdTN"="#117867","NAc"="#148e34")
lcb_affected = c("opc","microglia","astrocytes","vascular_cells","oligodendrocytes")
cell_types = c("GABAergic_neurons","glutamatergic_neurons","cerebellar_neurons",
               "medium_spiny_neurons","midbrain_neurons","basket_cells",
               "microglia","astrocytes","opc","oligodendrocytes","vascular_cells","ependymal_cells")
BASE_SIZE = 16

load_ct = function(ct, pipeline) {
  if(pipeline=="baseline") {
    base = if(ct %in% lcb_affected) "disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered" else "disp_age__FINAL_autosomeX_300c_QCfiltered"
  } else {
    base = if(ct %in% lcb_affected) "disp_age__mincells10_300c_noLcb_QCfiltered" else "disp_age__mincells10_300c_QCfiltered"
  }
  f = paste0("/scratch/easmit31/dispersion/dglm/",base,"/",ct,"/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv")
  if(!file.exists(f)) return(NULL)
  df = read.csv(f,sep="\t"); df$cell_type_broad=ct; df
}

cts = if(opt$ct=="all") cell_types else opt$ct
all_data = do.call(rbind, Filter(Negate(is.null), lapply(cts, load_ct, pipeline=opt$pipeline)))
sig = all_data[all_data$mash_lfsr<opt$lfsr,]
sig$subcluster_num = sub(".*_([0-9]+)$","\\1",sig$cell_type)

make_classic_panel = function(sig, groupby, direction, top_n, all_groups, dot_colors=NULL) {
  filt = if(direction=="inc") sig[sig$mash_beta>0,] else sig[sig$mash_beta<=0,]

  if(groupby=="region") {
    # columns = region combinations
    combos_per_sub = filt %>% group_by(cell_type) %>%
      summarise(combo=paste(sort(unique(region)),collapse=","), .groups="drop")
    combo_genes = filt %>% left_join(combos_per_sub, by="cell_type") %>%
      group_by(combo) %>% summarise(n_genes=n_distinct(ensembl_id), .groups="drop")
    combo_subs = combos_per_sub %>% count(combo, name="n_sub")
    combo_df = left_join(combo_genes, combo_subs, by="combo") %>%
      arrange(desc(n_genes)) %>% head(top_n)
    get_members = function(combo) strsplit(combo,",")[[1]]
    row_label = "Region"; bot_label = "N subclusters"
    dot_col_fn = function(g) if(!is.null(dot_colors)) dot_colors[g] else "grey30"
  } else {
    # columns = subcluster combinations
    combos_per_reg = filt %>% group_by(region) %>%
      summarise(combo=paste(sort(unique(subcluster_num)),collapse=","), .groups="drop")
    combo_genes = filt %>% left_join(combos_per_reg, by="region") %>%
      group_by(combo) %>% summarise(n_genes=n_distinct(ensembl_id), .groups="drop")
    combo_regs = combos_per_reg %>% count(combo, name="n_sub")
    combo_df = left_join(combo_genes, combo_regs, by="combo") %>%
      arrange(desc(n_genes)) %>% head(top_n)
    get_members = function(combo) strsplit(combo,",")[[1]]
    row_label = "Subcluster"; bot_label = "N regions"
    dot_col_fn = function(g) "grey30"
  }

  combos = combo_df$combo
  x_labels = as.character(seq_along(combos))
  combo_order = factor(x_labels, levels=x_labels)

  bar_df = data.frame(x=combo_order, genes=combo_df$n_genes, n_sub=combo_df$n_sub)

  p_top = ggplot(bar_df, aes(x=x, y=genes)) +
    geom_col(fill=if(direction=="inc") "firebrick" else "steelblue", width=0.7) +
    theme_minimal(base_size=BASE_SIZE) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), panel.grid.major.x=element_blank()) +
    labs(y=paste0("Unique genes\n(",if(direction=="inc") "increasing" else "decreasing",")"))

  dot_rows = list()
  for(i in seq_along(combos)) {
    members = get_members(as.character(combos[i]))
    for(g in all_groups) {
      dot_rows[[paste(i,g)]] = data.frame(x=combo_order[i], group=g, present=g %in% members)
    }
  }
  dot_df = do.call(rbind, dot_rows)
  dot_df$group = factor(dot_df$group, levels=rev(all_groups))

  if(groupby=="region" && !is.null(dot_colors)) {
    p_mid = ggplot(dot_df, aes(x=x, y=group, color=group)) +
      geom_point(aes(alpha=present, size=present)) +
      scale_color_manual(values=dot_colors, guide="none") +
      scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.08), guide="none") +
      scale_size_manual(values=c("TRUE"=3.5,"FALSE"=1.5), guide="none") +
      theme_minimal(base_size=BASE_SIZE) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
            panel.grid.major.x=element_line(color="grey90"), panel.grid.major.y=element_blank()) +
      labs(y=row_label)
  } else {
    p_mid = ggplot(dot_df, aes(x=x, y=group)) +
      geom_point(aes(alpha=present, size=present), color="grey30") +
      scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.08), guide="none") +
      scale_size_manual(values=c("TRUE"=3,"FALSE"=1.5), guide="none") +
      theme_minimal(base_size=BASE_SIZE) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
            panel.grid.major.x=element_line(color="grey90"), panel.grid.major.y=element_blank()) +
      labs(y=row_label)
  }

  p_bot = ggplot(bar_df, aes(x=x, y=n_sub)) +
    geom_col(fill="grey50", width=0.7) +
    theme_minimal(base_size=BASE_SIZE) +
    theme(panel.grid.major.x=element_blank()) +
    labs(y=bot_label, x="Combination #")

  p_top / p_mid / p_bot + plot_layout(heights=c(2,3,1))
}

if(opt$groupby=="region") {
  all_groups = sort(unique(sig$region))
  dot_colors = region_colors[all_groups]
} else {
  all_groups = sort(unique(sig$subcluster_num), decreasing=FALSE)
  dot_colors = NULL
}

p_inc = make_classic_panel(sig, opt$groupby, "inc", opt$top_n, all_groups, dot_colors)
p_dec = make_classic_panel(sig, opt$groupby, "dec", opt$top_n, all_groups, dot_colors)

ct_label = if(opt$ct=="all") "All cell types" else opt$ct
combined = (p_inc | p_dec) +
  plot_annotation(
    title=paste0(ct_label," (",opt$pipeline,", lfsr<",opt$lfsr,"): UpSet by ",opt$groupby," combination"),
    subtitle=paste0("Left=increasing dispersion, Right=decreasing. Top=unique genes, Middle=",opt$groupby,"s, Bottom=n ",if(opt$groupby=="region") "subclusters" else "regions","."),
    theme=theme(plot.title=element_text(size=BASE_SIZE+2), plot.subtitle=element_text(size=BASE_SIZE-2)))

if(is.null(opt$out)) {
  ct_file = if(opt$ct=="all") "ALL" else opt$ct
  opt$out = paste0("/scratch/easmit31/dispersion/dglm/figures/upset_classic_",opt$groupby,"_",ct_file,"_",opt$pipeline,".png")
}
h = max(10, length(all_groups)*0.4+5)
ggsave(opt$out, combined, width=20, height=h, dpi=150)
cat("Saved:", opt$out, "\n")
