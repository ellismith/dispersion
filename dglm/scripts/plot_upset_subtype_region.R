#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(patchwork); library(dplyr); library(optparse)})

option_list = list(
  make_option('--ct', type='character'),
  make_option('--pipeline', type='character', default='baseline'),
  make_option('--mode', type='character', default='subclusters'),
  make_option('--ymetric', type='character', default='pct'),
  make_option('--lfsr', type='double', default=0.05),
  make_option('--out', type='character', default=NULL)
)
opt = parse_args(OptionParser(option_list=option_list))

region_colors = c("ACC"="#5555b0","CN"="#2caa4d","dlPFC"="#345296","EC"="#532b70",
  "HIP"="#84db55","IPP"="#6c9cc0","lCb"="#e6cd69","M1"="#467aa1",
  "MB"="#fc7462","mdTN"="#117867","NAc"="#148e34")
lcb_affected = c("opc","microglia","astrocytes","vascular_cells","oligodendrocytes")
ct = opt$ct
BASE_SIZE = 18

if(opt$pipeline=="baseline") {
  base = if(ct %in% lcb_affected) "disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered" else "disp_age__FINAL_autosomeX_300c_QCfiltered"
} else {
  base = if(ct %in% lcb_affected) "disp_age__mincells10_300c_noLcb_QCfiltered" else "disp_age__mincells10_300c_QCfiltered"
}

f = paste0("/scratch/easmit31/dispersion/dglm/",base,"/",ct,"/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv")
df = read.csv(f,sep="\t")
sig = df[df$mash_lfsr<opt$lfsr,]
sig$subcluster_num = as.integer(sub(".*_([0-9]+)$","\\1",sig$cell_type))
df$subcluster_num = as.integer(sub(".*_([0-9]+)$","\\1",df$cell_type))
all_regions = sort(unique(sig$region))

stats = sig %>% group_by(subcluster_num, region) %>%
  summarise(n_inc=n_distinct(ensembl_id[mash_beta>0]), n_dec=n_distinct(ensembl_id[mash_beta<=0]), .groups="drop")
tested = df %>% group_by(subcluster_num, region) %>%
  summarise(n_tested=n_distinct(ensembl_id), .groups="drop")
stats = left_join(stats, tested, by=c("subcluster_num","region"))
stats$pct_inc = 100*stats$n_inc/stats$n_tested
stats$pct_dec = 100*stats$n_dec/stats$n_tested

if(opt$mode=="subclusters") {
  # use true unique genes across all regions (not sum of per-region counts)
  sub_stats = sig %>% group_by(subcluster_num) %>%
    summarise(
      n_inc_total=n_distinct(ensembl_id[mash_beta>0]),
      n_dec_total=n_distinct(ensembl_id[mash_beta<=0]),
      regions=list(sort(unique(region))),
      .groups="drop") %>% arrange(subcluster_num)
  tested_sub = df %>% group_by(subcluster_num) %>%
    summarise(n_tested_total=n_distinct(ensembl_id), .groups="drop")
  sub_stats = left_join(sub_stats, tested_sub, by="subcluster_num")
  sub_stats$pct_inc = 100*sub_stats$n_inc_total/sub_stats$n_tested_total
  sub_stats$pct_dec = 100*sub_stats$n_dec_total/sub_stats$n_tested_total
  sub_stats$row_label = factor(sub_stats$subcluster_num, levels=rev(sort(unique(sub_stats$subcluster_num))))

  if(opt$ymetric=="pct") {
    vals_inc=sub_stats$pct_inc; vals_dec=-sub_stats$pct_dec
    xlab="% of genes tested\n(right=increasing, left=decreasing)"
  } else {
    vals_inc=sub_stats$n_inc_total; vals_dec=-sub_stats$n_dec_total
    xlab="Unique significant genes\n(right=increasing, left=decreasing)"
  }

  bar_df = rbind(
    data.frame(row_label=sub_stats$row_label, val=vals_inc, direction="Increasing"),
    data.frame(row_label=sub_stats$row_label, val=vals_dec, direction="Decreasing"))

  p_bar = ggplot(bar_df, aes(x=val, y=row_label, fill=direction)) +
    geom_col(width=0.7) + geom_vline(xintercept=0, color="black", linewidth=0.5) +
    scale_fill_manual(values=c("Increasing"="firebrick","Decreasing"="steelblue"), name=NULL) +
    scale_x_continuous(labels=function(x) paste0(abs(round(x)), if(opt$ymetric=="pct") "%" else "")) +
    theme_minimal(base_size=BASE_SIZE) +
    theme(axis.text.y=element_blank(), axis.title.y=element_blank(),
          legend.position="top", panel.grid.major.y=element_blank()) +
    labs(x=xlab)

  dot_rows = list()
  for(i in seq_len(nrow(sub_stats))) {
    sn = sub_stats$subcluster_num[i]; sig_regs = sub_stats$regions[[i]]
    for(r in all_regions) dot_rows[[paste(sn,r)]] = data.frame(
      row_label=factor(sn, levels=levels(sub_stats$row_label)),
      region=factor(r, levels=all_regions), present=r %in% sig_regs)
  }
  dot_df = do.call(rbind, dot_rows)

  p_dot = ggplot(dot_df, aes(x=region, y=row_label, color=region)) +
    geom_point(aes(alpha=present, size=present)) +
    scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.1), guide="none") +
    scale_size_manual(values=c("TRUE"=4,"FALSE"=2), guide="none") +
    scale_color_manual(values=region_colors, name=NULL) +
    theme_minimal(base_size=BASE_SIZE) +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=BASE_SIZE-2),
          panel.grid.major.x=element_line(color="grey90"),
          panel.grid.major.y=element_blank(), legend.position="none") +
    labs(x=NULL, y="Subcluster")

  combined = (p_dot | p_bar) + plot_layout(widths=c(2,3)) +
    plot_annotation(title=paste0(ct," (",opt$pipeline,", lfsr<",opt$lfsr,")"),
      subtitle=paste0("Dots=regions where sig. Bars=",if(opt$ymetric=="pct") "% genes tested" else "unique sig genes"," pooled across sig regions."),
      theme=theme(plot.title=element_text(size=BASE_SIZE+2), plot.subtitle=element_text(size=BASE_SIZE-2)))

} else {
  stats$region = factor(stats$region, levels=all_regions)
  stats = stats %>% arrange(subcluster_num, region)
  stats$row_label = factor(paste0(stats$subcluster_num," × ",stats$region), levels=rev(unique(paste0(stats$subcluster_num," × ",stats$region))))

  bar_df = rbind(
    data.frame(row_label=stats$row_label, val=if(opt$ymetric=="pct") stats$pct_inc else stats$n_inc, direction="Increasing"),
    data.frame(row_label=stats$row_label, val=if(opt$ymetric=="pct") -stats$pct_dec else -stats$n_dec, direction="Decreasing"))

  combined = ggplot(bar_df, aes(x=val, y=row_label, fill=direction)) +
    geom_col(width=0.7) + geom_vline(xintercept=0, color="black", linewidth=0.5) +
    scale_fill_manual(values=c("Increasing"="firebrick","Decreasing"="steelblue"), name=NULL) +
    scale_x_continuous(labels=function(x) paste0(abs(round(x)), if(opt$ymetric=="pct") "%" else "")) +
    theme_minimal(base_size=BASE_SIZE) +
    theme(legend.position="top", panel.grid.major.y=element_blank()) +
    labs(x=if(opt$ymetric=="pct") "% of genes tested" else "Unique sig genes", y="Subcluster × region",
         title=paste0(ct," (",opt$pipeline,", lfsr<",opt$lfsr,")"))
}

if(is.null(opt$out)) opt$out = paste0("/scratch/easmit31/dispersion/dglm/",base,"/",ct,
  "/figures/upset_flexible_",ct,"_",opt$pipeline,"_",opt$mode,"_",opt$ymetric,".png")
h = if(opt$mode=="subclusters") max(8, length(unique(sig$subcluster_num))*0.5) else max(8, nrow(stats)*0.3)
ggsave(opt$out, combined, width=14, height=h, dpi=150)
cat("Saved:", opt$out, "\n")
