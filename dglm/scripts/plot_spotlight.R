#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(dplyr)})

make_spotlight = function(ct, reg, pipeline="baseline", out=NULL) {
  lcb = c("opc","microglia","astrocytes","vascular_cells","oligodendrocytes")
  if (pipeline=="baseline") {
    base = if(ct %in% lcb) "disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered" else "disp_age__FINAL_autosomeX_300c_QCfiltered"
  } else {
    base = if(ct %in% lcb) "disp_age__mincells10_300c_noLcb_QCfiltered" else "disp_age__mincells10_300c_QCfiltered"
  }
  f = paste0("/scratch/easmit31/dispersion/dglm/",base,"/",ct,"/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv")
  df = read.csv(f, sep="\t")
  sig = df[df$mash_lfsr<0.05 & df$region==reg,]
  all_tested = df[df$region==reg,]
  if(nrow(sig)==0) { cat("No sig results for",ct,reg,"\n"); return(NULL) }

  sig$subcluster_num = as.integer(sub(".*_([0-9]+)$","\\1",sig$cell_type))
  all_tested$subcluster_num = as.integer(sub(".*_([0-9]+)$","\\1",all_tested$cell_type))

  # per subcluster: n genes tested, n increasing sig, n decreasing sig
  tested_per_sub = all_tested %>%
    group_by(subcluster_num) %>%
    summarise(n_tested=n_distinct(ensembl_id), .groups="drop")

  inc_per_sub = sig %>% filter(mash_beta>0) %>%
    group_by(subcluster_num) %>%
    summarise(n_inc=n_distinct(ensembl_id), .groups="drop")

  dec_per_sub = sig %>% filter(mash_beta<=0) %>%
    group_by(subcluster_num) %>%
    summarise(n_dec=n_distinct(ensembl_id), .groups="drop")

  combined = tested_per_sub %>%
    left_join(inc_per_sub, by="subcluster_num") %>%
    left_join(dec_per_sub, by="subcluster_num") %>%
    mutate(n_inc=replace(n_inc, is.na(n_inc), 0),
           n_dec=replace(n_dec, is.na(n_dec), 0),
           pct_inc = 100*n_inc/n_tested,
           pct_dec = -100*n_dec/n_tested)  # negative so it goes below zero

  plot_df = rbind(
    data.frame(subcluster_num=combined$subcluster_num, pct=combined$pct_inc, direction="Increasing"),
    data.frame(subcluster_num=combined$subcluster_num, pct=combined$pct_dec, direction="Decreasing")
  )
  plot_df$direction = factor(plot_df$direction, levels=c("Increasing","Decreasing"))

  ct_short = gsub("_neurons","",gsub("_cells","",ct))
  title = paste0(ct_short," × ",reg," (lfsr<0.05, ",pipeline,")")
  subtitle = paste0(length(unique(sig$cell_type))," subclusters | y-axis = % of genes tested per subcluster")

  max_pct = max(abs(plot_df$pct), na.rm=TRUE)

  p = ggplot(plot_df, aes(x=factor(subcluster_num), y=pct, fill=direction)) +
    geom_col(width=0.7) +
    geom_hline(yintercept=0, color="black", linewidth=0.5) +
    scale_fill_manual(values=c("Increasing"="firebrick","Decreasing"="steelblue"), name=NULL) +
    scale_y_continuous(limits=c(-max_pct*1.1, max_pct*1.1),
                       labels=function(x) paste0(abs(x),"%")) +
    theme_minimal(base_size=12) +
    theme(legend.position="top", panel.grid.major.x=element_blank()) +
    labs(x="Louvain subcluster number",
         y="% of genes tested\n(above=increasing, below=decreasing dispersion with age)",
         title=title, subtitle=subtitle)

  if(is.null(out)) out = paste0("/scratch/easmit31/dispersion/dglm/figures/spotlight_",ct,"_",reg,".png")
  ggsave(out, p, width=12, height=6, dpi=150)
  cat("Saved:", out, "\n")
}

make_spotlight("oligodendrocytes", "dlPFC")
make_spotlight("astrocytes", "EC")
