#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(patchwork); library(dplyr)})

lcb = c("opc","microglia","astrocytes","vascular_cells","oligodendrocytes")
all_regions = c("ACC","CN","dlPFC","EC","HIP","IPP","lCb","M1","MB","mdTN","NAc")

make_upset = function(ct, pipeline="baseline", top_n=25, out=NULL) {
  if(pipeline=="baseline") {
    base = if(ct %in% lcb) "disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered" else "disp_age__FINAL_autosomeX_300c_QCfiltered"
  } else {
    base = if(ct %in% lcb) "disp_age__mincells10_300c_noLcb_QCfiltered" else "disp_age__mincells10_300c_QCfiltered"
  }
  f = paste0("/scratch/easmit31/dispersion/dglm/",base,"/",ct,"/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv")
  if(!file.exists(f)) { cat("MISSING:",f,"\n"); return(NULL) }
  df = read.csv(f,sep="\t")
  sig = df[df$mash_lfsr<0.05,]

  subclusters = unique(sig$cell_type)

  # for each subcluster, get its set of increasing regions and decreasing regions
  get_combo = function(s, dir) {
    regs = if(dir=="inc") unique(sig$region[sig$cell_type==s & sig$mash_beta>0]) else unique(sig$region[sig$cell_type==s & sig$mash_beta<=0])
    paste(sort(regs), collapse=",")
  }

  combo_inc = sapply(subclusters, get_combo, dir="inc")
  combo_dec = sapply(subclusters, get_combo, dir="dec")

  # count unique genes per combination
  count_genes = function(combos, dir) {
    tab = table(combos[combos!=""])
    sapply(names(tab), function(combo) {
      regs = strsplit(combo,",")[[1]]
      subs = subclusters[combos==combo]
      filt = sig[sig$cell_type %in% subs & sig$region %in% regs,]
      if(dir=="inc") filt = filt[filt$mash_beta>0,]
      else filt = filt[filt$mash_beta<=0,]
      length(unique(filt$ensembl_id))
    })
  }

  genes_inc = count_genes(combo_inc, "inc")
  genes_dec = count_genes(combo_dec, "dec")

  all_combos = union(names(genes_inc), names(genes_dec))
  total = setNames(numeric(length(all_combos)), all_combos)
  total[names(genes_inc)] = total[names(genes_inc)] + genes_inc
  total[names(genes_dec)] = total[names(genes_dec)] + genes_dec
  combo_order = names(sort(total, decreasing=TRUE))[1:min(top_n, length(total))]

  gi = setNames(rep(0L, length(combo_order)), combo_order)
  gd = setNames(rep(0L, length(combo_order)), combo_order)
  gi[intersect(names(genes_inc), combo_order)] = genes_inc[intersect(names(genes_inc), combo_order)]
  gd[intersect(names(genes_dec), combo_order)] = genes_dec[intersect(names(genes_dec), combo_order)]

  bar_df = data.frame(
    combo=factor(combo_order, levels=combo_order),
    n_inc=as.integer(gi), n_dec=as.integer(gd),
    total=as.integer(gi)+as.integer(gd)
  )
  bar_df$net = ifelse(bar_df$n_inc >= bar_df$n_dec, "Increasing", "Decreasing")

  # bar plot
  p_bar = ggplot(bar_df, aes(x=combo, y=total, fill=net)) +
    geom_col(width=0.7) +
    scale_fill_manual(values=c("Increasing"="firebrick","Decreasing"="steelblue"), name="Net direction") +
    theme_minimal(base_size=9) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
          panel.grid.major.x=element_blank(), legend.position="top") +
    labs(y="Unique sig genes")

  # split dot matrix
  dot_rows = list()
  for (cn in combo_order) {
    inc_regs = if(cn %in% names(genes_inc) && genes_inc[cn]>0) strsplit(cn,",")[[1]] else character(0)
    dec_regs = if(cn %in% names(genes_dec) && genes_dec[cn]>0) strsplit(cn,",")[[1]] else character(0)
    for (r in all_regions) {
      has_inc = r %in% inc_regs
      has_dec = r %in% dec_regs
      # left half dot (increasing)
      dot_rows[[paste(cn,r,"inc")]] = data.frame(combo=cn, region=r, xoff=-0.15, fill=ifelse(has_inc,"firebrick","grey90"), size=ifelse(has_inc,3,1.5), alpha=ifelse(has_inc,1,0.3))
      # right half dot (decreasing)
      dot_rows[[paste(cn,r,"dec")]] = data.frame(combo=cn, region=r, xoff=0.15, fill=ifelse(has_dec,"steelblue","grey90"), size=ifelse(has_dec,3,1.5), alpha=ifelse(has_dec,1,0.3))
    }
  }
  dot_df = do.call(rbind, dot_rows)
  dot_df$combo = factor(dot_df$combo, levels=combo_order)
  dot_df$region = factor(dot_df$region, levels=rev(all_regions))
  dot_df$x_num = as.numeric(dot_df$combo) + dot_df$xoff

  p_dot = ggplot(dot_df, aes(x=x_num, y=region)) +
    geom_point(aes(color=fill, size=size, alpha=alpha), shape=16) +
    scale_color_identity() +
    scale_size_identity() +
    scale_alpha_identity() +
    scale_x_continuous(breaks=seq_along(combo_order),
                       labels=levels(bar_df$combo),
                       limits=c(0.5, length(combo_order)+0.5)) +
    theme_minimal(base_size=8) +
    theme(axis.text.x=element_text(angle=90, hjust=1, size=6),
          panel.grid.major=element_line(color="grey94"),
          panel.grid.minor=element_blank()) +
    labs(x="Region combination", y=NULL)

  combined = p_bar / p_dot + plot_layout(heights=c(1.5,3)) +
    plot_annotation(
      title=paste0(ct," (",pipeline,", lfsr<0.05): UpSet by region combination"),
      subtitle="Bar color = net direction. Split dots: red left half = increasing, blue right half = decreasing dispersion."
    )

  if(is.null(out)) out = paste0("/scratch/easmit31/dispersion/dglm/figures/upset_splitdot_",ct,"_",pipeline,".png")
  ggsave(out, combined, width=16, height=10, dpi=150)
  cat("Saved:", out, "\n")
}

# run for GABAergic as test
make_upset("GABAergic_neurons", "baseline")
