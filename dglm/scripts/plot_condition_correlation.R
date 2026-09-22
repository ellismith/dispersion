#!/usr/bin/env Rscript
suppressMessages({library(optparse); library(ggplot2); library(reshape2)})

option_list = list(
  make_option('--master', type='character'),
  make_option('--all', action='store_true', default=FALSE),
  make_option('--pipeline', type='character', default='baseline'),
  make_option('--lfsr_thresh', type='double', default=0.05),
  make_option('--subset_by', type='character', default='none',
    help='none|top_n|direction'),
  make_option('--top_n', type='integer', default=30),
  make_option('--direction', type='character', default='inc',
    help='inc or dec (used when --subset_by direction)'),
  make_option('--title', type='character', default=NULL),
  make_option('--out', type='character', default='condition_correlation.png')
)
opt = parse_args(OptionParser(option_list=option_list))

lcb_affected = c("opc","microglia","astrocytes","vascular_cells","oligodendrocytes")

get_master = function(ct, pipeline) {
  if(pipeline=="baseline") {
    base = if(ct %in% lcb_affected) "disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered" else "disp_age__FINAL_autosomeX_300c_QCfiltered"
  } else {
    base = if(ct %in% lcb_affected) "disp_age__mincells10_300c_noLcb_QCfiltered" else "disp_age__mincells10_300c_QCfiltered"
  }
  file.path("/scratch/easmit31/dispersion/dglm",base,ct,"dglm_checkpoints_cutoff0.5","master_dglm_combined.tsv")
}

if(opt$all) {
  cts = c("GABAergic_neurons","glutamatergic_neurons","cerebellar_neurons",
          "medium_spiny_neurons","midbrain_neurons","basket_cells",
          "microglia","astrocytes","opc","oligodendrocytes","vascular_cells","ependymal_cells")
  dfs = Filter(Negate(is.null), lapply(cts, function(ct) {
    f = get_master(ct, opt$pipeline)
    if(!file.exists(f)) return(NULL)
    read.csv(f, sep="\t")
  }))
  df = do.call(rbind, dfs)
} else {
  df = read.csv(opt$master, sep="\t")
}

sig = df[df$mash_lfsr < opt$lfsr_thresh,]

# subset conditions
if(opt$subset_by == "top_n") {
  cond_counts = aggregate(ensembl_id ~ cell_type + region, data=sig, FUN=function(x) length(unique(x)))
  cond_counts = cond_counts[order(-cond_counts$ensembl_id),]
  top_conds = head(cond_counts, opt$top_n)
  sig = sig[paste(sig$cell_type, sig$region) %in% paste(top_conds$cell_type, top_conds$region),]
  cat("Subsetting to top", opt$top_n, "conditions by n significant genes\n")
} else if(opt$subset_by == "direction") {
  if(opt$direction == "inc") {
    sig = sig[sig$mash_beta > 0,]
    cat("Subsetting to increasing-dispersion genes only\n")
  } else {
    sig = sig[sig$mash_beta <= 0,]
    cat("Subsetting to decreasing-dispersion genes only\n")
  }
}

sig$subnum = sub(".*_([0-9]+)$","\\1",sig$cell_type)
sig$condition = paste0(sig$region,"_",sig$subnum)
sig_genes = unique(sig$ensembl_id)
cat("Genes significant in >=1 condition:", length(sig_genes), "\n")
cat("Unique conditions:", length(unique(sig$condition)), "\n")

sub = df[df$ensembl_id %in% sig_genes,]
sub$subnum = sub(".*_([0-9]+)$","\\1",sub$cell_type)
sub$condition = paste0(sub$region,"_",sub$subnum)
mat = reshape2::acast(sub, ensembl_id ~ condition, value.var="mash_beta", fun.aggregate=mean)
cat("Matrix:", nrow(mat), "genes x", ncol(mat), "conditions\n")

# pairwise correlation
conditions = colnames(mat)
cor_mat = matrix(NA, nrow=length(conditions), ncol=length(conditions),
                 dimnames=list(conditions,conditions))
for(i in seq_along(conditions)) {
  for(j in i:length(conditions)) {
    both = !is.na(mat[,i]) & !is.na(mat[,j])
    if(sum(both) >= 10) {
      r = cor(mat[both,i], mat[both,j], method="pearson")
      cor_mat[i,j] = r; cor_mat[j,i] = r
    }
  }
}
diag(cor_mat) = 1
cat("Correlation range:", round(range(cor_mat,na.rm=TRUE),3), "\n")

cor_for_clust = cor_mat; cor_for_clust[is.na(cor_for_clust)] = 0
ord = hclust(as.dist(1-cor_for_clust), method="ward.D2")$order
cor_ordered = cor_mat[ord,ord]
cor_melt = reshape2::melt(cor_ordered)
cor_melt$Var1 = factor(cor_melt$Var1, levels=rownames(cor_ordered))
cor_melt$Var2 = factor(cor_melt$Var2, levels=rownames(cor_ordered))

n_cond = ncol(cor_mat)
txt_size = max(5, min(10, 160/n_cond))
title = if(!is.null(opt$title)) opt$title else "Condition correlation of dispersion-age effect sizes"
subtitle = paste0(length(sig_genes)," genes (lfsr<",opt$lfsr_thresh,"), ",n_cond," conditions. Ward linkage. Grey=<10 shared genes.")

p = ggplot(cor_melt, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                       midpoint=0, limits=c(-1,1), name="Pearson r", na.value="grey90") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=90,hjust=1,size=txt_size),
        axis.text.y=element_text(size=txt_size)) +
  labs(x=NULL, y=NULL, title=title, subtitle=subtitle)

dim_in = max(8, min(22, n_cond*0.15))
ggsave(opt$out, p, width=dim_in, height=dim_in-1, dpi=150)
cat("Saved:", opt$out, "\n")
