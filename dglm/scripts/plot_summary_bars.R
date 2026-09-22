#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(patchwork)})
F <- "disp_age__mask1_915_300c"; LF <- 0.05
CTS = c(astrocytes='AST', basket_cells='BC', cerebellar_neurons='CER', ependymal_cells='EPEN',
        GABAergic_neurons='INH', glutamatergic_neurons='EXC', medium_spiny_neurons='MSN',
        microglia='MGL', midbrain_neurons='MBN', opc='OPC', oligodendrocytes='OLIG', vascular_cells='VASC')
CT_COLORS = setNames(c('#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#E07B39','#3B7DD8',
                       '#A6D854','#FFD92F','#B15928','#80B1D3','#BC80BD','#999999'),
                     sort(unname(CTS)))
REGIONS = c('ACC','CN','dlPFC','EC','HIP','IPP','lCb','M1','MB','mdTN','NAc')
REG_COLORS = setNames(c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#A65628',
                        '#F781BF','#1B9E77','#D95F02','#7570B3','#66A61E'), REGIONS)
alld <- list()
for (ct in list.dirs(F, recursive=FALSE)) {
  f <- file.path(ct, "dglm_checkpoints_cutoff0.5", "master_dglm_combined.tsv")
  if (!file.exists(f)) next
  d <- read.csv(f, sep="\t"); d$parent <- basename(ct); alld[[ct]] <- d
}
d <- do.call(rbind, alld)
s <- d[!is.na(d$mash_lfsr) & d$mash_lfsr < LF,]
mk <- function(key, pal, lab) {
  agg <- do.call(rbind, lapply(split(s, s[[key]]), function(g) data.frame(
    grp=g[[key]][1],
    inc=length(unique(g$ensembl_id[g$mash_beta>0])),
    dec=length(unique(g$ensembl_id[g$mash_beta<=0])))))
  long <- rbind(data.frame(grp=agg$grp, n=agg$inc, dir='Increase'),
                data.frame(grp=agg$grp, n=-agg$dec, dir='Decrease'))
  long$dir <- factor(long$dir, levels=c('Increase','Decrease'))
  ggplot(long, aes(x=grp, y=n, fill=grp, alpha=dir)) +
    geom_col(width=0.8, color='grey25', linewidth=0.25) +
    scale_fill_manual(values=pal, guide='none') +
    scale_alpha_manual(values=c(Increase=1, Decrease=0.42), name=NULL) +
    geom_hline(yintercept=0, color='black', linewidth=0.5) +
    scale_y_continuous(labels=function(x) abs(x)) +
    theme_minimal(base_size=24) +
    theme(axis.text.x=element_text(angle=45,hjust=1), legend.position='top',
          panel.grid.major.x=element_blank()) +
    labs(x=NULL, y='Unique significant genes', title=lab)
}
s$abbr <- CTS[s$parent]
p1 <- mk('abbr', CT_COLORS, 'By cell type')
p2 <- mk('region', REG_COLORS, 'By brain region')
p <- p1 / p2 + plot_annotation(
  title=paste0('Age-associated dispersion changes (mashr lfsr<', LF, ')'),
  subtitle='Solid = increased dispersion, transparent = decreased. Unique genes within each group.',
  theme=theme(plot.title=element_text(size=26, face='bold')))
ggsave('figures_915/summary_bars_celltype_region.png', p, width=13, height=13, dpi=150)
cat('Saved: figures_915/summary_bars_celltype_region.png\n')
