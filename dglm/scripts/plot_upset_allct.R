#!/usr/bin/env Rscript
suppressMessages({library(optparse); library(ggplot2); library(patchwork)})
opt = parse_args(OptionParser(option_list=list(
  make_option('--base_dir'), make_option('--lfsr', type='double', default=0.05),
  make_option('--out'), make_option('--min_pct', type='double', default=0))))
CTS = c(astrocytes='AST', basket_cells='BC', cerebellar_neurons='CER', ependymal_cells='EPEN',
        GABAergic_neurons='INH', glutamatergic_neurons='EXC', medium_spiny_neurons='MSN',
        microglia='MGL', midbrain_neurons='MBN', opc='OPC', oligodendrocytes='OLIG', vascular_cells='VASC')
CT_COLORS = setNames(c('#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#E07B39','#3B7DD8',
                       '#A6D854','#FFD92F','#B15928','#80B1D3','#BC80BD','#999999'),
                     sort(unname(CTS)))
REGIONS = c('ACC','CN','dlPFC','EC','HIP','IPP','lCb','M1','MB','mdTN','NAc')
rows = list(); tiles = list()
for (ct in names(CTS)) {
  f = file.path(opt$base_dir, ct, 'dglm_checkpoints_cutoff0.5', 'master_dglm_combined.tsv')
  if (!file.exists(f)) next
  d = read.csv(f, sep='\t'); s = d[!is.na(d$mash_lfsr) & d$mash_lfsr < opt$lfsr,]
  if (opt$min_pct > 0) {
    # pooled across ALL subclusters of this cell type, per region
    tst = tapply(!is.na(d$mash_lfsr), d$region, sum)
    sg  = tapply(!is.na(d$mash_lfsr) & d$mash_lfsr < opt$lfsr, d$region, sum)
    keep = names(tst)[100*sg/pmax(tst,1) >= opt$min_pct]
    s = s[s$region %in% keep,]
  }
  if (nrow(s)==0) next
  ab = CTS[[ct]]
  rows[[ab]] = data.frame(abbr=ab, n_inc=length(unique(s$ensembl_id[s$mash_beta>0])),
                          n_dec=length(unique(s$ensembl_id[s$mash_beta<=0])))
  for (r in unique(s$region)) {
    sr = s[s$region==r,]
    tiles[[paste(ab,r)]] = data.frame(abbr=ab, region=r,
      dircol=if (sum(sr$mash_beta>0) > sum(sr$mash_beta<=0)) 'inc' else 'dec')
  }
}
bar = do.call(rbind, rows); td = do.call(rbind, tiles)
ord = rev(sort(bar$abbr))          # alphabetical, AST at top
bar$abbr = factor(bar$abbr, levels=ord); td$abbr = factor(td$abbr, levels=ord)
td$region = factor(td$region, levels=REGIONS)
grid = expand.grid(abbr=factor(ord, levels=ord), region=factor(REGIONS, levels=REGIONS))
BS = 22
p_l = ggplot(bar, aes(x=-n_dec, y=abbr, fill=abbr)) + geom_col(width=0.7) +
  scale_fill_manual(values=CT_COLORS, guide='none') +
  scale_x_continuous(labels=function(x) abs(x)) + theme_minimal(base_size=BS) +
  theme(panel.grid.major.y=element_blank()) + labs(x='Decreased\n(unique genes)', y=NULL)
p_m = ggplot() +
  geom_tile(data=grid, aes(x=region, y=abbr), fill='grey93', width=0.7, height=0.7) +
  geom_tile(data=td, aes(x=region, y=abbr, fill=dircol), width=0.7, height=0.7) +
  scale_fill_manual(values=c(inc='firebrick', dec='steelblue'),
                    labels=c(inc='Mostly increasing', dec='Mostly decreasing'), name=NULL) +
  theme_minimal(base_size=BS) +
  theme(axis.text.x=element_text(angle=45,hjust=1), axis.text.y=element_blank(),
        panel.grid=element_blank(), legend.position='top') +
  labs(x='Brain region', y=NULL)
p_r = ggplot(bar, aes(x=n_inc, y=abbr, fill=abbr)) + geom_col(width=0.7) +
  scale_fill_manual(values=CT_COLORS, guide='none') + theme_minimal(base_size=BS) +
  theme(axis.text.y=element_blank(), panel.grid.major.y=element_blank()) +
  labs(x='Increased\n(unique genes)', y=NULL)
p = (p_l | p_m | p_r) + plot_layout(widths=c(1.4, 2.2, 1.4)) +
  plot_annotation(title=paste0('All cell types (lfsr<', opt$lfsr, '): dispersion changes by region'),
    subtitle='Tiles: region has significant genes; color = majority direction. Bars: unique gene counts per direction.',
    theme=theme(plot.title=element_text(size=BS+2, face='bold'), plot.subtitle=element_text(size=BS-6)))
ggsave(opt$out, p, width=13, height=10, dpi=150)
cat('Saved:', opt$out, '\n')
