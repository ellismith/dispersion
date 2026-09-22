#!/usr/bin/env Rscript
suppressMessages({library(optparse); library(ggplot2); library(patchwork)})
opt = parse_args(OptionParser(option_list=list(
  make_option('--base_dir'), make_option('--lfsr', type='double', default=0.05),
  make_option('--top_n', type='integer', default=8), make_option('--drop_lcb', action='store_true', default=FALSE), make_option('--font_scale', type='double', default=1), make_option('--out'))))
CTS = c(astrocytes='AST', ependymal_cells='EPEN', GABAergic_neurons='INH',
        glutamatergic_neurons='EXC', medium_spiny_neurons='MSN', microglia='MGL',
        midbrain_neurons='MBN', opc='OPC', oligodendrocytes='OLIG', vascular_cells='VASC')
CT_COLORS = c(AST='#66C2A5', EPEN='#E78AC3', INH='#3B7DD8', EXC='#E07B39', MSN='#A6D854',
              MGL='#FFD92F', MBN='#B15928', OPC='#80B1D3', OLIG='#BC80BD', VASC='#999999')
REGIONS = c('ACC','CN','dlPFC','EC','HIP','IPP','lCb','M1','MB','mdTN','NAc')
if (opt$drop_lcb) REGIONS = setdiff(REGIONS, 'lCb')   # no lCb data after condition filter
bars <- list(); tiles <- list(); order_rows <- c()
for (ct in names(CTS)) {
  f = file.path(opt$base_dir, ct, 'dglm_checkpoints_cutoff0.5', 'master_dglm_combined.tsv')
  if (!file.exists(f)) next
  d = read.csv(f, sep='\t'); s = d[!is.na(d$mash_lfsr) & d$mash_lfsr < opt$lfsr,]
  if (nrow(s)==0) next
  ab = CTS[[ct]]
  tot = sort(table(s$cell_type), decreasing=TRUE)
  keep = names(tot)[seq_len(min(opt$top_n, length(tot)))]
  keep = keep[order(as.integer(sub('.*_([0-9]+)$','\\1', keep)))]
  for (sc in keep) {
    ss = s[s$cell_type==sc,]
    lab = paste0(ab, "_", sub('.*_([0-9]+)$','\\1', sc))
    order_rows = c(order_rows, lab)
    bars[[lab]] = data.frame(row=lab, abbr=ab,
      inc=length(unique(ss$ensembl_id[ss$mash_beta>0])),
      dec=length(unique(ss$ensembl_id[ss$mash_beta<=0])))
    for (r in unique(ss$region)) {
      sr = ss[ss$region==r,]
      pi = mean(sr$mash_beta > 0)
      tiles[[paste(lab,r)]] = data.frame(row=lab, region=r,
        dircol=if (pi >= 0.6) 'inc' else if (pi <= 0.4) 'dec' else 'mixed')
    } } }
bar = do.call(rbind, bars); td = do.call(rbind, tiles)
lv = rev(order_rows)
bar$row = factor(bar$row, levels=lv); td$row = factor(td$row, levels=lv)
td$region = factor(td$region, levels=REGIONS)
grid = expand.grid(row=factor(lv, levels=lv), region=factor(REGIONS, levels=REGIONS))
BS = 15
p_l = ggplot(bar, aes(x=-dec, y=row, fill=abbr)) + geom_col(width=0.72) +
  scale_fill_manual(values=CT_COLORS, guide='none') +
  scale_x_continuous(labels=function(x) abs(x)) + theme_minimal(base_size=BS) +
  theme(panel.grid.major.y=element_blank(), axis.text.y=element_text(size=BS-5)) +
  labs(x='Decreased (unique genes)', y=NULL)
p_m = ggplot() +
  geom_tile(data=grid, aes(x=region, y=row), fill='grey94', width=0.72, height=0.72) +
  geom_tile(data=td, aes(x=region, y=row, fill=dircol), width=0.72, height=0.72) +
  scale_fill_manual(values=c(inc='firebrick', dec='steelblue', mixed='grey60'),
                    labels=c(inc='\u226560% increasing', dec='\u226560% decreasing', mixed='Mixed (40-60%)'), name=NULL) +
  theme_minimal(base_size=BS) +
  theme(axis.text.x=element_text(angle=45,hjust=1), axis.text.y=element_blank(),
        panel.grid=element_blank(), legend.position='top') + labs(x='Brain region', y=NULL)
p_r = ggplot(bar, aes(x=inc, y=row, fill=abbr)) + geom_col(width=0.72) +
  scale_fill_manual(values=CT_COLORS, guide='none') + theme_minimal(base_size=BS) +
  theme(axis.text.y=element_blank(), panel.grid.major.y=element_blank()) +
  labs(x='Increased (unique genes)', y=NULL)
p = (p_l | p_m | p_r) + plot_layout(widths=c(1.3, 2.1, 1.3)) +
  plot_annotation(title=paste0('Subcluster-level dispersion changes (lfsr<', opt$lfsr,
    '), top ', opt$top_n, '/cell type, non-cerebellar'),
    theme=theme(plot.title=element_text(size=BS+3, face='bold')))
update_geom_defaults('text', list(size = 3.88 * opt$font_scale))
.th = theme(text = element_text(size = 11 * opt$font_scale))
p = if (inherits(p, 'patchwork')) p & .th else p + .th
ggsave(opt$out, p, width=14, height=max(10, length(lv)*0.28), dpi=150)
cat('Saved:', opt$out, '| rows:', length(lv), '\n')
