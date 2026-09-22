#!/usr/bin/env Rscript
suppressMessages({library(optparse); library(ggplot2)})
opt = parse_args(OptionParser(option_list=list(
  make_option('--ct'), make_option('--base_dir'),
  make_option('--lfsr', type='double', default=0.05), make_option('--out'))))
f = file.path(opt$base_dir, opt$ct, 'dglm_checkpoints_cutoff0.5', 'master_dglm_combined.tsv')
d = read.csv(f, sep='\t')
sig = d[!is.na(d$mash_lfsr) & d$mash_lfsr < opt$lfsr,]
if (nrow(sig)==0) { cat('SKIP: no sig rows\n'); quit(save='no', status=0) }
sig$subnum = sub('.*_([0-9]+)$','\\1',sig$cell_type)
ag = aggregate(cbind(inc=sig$mash_beta>0, dec=sig$mash_beta<=0) ~ subnum + region, data=sig, FUN=sum)
ag$total = ag$inc + ag$dec
ag$prop_inc = ag$inc / ag$total
ag$subnum = factor(as.integer(ag$subnum), levels=sort(unique(as.integer(ag$subnum))))
p = ggplot(ag, aes(x=subnum, y=region, fill=prop_inc)) +
  geom_tile(color='white', linewidth=0.6) +
  geom_text(aes(label=total), size=2.8) +
  scale_fill_gradient2(low='#2166AC', mid='white', high='#B2182B', midpoint=0.5,
                       limits=c(0,1), name='Prop.\nincreasing') +
  theme_minimal(base_size=14) +
  labs(x='Subcluster #', y='Brain region',
       title=paste0(opt$ct, ' (lfsr<', opt$lfsr, '): direction of dispersion change'),
       subtitle='Red = mostly increasing dispersion, blue = mostly decreasing. Number = total significant genes.')
ggsave(opt$out, p, width=max(8, length(unique(ag$subnum))*0.55+3), height=6.5, dpi=150)
cat('Saved:', opt$out, '\n')
