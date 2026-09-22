#!/usr/bin/env Rscript
#
# 2D histogram: for each significant gene, how many distinct regions is it
# significant in (x-axis) vs. how many distinct subclusters (y-axis).
# Fill = number of genes at that (n_regions, n_subclusters) combination.
#
# Usage:
#   Rscript regions_vs_subclusters_hist.R --master <master_dglm_combined.tsv> \
#     --lfsr_thresh 0.2 --title "My Title" --out <output.png>

suppressMessages({library(optparse); library(ggplot2)})

option_list = list(
  make_option('--master', type='character'),
  make_option('--lfsr_thresh', type='double', default=0.2),
  make_option('--title', type='character', default=NULL),
  make_option('--out', type='character', default='regions_vs_subclusters.png')
)
opt = parse_args(OptionParser(option_list=option_list))
if (is.null(opt$master)) stop('--master is required')

df = read.csv(opt$master, sep='\t')
sig = df[df$mash_lfsr < opt$lfsr_thresh, ]

per_gene = aggregate(cbind(region, cell_type) ~ ensembl_id, data=sig, FUN=function(x) length(unique(x)))
colnames(per_gene) = c('ensembl_id', 'n_regions', 'n_subclusters')

cat('Genes summarized:', nrow(per_gene), '\n')
cat('n_regions range:', range(per_gene$n_regions), '| n_subclusters range:', range(per_gene$n_subclusters), '\n')

title = if (!is.null(opt$title)) opt$title else 'Significant genes: regions vs. subclusters spanned'

p = ggplot(per_gene, aes(x=n_regions, y=n_subclusters)) +
  geom_bin2d(binwidth=1) +
  scale_fill_gradient(low='#E6F1FB', high='#042C53', name='Gene count') +
  scale_x_continuous(breaks=seq(0, max(per_gene$n_regions), by=1)) +
  scale_y_continuous(breaks=seq(0, max(per_gene$n_subclusters), by=2)) +
  theme_minimal(base_size=12) +
  labs(x='Number of distinct regions gene is significant in',
       y='Number of distinct subclusters gene is significant in',
       title=title,
       subtitle=paste0(nrow(per_gene), ' significant genes (mash_lfsr < ', opt$lfsr_thresh, '), any direction'))

ggsave(opt$out, p, width=9, height=7, dpi=150)
cat('Saved:', opt$out, '\n')
