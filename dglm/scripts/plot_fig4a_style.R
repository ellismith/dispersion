#!/usr/bin/env Rscript
# plot_fig4a_style.R
# Diverging bar chart of significant dispersion-change genes per region,
# matching the reference Fig 4a style. Two independent toggles:
#   --count_mode rows  : gene x subcluster combinations (current pipeline's
#                        native counting -- a gene significant in multiple
#                        subclusters counts once per subcluster)
#   --count_mode genes : distinct genes only (deduped across subclusters --
#                        a gene counts as "increase" if it increases in
#                        >=1 subcluster, "decrease" if it decreases in
#                        >=1 subcluster; a gene can count in both if
#                        inconsistent across subclusters)
#   --display net      : one diverging bar per region (increase - decrease)
#   --display split     : two bars per region, positive up / negative down
#                        (matches the existing bar_sig_by_region style)

suppressMessages({library(optparse); library(ggplot2)})

option_list = list(
    make_option('--base_dir', type='character', default='/scratch/easmit31/dispersion/dglm/disp_age__min300cells_min30animals'),
    make_option('--cell_type', type='character', default='microglia'),
    make_option('--cutoff', type='character', default='0.5'),
    make_option('--lfsr_thresh', type='double', default=0.2),
    make_option('--count_mode', type='character', default='genes'),
    make_option('--display', type='character', default='net'),
    make_option('--figdir', type='character', default='/scratch/easmit31/dispersion/dglm/figures')
)
opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

region_colors = c("NAc"="#148e34","ACC"="#5555b0","IPP"="#6c9cc0","CN"="#2caa4d",
                   "HIP"="#84db55","lCb"="#e6cd69","mdTN"="#117867","dlPFC"="#345296",
                   "M1"="#467aa1","EC"="#532b70","MB"="#fc7462")

tsv = file.path(opt$base_dir, opt$cell_type, paste0('dglm_checkpoints_cutoff', opt$cutoff), 'master_dglm_combined.tsv')
df = read.table(tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
sig = df[!is.na(df$mash_lfsr) & df$mash_lfsr < opt$lfsr_thresh, ]

if (opt$count_mode == 'rows') {
    agg = aggregate(list(n=sig$ensembl_id), by=list(region=sig$region, dir=ifelse(sig$mash_beta>0,'increase','decrease')), FUN=length)
} else {
    inc = unique(sig[sig$mash_beta > 0, c('region','ensembl_id')]); inc$dir = 'increase'
    dec = unique(sig[sig$mash_beta <= 0, c('region','ensembl_id')]); dec$dir = 'decrease'
    both = rbind(inc, dec)
    agg = aggregate(list(n=both$ensembl_id), by=list(region=both$region, dir=both$dir), FUN=length)
}

wide = reshape(agg, idvar='region', timevar='dir', direction='wide')
colnames(wide) = sub('^n\\.', '', colnames(wide))
wide$increase[is.na(wide$increase)] = 0
wide$decrease[is.na(wide$decrease)] = 0
wide$net = wide$increase - wide$decrease
wide = wide[order(-wide$net), ]
wide$region = factor(wide$region, levels=wide$region)

title_suffix = paste0(opt$count_mode, '_', opt$display)
out = file.path(opt$figdir, paste0('fig4a_style_', opt$cell_type, '_', title_suffix, '.png'))

if (opt$display == 'net') {
    p = ggplot(wide, aes(x=region, y=net, fill=region)) +
        geom_col(color='black', linewidth=0.3) +
        scale_fill_manual(values=region_colors, guide='none') +
        theme_classic(base_size=18) +
        ylab('Number of genes') + xlab(NULL) +
        ggtitle(paste0(opt$cell_type, ': net dispersion change by region (', opt$count_mode, ' count)')) +
        theme(axis.text.x=element_text(angle=45, hjust=1, size=20, face='bold'),
              axis.text.y=element_text(size=16),
              axis.title.y=element_text(size=18),
              plot.title=element_text(size=16))
} else {
    long = rbind(
        data.frame(region=wide$region, dir='increase', n=wide$increase),
        data.frame(region=wide$region, dir='decrease', n=-wide$decrease)
    )
    p = ggplot(long, aes(x=region, y=n, fill=region, alpha=dir)) +
        geom_col(color='black', linewidth=0.3) +
        scale_fill_manual(values=region_colors, guide='none') +
        scale_alpha_manual(values=c(increase=1, decrease=0.55), name='direction') +
        theme_classic(base_size=18) +
        ylab('Number of genes') + xlab(NULL) +
        ggtitle(paste0(opt$cell_type, ': dispersion change by region (', opt$count_mode, ' count, split)')) +
        theme(axis.text.x=element_text(angle=45, hjust=1, size=20, face='bold'),
              axis.text.y=element_text(size=16),
              axis.title.y=element_text(size=18),
              plot.title=element_text(size=16))
}

ggsave(p, file=out, width=9, height=6, dpi=150)
message('Saved: ', out)
