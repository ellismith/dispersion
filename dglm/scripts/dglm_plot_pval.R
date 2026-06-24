#!/usr/bin/env Rscript
# dglm_plot_pval.R
# P-value histograms and QQ plots per region for DGLM results.
# Follows Chiou et al. 2022 Nat Neurosci dglm_visualize.R.
#
# Usage:
#   Rscript dglm_plot_pval.R --cell_type microglia

source('/scratch/easmit31/variability/dglm/scripts/_include_options.R')

library(optparse)
library(ggplot2)
library(reshape2)

option_list = list(
    make_option('--cell_type',   type='character', help='cell type to analyze'),
    make_option('--checkpoints', type='character', default='/scratch/easmit31/variability/dglm/checkpoints'),
    make_option('--figdir',      type='character', default='/scratch/easmit31/variability/dglm/figures'),
    make_option('--outfmt',      type='character', default='png', help='output format: png or pdf')
)
opt = parse_args(OptionParser(option_list=option_list))

cell.type = opt$cell_type
dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)
message('Cell type: ', cell.type)

save_fig = function(p, name, height=5, width=8) {
    out = file.path(opt$figdir, paste0(cell.type, '_', name, '.', opt$outfmt))
    if (opt$outfmt == 'png') {
        ggsave(p, file=out, height=height, width=width, dpi=150)
    } else {
        ggsave(p, file=out, height=height, width=width, useDingbats=FALSE)
    }
    message('Saved: ', out)
}

# ── load results ──────────────────────────────────────────────────────────
in.file = file.path(opt$checkpoints, paste0(cell.type, '_dglm_mashr_results.rds'))
message('Loading ', in.file)
obj          = readRDS(in.file)
dglm.results = obj$dglm_results
regions      = obj$regions

# ── extract pvals ─────────────────────────────────────────────────────────
dglm.pval = dglm.results[,'pval', regions, drop=FALSE][,1,]
if (is.null(dim(dglm.pval))) {
    dglm.pval = matrix(dglm.pval, ncol=1, dimnames=list(dimnames(dglm.results)[[1]], regions))
}

these.colors = region.colors[regions]

# ── p-value histogram ─────────────────────────────────────────────────────
pval.df = melt(dglm.pval, varnames=c('gene','region'), value.name='pval')
pval.df$region = factor(pval.df$region, levels=regions)

p = ggplot(pval.df, aes(pval, fill=region)) +
    geom_histogram(binwidth=0.005) +
    facet_wrap(~region, nrow=ceiling(sqrt(length(regions)))) +
    scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
    scale_fill_manual(values=these.colors) +
    xlab('p value') + ylab('Count') +
    ggtitle(paste0(cell.type, ' — DGLM p-value distribution')) +
    theme_classic() +
    theme(legend.position='none')
save_fig(p, 'dglm_pval_histogram')

# ── QQ plot ───────────────────────────────────────────────────────────────
qq.df = do.call(rbind, lapply(regions, function(r) {
    x = na.omit(dglm.pval[, r])
    n = length(x)
    data.frame(
        region   = factor(r, levels=regions),
        expected = -log10(seq(1/n, 1, 1/n)),
        observed = -log10(quantile(x, seq(1/n, 1, 1/n), na.rm=TRUE))
    )
}))

p = ggplot(qq.df, aes(expected, observed, color=region)) +
    geom_point(size=0.5) +
    geom_abline(slope=1, col='black', linewidth=0.2) +
    facet_wrap(~region, nrow=ceiling(sqrt(length(regions)))) +
    coord_cartesian(xlim=c(0,4), ylim=c(0,10)) +
    scale_color_manual(values=these.colors) +
    xlab(expression(-log[10]*('Expected'~italic(p)))) +
    ylab(expression(-log[10]*('Observed'~italic(p)))) +
    ggtitle(paste0(cell.type, ' — DGLM QQ plot')) +
    theme_classic() +
    theme(legend.position='none')
save_fig(p, 'dglm_pval_qqplot')

message('done.')
