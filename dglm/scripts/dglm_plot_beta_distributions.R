#!/usr/bin/env Rscript
# dglm_plot_beta_distributions.R
#
# Sanity-check plot: for each cell type, one PNG showing the distribution
# of raw DGLM beta values (dispersion-with-age coefficient), faceted by
# region. X-axis is scaled to a percentile window of the real data (not
# the true min/max, which gets stretched far out by a handful of outlier
# genes -- using true min/max squashes the actual distribution shape into
# one spike; using a fixed small clip instead cuts off cell types with
# genuinely wider real spread). True min/max is always reported in the
# console and plot subtitle regardless of what's plotted.
#
# Usage:
#   Rscript dglm_plot_beta_distributions.R --checkpoints /path/to/model_dir --figdir /path/to/output

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(ggplot2)

option_list = list(
    make_option('--checkpoints', type='character', default='/scratch/easmit31/dispersion/dglm/checkpoints'),
    make_option('--figdir',      type='character', default='/scratch/easmit31/dispersion/dglm/figures'),
    make_option('--outfmt',      type='character', default='png'),
    make_option('--percentile',  type='double', default=0.99,
                help='plot window = this percentile of |beta| per cell type (default 0.99 -- the middle 99%% of the real data), not the true min/max, so a few outlier genes don\'t stretch the axis and flatten the real distribution shape. Set to 1.0 to use the true full range instead.'),
    make_option('--clip',        type='double', default=NA,
                help='overrides --percentile: use this fixed +/- x-axis window for every cell type instead (useful for directly comparing axes across cell types).')
)
opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

rds.files = list.files(opt$checkpoints, pattern='_dglm_results\\.rds$')
cell.types = sub('_dglm_results\\.rds$', '', rds.files)
message('Found ', length(cell.types), ' unit(s): ', paste(cell.types, collapse=', '))

for (ct in cell.types) {
    obj = readRDS(file.path(opt$checkpoints, paste0(ct, '_dglm_results.rds')))
    arr = obj$array

    stat.names = dimnames(arr)[[2]]
    beta.col = if ('beta' %in% stat.names) 'beta' else 'beta_age'

    regions = dimnames(arr)[[3]]
    long.df = do.call(rbind, lapply(regions, function(r) {
        b = arr[, beta.col, r]
        b = b[!is.na(b)]
        if (length(b) == 0) return(NULL)
        data.frame(region = r, beta = b, stringsAsFactors = FALSE)
    }))
    if (is.null(long.df) || nrow(long.df) == 0) {
        message('  ', ct, ': no data, skipping')
        next
    }

    n.total = nrow(long.df)
    true.min = min(long.df$beta); true.max = max(long.df$beta)
    window = if (!is.na(opt$clip)) opt$clip else as.numeric(quantile(abs(long.df$beta), probs = opt$percentile))
    n.beyond = sum(abs(long.df$beta) > window)
    message('  ', ct, ': ', n.total, ' values, TRUE range [', round(true.min,2), ', ', round(true.max,2),
            '], plotted window +/-', round(window,2), ' (', round(100*opt$percentile,1), 'th pct), ',
            n.beyond, ' points (', round(100*n.beyond/n.total,2), '%) outside plotted window but still real/present in the data')

    long.df$region = factor(long.df$region, levels = region.levels[region.levels %in% unique(long.df$region)])
    these.colors = region.colors[levels(long.df$region)]

    p = ggplot(long.df, aes(beta, fill = region)) +
        geom_histogram(bins = 80, alpha = 0.8) +
        geom_vline(xintercept = 0, color = 'black', linewidth = 0.3, linetype = 2) +
        facet_wrap(~region, scales = 'free_y') +
        scale_fill_manual(values = these.colors) +
        coord_cartesian(xlim = c(-window, window)) +
        theme_classic(base_size = 12) +
        theme(legend.position = 'none', strip.background = element_blank()) +
        xlab(expression(italic(beta) ~ '(raw DGLM dispersion coefficient)')) +
        ylab('Number of genes') +
        ggtitle(paste0(ct, ' — raw beta distribution by region'),
                subtitle = paste0('n=', n.total, '; true range [', round(true.min,2), ', ', round(true.max,2),
                                   ']; showing ', round(100*opt$percentile,1), 'th pct window (+/-', round(window,2), '), ',
                                   n.beyond, ' pts (', round(100*n.beyond/n.total,2), '%) outside plotted window'))

    out.file = file.path(opt$figdir, paste0('beta_distribution_', ct, '.', opt$outfmt))
    ggsave(p, file = out.file, width = 12, height = 8, dpi = 150)
    message('  Saved: ', out.file)
}

message('done.')
