#!/usr/bin/env Rscript
# dglm_plot_summary_scatter.R
#
# One point per unit (cell type, or louvain subcluster -- whatever's in the
# master TSV's cell_type column) summarizing dispersion-with-age changes:
#   x = number of DISTINCT GENES significant in >= min_region_frac of that
#       unit's own tested regions (default 1/3 -- matches Chiou et al.'s
#       actual criterion: dglm_visualize.R's mash.sig requires LFSR<cutoff
#       AND consistent direction in >=5 of their 15 regions, i.e. 1/3; same
#       fraction as their topGO script's fraction.shared.cutoff override).
#       Each gene counted ONCE regardless of how many regions it hits --
#       NOT a raw sum of gene x region hits (that number conflates "many
#       genes each hitting once" with "few genes each hitting everywhere").
#       Uses each unit's OWN tested-region count as the denominator (not a
#       fixed 15 or 11), since cell types/subclusters vary in how many
#       regions actually have usable data.
#   y = mean |mash_beta| among that gene's significant regions, averaged
#       across those genes (one number per gene first, then averaged --
#       not averaged across raw hits)
#   color = % of those genes that are net-increasing (red) vs net-decreasing
#   label = unit name
#
# Answers "which units have the most, and biggest, dispersion changes" at a
# glance, across cell types or across louvain subclusters, without opening
# 12+ individual volcano plots.
#
# Usage:
#   Rscript dglm_plot_summary_scatter.R \
#       --master_tsv dglm/disp_age__mean_full/master_dglm_combined.tsv \
#       --figdir dglm/disp_age__mean_full/figures --outname summary_celltype
#   Rscript dglm_plot_summary_scatter.R \
#       --master_tsv dglm/disp_age__mean_full_louvain/master_dglm_combined.tsv \
#       --figdir dglm/disp_age__mean_full_louvain/figures --outname summary_louvain

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(ggplot2)
have.ggrepel = requireNamespace('ggrepel', quietly=TRUE)
if (!have.ggrepel) message('ggrepel not installed -- falling back to plain (non-repelled) labels')

option_list = list(
    make_option('--master_tsv',       type='character', help='master_dglm_combined.tsv path'),
    make_option('--figdir',           type='character', help='output directory for figures'),
    make_option('--outfmt',           type='character', default='png'),
    make_option('--outname',          type='character', default='summary_scatter'),
    make_option('--sig_col',          type='character', default='mash_lfsr'),
    make_option('--qthresh',          type='double',    default=0.2),
    make_option('--min_region_frac',  type='double',    default=1/3,
                help='a gene counts as significant for a unit only if it hits this fraction (or more) of that unit\'s own tested regions, consistently in one direction -- matches Chiou et al.\'s mash.sig criterion (>=5 of 15 regions = 1/3)'),
    make_option('--min_sig',          type='integer',   default=1,
                help='drop units with fewer than this many significant genes (avoids 0-count noise)'),
    make_option('--base_size',        type='double',    default=20,
                help='base font size for the whole plot (axis titles/text, legend, plot title all scale off this via theme_classic) -- bumped up from the original default of 12 for readability'),
    make_option('--point_size',       type='double',    default=5,
                help='point size, bumped up alongside base_size'),
    make_option('--label_size',       type='double',    default=6,
                help='ggrepel/geom_text label size, bumped up alongside base_size')
)
opt = parse_args(OptionParser(option_list=option_list))

dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

message('Loading: ', opt$master_tsv)
df = read.table(opt$master_tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)

beta.col = if ('mash_beta' %in% colnames(df)) 'mash_beta' else 'beta'
df$plot_beta = df[[beta.col]]
df = df[!is.na(df$plot_beta) & !is.na(df[[opt$sig_col]]) & abs(df$plot_beta) <= 100,]
df$sig = df[[opt$sig_col]] < opt$qthresh

summary.df = do.call(rbind, lapply(split(df, df$cell_type), function(x) {
    n.regions.tested = length(unique(x$region))
    min.regions.needed = max(1, ceiling(opt$min_region_frac * n.regions.tested))

    # per-gene: how many of ITS tested regions is it sig+increasing, or sig+decreasing
    gene.split = split(x, x$symbol)
    gene.summary = do.call(rbind, lapply(gene.split, function(g) {
        n.inc = sum(g$sig & g$plot_beta > 0)
        n.dec = sum(g$sig & g$plot_beta < 0)
        if (n.inc < min.regions.needed && n.dec < min.regions.needed) return(NULL)
        direction = if (n.inc >= min.regions.needed) 'Increase' else 'Decrease'
        hit.rows  = if (direction == 'Increase') g[g$sig & g$plot_beta > 0,] else g[g$sig & g$plot_beta < 0,]
        data.frame(symbol=unique(g$symbol), direction=direction,
                   mean_abs_beta=mean(abs(hit.rows$plot_beta)), stringsAsFactors=FALSE)
    }))

    n_sig = if (is.null(gene.summary)) 0 else nrow(gene.summary)
    if (n_sig == 0) return(NULL)

    data.frame(
        unit           = unique(x$cell_type),
        n_regions_tested = n.regions.tested,
        n_sig          = n_sig,
        mean_abs_beta  = mean(gene.summary$mean_abs_beta),
        pct_increase   = 100 * sum(gene.summary$direction == 'Increase') / n_sig,
        stringsAsFactors = FALSE
    )
}))
summary.df = summary.df[summary.df$n_sig >= opt$min_sig,]

message('Units plotted: ', nrow(summary.df))
message(paste(capture.output(print(summary.df[order(-summary.df$n_sig),])), collapse='\n'))

p = ggplot(summary.df, aes(n_sig, mean_abs_beta, color=pct_increase)) +
    geom_point(size=opt$point_size) +
    {if (have.ggrepel) ggrepel::geom_text_repel(aes(label=unit), size=opt$label_size, color='black',
                              max.overlaps=Inf, box.padding=0.6, point.padding=0.4,
                              force=3, force_pull=0.5, max.iter=10000, max.time=2,
                              min.segment.length=0, segment.size=0.3, seed=1)
     else geom_text(aes(label=unit), size=opt$label_size, color='black', vjust=-0.8, check_overlap=TRUE)} +
    scale_color_gradient2(low='#4575b4', mid='gray80', high='#d73027', midpoint=50,
                          name='% of sig genes\nthat increase', limits=c(0,100)) +
    xlab(paste0('Number of distinct genes sig. in >=', round(opt$min_region_frac*100), '% of tested regions\n(', opt$sig_col, '<', opt$qthresh, ')')) +
    ylab(expression(paste('Mean |', italic(beta), '| among significant genes'))) +
    theme_classic(base_size=opt$base_size) +
    theme(legend.text=element_text(size=opt$base_size*0.8),
          legend.title=element_text(size=opt$base_size*0.85),
          plot.title=element_text(size=opt$base_size*1.15, face='bold')) +
    ggtitle(paste0('Dispersion-with-age summary (', opt$sig_col, '<', opt$qthresh, ')'))

out = file.path(opt$figdir, paste0(opt$outname, '_', opt$sig_col, '_q', opt$qthresh, '.', opt$outfmt))
# Canvas bumped up from 7x9 -- more room is what actually lets ggrepel spread
# labels apart at a larger font size instead of cramming them into the same
# footprint as before
ggsave(p, file=out, height=10, width=13, dpi=150)
message('Saved: ', out)
message('done.')
