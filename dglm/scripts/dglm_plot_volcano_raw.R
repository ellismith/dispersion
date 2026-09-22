#!/usr/bin/env Rscript
# dglm_plot_volcano_raw.R
#
# EARLY-LOOK PLOT ONLY -- read this before trusting the significance axis.
#
# dglm()'s dispersion-submodel SE comes out IDENTICAL for every gene within
# a given cell_type x region (confirmed empirically: two totally unrelated
# synthetic response vectors fit against the same design matrix produce the
# exact same SE to full floating-point precision -- this is a property of
# how dglm computes dispersion-submodel SEs from the design alone, not a bug
# in our pipeline, and it's equally true of Chiou et al.'s own results since
# they use the identical package/method). Practical consequence: WITHIN one
# region, ranking genes by this raw q-value is mathematically equivalent to
# ranking them by |beta| alone -- it is not doing real per-gene significance
# discrimination, despite the name. Treat the y-axis here as "genes sorted
# by effect size," not "genes we're statistically confident about."
#
# mash_lfsr (see dglm_plot_volcano.R, post-mashr) does NOT have this problem
# -- mashr compares each gene's pattern ACROSS regions, where SE genuinely
# does vary (different animal counts per region), which is exactly what
# recovers real statistical discrimination. Use dglm_plot_volcano.R once
# mashr results exist; use this only for an early look before mashr
# finishes, and don't read significance claims out of it.
#
# Usage:
#   Rscript dglm_plot_volcano_raw.R --cell_type microglia_5 \
#       --checkpoints dglm/disp_age__mean_full_louvain \
#       --figdir dglm/disp_age__mean_full_louvain/figures

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(ggplot2)
have.ggrepel = requireNamespace('ggrepel', quietly=TRUE)
if (!have.ggrepel) message('ggrepel not installed -- falling back to plain (non-repelled) labels')

option_list = list(
    make_option('--cell_type',   type='character', help='cell type or subcluster label'),
    make_option('--checkpoints', type='character', help='dir containing {cell_type}_dglm_results.rds'),
    make_option('--figdir',      type='character', help='output directory for figures'),
    make_option('--outfmt',      type='character', default='png'),
    make_option('--qthresh',     type='double',    default=0.05,
                help='raw per-region BH-FDR threshold (0.05 is the usual plain-FDR convention, NOT the 0.2 LFSR one -- there is no mashr shrinkage here)'),
    make_option('--n_label',     type='integer',   default=15)
)
opt = parse_args(OptionParser(option_list=option_list))

dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

in.file = file.path(opt$checkpoints, paste0(opt$cell_type, '_dglm_results.rds'))
message('Loading: ', in.file)
obj = readRDS(in.file)
arr = obj$array
human.symbols = obj$human_symbols

regions = region.levels[region.levels %in% dimnames(arr)[[3]]]
regions = regions[apply(arr[,'beta',regions,drop=FALSE], 3, function(x) sum(!is.na(x) & x!=0) > 0)]
message('Regions with data: ', paste(regions, collapse=', '))

df = do.call(rbind, lapply(regions, function(r) {
    data.frame(
        ensembl_id = dimnames(arr)[[1]],
        symbol     = human.symbols[dimnames(arr)[[1]]],
        region     = r,
        beta       = arr[,'beta',r],
        qval       = arr[,'qval',r],
        stringsAsFactors = FALSE
    )
}))
df = df[!is.na(df$beta) & !is.na(df$qval) & abs(df$beta) <= 100,]
df$region = factor(df$region, levels=regions)

df$neglog10 = -log10(pmax(df$qval, 1e-300))
df$sig      = df$qval < opt$qthresh
df$direction = ifelse(df$sig & df$beta > 0, 'Increase',
                ifelse(df$sig & df$beta < 0, 'Decrease', 'Not significant'))
df$direction = factor(df$direction, levels=c('Increase','Decrease','Not significant'))

message('Total tests plotted: ', nrow(df))
message('Significant (raw q<', opt$qthresh, '): ', sum(df$sig, na.rm=TRUE))

label.df = do.call(rbind, lapply(split(df, df$region), function(x) {
    x.sig = x[x$sig,]
    if (nrow(x.sig) == 0) return(NULL)
    x.sig[order(-abs(x.sig$beta)),][1:min(opt$n_label, nrow(x.sig)),]
}))

p = ggplot(df, aes(beta, neglog10, color=direction)) +
    geom_point(data=subset(df, !sig), size=0.8, alpha=0.4) +
    geom_point(data=subset(df, sig),  size=1.3, alpha=0.8) +
    geom_hline(yintercept=-log10(opt$qthresh), linetype='dashed', color='gray50', linewidth=0.3) +
    scale_color_manual(values=c('Increase'='#d73027', 'Decrease'='#4575b4', 'Not significant'='gray70')) +
    {if (have.ggrepel) ggrepel::geom_text_repel(data=label.df, aes(label=symbol), size=2.5, max.overlaps=20,
                              segment.size=0.2, color='black')
     else geom_text(data=label.df, aes(label=symbol), size=2.5, color='black', vjust=-0.6, check_overlap=TRUE)} +
    facet_wrap(~region, scales='free_x') +
    xlab(expression(paste(italic(beta), ' (dispersion effect size, raw DGLM)'))) +
    ylab(expression(-log[10]~'(raw q -- effect-size proxy, see header)')) +
    theme_classic(base_size=11) +
    theme(legend.position='bottom', strip.background=element_blank()) +
    ggtitle(paste0(opt$cell_type, ' — genes ranked by effect size (pre-mashr, NOT a real significance test -- see header)'))

n.facets = length(regions)
out = file.path(opt$figdir, paste0('volcano_raw_', opt$cell_type, '_q', opt$qthresh, '.', opt$outfmt))
ggsave(p, file=out, height=3*ceiling(n.facets/4)+2, width=min(16, 4*n.facets+2), dpi=150)
message('Saved: ', out)
message('done.')
