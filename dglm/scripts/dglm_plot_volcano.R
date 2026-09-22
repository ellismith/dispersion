#!/usr/bin/env Rscript
# dglm_plot_volcano.R
#
# Volcano plot of dispersion effect size (mash_beta) vs significance
# (-log10(mash_lfsr)) for one cell type (or louvain subcluster -- anything
# that's a valid --cell_type value in the master TSV), faceted by region.
# Significant genes (lfsr<qthresh) colored by direction, top N by |beta|
# labeled with gene symbol.
#
# Works identically at cell-type or louvain resolution -- just point
# --master_tsv at the relevant checkpoints dir's master_dglm_combined.tsv.
#
# Usage:
#   Rscript dglm_plot_volcano.R --cell_type microglia \
#       --master_tsv dglm/disp_age__mean_full/master_dglm_combined.tsv \
#       --figdir dglm/disp_age__mean_full/figures
#   Rscript dglm_plot_volcano.R --cell_type microglia_5 \
#       --master_tsv dglm/disp_age__mean_full_louvain/master_dglm_combined.tsv \
#       --figdir dglm/disp_age__mean_full_louvain/figures

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(ggplot2)
have.ggrepel = requireNamespace('ggrepel', quietly=TRUE)
if (!have.ggrepel) message('ggrepel not installed -- falling back to plain (non-repelled) labels')

option_list = list(
    make_option('--cell_type',  type='character', help='cell type or subcluster label to plot'),
    make_option('--master_tsv', type='character', help='master_dglm_combined.tsv path'),
    make_option('--figdir',     type='character', help='output directory for figures'),
    make_option('--outfmt',     type='character', default='png'),
    make_option('--sig_col',    type='character', default='mash_lfsr'),
    make_option('--qthresh',    type='double',    default=0.2),
    make_option('--n_label',    type='integer',   default=15,
                help='label the top N significant genes by |beta| per facet'),
    make_option('--facet_ncol', type='integer',   default=3,
                help='number of facet columns -- lower means fewer, bigger panels per row (was an unset default that let ggplot pick ~4 columns regardless of facet count)'),
    make_option('--base_size',  type='double',    default=16,
                help='base font size (axis titles/text, legend, strip labels, title all scale off this) -- bumped up from the original default of 11'),
    make_option('--label_size', type='double',    default=4,
                help='gene-symbol label size, bumped up from the original 2.5'),
    make_option('--lfsr_floor', type='double',    default=1e-20,
                help='floor for mash_lfsr before -log10 (was hardcoded to Rs numerical limit 1e-300 -- a single near-zero-lfsr gene could blow the y-axis out to 300, squashing everything else)')
)
opt = parse_args(OptionParser(option_list=option_list))

dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

message('Loading: ', opt$master_tsv)
df = read.table(opt$master_tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
df = df[df$cell_type == opt$cell_type,]
if (nrow(df) == 0) stop('No rows for cell_type=', opt$cell_type, ' in ', opt$master_tsv)

beta.col = if ('mash_beta' %in% colnames(df)) 'mash_beta' else 'beta'
df$plot_beta = df[[beta.col]]
df = df[!is.na(df$plot_beta) & !is.na(df[[opt$sig_col]]) & abs(df$plot_beta) <= 100,]

# floor at --lfsr_floor (default 1e-20) rather than R's numerical limit
# (1e-300) -- a handful of genes can land at effectively-zero lfsr, and an
# uncapped floor turns -log10 into 300 for those, which (combined with a
# shared y-axis) drags every panel's scale out to 300 and squashes every
# other, still-meaningful point down near zero. 1e-20 is already far past
# "unambiguously significant" for a plot at this lfsr<0.2 threshold.
df$neglog10 = -log10(pmax(df[[opt$sig_col]], opt$lfsr_floor))
df$sig      = df[[opt$sig_col]] < opt$qthresh
df$direction = ifelse(df$sig & df$plot_beta > 0, 'Increase',
                ifelse(df$sig & df$plot_beta < 0, 'Decrease', 'Not significant'))
df$direction = factor(df$direction, levels=c('Increase','Decrease','Not significant'))
df$region    = factor(df$region, levels=region.levels[region.levels %in% unique(df$region)])

message('Total genes plotted: ', nrow(df))
message('Significant (', opt$sig_col, '<', opt$qthresh, '): ', sum(df$sig, na.rm=TRUE))

# top N by |beta| among significant, per region, for labeling
label.df = do.call(rbind, lapply(split(df, df$region), function(x) {
    x.sig = x[x$sig,]
    if (nrow(x.sig) == 0) return(NULL)
    x.sig[order(-abs(x.sig$plot_beta)),][1:min(opt$n_label, nrow(x.sig)),]
}))

p = ggplot(df, aes(plot_beta, neglog10, color=direction)) +
    geom_point(data=subset(df, !sig), size=1.2, alpha=0.4) +
    geom_point(data=subset(df, sig),  size=2.2, alpha=0.85) +
    geom_hline(yintercept=-log10(opt$qthresh), linetype='dashed', color='gray50', linewidth=0.3) +
    scale_color_manual(values=c('Increase'='#d73027', 'Decrease'='#4575b4', 'Not significant'='gray70')) +
    {if (have.ggrepel) ggrepel::geom_text_repel(data=label.df, aes(label=symbol), size=opt$label_size,
                              max.overlaps=Inf, box.padding=0.5, point.padding=0.3,
                              force=2, max.iter=8000, max.time=1.5, min.segment.length=0,
                              segment.size=0.25, color='black', seed=1)
     else geom_text(data=label.df, aes(label=symbol), size=opt$label_size, color='black', vjust=-0.6, check_overlap=TRUE)} +
    # extra x-margin so repelled labels near panel edges (e.g. long
    # ENSMMUG... IDs at the extremes) have room to sit inside the panel
    # instead of getting clipped -- this is what was cutting labels off at
    # the left/right edges before
    scale_x_continuous(expand=expansion(mult=0.18)) +
    facet_wrap(~region, scales='free', ncol=opt$facet_ncol) +
    xlab(expression(paste(italic(beta), ' (dispersion effect size)'))) +
    ylab(bquote(-log[10]~.(opt$sig_col))) +
    theme_classic(base_size=opt$base_size) +
    theme(legend.position='bottom', strip.background=element_blank(),
          strip.text=element_text(size=opt$base_size*0.9, face='bold'),
          plot.title=element_text(size=opt$base_size*1.2, face='bold')) +
    ggtitle(paste0(opt$cell_type, ' — dispersion volcano (', opt$sig_col, '<', opt$qthresh, ')'))

n.facets = length(unique(df$region))
ncol.actual = min(opt$facet_ncol, n.facets)
nrow.actual = ceiling(n.facets / ncol.actual)
out = file.path(opt$figdir, paste0('volcano_', opt$cell_type, '_', opt$sig_col, '_q', opt$qthresh, '.', opt$outfmt))
# Sized per-panel (5.5in wide x 4.5in tall) times actual grid dims, instead
# of the old min(16, 4*n.facets+2) formula -- that cap meant total width
# stayed ~constant no matter how many facets there were, so more regions
# just meant each panel got squeezed smaller. This scales properly: more
# facets -> bigger saved image, not more crowding.
ggsave(p, file=out, height=4.5*nrow.actual+2.5, width=5.5*ncol.actual+1.5, dpi=150, limitsize=FALSE)
message('Saved: ', out, ' (', ncol.actual, ' cols x ', nrow.actual, ' rows)')
message('done.')
