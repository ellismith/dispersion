#!/usr/bin/env Rscript
# dglm_plot_heatmap.R
# Summary heatmaps of sig genes across cell types x regions.
# Can read from master TSV (global FDR) or per-cell-type RDS (per-region FDR).
#
# Usage:
#   # global FDR (recommended)
#   Rscript dglm_plot_heatmap.R --master_tsv checkpoints/master_dglm_globalfdr.tsv
#   # per-region FDR
#   Rscript dglm_plot_heatmap.R

source('/scratch/easmit31/variability/dglm/scripts/_include_options.R')

library(optparse)
library(ggplot2)
library(reshape2)

option_list = list(
    make_option('--checkpoints', type='character', default='/scratch/easmit31/variability/dglm/checkpoints'),
    make_option('--master_tsv',  type='character', default=NULL,
                help='optional: master TSV with global FDR qvalues'),
    make_option('--figdir',      type='character', default='/scratch/easmit31/variability/dglm/figures'),
    make_option('--outfmt',      type='character', default='png'),
    make_option('--sig_col', type='character', default='qvalue', help='column to use for significance: qvalue or mash_lfsr'),
    make_option('--qthresh',     type='double',    default=0.05)
)
opt = parse_args(OptionParser(option_list=option_list))

dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

save_fig = function(p, name, height, width) {
    out = file.path(opt$figdir, paste0(name, '.', opt$outfmt))
    if (opt$outfmt == 'png') {
        ggsave(p, file=out, height=height, width=width, dpi=150)
    } else {
        ggsave(p, file=out, height=height, width=width, useDingbats=FALSE)
    }
    message('Saved: ', out)
}

# ── load results ──────────────────────────────────────────────────────────
if (!is.null(opt$master_tsv)) {
    message('Loading from master TSV (global FDR): ', opt$master_tsv)
    master = read.table(opt$master_tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    fdr_label = 'global FDR'
} else {
    message('Loading from per-cell-type RDS (per-region FDR)')
    all.results = list()
    for (ct in cell.type.levels) {
        rds.file = file.path(opt$checkpoints, paste0(ct, '_dglm_mashr_results.rds'))
        if (!file.exists(rds.file)) next
        obj          = readRDS(rds.file)
        dglm.results = obj$dglm_results
        regions      = obj$regions
        beta.mat     = dglm.results[,'beta', regions, drop=FALSE][,1,, drop=FALSE]
        extreme      = apply(beta.mat, 1, function(x) any(abs(x) > 100, na.rm=TRUE))
        if (sum(extreme) > 0) dglm.results = dglm.results[!extreme,,, drop=FALSE]
        for (r in regions) {
            pval = dglm.results[,'pval', r]
            beta = dglm.results[,'beta', r]
            qval = p.adjust(pval, method='fdr')
            all.results[[length(all.results)+1]] = data.frame(
                symbol=dimnames(dglm.results)[[1]], cell_type=ct, region=r,
                beta=beta, qvalue=qval, stringsAsFactors=FALSE)
        }
    }
    master    = do.call(rbind, all.results)
    fdr_label = 'per-region FDR'
}

master = master[!is.na(master$qvalue),]
message('Total tests: ', nrow(master))
message('Sig at q<', opt$qthresh, ': ', sum(master$qvalue < opt$qthresh))

# ── build count matrices ──────────────────────────────────────────────────
inc.mat = matrix(NA, nrow=length(cell.type.levels), ncol=length(region.levels),
                 dimnames=list(cell.type.levels, region.levels))
dec.mat = inc.mat

for (ct in cell.type.levels) {
    for (r in region.levels) {
        sub = master[master$cell_type == ct & master$region == r,]
        if (nrow(sub) == 0) next
        inc.mat[ct, r] = sum(sub[[opt$sig_col]] < opt$qthresh & sub$beta > 0, na.rm=TRUE)
        dec.mat[ct, r] = sum(sub[[opt$sig_col]] < opt$qthresh & sub$beta < 0, na.rm=TRUE)
    }
}

# ── plot helper ───────────────────────────────────────────────────────────
make_heatmap = function(mat, title, outpath, high_col) {
    cols  = region.levels[region.levels %in% colnames(mat)]
    mat   = mat[, cols, drop=FALSE]
    df    = melt(mat, varnames=c('cell_type','region'), value.name='n')
    df$cell_type = factor(df$cell_type, levels=rev(cell.type.levels))
    df$region    = factor(df$region,    levels=region.levels)
    vmax  = max(df$n, na.rm=TRUE)

    p = ggplot(df, aes(region, cell_type, fill=n)) +
        geom_tile(color='white', linewidth=0.5) +
        geom_text(aes(label=ifelse(is.na(n), '', as.integer(n))),
                  size=3, color='black') +
        scale_fill_gradient(low='white', high=high_col, na.value='gray90',
                            limits=c(0, vmax), name='n genes') +
        theme_classic(base_size=12) +
        theme(axis.text.x  = element_text(angle=45, hjust=1),
              axis.title    = element_blank(),
              panel.border  = element_rect(fill=NA, color='gray80')) +
        ggtitle(title)
    ggsave(p, file=outpath, height=0.4*length(cell.type.levels)+2,
           width=0.5*length(region.levels)+3, dpi=150)
    message('Saved: ', outpath)
}

suffix = ifelse(!is.null(opt$master_tsv), '_globalfdr', '_perregionfdr')

make_heatmap(inc.mat,
    title   = paste0('DGLM: variance increases with age (q<', opt$qthresh, ', ', fdr_label, ')'),
    outpath = file.path(opt$figdir, paste0('heatmap_dglm_increase', suffix, '.', opt$outfmt)),
    high_col = '#d73027')

make_heatmap(dec.mat,
    title   = paste0('DGLM: variance decreases with age (q<', opt$qthresh, ', ', fdr_label, ')'),
    outpath = file.path(opt$figdir, paste0('heatmap_dglm_decrease', suffix, '.', opt$outfmt)),
    high_col = '#4575b4')

message('done.')
