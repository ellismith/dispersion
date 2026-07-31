#!/usr/bin/env Rscript
# dglm_plot_heatmap.R
# Summary heatmaps of sig genes across cell types x regions.
# Default: increase/decrease/net heatmaps for age effects.
# --sex_mode: male/female separate heatmaps (no negatives).

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(ggplot2)
library(reshape2)

option_list = list(
    make_option('--master_tsv', type='character', required=TRUE),
    make_option('--figdir',     type='character', default='/scratch/easmit31/dispersion/dglm/figures'),
    make_option('--outfmt',     type='character', default='png'),
    make_option('--sig_col',    type='character', default='qvalue'),
    make_option('--qthresh',    type='double',    default=0.05),
    make_option('--sex_mode',   action='store_true', default=FALSE)
)
opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

message('Loading: ', opt$master_tsv)
master = read.table(opt$master_tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
if (!opt$sig_col %in% colnames(master))
    stop('Column not found: ', opt$sig_col)
master = master[!is.na(master[[opt$sig_col]]) & !is.na(master$beta),]
master$sig = master[[opt$sig_col]] < opt$qthresh
fdr_label  = ifelse(opt$sig_col == 'mash_lfsr', 'mashr lfsr', 'global FDR')

message('Total tests: ', nrow(master))
message('Sig at ', opt$sig_col, '<', opt$qthresh, ': ', sum(master$sig))

# ── count matrices ────────────────────────────────────────────────────────────
mat_a = matrix(NA, nrow=length(cell.type.levels), ncol=length(region.levels),
               dimnames=list(cell.type.levels, region.levels))
mat_b = mat_a

for (ct in cell.type.levels) {
    for (r in region.levels) {
        sub = master[master$cell_type == ct & master$region == r,]
        if (nrow(sub) == 0) next
        mat_a[ct, r] = sum(sub$sig & sub$beta > 0, na.rm=TRUE)
        mat_b[ct, r] = sum(sub$sig & sub$beta < 0, na.rm=TRUE)
    }
}

# ── plot helper ───────────────────────────────────────────────────────────────
save_heatmap = function(mat, title, fname, high_col, diverging=FALSE) {
    df           = melt(mat, varnames=c('cell_type','region'), value.name='n')
    # apply acronym labels and canonical order
    df$cell_type = factor(cell.type.labels[as.character(df$cell_type)],
                          levels=rev(cell.type.labels[cell.type.levels]))
    df$region    = factor(df$region, levels=region.levels)
    vmax = max(abs(mat), na.rm=TRUE)
    if (vmax == 0) vmax = 1

    base_theme = theme_classic(base_size=14) +
        theme(
            axis.text.x  = element_text(angle=45, hjust=1, size=13),
            axis.text.y  = element_text(size=13),
            axis.title   = element_blank(),
            plot.title   = element_text(size=14, hjust=0.5)
        )

    if (diverging) {
        p = ggplot(df, aes(region, cell_type, fill=n)) +
            geom_tile(color='white', linewidth=0.5) +
            geom_text(aes(label=ifelse(is.na(n), '', as.integer(n))), size=3, color='black') +
            scale_fill_gradient2(low='#4575b4', mid='white', high='#d73027',
                                 midpoint=0, limits=c(-vmax, vmax),
                                 na.value='gray90', name='n genes') +
            base_theme + ggtitle(title)
    } else {
        p = ggplot(df, aes(region, cell_type, fill=n)) +
            geom_tile(color='white', linewidth=0.5) +
            geom_text(aes(label=ifelse(is.na(n), '', as.integer(n))), size=3, color='black') +
            scale_fill_gradient(low='white', high=high_col, na.value='gray90',
                                limits=c(0, vmax), name='n genes') +
            base_theme + ggtitle(title)
    }

    out = file.path(opt$figdir, paste0(fname, '_final.', opt$outfmt))
    ggsave(p, file=out,
           height=0.45*length(cell.type.levels)+2,
           width=0.55*length(region.levels)+3, dpi=150)
    message('Saved: ', out)
}

if (opt$sex_mode) {
    save_heatmap(mat_a,
        title    = paste0('Increased dispersion in males (', fdr_label, '<', opt$qthresh, ')'),
        fname    = paste0('heatmap_sex_dispersion_male_', opt$sig_col, '_q', opt$qthresh),
        high_col = '#d73027')
    save_heatmap(mat_b,
        title    = paste0('Increased dispersion in females (', fdr_label, '<', opt$qthresh, ')'),
        fname    = paste0('heatmap_sex_dispersion_female_', opt$sig_col, '_q', opt$qthresh),
        high_col = '#4575b4')
} else {
    save_heatmap(mat_a,
        title    = paste0('DGLM: variance increases with age (', fdr_label, '<', opt$qthresh, ')'),
        fname    = paste0('heatmap_age_dispersion_increase_', opt$sig_col, '_q', opt$qthresh),
        high_col = '#d73027')
    save_heatmap(mat_b,
        title    = paste0('DGLM: variance decreases with age (', fdr_label, '<', opt$qthresh, ')'),
        fname    = paste0('heatmap_age_dispersion_decrease_', opt$sig_col, '_q', opt$qthresh),
        high_col = '#4575b4')
    net.mat = mat_a - mat_b
    save_heatmap(net.mat,
        title     = paste0('DGLM: net direction (', fdr_label, '<', opt$qthresh, ')'),
        fname     = paste0('heatmap_age_dispersion_net_', opt$sig_col, '_q', opt$qthresh),
        high_col  = '#d73027', diverging=TRUE)
}

message('done.')
