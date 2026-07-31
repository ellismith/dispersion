#!/usr/bin/env Rscript
# dglm_plot_heatmap_sex.R
# Heatmap for sex dispersion effects.
# Red = increased variability in males, Blue = increased variability in females.
# No negative values — uses absolute effect size with direction encoded by color.

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(ggplot2)
library(reshape2)

option_list = list(
    make_option('--master_tsv', type='character', required=TRUE),
    make_option('--figdir',     type='character', default='/scratch/easmit31/dispersion/dglm/figures'),
    make_option('--outfmt',     type='character', default='png'),
    make_option('--qthresh',    type='double',    default=0.05)
)
opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

message('Loading: ', opt$master_tsv)
master = read.table(opt$master_tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
master = master[!is.na(master$qvalue) & !is.na(master$beta),]
master$sig = master$qvalue < opt$qthresh

message('Total tests: ', nrow(master))
message('Sig at q<', opt$qthresh, ': ', sum(master$sig))
message('Sig male>female (beta>0): ', sum(master$sig & master$beta > 0))
message('Sig female>male (beta<0): ', sum(master$sig & master$beta < 0))

# ── count matrices ────────────────────────────────────────────────────────
male.mat   = matrix(NA, nrow=length(cell.type.levels), ncol=length(region.levels),
                    dimnames=list(cell.type.levels, region.levels))
female.mat = male.mat

for (ct in cell.type.levels) {
    for (r in region.levels) {
        sub = master[master$cell_type == ct & master$region == r,]
        if (nrow(sub) == 0) next
        male.mat[ct, r]   = sum(sub$sig & sub$beta > 0, na.rm=TRUE)
        female.mat[ct, r] = sum(sub$sig & sub$beta < 0, na.rm=TRUE)
    }
}

save_heatmap = function(mat, title, fname, high_col) {
    df           = melt(mat, varnames=c('cell_type','region'), value.name='n')
    df$cell_type = factor(df$cell_type, levels=rev(cell.type.levels))
    df$region    = factor(df$region,    levels=region.levels)
    vmax = max(mat, na.rm=TRUE)
    if (vmax == 0) vmax = 1

    p = ggplot(df, aes(region, cell_type, fill=n)) +
        geom_tile(color='white', linewidth=0.5) +
        geom_text(aes(label=ifelse(is.na(n), '', as.integer(n))),
                  size=2.8, color='black') +
        scale_fill_gradient(low='white', high=high_col, na.value='gray90',
                            limits=c(0, vmax), name='n genes') +
        theme_classic(base_size=11) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              axis.title=element_blank()) +
        ggtitle(title)

    out = file.path(opt$figdir, paste0(fname, '.', opt$outfmt))
    ggsave(p, file=out,
           height=0.4*length(cell.type.levels)+2,
           width=0.5*length(region.levels)+3, dpi=150)
    message('Saved: ', out)
}

save_heatmap(male.mat,
    title = paste0('Increased dispersion in males (q<', opt$qthresh, ')'),
    fname = 'heatmap_sex_disp_male',
    high_col = '#d73027')

save_heatmap(female.mat,
    title = paste0('Increased dispersion in females (q<', opt$qthresh, ')'),
    fname = 'heatmap_sex_disp_female',
    high_col = '#4575b4')

message('done.')
