#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(ggplot2)
  library(reshape2)
})

option_list = list(
    make_option('--base_dir',   type='character',
                default='/scratch/easmit31/dispersion/dglm/disp_age__raw_counts'),
    make_option('--cutoff',     type='character', default='0.5'),
    make_option('--dispersion', type='character', default='age'),
    make_option('--sig_col',    type='character', default='mash_lfsr'),
    make_option('--qthresh',    type='double', default=0.2),
    make_option('--figdir',     type='character',
                default='/scratch/easmit31/dispersion/dglm/figures'),
    make_option('--target_celltype', type='character', default=NULL,
                help='restrict to one cell type only (isolates it from the multi-cell-type figure)')
)
opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

ckpt_suffix = if (opt$dispersion == 'age_sex') paste0('_cutoff', opt$cutoff, '_agesex') else paste0('_cutoff', opt$cutoff)
tsv_name    = if (opt$dispersion == 'age_sex') 'master_dglm_age_disp_combined.tsv' else 'master_dglm_combined.tsv'

cell_type_dirs = list.dirs(opt$base_dir, recursive=FALSE)
if (!is.null(opt$target_celltype)) {
    cell_type_dirs = cell_type_dirs[basename(cell_type_dirs) == opt$target_celltype]
}
all_rows = list()
for (d in cell_type_dirs) {
    tsv = file.path(d, paste0('dglm_checkpoints', ckpt_suffix), tsv_name)
    if (!file.exists(tsv)) next
    df = read.table(tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    df$true_parent_cell_type = basename(d)
    all_rows[[basename(d)]] = df
}
message('Found master TSVs for ', length(all_rows), ' cell type(s)')
if (length(all_rows) == 0) stop('No master TSVs found')

master = do.call(rbind, all_rows)
rownames(master) = NULL
master = master[!is.na(master[[opt$sig_col]]) & !is.na(master$beta),]
master$sig = master[[opt$sig_col]] < opt$qthresh
master$parent_cell_type = master$true_parent_cell_type

opc_olig_idx = grepl('^opc-olig_', master$cell_type)
if (any(opc_olig_idx)) {
    suffix = sub('^opc-olig_', '', master$cell_type[opc_olig_idx])
    master$cell_type[opc_olig_idx] = paste0(master$true_parent_cell_type[opc_olig_idx], '_', suffix)
}

label_suffix = paste0(' (dispersion~', opt$dispersion, ', cutoff', opt$cutoff, ', ', opt$sig_col, '<', opt$qthresh, ')')
beta_col = if (opt$sig_col == 'mash_lfsr') 'mash_beta' else 'beta'

natural_sort = function(x) {
    prefix = sub('_[0-9]+$', '', x)
    suffix = suppressWarnings(as.numeric(sub('.*_([0-9]+)$', '\\1', x)))
    x[order(prefix, ifelse(is.na(suffix), 0, suffix))]
}

save_bar = function(df, group_col, title, fname) {
    bar_df = do.call(rbind, lapply(natural_sort(unique(df[[group_col]])), function(g) {
        sub_outer = df[df[[group_col]] == g,]
        do.call(rbind, lapply(sort(unique(sub_outer$region)), function(r) {
            sub = sub_outer[sub_outer$region == r,]
            data.frame(group=g, region=r,
                       increase = sum(sub$sig & sub[[beta_col]] > 0),
                       decrease = -sum(sub$sig & sub[[beta_col]] <= 0))
        }))
    }))
    bar_long = melt(bar_df, id.vars=c('group','region'), variable.name='direction', value.name='n')
    bar_long$group = factor(bar_long$group, levels=natural_sort(unique(bar_long$group)))
    bar_long$region = factor(bar_long$region, levels=sort(unique(bar_long$region)))
    p = ggplot(bar_long, aes(region, n, fill=direction)) +
        geom_col() +
        facet_wrap(~group, scales='free') +
        scale_fill_manual(values=c(increase='#d9534f', decrease='#428bca')) +
        theme_classic(base_size=11) +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        ylab('n significant genes') + xlab(NULL) + ggtitle(title)
    n_groups = length(unique(bar_long$group))
    ncol_facet = min(4, n_groups)
    nrow_facet = ceiling(n_groups / ncol_facet)
    out = file.path(opt$figdir, fname)
    ggsave(p, file=out, height=3.5*nrow_facet+1, width=4*ncol_facet+1, dpi=150, limitsize=FALSE)
    message('Saved: ', out)
}

save_bar(master, 'parent_cell_type', paste0('Sig genes by region, per cell type', label_suffix),
          paste0('bar_sig_by_region_celltype_', opt$dispersion, '_', opt$sig_col, opt$qthresh, '.png'))
save_bar(master, 'cell_type', paste0('Sig genes by region, per subcluster', label_suffix),
          paste0('bar_sig_by_region_subcluster_', opt$dispersion, '_', opt$sig_col, opt$qthresh, '.png'))
message('done.')
