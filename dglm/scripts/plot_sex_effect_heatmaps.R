#!/usr/bin/env Rscript
#
# plot_sex_effect_heatmaps.R
#
# Heatmaps of the SEX term from the age+sex dispersion model
# (master_dglm_sex_disp_combined.tsv -- never plotted before this script).
# Sex is coded as.factor(sex) with no relevel() in dglm_model.R; raw values
# are literally "F"/"M", so R's default alphabetical factor ordering makes
# F the reference level. The sex coefficient therefore represents M
# relative to F: positive beta = higher dispersion in males (male-biased),
# negative beta = higher dispersion in females (female-biased). Colorbars
# below are labeled accordingly, not as generic "increase/decrease".
#
# Same opc/oligodendrocytes relabel as the other plotting scripts (cosmetic
# label fix only, using the already-correct parent directory).
#
# Usage:
#   Rscript plot_sex_effect_heatmaps.R --metric n_net --per_celltype

suppressMessages({
  library(optparse)
  library(ggplot2)
  library(reshape2)
})

option_list = list(
    make_option('--base_dir',     type='character',
                default='/scratch/easmit31/dispersion/dglm/disp_age__raw_counts'),
    make_option('--cutoff',       type='character', default='0.5'),
    make_option('--figdir',       type='character',
                default='/scratch/easmit31/dispersion/dglm/figures'),
    make_option('--outfmt',       type='character', default='png'),
    make_option('--metric',       type='character', default='pct_sig',
                help='pct_sig | n_male_biased | n_female_biased | n_net | mean_beta_sig | mean_beta_mag_sig | genes_tested | genes_significant'),
    make_option('--sig_col',      type='character', default='qvalue'),
    make_option('--qthresh',      type='double',    default=0.05),
    make_option('--lfsr_thresh',  type='double',    default=0.2),
    make_option('--per_celltype', action='store_true', default=FALSE)
)
opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

ckpt_suffix = paste0('_cutoff', opt$cutoff, '_agesex')
tsv_name = 'master_dglm_sex_disp_combined.tsv'

cell_type_dirs = list.dirs(opt$base_dir, recursive=FALSE)
all_rows = list()
for (d in cell_type_dirs) {
    tsv = file.path(d, paste0('dglm_checkpoints', ckpt_suffix), tsv_name)
    if (!file.exists(tsv)) next
    df = read.table(tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    df$true_parent_cell_type = basename(d)
    all_rows[[basename(d)]] = df
}
message('Found sex-term master TSVs for ', length(all_rows), ' cell type(s): ', paste(names(all_rows), collapse=', '))
if (length(all_rows) == 0) stop('No master_dglm_sex_disp_combined.tsv found under ', opt$base_dir)

master = do.call(rbind, all_rows)
rownames(master) = NULL
master = master[!is.na(master$beta),]
master$parent_cell_type = master$true_parent_cell_type

opc_olig_idx = grepl('^opc-olig_', master$cell_type)
if (any(opc_olig_idx)) {
    suffix = sub('^opc-olig_', '', master$cell_type[opc_olig_idx])
    master$cell_type[opc_olig_idx] = paste0(master$true_parent_cell_type[opc_olig_idx], '_', suffix)
}

master$sig = if (opt$sig_col %in% colnames(master)) master[[opt$sig_col]] < opt$qthresh else NA
message('Total tested rows: ', nrow(master))

metric_fn = switch(opt$metric,
    genes_tested       = function(sub) nrow(sub),
    genes_significant  = function(sub) sum(sub$sig, na.rm=TRUE),
    pct_sig             = function(sub) { l = sub$mash_lfsr; 100 * sum(!is.na(l) & l < opt$lfsr_thresh) / nrow(sub) },
    n_male_biased        = function(sub) { l = sub$mash_lfsr; b = sub$mash_beta; sum(!is.na(l) & l < opt$lfsr_thresh & b > 0, na.rm=TRUE) },
    n_female_biased        = function(sub) { l = sub$mash_lfsr; b = sub$mash_beta; sum(!is.na(l) & l < opt$lfsr_thresh & b <= 0, na.rm=TRUE) },
    n_net                    = function(sub) { l = sub$mash_lfsr; b = sub$mash_beta; sum(!is.na(l) & l < opt$lfsr_thresh & b > 0, na.rm=TRUE) - sum(!is.na(l) & l < opt$lfsr_thresh & b <= 0, na.rm=TRUE) },
    mean_beta_sig             = function(sub) { l = sub$mash_lfsr; s = sub[!is.na(l) & l < opt$lfsr_thresh,]; if (nrow(s)==0) 0 else mean(s$mash_beta, na.rm=TRUE) },
    mean_beta_mag_sig          = function(sub) { l = sub$mash_lfsr; s = sub[!is.na(l) & l < opt$lfsr_thresh,]; if (nrow(s)==0) 0 else mean(abs(s$mash_beta), na.rm=TRUE) },
    stop('Unknown --metric: ', opt$metric)
)
is_diverging  = opt$metric %in% c('n_net', 'mean_beta_sig')
metric_label  = c(genes_tested='Genes tested', genes_significant=paste0('Genes significant (', opt$sig_col, '<', opt$qthresh, ')'),
                   pct_sig=paste0('% significant sex effect (lfsr<', opt$lfsr_thresh, ')'),
                   n_male_biased='# genes, male-biased', n_female_biased='# genes, female-biased',
                   n_net='Net (male-biased \u2212 female-biased)',
                   mean_beta_sig='Mean sex beta, signed (+ = male-biased, sig genes)',
                   mean_beta_mag_sig='Mean |sex beta| (sig genes)')[[opt$metric]]

natural_sort = function(x) {
    prefix = sub('_[0-9]+$', '', x)
    suffix = suppressWarnings(as.numeric(sub('.*_([0-9]+)$', '\\1', x)))
    x[order(prefix, ifelse(is.na(suffix), 0, suffix))]
}

build_matrix = function(df, unit_col) {
    units   = natural_sort(unique(df[[unit_col]]))
    regions = sort(unique(df$region))
    m = matrix(NA_real_, nrow=length(units), ncol=length(regions), dimnames=list(units, regions))
    for (u in units) for (r in regions) {
        sub = df[df[[unit_col]] == u & df$region == r,]
        if (nrow(sub) == 0) next
        m[u, r] = metric_fn(sub)
    }
    m
}

save_heatmap = function(mat, title, fname) {
    df = melt(mat, varnames=c('unit','region'), value.name='n')
    df$unit   = factor(df$unit, levels=rev(natural_sort(unique(df$unit))))
    df$region = factor(df$region, levels=sort(unique(df$region)))

    p = ggplot(df, aes(region, unit, fill=n)) +
        geom_tile(color='white', linewidth=0.4) +
        geom_text(aes(label=ifelse(is.na(n), '', round(n, 1))), size=2.6, color='black') +
        theme_classic(base_size=13) +
        theme(axis.text.x=element_text(angle=45, hjust=1), axis.title=element_blank(),
              plot.title=element_text(hjust=0.5, size=10)) +
        ggtitle(paste(strwrap(title, width=70), collapse='\n'))

    if (is_diverging) {
        bound = max(abs(mat), na.rm=TRUE); if (!is.finite(bound) || bound == 0) bound = 1
        p = p + scale_fill_gradient2(low='#2166ac', mid='white', high='#b2182b', midpoint=0,
                                      limits=c(-bound, bound), na.value='gray90',
                                      name=paste0(metric_label, '\n\u2190 female-biased   male-biased \u2192'))
    } else if (opt$metric == 'n_female_biased') {
        vmax = max(mat, na.rm=TRUE); if (!is.finite(vmax) || vmax == 0) vmax = 1
        p = p + scale_fill_gradient(low='white', high='#2166ac', na.value='gray90', limits=c(0, vmax), name=metric_label)
    } else if (opt$metric == 'n_male_biased') {
        vmax = max(mat, na.rm=TRUE); if (!is.finite(vmax) || vmax == 0) vmax = 1
        p = p + scale_fill_gradient(low='white', high='#b2182b', na.value='gray90', limits=c(0, vmax), name=metric_label)
    } else {
        vmax = max(mat, na.rm=TRUE); if (!is.finite(vmax) || vmax == 0) vmax = 1
        p = p + scale_fill_gradient(low='white', high='#2166ac', na.value='gray90', limits=c(0, vmax), name=metric_label)
    }

    n_units = length(unique(df$unit))
    out = file.path(opt$figdir, paste0(fname, '.', opt$outfmt))
    ggsave(p, file=out, height=0.35*n_units+2, width=0.55*length(unique(df$region))+3.5, dpi=150, limitsize=FALSE)
    message('Saved: ', out)
}

label_suffix = paste0(' -- sex effect (from dispersion~age+sex), cutoff', opt$cutoff)

mat_sub = build_matrix(master, 'cell_type')
save_heatmap(mat_sub, paste0(metric_label, ', subcluster x region', label_suffix),
             paste0('heatmap_sexeffect_', opt$metric, '_subcluster'))

mat_ct = build_matrix(master, 'parent_cell_type')
save_heatmap(mat_ct, paste0(metric_label, ', cell type x region', label_suffix),
             paste0('heatmap_sexeffect_', opt$metric, '_celltype'))

if (opt$per_celltype) {
    for (ct in natural_sort(unique(master$parent_cell_type))) {
        sub_master = master[master$parent_cell_type == ct,]
        mat_one = build_matrix(sub_master, 'cell_type')
        save_heatmap(mat_one, paste0(metric_label, ', ', ct, ' subclusters x region', label_suffix),
                     paste0('heatmap_sexeffect_', opt$metric, '_subcluster_', ct))
    }
}

message('done.')
