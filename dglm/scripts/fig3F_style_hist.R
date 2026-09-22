#!/usr/bin/env Rscript
#
# Fig 3F-style histogram: x-axis = number of regions a subcluster shows
# significance in, y-axis = number of subclusters. Colored by
# neuronal/non-neuronal. Optionally split by NET direction per region
# (majority of unique significant genes increasing vs decreasing) --
# NOT "any gene present", which is too permissive to distinguish directions
# when regions have thousands of significant genes in both directions.
#
# Usage (single cell type):
#   Rscript fig3F_style_hist.R --master <master_dglm_combined.tsv> \
#     [--split_direction] --out <output.png>
#
# Usage (all cell types combined):
#   Rscript fig3F_style_hist.R --all [--split_direction] --out <output.png>

suppressMessages({library(optparse); library(ggplot2); library(stringr)})

integer_breaks = function(n = 5) {
  function(x) unique(floor(pretty(x, n)))
}

option_list = list(
  make_option('--master', type='character', default=NULL),
  make_option('--all', action='store_true', default=FALSE),
  make_option('--split_direction', action='store_true', default=FALSE),
  make_option('--lfsr_thresh', type='double', default=0.05),
  make_option('--base_dir', type='character', default=NULL, help='override pipeline base dir for non-lCb cell types (lCb types get _noLcb_ version automatically)'),
  make_option('--min_genes', type='integer', default=1, help='minimum unique sig genes for a region to count'),
  make_option('--title', type='character', default=NULL),
  make_option('--fixed_y', action='store_true', default=FALSE, help='same y-axis for Increasing/Decreasing facets'),
  make_option('--dir_thresh', type='double', default=0.6, help='inc if p_inc>=t, dec if p_inc<=1-t; 0.5 = strict majority, ties dropped'),
  make_option('--out', type='character', default='fig3F_style.png')
)
opt = parse_args(OptionParser(option_list=option_list))
if (is.null(opt$master) && !opt$all) stop('Provide either --master <file> or --all')

neuronal_types = c('GABAergic_neurons','glutamatergic_neurons','cerebellar_neurons',
                    'medium_spiny_neurons','midbrain_neurons','basket_cells')

lcb_affected = c('opc','microglia','astrocytes','vascular_cells','oligodendrocytes')
if (!is.null(opt$base_dir)) {
  base_noLcb = opt$base_dir  # lCb already excluded at condition-filter time; same dir
  all_cts = c('GABAergic_neurons','glutamatergic_neurons','medium_spiny_neurons','midbrain_neurons','basket_cells','ependymal_cells','cerebellar_neurons','opc','microglia','astrocytes','vascular_cells','oligodendrocytes')
  cell_dirs = lapply(setNames(all_cts, all_cts),
    function(ct) if (ct %in% lcb_affected) base_noLcb else opt$base_dir)
} else {
  cell_dirs = list(
    GABAergic_neurons='disp_age__FINAL_autosomeX_300c_QCfiltered',
    glutamatergic_neurons='disp_age__FINAL_autosomeX_300c_QCfiltered',
    medium_spiny_neurons='disp_age__FINAL_autosomeX_300c_QCfiltered',
    midbrain_neurons='disp_age__FINAL_autosomeX_300c_QCfiltered',
    basket_cells='disp_age__FINAL_autosomeX_300c_QCfiltered',
    ependymal_cells='disp_age__FINAL_autosomeX_300c_QCfiltered',
    cerebellar_neurons='disp_age__FINAL_autosomeX_300c_QCfiltered',
    opc='disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered',
    microglia='disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered',
    astrocytes='disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered',
    vascular_cells='disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered',
    oligodendrocytes='disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered'
  )
}

# Simple version (no split): per subcluster, count DISTINCT regions with
# >=1 significant gene (unique by ensembl_id, though presence doesn't need dedup).
compute_simple = function(df, ct_label) {
  sig = df[df$mash_lfsr < opt$lfsr_thresh, ]
  if (nrow(sig) == 0) return(NULL)
  # per subcluster: count regions with >= min_genes unique sig genes
  per_sub_rows = list()
  for (subcl in unique(sig$cell_type)) {
    sub_sig = sig[sig$cell_type == subcl, ]
    regions_passing = sapply(unique(sub_sig$region), function(r) {
      length(unique(sub_sig$ensembl_id[sub_sig$region == r])) >= opt$min_genes
    })
    n_regions = sum(regions_passing)
    if (n_regions > 0) per_sub_rows[[subcl]] = data.frame(cell_type=subcl, region=n_regions)
  }
  per_sub = if(length(per_sub_rows)>0) do.call(rbind, per_sub_rows) else data.frame(cell_type=character(), region=integer())
  colnames(per_sub) = c("subcluster","n_regions")
  colnames(per_sub) = c('subcluster', 'n_regions')
  per_sub$parent_celltype = ct_label
  per_sub$class = ifelse(ct_label %in% neuronal_types, 'Neuronal', 'Non-neuronal')
  per_sub
}

# Net-direction version: per (subcluster, region), count UNIQUE genes
# increasing vs decreasing, region's net direction = whichever has more
# unique genes. Then per subcluster, count regions net-Increasing vs
# net-Decreasing.
compute_net_direction = function(df, ct_label) {
  sig = df[df$mash_lfsr < opt$lfsr_thresh, ]
  if (nrow(sig) == 0) return(NULL)
  sub_region = unique(sig[, c('cell_type','region')])
  results = list()
  for (i in seq_len(nrow(sub_region))) {
    s = sub_region$cell_type[i]; r = sub_region$region[i]
    rows = sig[sig$cell_type==s & sig$region==r, ]
    n_inc = length(unique(rows$ensembl_id[rows$mash_beta > 0]))
    n_dec = length(unique(rows$ensembl_id[rows$mash_beta <= 0]))
    if (max(n_inc, n_dec) < opt$min_genes) next  # region ignored below gene floor
    p_inc = n_inc / (n_inc + n_dec)
    if ((p_inc > 1 - opt$dir_thresh && p_inc < opt$dir_thresh) || p_inc == 0.5) next  # mixed (40-60%): no net claim, region dropped
    net_dir = if (p_inc >= opt$dir_thresh && p_inc > 0.5) 'Increasing' else 'Decreasing'
    results[[paste(s,r)]] = data.frame(subcluster=s, region=r, n_inc=n_inc, n_dec=n_dec, net_direction=net_dir)
  }
  net_df = do.call(rbind, results)
  per_sub = aggregate(region ~ subcluster + net_direction, data=net_df, FUN=length)
  colnames(per_sub) = c('subcluster', 'direction', 'n_regions')
  per_sub$parent_celltype = ct_label
  per_sub$class = ifelse(ct_label %in% neuronal_types, 'Neuronal', 'Non-neuronal')
  per_sub
}

all_data = list()
if (opt$all) {
  for (ct in names(cell_dirs)) {
    bd = cell_dirs[[ct]]
    f = if (startsWith(bd, '/')) file.path(bd, ct, 'dglm_checkpoints_cutoff0.5', 'master_dglm_combined.tsv') else file.path('/scratch/easmit31/dispersion/dglm', bd, ct, 'dglm_checkpoints_cutoff0.5', 'master_dglm_combined.tsv')
    if (!file.exists(f)) { cat(ct, ': MISSING, skipping\n'); next }
    df = read.csv(f, sep='\t')
    all_data[[ct]] = if (opt$split_direction) compute_net_direction(df, ct) else compute_simple(df, ct)
  }
} else {
  df = read.csv(opt$master, sep='\t')
  ct_label = sub('_[0-9]+$', '', df$cell_type[1])
  all_data[['single']] = if (opt$split_direction) compute_net_direction(df, ct_label) else compute_simple(df, ct_label)
}

per_sub = do.call(rbind, all_data)
cat('Total subcluster', if(opt$split_direction) '-direction' else '', 'rows summarized:', nrow(per_sub), '\n')

title = if (!is.null(opt$title)) opt$title else if (opt$all) 'All cell types: subclusters by number of regions showing significant effects' else paste(per_sub$parent_celltype[1], ': subclusters by number of regions showing significant effects')

neuronal_present = intersect(neuronal_types, unique(per_sub$parent_celltype))
nonneuronal_present = setdiff(unique(per_sub$parent_celltype), neuronal_types)
subtitle = paste0('Neuronal: ', paste(neuronal_present, collapse=', '),
                   '  |  Non-neuronal: ', paste(nonneuronal_present, collapse=', '))

p = ggplot(per_sub, aes(x=n_regions, fill=class)) +
  geom_histogram(binwidth=1, color='white', position='stack') +
  scale_fill_manual(values=c('Neuronal'='#E07B39','Non-neuronal'='#3B7DD8'), name=NULL) +
  scale_x_continuous(breaks=1:max(per_sub$n_regions)) +
  theme_minimal(base_size=24) +
  theme(legend.position='top') +
  labs(x=if(opt$split_direction) 'Number of regions with >=60% of significant genes in this direction' else paste0('Number of regions with >=', opt$min_genes, ' significant gene(s)'),
       y='Number of subclusters', title=title, subtitle=subtitle)

if (opt$split_direction) {
  per_sub$direction = factor(per_sub$direction, levels=c('Increasing','Decreasing'))
  p = p %+% per_sub + facet_wrap(~direction, ncol=1, scales = if (opt$fixed_y) 'fixed' else 'free_y')
}

p = p + labs(title=stringr::str_wrap(title, 50), subtitle=stringr::str_wrap(subtitle, 70)) + scale_y_continuous(breaks=integer_breaks())
h = if (opt$split_direction) 13 else 8
dirlab = if (opt$dir_thresh == 0.5) '>50%' else paste0('>=', round(100*opt$dir_thresh), '%')
if (!is.null(p$labels$x)) p$labels$x = gsub('>=60%', dirlab, p$labels$x, fixed=TRUE)
ggsave(opt$out, p, width=15, height=h, dpi=150)
cat('Saved:', opt$out, '\n')
# note: min_genes option added below -- need to add to option_list
