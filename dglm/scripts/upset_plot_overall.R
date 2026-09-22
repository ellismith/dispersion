#!/usr/bin/env Rscript
#
# Overall UpSet across ALL cell types: for each gene, is it significant
# (any direction, any condition) within each of the 12 broad cell types?
# Groups = cell types. Same bar+dot structure as upset_plot.R.
#
# Usage:
#   Rscript upset_plot_overall.R --top_n 30 --lfsr_thresh 0.2 --out <output.png>

suppressMessages({library(optparse); library(ggplot2); library(patchwork)})

option_list = list(
  make_option('--top_n', type='integer', default=30),
  make_option('--lfsr_thresh', type='double', default=0.2),
  make_option('--out', type='character', default='upset_overall.png')
)
opt = parse_args(OptionParser(option_list=option_list))

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

all_inc = list(); all_dec = list()
for (ct in names(cell_dirs)) {
  f = file.path('/scratch/easmit31/dispersion/dglm', cell_dirs[[ct]], ct, 'dglm_checkpoints_cutoff0.5', 'master_dglm_combined.tsv')
  if (!file.exists(f)) { cat(ct, ': MISSING, skipping\n'); next }
  df = read.csv(f, sep='\t')
  sig = df[df$mash_lfsr < opt$lfsr_thresh, ]
  inc_genes = unique(sig$ensembl_id[sig$mash_beta > 0])
  dec_genes = unique(sig$ensembl_id[sig$mash_beta <= 0])
  all_inc[[ct]] = inc_genes
  all_dec[[ct]] = dec_genes
}

build_membership_from_list = function(gene_lists) {
  all_genes = unique(unlist(gene_lists))
  membership = sapply(all_genes, function(g) {
    paste(sort(names(gene_lists)[sapply(gene_lists, function(v) g %in% v)]), collapse=',')
  })
  table(membership)
}

inc_combos = build_membership_from_list(all_inc)
dec_combos = build_membership_from_list(all_dec)

top_inc_names = names(sort(inc_combos, decreasing=TRUE))[1:min(opt$top_n, length(inc_combos))]
top_dec_names = names(sort(dec_combos, decreasing=TRUE))[1:min(opt$top_n, length(dec_combos))]
all_combo_names = union(top_inc_names, top_dec_names)

inc_n = setNames(as.numeric(inc_combos[all_combo_names]), all_combo_names); inc_n[is.na(inc_n)] = 0
dec_n = setNames(as.numeric(dec_combos[all_combo_names]), all_combo_names); dec_n[is.na(dec_n)] = 0

combo_order = all_combo_names[order(-(inc_n + dec_n))]
combo_order = combo_order[1:min(opt$top_n, length(combo_order))]

bar_df = data.frame(combo=factor(combo_order, levels=combo_order),
                     n_increasing=inc_n[combo_order], n_decreasing=dec_n[combo_order])

p_top = ggplot(bar_df, aes(x=combo, y=n_increasing)) +
  geom_col(fill='firebrick') +
  theme_minimal(base_size=9) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), panel.grid.major.x=element_blank()) +
  labs(y='Gene count\n(increasing)')

p_bottom = ggplot(bar_df, aes(x=combo, y=-n_decreasing)) +
  geom_col(fill='steelblue') +
  theme_minimal(base_size=9) +
  theme(axis.text.x=element_blank(), panel.grid.major.x=element_blank()) +
  labs(y='Gene count\n(decreasing)', x='Cell type combination')

all_celltypes = names(cell_dirs)
dot_rows = list()
for (cn in combo_order) {
  members = strsplit(cn, ',')[[1]]
  for (g in all_celltypes) dot_rows[[paste(cn,g)]] = data.frame(combo=cn, grp=g, is_member=g %in% members)
}
dot_data = do.call(rbind, dot_rows)
dot_data$combo = factor(dot_data$combo, levels=combo_order)
dot_data$grp = factor(dot_data$grp, levels=rev(all_celltypes))

p_mid = ggplot(dot_data, aes(x=combo, y=grp)) +
  geom_point(aes(alpha=is_member), color='black', size=2.5) +
  scale_alpha_manual(values=c('TRUE'=1, 'FALSE'=0.08), guide='none') +
  theme_minimal(base_size=9) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
        panel.grid.major=element_line(color='grey92')) +
  labs(y='Cell type')

combined = p_top / p_mid / p_bottom + plot_layout(heights=c(1.5, 2, 1.5)) +
  plot_annotation(title='Overall UpSet: significant genes shared across cell types',
                   subtitle=paste0('A gene counts as significant in a cell type if significant (any direction, any condition) anywhere within it. Top ', length(combo_order), ' combinations shown.'))

ggsave(opt$out, combined, width=15, height=9, dpi=150)
cat('Saved:', opt$out, '\n')
