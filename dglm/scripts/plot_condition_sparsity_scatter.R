#!/usr/bin/env Rscript
#
# plot_condition_sparsity_scatter.R
#
# Scatter of n_animals (x) vs n_cells (y, log10) per subcluster x region
# combination, colored by region, labeled with the louvain subcluster
# number. Dashed reference lines show candidate cell/animal thresholds so
# pass/fail is visible by which quadrant a point falls in.
#
# Usage:
#   Rscript plot_condition_sparsity_scatter.R --cell_type GABAergic_neurons

suppressMessages({
  library(optparse)
  library(ggplot2)
})

option_list = list(
    make_option('--cell_type', type='character', default=NULL),
    make_option('--sparsity_csv', type='character', default=NULL),
    make_option('--out_base', type='character',
                default='/scratch/easmit31/dispersion/dglm/disp_age__min300cells_min30animals'),
    make_option('--cell_thresholds', type='character', default='300,100',
                help='comma-separated horizontal reference lines (cells)'),
    make_option('--animal_threshold', type='double', default=30,
                help='vertical reference line (animals)'),
    make_option('--figdir', type='character',
                default='/scratch/easmit31/dispersion/dglm/disp_age__min300cells_min30animals/figures')
)
opt = parse_args(OptionParser(option_list=option_list))
if (is.null(opt$cell_type)) stop('--cell_type is required')
dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

sparsity_csv = if (!is.null(opt$sparsity_csv)) opt$sparsity_csv else
    file.path(opt$out_base, paste0(opt$cell_type, '_condition_sparsity.csv'))
if (!file.exists(sparsity_csv)) stop('Sparsity CSV not found: ', sparsity_csv)

df = read.csv(sparsity_csv, stringsAsFactors=FALSE)
df$louvain_num = sub(paste0('^', opt$cell_type, '_'), '', df$subcluster)
cell_lines = as.numeric(strsplit(opt$cell_thresholds, ',')[[1]])

p = ggplot(df, aes(x=n_animals, y=n_cells, color=region, label=louvain_num)) +
    geom_hline(yintercept=cell_lines, linetype='dashed', color='gray50') +
    geom_vline(xintercept=opt$animal_threshold, linetype='dashed', color='gray50') +
    geom_point(size=2, alpha=0.8) +
    geom_text(size=2.6, hjust=-0.3, vjust=-0.3, show.legend=FALSE) +
    scale_y_log10() +
    theme_classic(base_size=13) +
    xlab('n animals') + ylab('n cells (log10 scale)') +
    ggtitle(paste0(opt$cell_type, ': cells vs animals per subcluster x region\n',
                    'dashed lines: cells = ', paste(cell_lines, collapse=' / '),
                    ', animals = ', opt$animal_threshold))

out = file.path(opt$figdir, paste0(opt$cell_type, '_condition_sparsity_scatter.png'))
ggsave(p, file=out, height=8, width=11, dpi=150)
message('Saved: ', out)
