#!/usr/bin/env Rscript
#
# Expression-vs-age scatter, one panel per gene, using the CORRECTED pipeline's
# actual response values (TMM-normalized log2-CPM straight from dglm_model.R's
# own transformation -- not a separately-recomputed log1p(CPM)).
#
# Three fixes from the earlier version of this plot:
#   1. Each panel's y-axis is scaled independently (facet_wrap(scales="free_y"))
#      instead of one shared range across all genes.
#   2. Confidence-interval shading is back on (geom_smooth(se=TRUE)) -- an
#      earlier version of this script had it stripped out; this restores it.
#   3. Y-values come directly from the same TMM+logCPM matrix DGLM was fit on,
#      not a separate log(1+CPM) recomputation -- these are NOT identical
#      transforms (TMM adjusts for composition bias beyond plain CPM, and
#      prior.count=0.5 differs from log1p's implicit +1), so this makes the
#      plot consistent with what the model actually saw.
#
# Usage:
#   Rscript plot_expr_vs_age_fixed.R \
#     --master <master_dglm_combined.tsv> \
#     --raw_counts_dir <dir with *_pseudobulk.csv and *_metadata.csv> \
#     --cell_type <e.g. GABAergic_neurons_0> --region <e.g. ACC> \
#     --lfsr_thresh 0.2 --n_examples 6 \
#     --out <output.png>

suppressMessages({
  library(optparse)
  library(edgeR)
  library(ggplot2)
})

option_list = list(
  make_option('--master', type='character'),
  make_option('--raw_counts_dir', type='character'),
  make_option('--cell_type', type='character'),
  make_option('--region', type='character'),
  make_option('--lfsr_thresh', type='double', default=0.2),
  make_option('--n_examples', type='integer', default=6,
              help='genes per direction (increasing/decreasing) to plot'),
  make_option('--out', type='character', default='expr_vs_age.png')
)
opt = parse_args(OptionParser(option_list=option_list))

for (req in c('master','raw_counts_dir','cell_type','region')) {
  if (is.null(opt[[req]])) stop('--', req, ' is required')
}

master = read.csv(opt$master, sep='\t', stringsAsFactors=FALSE)
sub = master[master$cell_type == opt$cell_type & master$region == opt$region &
             master$mash_lfsr < opt$lfsr_thresh, ]
if (nrow(sub) == 0) stop('No significant genes found for ', opt$cell_type, ' / ', opt$region,
                          ' at mash_lfsr < ', opt$lfsr_thresh)

inc = sub[sub$mash_beta > 0, ]
dec = sub[sub$mash_beta <= 0, ]
inc = inc[order(-inc$mash_beta), ][seq_len(min(opt$n_examples, nrow(inc))), ]
dec = dec[order(dec$mash_beta), ][seq_len(min(opt$n_examples, nrow(dec))), ]

genes = rbind(
  data.frame(ensembl_id=inc$ensembl_id, symbol=inc$symbol, mash_beta=inc$mash_beta, direction='Increasing dispersion'),
  data.frame(ensembl_id=dec$ensembl_id, symbol=dec$symbol, mash_beta=dec$mash_beta, direction='Decreasing dispersion')
)

# --- Recompute the SAME TMM-normalized log2-CPM values DGLM was fit on ---
pb_file = file.path(opt$raw_counts_dir, paste0(opt$cell_type, '_', opt$region, '_pseudobulk.csv'))
filtered_file = file.path(opt$raw_counts_dir, 'filtered',
                           paste0(opt$cell_type, '_', opt$region, '_filtered_cutoff0.5.csv'))
counts_file = if (file.exists(filtered_file)) filtered_file else pb_file
if (!file.exists(counts_file)) stop('Counts file not found: ', counts_file)

counts = read.csv(counts_file, row.names=1, check.names=FALSE)
counts = as.matrix(counts); storage.mode(counts) = 'numeric'

meta_file = file.path(opt$raw_counts_dir, paste0(sub('_[0-9]+$', '', opt$cell_type), '_metadata.csv'))
if (!file.exists(meta_file)) {
  # fall back to exact cell_type name (covers subcluster-suffixed metadata files)
  meta_file = file.path(opt$raw_counts_dir, paste0(opt$cell_type, '_metadata.csv'))
}
meta = read.csv(meta_file, stringsAsFactors=FALSE)
meta = meta[meta$region == opt$region, ]
meta = meta[match(colnames(counts), meta$animal_id), ]

y = DGEList(counts=counts)
y = calcNormFactors(y, method='TMM')
logCPM = cpm(y, log=TRUE, prior.count=0.5)

# --- Assemble long-format plotting data ---
plot_df = do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
  g = genes$ensembl_id[i]
  if (!g %in% rownames(logCPM)) return(NULL)
  data.frame(
    animal_id = colnames(logCPM),
    age = meta$age,
    expr = as.numeric(logCPM[g, ]),
    gene_label = sprintf('%s (beta=%.3f)', ifelse(is.na(genes$symbol[i]) | genes$symbol[i]=='', g, genes$symbol[i]), genes$mash_beta[i]),
    direction = genes$direction[i]
  )
}))
plot_df = plot_df[!is.na(plot_df$age) & !is.na(plot_df$expr), ]

p = ggplot(plot_df, aes(x=age, y=expr)) +
  geom_point(aes(color=direction), size=2, alpha=0.75) +
  geom_smooth(method='loess', se=TRUE, color='black', linewidth=0.6, fill='grey70') +
  facet_wrap(~gene_label, scales='free_y') +
  scale_color_manual(values=c('Increasing dispersion'='firebrick', 'Decreasing dispersion'='steelblue')) +
  theme_classic(base_size=13) +
  theme(legend.position='top', strip.text=element_text(face='bold', size=9)) +
  labs(x='Age', y='TMM-normalized log2(CPM)', color=NULL,
       title=paste0(opt$cell_type, ', ', opt$region, ' -- expression vs. age'))

ggsave(opt$out, p, width=4*min(3, ceiling(sqrt(nrow(genes)))), height=4*ceiling(nrow(genes)/min(3, ceiling(sqrt(nrow(genes))))), dpi=150, limitsize=FALSE)
cat('Saved:', opt$out, '\n')
