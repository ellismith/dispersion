#!/usr/bin/env Rscript

suppressMessages({
  library(optparse)
  library(edgeR)
  library(ggplot2)
})

option_list = list(
  make_option('--master', type='character'),
  make_option('--raw_counts_dir', type='character'),
  make_option('--n_top', type='integer', default=6),
  make_option('--lfsr_thresh', type='double', default=0.05),
  make_option('--title', type='character', default=NULL),
  make_option('--out', type='character', default='top_genes_expr_vs_age.png'),
  make_option('--ncol', type='integer', default=4,
              help='Number of facet columns [default: %default]'),
  make_option('--width', type='double', default=30,
              help='Figure width in inches [default: %default]'),
  make_option('--height', type='double', default=16,
              help='Figure height in inches [default: %default]'),
  make_option('--dpi', type='integer', default=180,
              help='Output DPI [default: %default]')
)

opt = parse_args(OptionParser(option_list=option_list))

if (is.null(opt$master) || !file.exists(opt$master)) {
  stop('Missing or invalid --master file: ', opt$master)
}

if (is.null(opt$raw_counts_dir) || !dir.exists(opt$raw_counts_dir)) {
  stop('Missing or invalid --raw_counts_dir: ', opt$raw_counts_dir)
}

df = read.csv(opt$master, sep='\t', stringsAsFactors=FALSE)

required_cols = c(
  'mash_lfsr', 'mash_beta', 'ensembl_id',
  'symbol', 'cell_type', 'region'
)

missing_cols = setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop('Missing required columns in master file: ',
       paste(missing_cols, collapse=', '))
}

sig = df[!is.na(df$mash_lfsr) & df$mash_lfsr < opt$lfsr_thresh, ]

if (nrow(sig) == 0) {
  stop('No rows pass mash_lfsr < ', opt$lfsr_thresh)
}

inc = sig[!is.na(sig$mash_beta) & sig$mash_beta > 0, ]
dec = sig[!is.na(sig$mash_beta) & sig$mash_beta <= 0, ]

top_inc = inc[order(-inc$mash_beta), , drop=FALSE]
top_inc = top_inc[seq_len(min(opt$n_top, nrow(top_inc))), , drop=FALSE]

top_dec = dec[order(dec$mash_beta), , drop=FALSE]
top_dec = top_dec[seq_len(min(opt$n_top, nrow(top_dec))), , drop=FALSE]

selected = rbind(
  data.frame(
    ensembl_id=top_inc$ensembl_id,
    symbol=top_inc$symbol,
    subcluster=top_inc$cell_type,
    region=top_inc$region,
    mash_beta=top_inc$mash_beta,
    direction='Increasing',
    stringsAsFactors=FALSE
  ),
  data.frame(
    ensembl_id=top_dec$ensembl_id,
    symbol=top_dec$symbol,
    subcluster=top_dec$cell_type,
    region=top_dec$region,
    mash_beta=top_dec$mash_beta,
    direction='Decreasing',
    stringsAsFactors=FALSE
  )
)

if (nrow(selected) == 0) {
  stop('No increasing or decreasing genes available after filtering.')
}

cat('Selected genes:\n')
print(selected[, c('symbol', 'subcluster', 'region', 'mash_beta', 'direction')])

plot_rows = list()
row_i = 1

for (i in seq_len(nrow(selected))) {
  s = selected[i, ]

  filt_file = file.path(
    opt$raw_counts_dir,
    'filtered',
    paste0(s$subcluster, '_', s$region, '_filtered_cutoff0.5.csv')
  )

  pb_file = file.path(
    opt$raw_counts_dir,
    paste0(s$subcluster, '_', s$region, '_pseudobulk.csv')
  )

  counts_file = if (file.exists(filt_file)) filt_file else pb_file

  if (!file.exists(counts_file)) {
    cat('MISSING counts:', s$symbol, '|', s$subcluster, '|', s$region, '\n')
    next
  }

  counts = as.matrix(read.csv(
    counts_file,
    row.names=1,
    check.names=FALSE,
    stringsAsFactors=FALSE
  ))

  storage.mode(counts) = 'numeric'

  if (!(s$ensembl_id %in% rownames(counts))) {
    cat('Gene not found:', s$symbol, '|', s$ensembl_id, '\n')
    next
  }

  meta_file = list.files(
    opt$raw_counts_dir,
    pattern=paste0('^', s$subcluster, '_metadata\\.csv$'),
    full.names=TRUE
  )

  if (length(meta_file) == 0) {
    cat('MISSING metadata:', s$subcluster, '\n')
    next
  }

  meta = read.csv(meta_file[1], stringsAsFactors=FALSE)

  if (!all(c('region', 'animal_id', 'age') %in% colnames(meta))) {
    cat('Invalid metadata columns:', meta_file[1], '\n')
    next
  }

  meta = meta[meta$region == s$region, , drop=FALSE]
  meta = meta[match(colnames(counts), meta$animal_id), , drop=FALSE]

  if (all(is.na(meta$animal_id))) {
    cat('No animal-ID match:', s$symbol, '|', s$subcluster, '|', s$region, '\n')
    next
  }

  y = DGEList(counts=counts)
  y = calcNormFactors(y, method='TMM')
  logCPM = cpm(y, log=TRUE, prior.count=0.5)

  sym = if (is.na(s$symbol) || s$symbol == '') s$ensembl_id else s$symbol

  label = sprintf(
    '%s\n%s | %s\nbeta = %.3f',
    sym,
    s$subcluster,
    s$region,
    s$mash_beta
  )

  plot_rows[[row_i]] = data.frame(
    age=as.numeric(meta$age),
    expr=as.numeric(logCPM[s$ensembl_id, ]),
    gene_label=label,
    direction=s$direction,
    stringsAsFactors=FALSE
  )

  row_i = row_i + 1
}

if (length(plot_rows) == 0) {
  stop('No plottable gene-expression data were found.')
}

plot_df = do.call(rbind, plot_rows)
plot_df = plot_df[!is.na(plot_df$age) & !is.na(plot_df$expr), , drop=FALSE]

if (nrow(plot_df) == 0) {
  stop('All expression or age values were missing after matching metadata.')
}

label_order = unique(
  plot_df$gene_label[
    order(match(plot_df$direction, c('Increasing', 'Decreasing')))
  ]
)

plot_df$gene_label = factor(
  plot_df$gene_label,
  levels=label_order
)

plot_df$direction = factor(
  plot_df$direction,
  levels=c('Increasing', 'Decreasing')
)

title = if (!is.null(opt$title)) {
  opt$title
} else {
  'Top increasing and decreasing dispersion-age genes'
}

p = ggplot(plot_df, aes(x=age, y=expr)) +
  geom_point(
    aes(color=direction),
    size=2.8,
    alpha=0.80
  ) +
  geom_smooth(
    method='lm',
    se=TRUE,
    color='black',
    linewidth=0.9,
    fill='grey70',
    alpha=0.35
  ) +
  facet_wrap(
    ~gene_label,
    scales='free_y',
    ncol=opt$ncol,
    strip.position='top'
  ) +
  scale_x_continuous(
    position='top',
    breaks=pretty(plot_df$age, n=5)
  ) +
  scale_color_manual(
    values=c(
      'Increasing'='firebrick',
      'Decreasing'='steelblue'
    )
  ) +
  theme_classic(base_size=16) +
  theme(
    legend.position='top',
    legend.text=element_text(size=15),
    legend.key.width=grid::unit(1.2, 'cm'),
    strip.placement='outside',
    strip.background=element_rect(
      fill='grey95',
      color='grey35',
      linewidth=0.5
    ),
    strip.text=element_text(
      size=12,
      face='bold',
      lineheight=1.05,
      margin=margin(t=7, r=6, b=7, l=6)
    ),
    axis.title.x.top=element_text(
      size=17,
      face='bold',
      margin=margin(b=10)
    ),
    axis.title.y=element_text(
      size=17,
      face='bold',
      margin=margin(r=10)
    ),
    axis.text.x.top=element_text(
      size=14,
      margin=margin(b=3)
    ),
    axis.text.y=element_text(size=14),
    axis.ticks.length=grid::unit(0.20, 'cm'),
    plot.title=element_text(
      size=23,
      face='bold',
      margin=margin(b=12)
    ),
    plot.margin=margin(t=10, r=20, b=15, l=15),
    panel.spacing=grid::unit(1.1, 'lines')
  ) +
  labs(
    x='Age',
    y='TMM-normalized log2(CPM)',
    color=NULL,
    title=title
  )

dir.create(dirname(opt$out), recursive=TRUE, showWarnings=FALSE)

ggsave(
  opt$out,
  p,
  width=opt$width,
  height=opt$height,
  dpi=opt$dpi,
  limitsize=FALSE
)

cat('\nSaved:', opt$out, '\n')
