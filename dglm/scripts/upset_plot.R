#!/usr/bin/env Rscript
#
# Real UpSet plot: bars = count of genes significant in EXACTLY that combination
# of groups, dots = which groups make up each combination. No gene identities
# shown anywhere -- purely counts and group membership.
#
# Usage:
#   Rscript upset_plot.R --master <master_dglm_combined.tsv> \
#     --group_by region|subcluster|combo --top_n 30 --lfsr_thresh 0.2 \
#     --title "My Title" --out <output.png>

suppressMessages({library(optparse); library(ggplot2); library(patchwork)})

option_list = list(
  make_option('--master', type='character'),
  make_option('--group_by', type='character', default='region',
              help='region | subcluster | combo'),
  make_option('--top_n', type='integer', default=30),
  make_option('--lfsr_thresh', type='double', default=0.2),
  make_option('--title', type='character', default=NULL),
  make_option('--out', type='character', default='upset.png')
)
opt = parse_args(OptionParser(option_list=option_list))
if (is.null(opt$master)) stop('--master is required')
if (!opt$group_by %in% c('region','subcluster','combo')) stop('--group_by must be region, subcluster, or combo')

df = read.csv(opt$master, sep='\t')
sig = df[df$mash_lfsr < opt$lfsr_thresh, ]

if (opt$group_by == 'region') {
  sig$grp = sig$region
} else if (opt$group_by == 'subcluster') {
  sig$grp = sig$cell_type
} else {
  sig$grp = paste(sig$cell_type, sig$region, sep=' | ')
}

all_groups = sort(unique(sig$grp))
cat('Number of distinct groups (', opt$group_by, '):', length(all_groups), '\n')

build_membership = function(d) {
  agg = aggregate(d$grp, by=list(ensembl_id=d$ensembl_id), FUN=function(x) paste(sort(unique(x)), collapse=','))
  table(agg$x)
}

inc_combos = build_membership(sig[sig$mash_beta > 0, ])
dec_combos = build_membership(sig[sig$mash_beta <= 0, ])

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
  labs(y='Gene count\n(decreasing)', x=paste(opt$group_by, 'combination'))

dot_rows = list()
for (cn in combo_order) {
  members = strsplit(cn, ',')[[1]]
  for (g in all_groups) dot_rows[[paste(cn,g)]] = data.frame(combo=cn, grp=g, is_member=g %in% members)
}
dot_data = do.call(rbind, dot_rows)
dot_data$combo = factor(dot_data$combo, levels=combo_order)
dot_data$grp = factor(dot_data$grp, levels=rev(all_groups))

p_mid = ggplot(dot_data, aes(x=combo, y=grp)) +
  geom_point(aes(alpha=is_member), color='black', size=2.5) +
  scale_alpha_manual(values=c('TRUE'=1, 'FALSE'=0.08), guide='none') +
  theme_minimal(base_size=ifelse(length(all_groups)>15, 5, 10)) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
        panel.grid.major=element_line(color='grey92')) +
  labs(y=opt$group_by)

title = if (!is.null(opt$title)) opt$title else paste('UpSet of significant genes by', opt$group_by)
combined = p_top / p_mid / p_bottom + plot_layout(heights=c(1.5, 2, 1.5)) +
  plot_annotation(title=title,
                   subtitle=paste0('Bars = gene count per exact combination. Dots = group membership. Top ', length(combo_order), ' combinations shown.'))

h = if (length(all_groups) > 15) 14 else 9
ggsave(opt$out, combined, width=14, height=h, dpi=150)
cat('Saved:', opt$out, '\n')
