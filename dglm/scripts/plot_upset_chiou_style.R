#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(patchwork); library(dplyr); library(optparse)})

option_list = list(
  make_option('--ct', type='character', default='all'),
  make_option('--pipeline', type='character', default='baseline'),
  make_option('--lfsr', type='double', default=0.05),
  make_option('--top_n', type='integer', default=20),
  make_option('--groupby', type='character', default='region'),
  make_option('--base_dir', type='character', default=NULL),
  make_option('--font_scale', type='double', default=1),
  make_option('--out', type='character', default=NULL),
  make_option('--min_pct', type='double', default=0)
)
opt = parse_args(OptionParser(option_list=option_list))

region_colors = c("ACC"="#5555b0","CN"="#2caa4d","dlPFC"="#345296","EC"="#532b70",
  "HIP"="#84db55","IPP"="#6c9cc0","lCb"="#e6cd69","M1"="#467aa1",
  "MB"="#fc7462","mdTN"="#117867","NAc"="#148e34")
lcb_affected = c("opc","microglia","astrocytes","vascular_cells","oligodendrocytes")
cell_types = c("GABAergic_neurons","glutamatergic_neurons","cerebellar_neurons",
               "medium_spiny_neurons","midbrain_neurons","basket_cells",
               "microglia","astrocytes","opc","oligodendrocytes","vascular_cells","ependymal_cells")
BASE_SIZE = 14

load_ct = function(ct, pipeline) {
  if(!is.null(opt$base_dir)) {
    base = opt$base_dir  # lCb excluded at condition-filter time; same dir for all
  } else if(pipeline=="baseline") {
    base = if(ct %in% lcb_affected) "disp_age__FINAL_autosomeX_300c_noLcb_QCfiltered" else "disp_age__FINAL_autosomeX_300c_QCfiltered"
  } else {
    base = if(ct %in% lcb_affected) "disp_age__mincells10_300c_noLcb_QCfiltered" else "disp_age__mincells10_300c_QCfiltered"
  }
  f = if (startsWith(base, "/")) paste0(base,"/",ct,"/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv") else paste0("/scratch/easmit31/dispersion/dglm/",base,"/",ct,"/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv")
  if(!file.exists(f)) return(NULL)
  df = read.csv(f,sep="\t"); df$cell_type_broad=ct; df
}

cts = if(opt$ct=="all") cell_types else opt$ct
loaded = Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0,
                lapply(cts, load_ct, pipeline=opt$pipeline))
if (length(loaded) == 0) { cat("SKIP: no master tables loaded for", opt$ct, "\n"); quit(save="no", status=0) }
all_data = do.call(rbind, loaded)
sig = all_data[!is.na(all_data$mash_lfsr) & all_data$mash_lfsr < opt$lfsr, , drop=FALSE]
if (opt$min_pct > 0) {
  all_data$cond = paste(all_data$cell_type, all_data$region)
  tst = tapply(!is.na(all_data$mash_lfsr), all_data$cond, sum)
  sg  = tapply(!is.na(all_data$mash_lfsr) & all_data$mash_lfsr < opt$lfsr, all_data$cond, sum)
  keep = names(tst)[100*sg/pmax(tst,1) >= opt$min_pct]
  sig = sig[paste(sig$cell_type, sig$region) %in% keep, , drop=FALSE]
  cat("min_pct", opt$min_pct, ": kept", length(keep), "of", length(tst), "conditions\n")
  if (nrow(sig) == 0) { cat("SKIP: nothing above min_pct\n"); quit(save="no", status=0) }
}
if (nrow(sig) == 0) { cat("SKIP: no significant rows at lfsr<", opt$lfsr, "for", opt$ct, "\n"); quit(save="no", status=0) }
sig$subcluster_num = sub(".*_([0-9]+)$","\\1",sig$cell_type)
ct_label = if(opt$ct=="all") "All cell types" else opt$ct

if(opt$groupby == "region") {
  # each row = unique set of regions a subcluster shows effects in
  group_col = "region"; member_col = "subcluster"
  inc_combos = sig %>% filter(mash_beta>0) %>%
    group_by(cell_type) %>% summarise(combo=paste(sort(unique(region)),collapse=","), .groups="drop")
  dec_combos = sig %>% filter(mash_beta<=0) %>%
    group_by(cell_type) %>% summarise(combo=paste(sort(unique(region)),collapse=","), .groups="drop")
  all_members = sort(unique(sig$region))
  member_colors = region_colors[all_members]
  use_color = TRUE
  title = paste0(ct_label," (",opt$pipeline,", lfsr<",opt$lfsr,"): transcriptional dispersion changes by brain region combination")
  subtitle = paste0("Each row = unique set of brain regions in which one or more subclusters show significant age-associated dispersion changes. ",
                    "Colored tiles = regions in that combination. Gene counts are unique across all subclusters sharing that combination. ",
                    "Sorted by total gene count.")
  bot_label = "N subclusters"
  xmid_label = "Brain region"
} else {
  # each row = unique set of subclusters a region shows effects in
  inc_combos = sig %>% filter(mash_beta>0) %>%
    group_by(region) %>% summarise(combo=paste(sort(unique(subcluster_num)),collapse=","), .groups="drop")
  dec_combos = sig %>% filter(mash_beta<=0) %>%
    group_by(region) %>% summarise(combo=paste(sort(unique(subcluster_num)),collapse=","), .groups="drop")
  all_members = as.character(sort(as.integer(unique(sig$subcluster_num))))
  member_colors = NULL
  use_color = FALSE
  title = paste0(ct_label," (",opt$pipeline,", lfsr<",opt$lfsr,"): transcriptional dispersion changes by subcluster combination")
  subtitle = paste0("Each row = unique set of subclusters in which one or more brain regions show significant age-associated dispersion changes. ",
                    "Filled tiles = subclusters present in that combination. Gene counts are unique across all regions sharing that combination. ",
                    "Sorted by total gene count.")
  bot_label = "N regions"
  xmid_label = "Subcluster #"
}

# unique genes per combo
get_combo_df = function(combos_df, join_col, direction_filter) {
  genes = sig %>% filter(if(direction_filter=="inc") mash_beta>0 else mash_beta<=0) %>%
    left_join(combos_df, by=join_col) %>%
    group_by(combo) %>% summarise(n_genes=n_distinct(ensembl_id), .groups="drop")
  subs = combos_df %>% count(combo, name="n_sub")
  left_join(genes, subs, by="combo") %>% arrange(desc(n_genes)) %>% head(opt$top_n)
}

join_col = if(opt$groupby=="region") "cell_type" else "region"
inc_df = get_combo_df(inc_combos, join_col, "inc")
dec_df = get_combo_df(dec_combos, join_col, "dec")

all_combos = union(inc_df$combo, dec_df$combo)
combo_totals = data.frame(combo=all_combos) %>%
  left_join(inc_df %>% select(combo, n_inc=n_genes, n_sub_inc=n_sub), by="combo") %>%
  left_join(dec_df %>% select(combo, n_dec=n_genes, n_sub_dec=n_sub), by="combo") %>%
  mutate(n_inc=ifelse(is.na(n_inc),0,n_inc), n_dec=ifelse(is.na(n_dec),0,n_dec),
         n_sub=pmax(ifelse(is.na(n_sub_inc),0,n_sub_inc), ifelse(is.na(n_sub_dec),0,n_sub_dec))) %>%
  arrange(n_inc+n_dec)
if (nrow(combo_totals) == 0) { cat("SKIP: no combos to plot for", opt$ct, "\n"); quit(save="no", status=0) }
combo_totals$row_label = factor(seq_len(nrow(combo_totals)), levels=seq_len(nrow(combo_totals)))

# bars
p_left = ggplot(combo_totals, aes(x=-n_dec, y=row_label)) +
  geom_col(fill="steelblue", width=0.7) + geom_vline(xintercept=0,color="black",linewidth=0.4) +
  scale_x_continuous(labels=function(x) paste0(abs(x))) +
  theme_minimal(base_size=BASE_SIZE) +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), panel.grid.major.y=element_blank()) +
  labs(x="Decreased dispersion\n(unique genes)")

p_right = ggplot(combo_totals, aes(x=n_inc, y=row_label)) +
  geom_col(fill="firebrick", width=0.7) + geom_vline(xintercept=0,color="black",linewidth=0.4) +
  theme_minimal(base_size=BASE_SIZE) +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), panel.grid.major.y=element_blank()) +
  labs(x="Increased dispersion\n(unique genes)")

# dot matrix: one box per member, colored by that member's MAJORITY direction in this combo
mem_col = if(opt$groupby == "region") "region" else "subcluster_num"
grp_col = if(opt$groupby == "region") "cell_type" else "region"
dot_rows = list()
for(i in seq_len(nrow(combo_totals))) {
  cmb = as.character(combo_totals$combo[i])
  members = strsplit(cmb, ",")[[1]]
  grps = unique(c(as.character(inc_combos[[grp_col]][inc_combos$combo == cmb]),
                  as.character(dec_combos[[grp_col]][dec_combos$combo == cmb])))
  sub = sig[as.character(sig[[grp_col]]) %in% grps,]
  for(m in all_members) {
    ms = sub[sub[[mem_col]] == m,]
    n_i = sum(ms$mash_beta > 0); n_d = sum(ms$mash_beta <= 0)
    dot_rows[[paste(i,m)]] = data.frame(row_label=combo_totals$row_label[i],
      member=factor(m, levels=all_members), present=m %in% members,
      dircol=if (n_i + n_d == 0) "none" else if (n_i > n_d) "inc" else "dec")
  }
}
dot_df = do.call(rbind, dot_rows)
pres = dot_df[dot_df$present,]
abs_df = dot_df[!dot_df$present,]
p_mid = ggplot() +
  geom_tile(data=abs_df, aes(x=member, y=row_label), fill="grey93", width=0.62, height=0.62) +
  geom_tile(data=pres, aes(x=member, y=row_label, fill=dircol), width=0.62, height=0.62) +
  scale_fill_manual(values=c(inc="firebrick", dec="steelblue", none="grey80"), guide="none")

p_mid = p_mid + theme_minimal(base_size=BASE_SIZE) +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=BASE_SIZE-4),
        axis.text.y=element_blank(), panel.grid=element_blank()) +
  labs(x=xmid_label, y=NULL)

combined = (p_left | p_mid | p_right) +
  plot_layout(widths=c(2,if(opt$groupby=="region") 1.2 else 2.5, 2)) +
  plot_annotation(title=title, subtitle=subtitle,
    theme=theme(plot.title=element_text(size=BASE_SIZE+2,face="bold"),
                plot.subtitle=element_text(size=BASE_SIZE-3)))

if(is.null(opt$out)) {
  ct_file = if(opt$ct=="all") "ALL" else opt$ct
  opt$out = paste0("/scratch/easmit31/dispersion/dglm/figures/upset_chiou_",opt$groupby,"_",ct_file,"_",opt$pipeline,".png")
}
h = max(8, nrow(combo_totals)*0.45)
update_geom_defaults('text', list(size = 3.88 * opt$font_scale))
combined = combined & theme(text = element_text(size = 11 * opt$font_scale))
ggsave(opt$out, combined, width=18, height=h, dpi=300)
cat("Saved:", opt$out, "\n")
