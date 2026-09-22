#!/usr/bin/env Rscript
#
# plot_gene_count_heatmaps.R
#
# ONE flexible heatmap script: unit x region, for either gene-count metrics
# (genes tested / genes significant, from raw qvalue) OR mash-based effect
# metrics (pct significant, direction counts, mean mash beta) -- selected via
# --metric. Same formulas as the explorer's build_heatmap(), so static PNGs
# and the interactive GUI always agree.
#
# Fixes opc/oligodendrocytes subcluster labels: both share a raw "opc-olig_N"
# cell_type string from their shared source h5ad (pre-split louvain naming),
# but pseudobulk already correctly separates them into distinct files per
# OPC_LOUVAIN={12,13} -- confirmed by opc/ containing ONLY 12,13 and
# oligodendrocytes/ containing ONLY 0-11, no overlap. This script relabels
# using the (already-correct) parent directory name -- cosmetic only, no
# data touched.
#
# Metrics:
#   genes_tested        - n rows tested, unit x region
#   genes_significant   - n rows with sig_col < qthresh (raw q-value based)
#   raw_beta_all         - mean |raw beta|, all genes
#   raw_beta_sig          - mean |raw beta|, raw q<qthresh genes only
#   pct_sig               - % of genes with mash_lfsr < lfsr_thresh
#   n_increase             - n genes with mash_lfsr<lfsr_thresh & mash_beta>0
#   n_decrease              - n genes with mash_lfsr<lfsr_thresh & mash_beta<=0
#   n_net                    - n_increase - n_decrease (diverging)
#   mash_beta_sig             - mean mash_beta among sig genes (diverging)
#   mash_beta_mag_sig          - mean |mash_beta| among sig genes
#
# Usage:
#   Rscript plot_gene_count_heatmaps.R --metric genes_tested --dispersion age
#   Rscript plot_gene_count_heatmaps.R --metric n_net --dispersion age_sex --lfsr_thresh 0.2
#   Rscript plot_gene_count_heatmaps.R --metric genes_tested --dispersion age --per_celltype

suppressMessages({
  library(optparse)
  library(ggplot2)
  library(reshape2)
})

option_list = list(
    make_option('--base_dir',     type='character',
                default='/scratch/easmit31/dispersion/dglm/disp_age__raw_counts'),
    make_option('--cutoff',       type='character', default='0.5'),
    make_option('--dispersion',   type='character', default='age',
                help="'age' (V2) or 'age_sex' (V3)"),
    make_option('--figdir',       type='character',
                default='/scratch/easmit31/dispersion/dglm/figures'),
    make_option('--outfmt',       type='character', default='png'),
    make_option('--metric',       type='character', default='genes_tested',
                help='genes_tested | genes_significant | raw_beta_all | raw_beta_sig | pct_sig | n_increase | n_decrease | n_net | mash_beta_sig | mash_beta_mag_sig'),
    make_option('--sig_col',      type='character', default='qvalue'),
    make_option('--qthresh',      type='double',    default=0.05),
    make_option('--lfsr_thresh',  type='double',    default=0.2),
    make_option('--per_celltype', action='store_true', default=FALSE,
                help='Also save one PNG per parent cell type (subcluster-level only)'),
    make_option('--target_celltype', type='character', default=NULL,
                help='restrict to one cell type only')
)
opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

ckpt_suffix = if (opt$dispersion == 'age_sex') paste0('_cutoff', opt$cutoff, '_agesex') else paste0('_cutoff', opt$cutoff)

# ── find and load every cell type's master TSV ─────────────────────────────
cell_type_dirs = list.dirs(opt$base_dir, recursive=FALSE)
if (!is.null(opt$target_celltype)) {
    cell_type_dirs = cell_type_dirs[basename(cell_type_dirs) == opt$target_celltype]
}
all_rows = list()
for (d in cell_type_dirs) {
    tsv_name = if (opt$dispersion == 'age_sex') 'master_dglm_age_disp_combined.tsv' else 'master_dglm_combined.tsv'
    tsv = file.path(d, paste0('dglm_checkpoints', ckpt_suffix), tsv_name)
    if (!file.exists(tsv)) next
    df = read.table(tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    df$true_parent_cell_type = basename(d)
    all_rows[[basename(d)]] = df
}
message('Found master TSVs for ', length(all_rows), ' cell type(s): ', paste(names(all_rows), collapse=', '))
if (length(all_rows) == 0) stop('No master TSV found under ', opt$base_dir, ' with suffix ', ckpt_suffix)

master = do.call(rbind, all_rows)
rownames(master) = NULL
master = master[!is.na(master$beta),]
master$parent_cell_type = master$true_parent_cell_type

# ── fix opc/oligodendrocytes shared "opc-olig_N" label using the correct,
#    already-split parent directory (cosmetic relabel only) ────────────────
opc_olig_idx = grepl('^opc-olig_', master$cell_type)
if (any(opc_olig_idx)) {
    suffix = sub('^opc-olig_', '', master$cell_type[opc_olig_idx])
    master$cell_type[opc_olig_idx] = paste0(master$true_parent_cell_type[opc_olig_idx], '_', suffix)
}

master$sig = if (opt$sig_col %in% colnames(master)) master[[opt$sig_col]] < opt$qthresh else NA
message('Total tested rows: ', nrow(master))

# ── metric value function: given a subset df (one unit x region), returns a scalar ──
metric_fn = switch(opt$metric,
    genes_tested       = function(sub) nrow(sub),
    genes_significant  = function(sub) sum(sub$sig, na.rm=TRUE),
    raw_beta_all        = function(sub) mean(abs(sub$beta), na.rm=TRUE),
    raw_beta_sig          = function(sub) { s = sub[sub$sig,]; if (nrow(s)==0) 0 else mean(abs(s$beta), na.rm=TRUE) },
    pct_sig                = function(sub) { l = sub$mash_lfsr; 100 * sum(!is.na(l) & l < opt$lfsr_thresh) / nrow(sub) },
    n_increase               = function(sub) { l = sub$mash_lfsr; b = sub$mash_beta; sum(!is.na(l) & l < opt$lfsr_thresh & b > 0, na.rm=TRUE) },
    n_decrease                = function(sub) { l = sub$mash_lfsr; b = sub$mash_beta; sum(!is.na(l) & l < opt$lfsr_thresh & b <= 0, na.rm=TRUE) },
    n_net                      = function(sub) { l = sub$mash_lfsr; b = sub$mash_beta; sum(!is.na(l) & l < opt$lfsr_thresh & b > 0, na.rm=TRUE) - sum(!is.na(l) & l < opt$lfsr_thresh & b <= 0, na.rm=TRUE) },
    mash_beta_sig               = function(sub) { l = sub$mash_lfsr; s = sub[!is.na(l) & l < opt$lfsr_thresh,]; if (nrow(s)==0) 0 else mean(s$mash_beta, na.rm=TRUE) },
    mash_beta_mag_sig            = function(sub) { l = sub$mash_lfsr; s = sub[!is.na(l) & l < opt$lfsr_thresh,]; if (nrow(s)==0) 0 else mean(abs(s$mash_beta), na.rm=TRUE) },
    stop('Unknown --metric: ', opt$metric)
)
is_diverging  = opt$metric %in% c('n_net', 'mash_beta_sig')
high_col      = if (is_diverging) NULL else if (opt$metric == 'n_decrease') '#2166ac' else if (opt$metric %in% c('genes_significant')) '#d73027' else '#2166ac'
metric_label  = c(genes_tested='Genes tested', genes_significant=paste0('Genes significant (', opt$sig_col, '<', opt$qthresh, ')'),
                   raw_beta_all='Mean |raw beta| (all genes)', raw_beta_sig=paste0('Mean |raw beta| (raw q<', opt$qthresh, ')'),
                   pct_sig=paste0('% significant (lfsr<', opt$lfsr_thresh, ')'),
                   n_increase='# sig genes, increase', n_decrease='# sig genes, decrease',
                   n_net='Net (increase - decrease)', mash_beta_sig='Mean mash beta (sig genes)',
                   mash_beta_mag_sig='Mean |mash beta| (sig genes)')[[opt$metric]]

# ── generic matrix builder + heatmap plotter ────────────────────────────────
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
                                      limits=c(-bound, bound), na.value='gray90', name=metric_label)
    } else {
        vmax = max(mat, na.rm=TRUE); if (!is.finite(vmax) || vmax == 0) vmax = 1
        p = p + scale_fill_gradient(low='white', high=high_col, na.value='gray90', limits=c(0, vmax), name=metric_label)
    }

    n_units = length(unique(df$unit))
    out = file.path(opt$figdir, paste0(fname, '.', opt$outfmt))
    ggsave(p, file=out, height=0.35*n_units+2, width=0.55*length(unique(df$region))+3, dpi=150, limitsize=FALSE)
    message('Saved: ', out)
}

label_suffix = paste0(' (dispersion~', opt$dispersion, ', cutoff', opt$cutoff, ')')

# ── subcluster level (all cell types combined) ──────────────────────────────
mat_sub = build_matrix(master, 'cell_type')
save_heatmap(mat_sub, paste0(metric_label, ', subcluster x region', label_suffix),
             paste0('heatmap_', opt$metric, '_subcluster_', opt$dispersion))

# ── cell-type level (summed/aggregated across that cell type's subclusters) ─
mat_ct = build_matrix(master, 'parent_cell_type')
save_heatmap(mat_ct, paste0(metric_label, ', cell type x region', label_suffix),
             paste0('heatmap_', opt$metric, '_celltype_', opt$dispersion))

# ── optional: one PNG per parent cell type (subcluster-level only) ─────────
if (opt$per_celltype) {
    for (ct in natural_sort(unique(master$parent_cell_type))) {
        sub_master = master[master$parent_cell_type == ct,]
        mat_one = build_matrix(sub_master, 'cell_type')
        save_heatmap(mat_one, paste0(metric_label, ', ', ct, ' subclusters x region', label_suffix),
                     paste0('heatmap_', opt$metric, '_subcluster_', opt$dispersion, '_', ct))
    }
}

message('done.')
