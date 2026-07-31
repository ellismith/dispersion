#!/usr/bin/env Rscript
# plot_overlap.R
# Venn diagrams and direction concordance between gene_variance and DGLM.
# Accepts either --dglm_master (TSV) or --dglm_dir (checkpoints with RDS files).

library(optparse)

option_list = list(
    make_option('--gv_master',   type='character',
                default='/scratch/easmit31/dispersion/gene_variance/results_log/master_between.tsv'),
    make_option('--dglm_master', type='character', default=NULL,
                help='DGLM master TSV (preferred)'),
    make_option('--dglm_dir',    type='character', default=NULL,
                help='DGLM checkpoints dir with RDS files (legacy)'),
    make_option('--outdir',      type='character',
                default='/scratch/easmit31/dispersion/plotting_scripts/figures'),
    make_option('--q_gv',         type='double',  default=0.05),
    make_option('--q_dglm',       type='double',  default=0.05),
    make_option('--outfmt',       type='character', default='png'),
    make_option('--by_ct_region', action='store_true', default=FALSE),
    make_option('--cell_type',    type='character', default=NULL)
)
opt = parse_args(OptionParser(option_list=option_list))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

library(ggplot2)
library(ggvenn)

CELL_TYPES = c('astrocytes','basket_cells','cerebellar_neurons','ependymal_cells',
               'GABAergic_neurons','glutamatergic_neurons','medium_spiny_neurons',
               'microglia','midbrain_neurons','opc','oligodendrocytes','vascular_cells')

save_fig = function(p, name, height=5, width=6) {
    out = file.path(opt$outdir, paste0(name, '.', opt$outfmt))
    if (opt$outfmt == 'png') ggsave(p, file=out, height=height, width=width, dpi=150)
    else ggsave(p, file=out, height=height, width=width, useDingbats=FALSE)
    message('Saved: ', out)
}

# ── load ortholog map ─────────────────────────────────────────────────────
orth = read.csv('/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv')
orth = orth[orth$Human.homology.type == 'ortholog_one2one' & orth$Human.gene.name != '',
            c('Gene.stable.ID','Human.gene.name')]
orth = orth[!duplicated(orth$Gene.stable.ID),]
rownames(orth) = orth$Gene.stable.ID

# ── load gene_variance ────────────────────────────────────────────────────
message('Loading gene_variance results')
gv = read.table(opt$gv_master, sep='\t', header=TRUE, stringsAsFactors=FALSE)
if (!is.null(opt$cell_type)) gv = gv[gv$cell_type == opt$cell_type,]
gv$symbol = ifelse(!is.na(gv$human_symbol) & gv$human_symbol != '', gv$human_symbol, gv$ensembl_id)
gv$sig    = !is.na(gv$qvalue) & gv$qvalue < opt$q_gv
gv$up     = gv$sig & gv$age_slope > 0
gv$key    = if (opt$by_ct_region) paste(gv$symbol, gv$cell_type, gv$region, sep='|') else gv$symbol
message('  sig: ', sum(gv$sig))

# ── load DGLM ─────────────────────────────────────────────────────────────
message('Loading DGLM results')
if (!is.null(opt$dglm_master)) {
    # load from master TSV
    dglm = read.table(opt$dglm_master, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    if (!is.null(opt$cell_type)) dglm = dglm[dglm$cell_type == opt$cell_type,]
    dglm$sig = !is.na(dglm$qvalue) & dglm$qvalue < opt$q_dglm
    dglm$up  = dglm$sig & dglm$beta > 0
    dglm$key = if (opt$by_ct_region) paste(dglm$symbol, dglm$cell_type, dglm$region, sep='|') else dglm$symbol
} else if (!is.null(opt$dglm_dir)) {
    # load from RDS files
    dglm.list = list()
    cts = if (!is.null(opt$cell_type)) opt$cell_type else CELL_TYPES
    for (ct in cts) {
        # try mashr RDS first, then raw DGLM RDS
        rds.file = file.path(opt$dglm_dir, paste0(ct, '_dglm_mashr_results.rds'))
        if (!file.exists(rds.file))
            rds.file = file.path(opt$dglm_dir, paste0(ct, '_dglm_results.rds'))
        if (!file.exists(rds.file)) next
        obj          = readRDS(rds.file)
        dglm.results = if (!is.null(obj$dglm_results)) obj$dglm_results else obj$array
        human.sym    = obj$human_symbols
        regions      = if (!is.null(obj$regions)) obj$regions else
            region.levels[apply(dglm.results[,'beta',, drop=FALSE], 3,
                                function(x) sum(!is.na(x) & x != 0) > 0)]
        for (r in regions) {
            pval = dglm.results[,'pval', r]
            beta = dglm.results[,'beta', r]
            qval = p.adjust(pval, method='fdr')
            syms = sapply(dimnames(dglm.results)[[1]], function(g) {
                if (g %in% rownames(orth)) return(orth[g,'Human.gene.name'])
                s = human.sym[g]; if (!is.na(s) && s != g) return(s); return(g)
            })
            dglm.list[[length(dglm.list)+1]] = data.frame(
                ensembl_id=dimnames(dglm.results)[[1]], symbol=syms,
                cell_type=ct, region=r, beta=beta, qvalue=qval,
                stringsAsFactors=FALSE)
        }
    }
    dglm      = do.call(rbind, dglm.list)
    dglm      = dglm[abs(dglm$beta) <= 100,]
    dglm$sig  = !is.na(dglm$qvalue) & dglm$qvalue < opt$q_dglm
    dglm$up   = dglm$sig & dglm$beta > 0
    dglm$key  = if (opt$by_ct_region) paste(dglm$symbol, dglm$cell_type, dglm$region, sep='|') else dglm$symbol
} else {
    stop('Must provide --dglm_master or --dglm_dir')
}
message('  sig: ', sum(dglm$sig))

id_label = if (opt$by_ct_region) 'gene x cell type x region' else 'gene symbol'
suffix   = if (opt$by_ct_region) '_ctregion' else ''

# ── venn helper ───────────────────────────────────────────────────────────
make_venn = function(gv.sub, dglm.sub, title, fname) {
    gv.sig   = unique(gv.sub$key[gv.sub$sig])
    dglm.sig = unique(dglm.sub$key[dglm.sub$sig])
    if (length(gv.sig) == 0 && length(dglm.sig) == 0) return(invisible(NULL))

    overlap = intersect(gv.sig, dglm.sig)
    if (length(overlap) > 0) {
        gv.dir   = gv.sub[gv.sub$sig & gv.sub$key %in% overlap, c('key','up')]
        gv.dir   = gv.dir[!duplicated(gv.dir$key),]
        dglm.dir = dglm.sub[dglm.sub$sig & dglm.sub$key %in% overlap, c('key','up')]
        dglm.dir = dglm.dir[!duplicated(dglm.dir$key),]
        merged   = merge(gv.dir, dglm.dir, by='key', suffixes=c('.gv','.dglm'))
        if (nrow(merged) > 0) {
            pct = round(100*sum(merged$up.gv == merged$up.dglm)/nrow(merged), 1)
            subtitle = paste0('overlap n=', length(overlap), ' (', id_label, ')',
                              ', concordant: ', pct, '%')
        } else subtitle = paste0('overlap n=', length(overlap))
    } else subtitle = paste0('no overlap (', id_label, ')')

    venn.list = list(a=gv.sig, b=dglm.sig)
    names(venn.list) = c(paste0('gene_variance\n(q<',opt$q_gv,')'),
                         paste0('DGLM\n(q<',opt$q_dglm,')'))
    p = ggvenn(venn.list, fill_color=c('#4575b4','#d73027'), fill_alpha=0.4,
               stroke_size=0.5, set_name_size=3.5, text_size=3.5) +
        labs(title=title, subtitle=subtitle) +
        theme(plot.title=element_text(size=11, face='bold'),
              plot.subtitle=element_text(size=9))
    save_fig(p, fname)
}

# ── per cell type venns ───────────────────────────────────────────────────
cts = if (!is.null(opt$cell_type)) opt$cell_type else CELL_TYPES
for (ct in cts) {
    make_venn(gv[gv$cell_type==ct,], dglm[dglm$cell_type==ct,],
              paste0(ct, ' — between-individual dispersion'),
              paste0('venn_', ct, suffix))
}

# ── global venn ───────────────────────────────────────────────────────────
make_venn(gv, dglm, 'All cell types — between-individual dispersion',
          paste0('venn_global', suffix))

# ── concordance bar chart ─────────────────────────────────────────────────
concordance.df = do.call(rbind, lapply(cts, function(ct) {
    gv.ct   = gv[gv$cell_type==ct & gv$sig,]
    dglm.ct = dglm[dglm$cell_type==ct & dglm$sig,]
    overlap = intersect(unique(gv.ct$key), unique(dglm.ct$key))
    if (length(overlap) == 0) return(NULL)
    gv.dir   = gv.ct[gv.ct$key %in% overlap, c('key','up')]
    gv.dir   = gv.dir[!duplicated(gv.dir$key),]
    dglm.dir = dglm.ct[dglm.ct$key %in% overlap, c('key','up')]
    dglm.dir = dglm.dir[!duplicated(dglm.dir$key),]
    merged   = merge(gv.dir, dglm.dir, by='key', suffixes=c('.gv','.dglm'))
    if (nrow(merged) == 0) return(NULL)
    data.frame(cell_type=ct, n_overlap=nrow(merged),
               pct_concordant=100*sum(merged$up.gv==merged$up.dglm)/nrow(merged))
}))

if (!is.null(concordance.df) && nrow(concordance.df) > 0) {
    concordance.df$cell_type = factor(concordance.df$cell_type, levels=CELL_TYPES)
    p = ggplot(concordance.df, aes(cell_type, pct_concordant, fill=pct_concordant)) +
        geom_bar(stat='identity') +
        geom_text(aes(label=paste0('n=',n_overlap)), vjust=-0.3, size=3) +
        scale_fill_gradient(low='#f7f7f7', high='#2166ac', limits=c(0,100)) +
        scale_y_continuous(limits=c(0,115), breaks=seq(0,100,25)) +
        geom_hline(yintercept=50, linetype='dashed', color='gray50') +
        theme_classic(base_size=12) +
        theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank()) +
        ylab('% directional concordance') +
        ggtitle(paste0('Direction concordance (', id_label, ')'))
    save_fig(p, paste0('concordance_by_celltype', suffix), height=5, width=8)
}

message('done.')
