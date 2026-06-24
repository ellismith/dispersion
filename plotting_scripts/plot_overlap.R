#!/usr/bin/env Rscript
# plot_overlap.R
#
# Compares gene-level significance between gene_variance (OLS) and DGLM.
# Uses human_symbol as primary gene identifier throughout.
#
# Usage:
#   Rscript plot_overlap.R --q_gv 0.05 --q_dglm 0.05
#   Rscript plot_overlap.R --by_ct_region

library(optparse)

option_list = list(
    make_option('--gv_master',    type='character',
                default='/scratch/easmit31/variability/gene_variance/results_log/master_between.tsv'),
    make_option('--dglm_dir',     type='character',
                default='/scratch/easmit31/variability/dglm/checkpoints'),
    make_option('--outdir',       type='character',
                default='/scratch/easmit31/variability/plotting_scripts/figures'),
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
    if (opt$outfmt == 'png') {
        ggsave(p, file=out, height=height, width=width, dpi=150)
    } else {
        ggsave(p, file=out, height=height, width=width, useDingbats=FALSE)
    }
    message('Saved: ', out)
}

# ── load ortholog map for DGLM results ────────────────────────────────────
orth = read.csv('/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv')
orth = orth[orth$Human.homology.type == 'ortholog_one2one' & orth$Human.gene.name != '',
            c('Gene.stable.ID','Human.gene.name')]
orth = orth[!duplicated(orth$Gene.stable.ID),]
rownames(orth) = orth$Gene.stable.ID

get_symbol = function(ensembl.ids, gene.names.df=NULL) {
    sapply(ensembl.ids, function(g) {
        if (g %in% rownames(orth)) return(orth[g,'Human.gene.name'])
        if (!is.null(gene.names.df) && g %in% rownames(gene.names.df)) {
            nm = gene.names.df[g,'external_gene_name']
            if (!is.na(nm) && nm != '' && nm != g) return(nm)
        }
        return(g)
    })
}

# ── load gene_variance results ─────────────────────────────────────────────
message('Loading gene_variance results')
gv = read.table(opt$gv_master, sep='\t', header=TRUE, stringsAsFactors=FALSE)
if (!is.null(opt$cell_type)) gv = gv[gv$cell_type == opt$cell_type,]

# use human_symbol as primary ID, fall back to ensembl_id
gv$symbol = ifelse(!is.na(gv$human_symbol) & gv$human_symbol != '',
                   gv$human_symbol, gv$ensembl_id)
gv$sig  = !is.na(gv$qvalue) & gv$qvalue < opt$q_gv
gv$up   = gv$sig & gv$age_slope > 0
gv$down = gv$sig & gv$age_slope < 0
gv$key  = if (opt$by_ct_region) paste(gv$symbol, gv$cell_type, gv$region, sep='|') else gv$symbol
message('  gene_variance sig: ', sum(gv$sig))

# ── load DGLM results ─────────────────────────────────────────────────────
message('Loading DGLM results')
dglm.list = list()
cts = if (!is.null(opt$cell_type)) opt$cell_type else CELL_TYPES

for (ct in cts) {
    rds.file = file.path(opt$dglm_dir, paste0(ct, '_dglm_mashr_results.rds'))
    if (!file.exists(rds.file)) next
    obj          = readRDS(rds.file)
    dglm.results = obj$dglm_results
    human.sym    = obj$human_symbols
    regions      = obj$regions

    beta.mat = dglm.results[,'beta', regions, drop=FALSE][,1,, drop=FALSE]
    extreme  = apply(beta.mat, 1, function(x) any(abs(x) > 100, na.rm=TRUE))
    if (sum(extreme) > 0) {
        dglm.results = dglm.results[!extreme,,, drop=FALSE]
        human.sym    = human.sym[!extreme]
    }

    gene.names.file = file.path(opt$dglm_dir, paste0(ct, '_gene_names.csv'))
    gene.names.df   = if (file.exists(gene.names.file)) {
        read.csv(gene.names.file, row.names=1, stringsAsFactors=FALSE)
    } else NULL

    for (r in regions) {
        ensembl.ids = dimnames(dglm.results)[[1]]
        pval        = dglm.results[,'pval', r]
        beta        = dglm.results[,'beta', r]
        qval        = p.adjust(pval, method='fdr')
        symbols     = sapply(ensembl.ids, function(g) {
            if (!is.na(human.sym[g]) && human.sym[g] != g && human.sym[g] != '')
                return(human.sym[g])
            if (!is.null(gene.names.df) && g %in% rownames(gene.names.df)) {
                nm = gene.names.df[g, 'external_gene_name']
                if (!is.na(nm) && nm != g && nm != '') return(nm)
            }
            return(g)
        })

        dglm.list[[length(dglm.list)+1]] = data.frame(
            ensembl_id = ensembl.ids,
            symbol     = symbols,
            cell_type  = ct,
            region     = r,
            beta       = beta,
            qvalue     = qval,
            stringsAsFactors = FALSE
        )
    }
}

dglm = do.call(rbind, dglm.list)
dglm$sig  = !is.na(dglm$qvalue) & dglm$qvalue < opt$q_dglm
dglm$up   = dglm$sig & dglm$beta > 0
dglm$down = dglm$sig & dglm$beta < 0
dglm$key  = if (opt$by_ct_region) paste(dglm$symbol, dglm$cell_type, dglm$region, sep='|') else dglm$symbol
message('  DGLM sig: ', sum(dglm$sig))

# ── helper: make venn ─────────────────────────────────────────────────────
id_label = if (opt$by_ct_region) 'gene x cell type x region' else 'gene symbol'

make_venn = function(gv.sub, dglm.sub, title, fname) {
    gv.sig   = unique(gv.sub$key[gv.sub$sig])
    dglm.sig = unique(dglm.sub$key[dglm.sub$sig])

    if (length(gv.sig) == 0 && length(dglm.sig) == 0) {
        message('  no sig genes for: ', title)
        return(invisible(NULL))
    }

    overlap   = intersect(gv.sig, dglm.sig)
    n.overlap = length(overlap)

    if (n.overlap > 0) {
        gv.dir   = gv.sub[gv.sub$sig & gv.sub$key %in% overlap, c('key','up')]
        gv.dir   = gv.dir[!duplicated(gv.dir$key),]
        dglm.dir = dglm.sub[dglm.sub$sig & dglm.sub$key %in% overlap, c('key','up')]
        dglm.dir = dglm.dir[!duplicated(dglm.dir$key),]
        merged   = merge(gv.dir, dglm.dir, by='key', suffixes=c('.gv','.dglm'))
        if (nrow(merged) > 0) {
            n.concordant   = sum(merged$up.gv == merged$up.dglm, na.rm=TRUE)
            pct.concordant = round(100 * n.concordant / nrow(merged), 1)
            subtitle = paste0('overlap n=', n.overlap, ' (', id_label, ')',
                              ', concordant direction: ', pct.concordant, '%')
        } else {
            subtitle = paste0('overlap n=', n.overlap, ' (key mismatch in merge)')
        }
    } else {
        subtitle = paste0('no overlap (', id_label, ')')
    }

    venn.list = list(a=gv.sig, b=dglm.sig)
    names(venn.list) = c(
        paste0('gene_variance\n(q<', opt$q_gv, ')'),
        paste0('DGLM\n(q<', opt$q_dglm, ')')
    )

    p = ggvenn(venn.list,
               fill_color    = c('#4575b4','#d73027'),
               fill_alpha    = 0.4,
               stroke_size   = 0.5,
               set_name_size = 3.5,
               text_size     = 3.5) +
        labs(title=title, subtitle=subtitle) +
        theme(plot.title    = element_text(size=11, face='bold'),
              plot.subtitle = element_text(size=9))

    save_fig(p, fname, height=5, width=6)
}

# ── per cell type venns ───────────────────────────────────────────────────
message('\nMaking per cell type Venn diagrams')
suffix = if (opt$by_ct_region) '_ctregion' else ''
for (ct in cts) {
    gv.ct   = gv[gv$cell_type == ct,]
    dglm.ct = dglm[dglm$cell_type == ct,]
    if (nrow(gv.ct) == 0 && nrow(dglm.ct) == 0) next
    make_venn(gv.ct, dglm.ct,
              title = paste0(ct, ' — between-individual dispersion'),
              fname = paste0('venn_', ct, suffix))
}

# ── global venn ───────────────────────────────────────────────────────────
message('\nMaking global Venn diagram')
make_venn(gv, dglm,
          title = 'All cell types — between-individual dispersion',
          fname = paste0('venn_global', suffix))

# ── direction concordance summary ─────────────────────────────────────────
message('\nMaking direction concordance summary')
concordance.df = do.call(rbind, lapply(cts, function(ct) {
    gv.ct   = gv[gv$cell_type == ct & gv$sig,]
    dglm.ct = dglm[dglm$cell_type == ct & dglm$sig,]
    overlap = intersect(unique(gv.ct$key), unique(dglm.ct$key))
    if (length(overlap) == 0) return(NULL)

    gv.dir   = gv.ct[gv.ct$key %in% overlap, c('key','up')]
    gv.dir   = gv.dir[!duplicated(gv.dir$key),]
    dglm.dir = dglm.ct[dglm.ct$key %in% overlap, c('key','up')]
    dglm.dir = dglm.dir[!duplicated(dglm.dir$key),]
    merged   = merge(gv.dir, dglm.dir, by='key', suffixes=c('.gv','.dglm'))
    if (nrow(merged) == 0) return(NULL)

    data.frame(
        cell_type      = ct,
        n_overlap      = nrow(merged),
        n_concordant   = sum(merged$up.gv == merged$up.dglm, na.rm=TRUE),
        pct_concordant = 100 * sum(merged$up.gv == merged$up.dglm, na.rm=TRUE) / nrow(merged),
        stringsAsFactors = FALSE
    )
}))

if (!is.null(concordance.df) && nrow(concordance.df) > 0) {
    concordance.df$cell_type = factor(concordance.df$cell_type, levels=CELL_TYPES)

    p = ggplot(concordance.df, aes(cell_type, pct_concordant, fill=pct_concordant)) +
        geom_bar(stat='identity') +
        geom_text(aes(label=paste0('n=', n_overlap)), vjust=-0.3, size=3) +
        scale_fill_gradient(low='#f7f7f7', high='#2166ac', limits=c(0,100),
                            name='% concordant') +
        scale_y_continuous(limits=c(0,115), breaks=seq(0,100,25)) +
        geom_hline(yintercept=50, linetype='dashed', color='gray50', linewidth=0.5) +
        theme_classic(base_size=12) +
        theme(axis.text.x  = element_text(angle=45, hjust=1),
              axis.title.x = element_blank()) +
        ylab('% directional concordance') +
        ggtitle(paste0('Direction concordance in overlapping sig genes (', id_label, ')'))

    save_fig(p, paste0('concordance_by_celltype', suffix), height=5, width=8)
}

message('done.')
