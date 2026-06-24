#!/usr/bin/env Rscript
# dglm_plot_fig4.R
# Recreates Chiou et al. Fig 4 panels a and b using DGLM results.
#
# Panel a: bar chart per region, single dominant direction bar
# Panel b: per-animal z-scores for one example gene across regions as facets
#
# Usage:
#   Rscript dglm_plot_fig4.R
#   Rscript dglm_plot_fig4.R --example_symbol SERPINA3 --example_ct astrocytes

source('/scratch/easmit31/variability/dglm/scripts/_include_options.R')

library(optparse)
library(ggplot2)

option_list = list(
    make_option('--checkpoints',    type='character', default='/scratch/easmit31/variability/dglm/checkpoints'),
    make_option('--master_tsv',     type='character', default='/scratch/easmit31/variability/dglm/checkpoints/master_dglm_globalfdr.tsv'),
    make_option('--figdir',         type='character', default='/scratch/easmit31/variability/dglm/figures'),
    make_option('--outfmt',         type='character', default='png'),
    make_option('--by_celltype',    action='store_true', default=FALSE),
    make_option('--min_age',        type='double',    default=1.0),
    make_option('--qthresh',        type='double',    default=0.05),
    make_option('--example_symbol', type='character', default='SERPINA3'),
    make_option('--example_ct',     type='character', default='astrocytes'),
    make_option('--panel_a_only',   action='store_true', default=FALSE),
    make_option('--panel_b_only',   action='store_true', default=FALSE)
)
opt = parse_args(OptionParser(option_list=option_list))

dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)

save_fig = function(p, name, height=5, width=8) {
    out = file.path(opt$figdir, paste0(name, '.', opt$outfmt))
    if (opt$outfmt == 'png') {
        ggsave(p, file=out, height=height, width=width, dpi=150)
    } else {
        ggsave(p, file=out, height=height, width=width, useDingbats=FALSE)
    }
    message('Saved: ', out)
}

# ── load master TSV ───────────────────────────────────────────────────────
if (!opt$panel_b_only) {
    message('Loading master TSV: ', opt$master_tsv)
    all.df = read.table(opt$master_tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    all.df = all.df[!is.na(all.df$qvalue) & abs(all.df$beta) <= 100,]
    all.df$region    = factor(all.df$region,    levels=region.levels)
    all.df$cell_type = factor(all.df$cell_type, levels=cell.type.levels)
    message('Total tests: ', nrow(all.df))
    message('Sig at q<', opt$qthresh, ': ', sum(all.df$qvalue < opt$qthresh, na.rm=TRUE))
}

# ═══════════════════════════════════════════════════════════════════════════
# PANEL A
# ═══════════════════════════════════════════════════════════════════════════

if (!opt$panel_b_only) {
    make_panel_a = function(df, title, fname, width=8, height=5) {
        group.vars = if (opt$by_celltype) c('region','cell_type') else 'region'
        groups     = unique(df[, group.vars, drop=FALSE])

        counts = do.call(rbind, lapply(1:nrow(groups), function(i) {
            mask = rep(TRUE, nrow(df))
            for (v in group.vars) mask = mask & df[[v]] == groups[i, v]
            idx   = which(mask)
            n_inc = sum(df$qvalue[idx] < opt$qthresh & df$beta[idx] > 0, na.rm=TRUE)
            n_dec = sum(df$qvalue[idx] < opt$qthresh & df$beta[idx] < 0, na.rm=TRUE)
            net   = ifelse(n_inc >= n_dec, n_inc, -n_dec)
            data.frame(groups[i,, drop=FALSE], n_inc=n_inc, n_dec=n_dec, net=net,
                       stringsAsFactors=FALSE)
        }))

        counts$region = factor(counts$region, levels=region.levels)
        if (opt$by_celltype) counts$cell_type = factor(counts$cell_type, levels=cell.type.levels)

        ylimit = ceiling(max(abs(counts$net), na.rm=TRUE) / 100) * 100
        if (ylimit == 0) ylimit = 10

        p = ggplot(counts, aes(region, net, fill=region)) +
            geom_bar(stat='identity') +
            scale_fill_manual(values=region.colors) +
            scale_y_continuous(
                limits = c(-ylimit, ylimit),
                breaks = c(-ylimit, -ylimit/2, 0, ylimit/2, ylimit)
            ) +
            geom_hline(yintercept=0, color='black', linewidth=0.3) +
            annotate('text', x=length(region.levels)*0.75, y=ylimit*0.85,
                     label='Increased dispersion with age', size=3.5, fontface='italic') +
            annotate('text', x=length(region.levels)*0.75, y=-ylimit*0.85,
                     label='Decreased dispersion with age', size=3.5, fontface='italic') +
            theme_classic(base_size=12) +
            theme(legend.position='none',
                  axis.text.x=element_text(angle=-45, hjust=0, vjust=1),
                  axis.title.x=element_blank()) +
            ylab('Number of genes') +
            ggtitle(title)

        if (opt$by_celltype) {
            p     = p + facet_wrap(~cell_type, scales='free_y',
                                   nrow=ceiling(sqrt(length(cell.type.levels))))
            width = 14; height = 10
        }
        save_fig(p, fname, height=height, width=width)
    }

    suffix = if (opt$by_celltype) 'fig4a_by_celltype' else 'fig4a'
    make_panel_a(all.df,
        title = paste0('DGLM: variance changes with age (q<', opt$qthresh, ')'),
        fname = suffix)
}

# ═══════════════════════════════════════════════════════════════════════════
# PANEL B
# ═══════════════════════════════════════════════════════════════════════════

if (!opt$panel_a_only) {
    gene.symbol = opt$example_symbol
    top.ct      = opt$example_ct
    message('Panel b: ', gene.symbol, ' in ', top.ct)

    # find ensembl id from master TSV
    if (!opt$panel_b_only) {
        gene.row = all.df[all.df$symbol == gene.symbol,]
    } else {
        all.df   = read.table(opt$master_tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
        gene.row = all.df[all.df$symbol == gene.symbol,]
    }
    if (nrow(gene.row) == 0) {
        message('Gene ', gene.symbol, ' not found in master TSV — skipping panel b')
        quit(save='no', status=0)
    }
    gene.id = unique(gene.row$ensembl_id)[1]
    message('Ensembl ID: ', gene.id)

    # load pseudobulk
    pb.data = list()
    for (r in region.levels) {
        pb.file = file.path(opt$checkpoints, paste0(top.ct, '_', r, '_pseudobulk.csv'))
        if (!file.exists(pb.file)) next
        pb = tryCatch(read.csv(pb.file, row.names=1, check.names=FALSE), error=function(e) NULL)
        if (is.null(pb) || nrow(pb) == 0 || !gene.id %in% rownames(pb)) next
        pb.data[[r]] = pb
    }

    if (length(pb.data) == 0) {
        message('Gene not found in any pseudobulk file — skipping panel b')
        quit(save='no', status=0)
    }

    meta.all           = read.csv(file.path(opt$checkpoints, paste0(top.ct, '_metadata.csv')),
                                   stringsAsFactors=FALSE)
    meta.all           = meta.all[meta.all$age >= opt$min_age,]
    median.age         = median(unique(meta.all$age))
    meta.all$age_group = factor(ifelse(meta.all$age <= median.age, 'Young', 'Old'),
                                levels=c('Young','Old'))

    panel.b.data = do.call(rbind, lapply(names(pb.data), function(r) {
        pb     = pb.data[[r]]
        expr   = as.numeric(pb[gene.id,])
        names(expr) = colnames(pb)
        z      = (expr - mean(expr, na.rm=TRUE)) / sd(expr, na.rm=TRUE)
        meta.r = meta.all[meta.all$region == r & meta.all$animal_id %in% names(z),]
        if (nrow(meta.r) == 0) return(NULL)
        data.frame(region=r, animal_id=meta.r$animal_id, age=meta.r$age,
                   age_group=meta.r$age_group, zscore=z[meta.r$animal_id],
                   stringsAsFactors=FALSE)
    }))

    panel.b.data$region    = factor(panel.b.data$region,    levels=region.levels)
    panel.b.data$age_group = factor(panel.b.data$age_group, levels=c('Young','Old'))

    panel.b.summary = do.call(rbind, lapply(
        split(panel.b.data, list(panel.b.data$region, panel.b.data$age_group)),
        function(x) {
            if (nrow(x) == 0) return(NULL)
            data.frame(region=unique(x$region), age_group=unique(x$age_group),
                       mean_z=mean(x$zscore, na.rm=TRUE), sd_z=sd(x$zscore, na.rm=TRUE),
                       stringsAsFactors=FALSE)
        }
    ))
    panel.b.summary$region    = factor(panel.b.summary$region,    levels=region.levels)
    panel.b.summary$age_group = factor(panel.b.summary$age_group, levels=c('Young','Old'))

    zlim = max(abs(panel.b.data$zscore), na.rm=TRUE)

    p = ggplot(panel.b.data, aes(age_group, zscore)) +
        geom_jitter(width=0.2, height=0, size=1, alpha=0.6, color='#555555') +
        geom_point(data=panel.b.summary, aes(x=age_group, y=mean_z),
                   inherit.aes=FALSE, size=3, color='#222222') +
        geom_errorbar(data=panel.b.summary,
                      aes(x=age_group, ymin=mean_z-sd_z, ymax=mean_z+sd_z),
                      inherit.aes=FALSE, width=0.2, linewidth=0.5, color='#222222') +
        facet_wrap(~region, nrow=2) +
        scale_y_continuous(limits=c(-zlim, zlim),
                           breaks=seq(-floor(zlim), floor(zlim), 2)) +
        xlab(predictor.label) +
        ylab(expression(italic(Z)*'-score')) +
        ggtitle(paste0(gene.symbol, ' — ', top.ct)) +
        theme_classic(base_size=12) +
        theme(strip.background=element_blank())

    # output named by gene and cell type
    fname = paste0('fig4b_', gene.symbol, '_', top.ct)
    save_fig(p, fname, height=5, width=10)
}

message('done.')
