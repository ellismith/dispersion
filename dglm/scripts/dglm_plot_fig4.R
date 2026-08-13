#!/usr/bin/env Rscript
# dglm_plot_fig4.R
# Recreates Chiou et al. Fig 4 panels a and b using DGLM results.
#
# Panel a: bar chart per region, single dominant direction bar
# Panel b: per-animal z-scores for one example gene across regions as facets
#
# Panel b's example gene/cell type auto-selects (top |mash_beta| among genes
# significant at --sig_col<--qthresh) when --example_symbol isn't passed --
# no separate manual awk/grep step needed before running this script.
#
# Usage:
#   Rscript dglm_plot_fig4.R                                   # auto-picks gene+cell type
#   Rscript dglm_plot_fig4.R --example_ct astrocytes            # auto-picks gene within astrocytes
#   Rscript dglm_plot_fig4.R --example_symbol SERPINA3          # auto-picks best cell type for this gene
#   Rscript dglm_plot_fig4.R --example_symbol SERPINA3 --example_ct astrocytes  # fully manual

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(ggplot2)

option_list = list(
    make_option('--checkpoints',    type='character', default='/scratch/easmit31/dispersion/dglm/checkpoints'),
    make_option('--master_tsv',     type='character', default='/scratch/easmit31/dispersion/dglm/checkpoints/master_dglm_combined.tsv'),
    make_option('--figdir',         type='character', default='/scratch/easmit31/dispersion/dglm/figures'),
    make_option('--outfmt',         type='character', default='png'),
    make_option('--by_celltype',    action='store_true', default=FALSE),
    make_option('--min_age',        type='double',    default=1.0),
    make_option('--sig_col',        type='character', default='mash_lfsr',
                help='column to threshold significance on: mash_lfsr or qvalue'),
    make_option('--qthresh',        type='double',    default=0.2,
                help='significance threshold — 0.2 for mash_lfsr (Chiou et al.), 0.05 typical for qvalue'),
    make_option('--example_symbol', type='character', default=NULL,
                help='panel b example gene symbol. If unset, auto-picked as the gene with the largest |mash_beta| among genes significant at --sig_col<--qthresh (within --example_ct if that is also set, else across the whole dataset).'),
    make_option('--example_ct',     type='character', default=NULL,
                help='panel b example cell type. If unset, auto-picked alongside --example_symbol -- or, if --example_symbol IS set but this is not, picked as whichever cell type that gene is most significant in.'),
    make_option('--panel_a_only',   action='store_true', default=FALSE),
    make_option('--panel_b_only',   action='store_true', default=FALSE),
    make_option('--regions',        type='character', default=NULL,
                help='comma-separated region list for panel b (e.g. "CN,dlPFC,HIP"); overrides --sig_only_regions'),
    make_option('--all_regions',    action='store_true', default=FALSE,
                help='panel b: show all regions instead of only the ones significant for this gene+cell_type')
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
# Always loaded now (previously skipped under --panel_b_only) -- panel b's
# auto-selected example gene needs the full significance table regardless
# of whether panel a is being rendered, and this avoids a second, separate
# reload later that duplicated this same logic.
message('Loading master TSV: ', opt$master_tsv)
all.df = read.table(opt$master_tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
if (!opt$sig_col %in% colnames(all.df))
    stop('Column not found: ', opt$sig_col, ' — check --master_tsv has this column (mash_lfsr requires the combined mashr master TSV)')
if (!'mash_beta' %in% colnames(all.df))
    stop('Column not found: mash_beta — required for auto-selecting the panel b example gene')
all.df = all.df[!is.na(all.df[[opt$sig_col]]) & abs(all.df$beta) <= 100,]
all.df$region    = factor(all.df$region,    levels=region.levels)
# cell.type.levels only covers the 12 canonical cell types -- louvain
# subcluster labels (e.g. microglia_0..microglia_16) aren't in it, so
# factor() would silently turn them all into NA and collapse every
# subcluster into one empty facet. Fall back to the data's own sorted
# values whenever it isn't a pure subset of the canonical list.
effective.ct.levels = if (all(unique(all.df$cell_type) %in% cell.type.levels)) {
    cell.type.levels
} else {
    message('cell_type values not fully covered by the canonical 12 -- ',
            'using data-derived levels instead (likely louvain-resolution input)')
    sort(unique(all.df$cell_type))
}
all.df$cell_type = factor(all.df$cell_type, levels=effective.ct.levels)
message('Total tests: ', nrow(all.df))
message('Sig at ', opt$sig_col, '<', opt$qthresh, ': ',
        sum(all.df[[opt$sig_col]] < opt$qthresh, na.rm=TRUE))

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
            n_inc = sum(df[[opt$sig_col]][idx] < opt$qthresh & df$beta[idx] > 0, na.rm=TRUE)
            n_dec = sum(df[[opt$sig_col]][idx] < opt$qthresh & df$beta[idx] < 0, na.rm=TRUE)
            net   = ifelse(n_inc >= n_dec, n_inc, -n_dec)
            data.frame(groups[i,, drop=FALSE], n_inc=n_inc, n_dec=n_dec, net=net,
                       stringsAsFactors=FALSE)
        }))

        counts$region = factor(counts$region, levels=region.levels)
        if (opt$by_celltype) counts$cell_type = factor(counts$cell_type, levels=effective.ct.levels)

        # round up to a "nice" step scaled to the data's own magnitude, not a
        # fixed step of 100 -- with counts in the tens, a fixed-100 step makes
        # every panel look nearly empty regardless of by_celltype
        nice_ylimit = function(x) {
            m = max(abs(x), na.rm=TRUE)
            if (!is.finite(m) || m == 0) return(10)
            step = 10 ^ floor(log10(m))
            if (m / step > 5) step = step * 2
            ceiling(m / step) * step
        }

        p = ggplot(counts, aes(region, net, fill=region)) +
            geom_bar(stat='identity') +
            scale_fill_manual(values=region.colors) +
            geom_hline(yintercept=0, color='black', linewidth=0.3) +
            theme_classic(base_size=12) +
            theme(legend.position='none',
                  axis.text.x=element_text(angle=-45, hjust=0, vjust=1),
                  axis.title.x=element_blank()) +
            ylab('Number of genes') +
            ggtitle(title)

        if (opt$by_celltype) {
            # free_y: each cell type's panel scales to its own data, with the
            # "increased/decreased" labels placed relative to each panel's own
            # range (Inf/-Inf + vjust) so they sit correctly under free scales
            p = p +
                geom_text(data=unique(counts['cell_type']), x=Inf, y=Inf,
                          label='Increased dispersion with age',
                          size=2.8, fontface='italic', hjust=1.05, vjust=1.5,
                          inherit.aes=FALSE) +
                geom_text(data=unique(counts['cell_type']), x=Inf, y=-Inf,
                          label='Decreased dispersion with age',
                          size=2.8, fontface='italic', hjust=1.05, vjust=-0.8,
                          inherit.aes=FALSE) +
                facet_wrap(~cell_type, scales='free_y',
                           nrow=ceiling(sqrt(length(effective.ct.levels))))
            width = 14; height = 10
        } else {
            ylimit = nice_ylimit(counts$net)
            p = p +
                scale_y_continuous(
                    limits = c(-ylimit, ylimit),
                    breaks = c(-ylimit, -ylimit/2, 0, ylimit/2, ylimit)
                ) +
                annotate('text', x=length(region.levels)*0.75, y=ylimit*0.85,
                         label='Increased dispersion with age', size=3.5, fontface='italic') +
                annotate('text', x=length(region.levels)*0.75, y=-ylimit*0.85,
                         label='Decreased dispersion with age', size=3.5, fontface='italic')
        }
        save_fig(p, fname, height=height, width=width)
    }

    suffix = if (opt$by_celltype) 'fig4a_by_celltype' else 'fig4a'
    sig_label = if (opt$sig_col == 'mash_lfsr') 'mashr lfsr' else 'global FDR'
    make_panel_a(all.df,
        title = paste0('DGLM: variance changes with age (', sig_label, '<', opt$qthresh, ')'),
        fname = paste0(suffix, '_', opt$sig_col, '_q', opt$qthresh))
}

# ═══════════════════════════════════════════════════════════════════════════
# PANEL B
# ═══════════════════════════════════════════════════════════════════════════

if (!opt$panel_a_only) {
    # ── auto-select example gene/cell type if not explicitly given ────────
    if (is.null(opt$example_symbol)) {
        sig.rows = all.df[!is.na(all.df[[opt$sig_col]]) & all.df[[opt$sig_col]] < opt$qthresh, ]
        if (!is.null(opt$example_ct)) sig.rows = sig.rows[sig.rows$cell_type == opt$example_ct, ]
        if (nrow(sig.rows) == 0) {
            stop('No genes significant at ', opt$sig_col, '<', opt$qthresh,
                 if (!is.null(opt$example_ct)) paste0(' in cell_type=', opt$example_ct) else '',
                 ' -- cannot auto-select a panel b example gene. Pass --example_symbol manually, or ',
                 'loosen --qthresh.')
        }
        top.row     = sig.rows[which.max(abs(sig.rows$mash_beta)), ]
        gene.symbol = as.character(top.row$symbol)
        top.ct      = if (!is.null(opt$example_ct)) opt$example_ct else as.character(top.row$cell_type)
        message('Auto-selected panel b example: ', gene.symbol, ' (', top.ct, '), |mash_beta|=',
                round(abs(top.row$mash_beta), 3), ', ', opt$sig_col, '=', signif(top.row[[opt$sig_col]], 3))
    } else {
        gene.symbol = opt$example_symbol
        if (!is.null(opt$example_ct)) {
            top.ct = opt$example_ct
        } else {
            g.rows = all.df[all.df$symbol == gene.symbol & !is.na(all.df[[opt$sig_col]]) &
                             all.df[[opt$sig_col]] < opt$qthresh, ]
            if (nrow(g.rows) == 0) {
                stop('Gene ', gene.symbol, ' is not significant at ', opt$sig_col, '<', opt$qthresh,
                     ' in any cell type -- pass --example_ct manually, or loosen --qthresh.')
            }
            top.ct = as.character(g.rows$cell_type[which.max(abs(g.rows$mash_beta))])
            message('Auto-selected cell type for ', gene.symbol, ': ', top.ct)
        }
    }
    message('Panel b: ', gene.symbol, ' in ', top.ct)

    # find ensembl id from master TSV
    gene.row = all.df[all.df$symbol == gene.symbol,]
    if (nrow(gene.row) == 0) {
        message('Gene ', gene.symbol, ' not found in master TSV — skipping panel b')
        quit(save='no', status=0)
    }
    gene.id = unique(gene.row$ensembl_id)[1]
    message('Ensembl ID: ', gene.id)

    # decide which regions to show:
    #   --regions overrides everything (manual pick)
    #   --all_regions shows every region in region.levels
    #   default: only regions significant for this gene x cell_type, matching
    #     how the paper's own Fig 4b only shows a handful of regions, not all
    if (!is.null(opt$regions)) {
        plot.regions = trimws(strsplit(opt$regions, ',')[[1]])
        message('Regions (manual): ', paste(plot.regions, collapse=', '))
    } else if (opt$all_regions) {
        plot.regions = region.levels
        message('Regions (all): ', paste(plot.regions, collapse=', '))
    } else {
        ct.rows = gene.row[gene.row$cell_type == top.ct & !is.na(gene.row[[opt$sig_col]]), ]
        plot.regions = unique(ct.rows$region[ct.rows[[opt$sig_col]] < opt$qthresh])
        if (length(plot.regions) == 0) {
            message('No significant regions for ', gene.symbol, ' in ', top.ct,
                    ' at ', opt$sig_col, '<', opt$qthresh,
                    ' -- falling back to all regions (pass --regions to pick manually)')
            plot.regions = region.levels
        } else {
            message('Regions (sig at ', opt$sig_col, '<', opt$qthresh, '): ',
                    paste(plot.regions, collapse=', '))
        }
    }

    # load pseudobulk
    pb.data = list()
    for (r in plot.regions) {
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

    panel.b.data$region    = factor(panel.b.data$region,    levels=plot.regions)
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
    panel.b.summary$region    = factor(panel.b.summary$region,    levels=plot.regions)
    panel.b.summary$age_group = factor(panel.b.summary$age_group, levels=c('Young','Old'))

    zlim = max(abs(panel.b.data$zscore), na.rm=TRUE)

    p = ggplot(panel.b.data, aes(age_group, zscore)) +
        geom_jitter(width=0.2, height=0, size=1, alpha=0.6, color='#555555') +
        geom_point(data=panel.b.summary, aes(x=age_group, y=mean_z),
                   inherit.aes=FALSE, size=3, color='#222222') +
        geom_errorbar(data=panel.b.summary,
                      aes(x=age_group, ymin=mean_z-sd_z, ymax=mean_z+sd_z),
                      inherit.aes=FALSE, width=0.2, linewidth=0.5, color='#222222') +
        facet_wrap(~region, nrow=ceiling(sqrt(length(plot.regions)))) +
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
