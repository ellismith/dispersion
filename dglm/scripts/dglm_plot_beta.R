#!/usr/bin/env Rscript
# dglm_plot_beta.R
# Beta density plots and bar charts of sig genes per region (up vs down).
# Follows Chiou et al. 2022 Nat Neurosci dglm_visualize.R.
#
# DGLM-only plots read straight from {cell_type}_dglm_results.rds (raw
# dglm_model.R output) -- no mashr dependency at all.
# Mashr plots read from the COMBINED master TSV (master_dglm_combined.tsv),
# filtered to this cell type -- no per-cell-type mashr rerun required.
#
# Usage:
#   Rscript dglm_plot_beta.R --cell_type microglia
#   Rscript dglm_plot_beta.R --cell_type microglia --term age   # two-term (age_sex) RDS

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(ggplot2)
library(reshape2)

option_list = list(
    make_option('--cell_type',    type='character', help='cell type to analyze'),
    make_option('--checkpoints',  type='character', default='/scratch/easmit31/dispersion/dglm/checkpoints'),
    make_option('--master_tsv',   type='character', default='/scratch/easmit31/dispersion/dglm/checkpoints/master_dglm_combined.tsv',
                help='combined mashr master TSV -- source for the mashr plots'),
    make_option('--figdir',       type='character', default='/scratch/easmit31/dispersion/dglm/figures'),
    make_option('--outfmt',       type='character', default='png'),
    make_option('--dglm_qthresh', type='double',    default=0.05,
                help='plain per-region BH-FDR threshold for the DGLM-only plots (NOT fsr.cutoff -- that is mashr LFSR-only, matches the 0.05 convention used elsewhere, e.g. dglm_plot_heatmap.R/gene_variance)'),
    make_option('--term',         type='character', default=NULL,
                help='For two-term (age_sex) DGLM results: which term to plot (age|sex). Selects beta_<term>/pval_<term> columns instead of bare beta/pval. Ignored for single-term results.')
)
opt = parse_args(OptionParser(option_list=option_list))

cell.type = opt$cell_type
dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)
message('Cell type: ', cell.type)

beta.col = if (!is.null(opt$term)) paste0('beta_', opt$term) else 'beta'
pval.col = if (!is.null(opt$term)) paste0('pval_', opt$term) else 'pval'
if (!is.null(opt$term)) message('Term: ', opt$term, ' (columns: ', beta.col, ', ', pval.col, ')')

save_fig = function(p, name, height=5, width=7) {
    out = file.path(opt$figdir, paste0(cell.type, '_', name, '.', opt$outfmt))
    if (opt$outfmt == 'png') {
        ggsave(p, file=out, height=height, width=width, dpi=150)
    } else {
        ggsave(p, file=out, height=height, width=width, useDingbats=FALSE)
    }
    message('Saved: ', out)
}

# ═══════════════════════════════════════════════════════════════════════════
# DGLM-ONLY PLOTS -- straight from raw dglm_model.R output, no mashr needed
# ═══════════════════════════════════════════════════════════════════════════

in.file = file.path(opt$checkpoints, paste0(cell.type, '_dglm_results.rds'))
message('Loading ', in.file)
obj          = readRDS(in.file)
dglm.results = obj$array

# regions with actual data for this cell type (raw RDS has no separate
# 'regions' field -- derive it the same way dglm_fdr.R's raw-RDS fallback does)
regions = region.levels[region.levels %in% dimnames(dglm.results)[[3]]]
regions = regions[apply(dglm.results[, beta.col, regions, drop=FALSE], 3,
                        function(x) sum(!is.na(x) & x != 0) > 0)]
message('Regions with data: ', paste(regions, collapse=', '))

# ── filter extreme betas ──────────────────────────────────────────────────
beta.mat = dglm.results[, beta.col, regions, drop=FALSE][,1,, drop=FALSE]
if (length(regions) == 1) beta.mat = array(beta.mat, dim=c(dim(beta.mat)[1],1,1),
    dimnames=list(dimnames(dglm.results)[[1]], beta.col, regions))
extreme  = apply(beta.mat[,1,, drop=FALSE], 1, function(x) any(abs(x) > 100, na.rm=TRUE))
if (sum(extreme) > 0) {
    message('Filtering ', sum(extreme), ' genes with extreme beta (|beta|>100)')
    dglm.results = dglm.results[!extreme,,, drop=FALSE]
}

these.colors = region.colors[regions]

# ── extract beta and qval matrices ────────────────────────────────────────
dglm.beta = dglm.results[,beta.col, regions, drop=FALSE][,1,]
dglm.qval = apply(
    matrix(dglm.results[,pval.col, regions, drop=FALSE][,1,],
           nrow=nrow(dglm.results), ncol=length(regions),
           dimnames=list(dimnames(dglm.results)[[1]], regions)),
    2, function(x) p.adjust(x, method='fdr')
)
if (is.null(dim(dglm.beta))) {
    dglm.beta = matrix(dglm.beta, ncol=1, dimnames=list(dimnames(dglm.results)[[1]], regions))
    dglm.qval = matrix(dglm.qval, ncol=1, dimnames=list(dimnames(dglm.results)[[1]], regions))
}

# ── longform DGLM data frame ──────────────────────────────────────────────
b.dglm = data.frame(
    expand.grid(gene=rownames(dglm.beta), region=regions),
    beta = as.numeric(dglm.beta),
    qval = as.numeric(dglm.qval),
    stringsAsFactors=FALSE
)
b.dglm$region = factor(b.dglm$region, levels=regions)

# ── DGLM beta density ─────────────────────────────────────────────────────
p = ggplot(subset(b.dglm, qval < opt$dglm_qthresh), aes(beta, color=region)) +
    scale_color_manual(name='Region', values=these.colors) +
    geom_density() +
    theme_classic(base_size=12) +
    guides(color=guide_legend(ncol=2)) +
    xlab(expression(italic(beta))) +
    ylab('Density') +
    ggtitle(paste0(cell.type, ' — DGLM beta density (q<', opt$dglm_qthresh, ')'))
save_fig(p, 'dglm_beta_density_dglm')

# ── DGLM sig gene counts per region ──────────────────────────────────────
b.dglm$qval.signed = with(b.dglm, ifelse(beta > 0, 1, -1) * qval)

counts.df = rbind(
    within(melt(tapply(b.dglm$qval.signed, b.dglm$region,
                       function(x) -sum(abs(x) < opt$dglm_qthresh & x < 0, na.rm=TRUE))),
           { direction='Decrease' }),
    within(melt(tapply(b.dglm$qval.signed, b.dglm$region,
                       function(x) sum(abs(x) < opt$dglm_qthresh & x >= 0, na.rm=TRUE))),
           { direction='Increase' })
)
counts.df$Var1      = factor(counts.df$Var1, levels=regions)
counts.df$direction = factor(counts.df$direction, levels=c('Increase','Decrease'))
ylimit = ceiling(max(abs(counts.df$value)) / 100) * 100
if (ylimit == 0) ylimit = 10

p = ggplot(counts.df, aes(Var1, value, fill=Var1, alpha=direction)) +
    geom_bar(stat='identity') +
    scale_fill_manual(name='Region', values=these.colors) +
    scale_alpha_manual(values=c(1, 0.75)) +
    scale_y_continuous(
        limits = c(-ylimit, ylimit),
        breaks = c(-ylimit, -ylimit*0.5, 0, ylimit*0.5, ylimit),
        labels = c(formatC(ylimit, width=5, flag=' '), 'Decrease',
                   formatC(0, width=5, flag=' '), 'Increase',
                   formatC(ylimit, width=5, flag=' '))
    ) +
    theme_classic(base_size=12) +
    theme(
        legend.position = 'none',
        axis.text.x     = element_text(angle=-45, hjust=0, vjust=1),
        axis.title.x    = element_blank()
    ) +
    ylab('Number of genes') +
    ggtitle(paste0(cell.type, ' — DGLM sig genes per region (q<', opt$dglm_qthresh, ')'))
save_fig(p, 'dglm_beta_count_dglm')

# ═══════════════════════════════════════════════════════════════════════════
# MASHR PLOTS -- from the combined master TSV, filtered to this cell type
# ═══════════════════════════════════════════════════════════════════════════

if (file.exists(opt$master_tsv)) {
    message('Loading master TSV: ', opt$master_tsv)
    master.all = read.table(opt$master_tsv, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    b.mash     = master.all[master.all$cell_type == cell.type & !is.na(master.all$mash_lfsr),]

    if (nrow(b.mash) == 0) {
        message('No combined mashr rows for cell_type=', cell.type, ' — skipping mashr plots')
    } else {
        mash.regions = region.levels[region.levels %in% unique(b.mash$region)]
        b.mash$region = factor(b.mash$region, levels=mash.regions)
        b.mash$qval   = b.mash$mash_lfsr
        b.mash$beta   = b.mash$mash_beta
        these.colors.mash = region.colors[mash.regions]

        # mashr beta density
        p = ggplot(subset(b.mash, qval < fsr.cutoff), aes(beta, color=region)) +
            scale_color_manual(name='Region', values=these.colors.mash) +
            geom_density() +
            theme_classic(base_size=12) +
            guides(color=guide_legend(ncol=2)) +
            xlab(expression(italic(beta))) +
            ylab('Density') +
            ggtitle(paste0(cell.type, ' — mashr beta density (lfsr<', fsr.cutoff, ')'))
        save_fig(p, 'dglm_beta_density_mash')

        # mashr sig gene counts
        b.mash$qval.signed = with(b.mash, ifelse(beta > 0, 1, -1) * qval)
        mash.counts.df = rbind(
            within(melt(tapply(b.mash$qval.signed, b.mash$region,
                               function(x) -sum(abs(x) < fsr.cutoff & x < 0, na.rm=TRUE))),
                   { direction='Decrease' }),
            within(melt(tapply(b.mash$qval.signed, b.mash$region,
                               function(x) sum(abs(x) < fsr.cutoff & x >= 0, na.rm=TRUE))),
                   { direction='Increase' })
        )
        mash.counts.df$Var1      = factor(mash.counts.df$Var1, levels=mash.regions)
        mash.counts.df$direction = factor(mash.counts.df$direction, levels=c('Increase','Decrease'))
        ylimit = ceiling(max(abs(mash.counts.df$value)) / 100) * 100
        if (ylimit == 0) ylimit = 10

        p = ggplot(mash.counts.df, aes(Var1, value, fill=Var1, alpha=direction)) +
            geom_bar(stat='identity') +
            scale_fill_manual(name='Region', values=these.colors.mash) +
            scale_alpha_manual(values=c(1, 0.75)) +
            scale_y_continuous(
                limits = c(-ylimit, ylimit),
                breaks = c(-ylimit, -ylimit*0.5, 0, ylimit*0.5, ylimit),
                labels = c(formatC(ylimit, width=5, flag=' '), 'Decrease',
                           formatC(0, width=5, flag=' '), 'Increase',
                           formatC(ylimit, width=5, flag=' '))
            ) +
            theme_classic(base_size=12) +
            theme(
                legend.position = 'none',
                axis.text.x     = element_text(angle=-45, hjust=0, vjust=1),
                axis.title.x    = element_blank()
            ) +
            ylab('Number of genes') +
            ggtitle(paste0(cell.type, ' — mashr sig genes per region (lfsr<', fsr.cutoff, ')'))
        save_fig(p, 'dglm_beta_count_mash')

        message('mashr sig at lfsr<', fsr.cutoff, ': ', sum(b.mash$mash_lfsr < fsr.cutoff, na.rm=TRUE))
    }
} else {
    message('Master TSV not found (', opt$master_tsv, ') — skipping mashr plots')
}

message('DGLM sig at q<', opt$dglm_qthresh, ': ', sum(dglm.qval < opt$dglm_qthresh, na.rm=TRUE))
message('done.')
