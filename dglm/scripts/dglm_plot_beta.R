#!/usr/bin/env Rscript
# dglm_plot_beta.R
# Beta density plots and bar charts of sig genes per region (up vs down).
# Follows Chiou et al. 2022 Nat Neurosci dglm_visualize.R.
#
# Usage:
#   Rscript dglm_plot_beta.R --cell_type microglia

source('/scratch/easmit31/variability/dglm/scripts/_include_options.R')

library(optparse)
library(ggplot2)
library(reshape2)
library(mashr)

option_list = list(
    make_option('--cell_type',   type='character', help='cell type to analyze'),
    make_option('--checkpoints', type='character', default='/scratch/easmit31/variability/dglm/checkpoints'),
    make_option('--figdir',      type='character', default='/scratch/easmit31/variability/dglm/figures'),
    make_option('--outfmt',      type='character', default='png')
)
opt = parse_args(OptionParser(option_list=option_list))

cell.type = opt$cell_type
dir.create(opt$figdir, showWarnings=FALSE, recursive=TRUE)
message('Cell type: ', cell.type)

save_fig = function(p, name, height=5, width=7) {
    out = file.path(opt$figdir, paste0(cell.type, '_', name, '.', opt$outfmt))
    if (opt$outfmt == 'png') {
        ggsave(p, file=out, height=height, width=width, dpi=150)
    } else {
        ggsave(p, file=out, height=height, width=width, useDingbats=FALSE)
    }
    message('Saved: ', out)
}

# ── load results ──────────────────────────────────────────────────────────
in.file = file.path(opt$checkpoints, paste0(cell.type, '_dglm_mashr_results.rds'))
message('Loading ', in.file)
obj          = readRDS(in.file)
dglm.results = obj$dglm_results
regions      = obj$regions
m            = obj$mash

# ── filter extreme betas from dglm.results ────────────────────────────────
beta.mat = dglm.results[,'beta', regions, drop=FALSE][,1,, drop=FALSE]
if (length(regions) == 1) beta.mat = array(beta.mat, dim=c(dim(beta.mat)[1],1,1),
    dimnames=list(dimnames(dglm.results)[[1]], 'beta', regions))
extreme  = apply(beta.mat[,1,, drop=FALSE], 1, function(x) any(abs(x) > 100, na.rm=TRUE))
if (sum(extreme) > 0) {
    message('Filtering ', sum(extreme), ' genes with extreme beta (|beta|>100)')
    dglm.results = dglm.results[!extreme,,, drop=FALSE]
}

these.colors = region.colors[regions]

# ── extract beta and qval matrices ────────────────────────────────────────
dglm.beta = dglm.results[,'beta', regions, drop=FALSE][,1,]
dglm.qval = apply(
    matrix(dglm.results[,'pval', regions, drop=FALSE][,1,],
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
p = ggplot(subset(b.dglm, qval < fsr.cutoff), aes(beta, color=region)) +
    scale_color_manual(name='Region', values=these.colors) +
    geom_density() +
    theme_classic(base_size=12) +
    guides(color=guide_legend(ncol=2)) +
    xlab(expression(italic(beta))) +
    ylab('Density') +
    ggtitle(paste0(cell.type, ' — DGLM beta density (q<', fsr.cutoff, ')'))
save_fig(p, 'dglm_beta_density_dglm')

# ── DGLM sig gene counts per region ──────────────────────────────────────
b.dglm$qval.signed = with(b.dglm, ifelse(beta > 0, 1, -1) * qval)

counts.df = rbind(
    within(melt(tapply(b.dglm$qval.signed, b.dglm$region,
                       function(x) -sum(abs(x) < fsr.cutoff & x < 0, na.rm=TRUE))),
           { direction='Decrease' }),
    within(melt(tapply(b.dglm$qval.signed, b.dglm$region,
                       function(x) sum(abs(x) < fsr.cutoff & x >= 0, na.rm=TRUE))),
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
    ggtitle(paste0(cell.type, ' — DGLM sig genes per region (q<', fsr.cutoff, ')'))
save_fig(p, 'dglm_beta_count_dglm')

# ── mashr results (if available) ──────────────────────────────────────────
if (!is.null(m)) {
    mash.beta = get_pm(m)
    mash.lfsr = get_lfsr(m)
    mash.sbet = mash.beta / get_psd(m)

    b.mash = data.frame(
        expand.grid(gene=rownames(mash.beta), region=colnames(mash.beta)),
        qval = as.numeric(mash.lfsr),
        beta = as.numeric(mash.beta),
        sbet = as.numeric(mash.sbet),
        stringsAsFactors=FALSE
    )
    b.mash$region = factor(b.mash$region, levels=regions)

    # mashr beta density
    p = ggplot(subset(b.mash, qval < fsr.cutoff), aes(beta, color=region)) +
        scale_color_manual(name='Region', values=these.colors) +
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
    mash.counts.df$Var1      = factor(mash.counts.df$Var1, levels=regions)
    mash.counts.df$direction = factor(mash.counts.df$direction, levels=c('Increase','Decrease'))
    ylimit = ceiling(max(abs(mash.counts.df$value)) / 100) * 100
    if (ylimit == 0) ylimit = 10

    p = ggplot(mash.counts.df, aes(Var1, value, fill=Var1, alpha=direction)) +
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
        ggtitle(paste0(cell.type, ' — mashr sig genes per region (lfsr<', fsr.cutoff, ')'))
    save_fig(p, 'dglm_beta_count_mash')

    message('mashr sig at lfsr<', fsr.cutoff, ': ', sum(mash.lfsr < fsr.cutoff))
} else {
    message('No mashr results — skipping mashr plots')
}

message('DGLM sig at q<', fsr.cutoff, ': ', sum(dglm.qval < fsr.cutoff, na.rm=TRUE))
message('done.')
