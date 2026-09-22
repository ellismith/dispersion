#!/usr/bin/env Rscript
source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')
library(optparse)
library(mashr)
option_list = list(
    make_option('--mode', type='character', default='combined'),
    make_option('--checkpoints', type='character',
                default='/scratch/easmit31/dispersion/dglm/checkpoints'),
    make_option('--shat_mode', type='character', default='sqrt'),
    make_option('--na_shat_fill', type='double', default=1000,
                help='Shat value used for gene x condition combos with no real data (default 1000, deliberately uninformative)'),
    make_option('--strong_q', type='double', default=0.05,
                help='q threshold for qval_single strong subset'),
    make_option('--approach', type='character', default='qval_recurrence',
                help='strong-subset selection + fit structure. One of:
  qval_third      : DGLM q<0.05 in >=1/3 of all conditions; two-stage fit
                     (random-subset null-correlation estimate, canonical +
                     data-driven covariances fit together, applied to all genes)
  qval_recurrence : DGLM q<0.01 in >=10 of all conditions; two-stage fit (as above)
  qval_tenth      : DGLM q<0.05 in >=1/10 of all conditions; single-stage fit
                     (data-driven covariances only, fit directly on all genes,
                     no random subset / no null-correlation step -- avoids the
                     null-correlation crash seen with --shat_mode raw)
  lfsr_single     : mash_1by1 LFSR<0.05 in >=1 condition; single-stage fit
                     (data-driven covariances only, fit directly on all genes,
                     no random subset / no null-correlation step)')
)
opt = parse_args(OptionParser(option_list=option_list))
if (opt$mode != 'combined') stop('This script is configured for --mode combined only.')
if (!opt$shat_mode %in% c('sqrt', 'raw')) stop('--shat_mode must be sqrt or raw')
if (!opt$approach %in% c('qval_third', 'qval_recurrence', 'qval_tenth', 'qval_single', 'pval_single', 'lfsr_single'))
    stop('--approach must be qval_third, qval_recurrence, qval_tenth, qval_single, or lfsr_single')
apply_shat_mode = function(x) {
    if (opt$shat_mode == 'sqrt') sqrt(x) else x
}
files = list.files(opt$checkpoints, pattern='_dglm_results\\.rds$', full.names=TRUE)
if (length(files) == 0) stop('No *_dglm_results.rds files found.')
objs = lapply(files, readRDS)
names(objs) = sub('_dglm_results\\.rds$', '', basename(files))
all.genes = Reduce(union, lapply(objs, function(x) dimnames(x$array)[[1]]))
conditions = character(0)
for (nm in names(objs)) {
    a = objs[[nm]]$array
    stat_names = dimnames(a)[[2]]
    beta_col = if ('beta' %in% stat_names) 'beta' else 'beta_age'
    for (r in dimnames(a)[[3]]) {
        b = a[, beta_col, r]
        if (sum(!is.na(b) & b != 0) > 0) {
            conditions = c(conditions, paste(nm, r, sep='|'))
        }
    }
}
message('Units: ', length(objs))
message('Conditions: ', length(conditions))
message('Genes: ', length(all.genes))
message('NA Shat fill value: ', opt$na_shat_fill)
Bhat = matrix(0, nrow=length(all.genes), ncol=length(conditions),
              dimnames=list(all.genes, conditions))
Shat = matrix(opt$na_shat_fill, nrow=length(all.genes), ncol=length(conditions),
              dimnames=list(all.genes, conditions))
pval = matrix(1, nrow=length(all.genes), ncol=length(conditions),
    dimnames=list(all.genes, conditions))
qval = matrix(1, nrow=length(all.genes), ncol=length(conditions),
              dimnames=list(all.genes, conditions))

# kept alongside Bhat/Shat/qval so downstream (dglm_fdr_combined.R) can
# rebuild per-condition raw beta/bvar/pval and gene symbols -- these are
# NOT derived from Bhat/Shat (which are autosome+X filtered in some modes;
# these are the full, unfiltered raw arrays as loaded)
dglm.results.list  = list()
human.symbols.list = list()

for (nm in names(objs)) {
    a = objs[[nm]]$array
    genes = dimnames(a)[[1]]
    stat_names = dimnames(a)[[2]]
    beta_col = if ('beta' %in% stat_names) 'beta' else 'beta_age'
    bvar_col = if ('bvar' %in% stat_names) 'bvar' else 'bvar_age'
    qval_col = if ('qval' %in% stat_names) 'qval' else 'qval_age'
    dglm.results.list[[nm]]  = a
    human.symbols.list[[nm]] = objs[[nm]]$human_symbols
    for (r in dimnames(a)[[3]]) {
        cond = paste(nm, r, sep='|')
        if (!cond %in% conditions) next
        beta = a[, beta_col, r]
        bvar = a[, bvar_col, r]
        q = a[, qval_col, r]
        idx = genes[genes %in% all.genes]
        no.info = is.na(bvar[idx]) | bvar[idx] <= 0
        Bhat[idx, cond] = ifelse(is.na(beta[idx]) | no.info, 0, beta[idx])
        Shat[idx, cond] = ifelse(no.info, opt$na_shat_fill, apply_shat_mode(bvar[idx]))
        pval_col = if ('pval' %in% stat_names) 'pval' else 'pval_age'
        pv = a[, pval_col, r]
        pval[idx, cond] = ifelse(is.na(pv[idx]) | no.info, 1, pv[idx])
        qval[idx, cond] = ifelse(is.na(q[idx]) | no.info, 1, q[idx])
    }
}

human.symbols = setNames(all.genes, all.genes)
for (nm in names(human.symbols.list)) {
    sym = human.symbols.list[[nm]]
    for (g in names(sym)) {
        if (g %in% names(human.symbols) && human.symbols[g] == g &&
            !is.na(sym[g]) && sym[g] != g && sym[g] != '') {
            human.symbols[g] = sym[g]
        }
    }
}

mash.data = mash_set_data(Bhat, Shat)

message('Approach: ', opt$approach)
if (opt$approach == 'qval_third') {
    strong.subset = which(apply(qval, 1, function(x) sum(x < 0.05, na.rm=TRUE)) >= (length(conditions)/3))
} else if (opt$approach == 'qval_recurrence') {
    strong.subset = which(apply(qval, 1, function(x) sum(x < 0.01, na.rm=TRUE)) >= 10)
} else if (opt$approach == 'qval_tenth') {
    strong.subset = which(apply(qval, 1, function(x) sum(x < 0.05, na.rm=TRUE)) >= ceiling(length(conditions)/10))
} else if (opt$approach == 'qval_single') {
    strong.subset = which(apply(qval, 1, function(x) sum(x < opt$strong_q, na.rm=TRUE)) >= 1)
} else if (opt$approach == 'pval_single') {
    strong.subset = which(apply(pval, 1, function(x) sum(x < 0.01, na.rm=TRUE)) >= 1)
} else {
    m.1by1 = mash_1by1(mash.data)
    strong.subset = get_significant_results(m.1by1, thresh=0.05)
}
message('Strong subset: ', length(strong.subset))

if (opt$approach %in% c('lfsr_single', 'qval_tenth', 'qval_single')) {
    n_pca = min(5, ncol(mash.data$Bhat) - 1)
    U.pca = cov_pca(mash.data, n_pca, subset=strong.subset)
    U.ed  = cov_ed(mash.data, U.pca, subset=strong.subset)
    m = mash(mash.data, Ulist=U.ed)
} else {
    set.seed(seed)
    random.subset = sample(1:nrow(Bhat), ceiling(nrow(Bhat) / 2))
    message('Random subset: ', length(random.subset), ' genes (half of all genes)')
    temp = mash_set_data(Bhat[random.subset, ], Shat[random.subset, ])
    Vhat = estimate_null_correlation_simple(temp)
    rm(temp)
    mash.random = mash_set_data(
        Bhat[random.subset, ],
        Shat[random.subset, ],
        V=Vhat
    )
    mash.strong = mash_set_data(
        Bhat[strong.subset, ],
        Shat[strong.subset, ],
        V=Vhat
    )
    n_pca2 = min(5, ncol(mash.strong$Bhat) - 1)
    U.pca = cov_pca(mash.strong, n_pca2)
    U.ed = cov_ed(mash.strong, U.pca)
    U.c = cov_canonical(mash.random)
    m.r = mash(
        mash.random,
        Ulist=c(U.ed, U.c),
        outputlevel=1
    )
    m = mash(
        mash.data,
        g=get_fitted_g(m.r),
        fixg=TRUE
    )
}
out.file = file.path(
    opt$checkpoints,
    paste0('combined_dglm_mashr_results_', opt$approach, '_', opt$shat_mode, '_fill', opt$na_shat_fill, '.rds')
)
saveRDS(
    list(
        mash=m,
        Bhat=Bhat,
        Shat=Shat,
        dglm_qval=qval,
        human_symbols=human.symbols,
        dglm_results=dglm.results.list,
        conditions=conditions
    ),
    out.file
)
message('Saved: ', out.file)
