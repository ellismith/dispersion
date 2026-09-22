#!/usr/bin/env Rscript
#
# dglm_mashr_datadriven_v2.R
#
# Full data-driven mashr fit, testing all three strong-subset definitions,
# always using COMBINED canonical + data-driven covariances (Ulist =
# c(U.c, U.ed)) -- matching every variant in the reference script, unlike
# our earlier data-driven-only (U.ed alone) attempt.
#
# --strong_method:
#   qval_recurrence : DGLM q<0.05 in >= --qval_min_conditions conditions
#                     (Option 1; default recreates the original ~1/3 rule)
#   mash1by1        : mash_1by1 + get_significant_results (Option 2)
#   canonical_lfdr  : fit a preliminary canonical-only joint mash, then
#                     get_significant_results(..., sig_fn=ashr::get_lfdr)
#                     (Option 3 -- the expensive one, ~35+ min just for
#                     this preliminary step)
#
# Usage:
#   Rscript dglm_mashr_datadriven_v2.R --checkpoints <dir> --out_checkpoints <new dir> --strong_method mash1by1

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')
suppressMessages({library(optparse); library(mashr)})

option_list = list(
    make_option('--checkpoints', type='character', required=TRUE),
    make_option('--out_checkpoints', type='character', required=TRUE),
    make_option('--shat_mode', type='character', default='sqrt'),
    make_option('--gene_keep_list', type='character',
                default='/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv'),
    make_option('--strong_method', type='character', default='mash1by1',
                help='qval_recurrence | mash1by1 | canonical_lfdr'),
    make_option('--strong_thresh', type='double', default=0.1,
                help='threshold used by mash1by1 or canonical_lfdr'),
    make_option('--qval_min_conditions', type='integer', default=NULL,
                help='for qval_recurrence: min conditions with q<0.05. Default: ceiling(n_conditions/3), matching the original rule')
)
opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$out_checkpoints, showWarnings=FALSE, recursive=TRUE)

autosome.x.genes = read.csv(opt$gene_keep_list, stringsAsFactors=FALSE)$ensembl_gene_id
files = list.files(opt$checkpoints, pattern='_dglm_results\\.rds$', full.names=TRUE)
message('Found ', length(files), ' subcluster RDS file(s)')
objs = lapply(files, readRDS)
names(objs) = sub('_dglm_results\\.rds$', '', basename(files))

all.genes = Reduce(union, lapply(objs, function(x) dimnames(x$array)[[1]]))
n.before = length(all.genes)
all.genes = intersect(all.genes, autosome.x.genes)
message('Autosome+X filter: ', n.before, ' -> ', length(all.genes), ' genes')

conditions = character(0)
for (nm in names(objs)) {
    a = objs[[nm]]$array
    for (r in dimnames(a)[[3]]) {
        if (sum(!is.na(a[, 'beta', r]) & a[, 'beta', r] != 0) > 0) conditions = c(conditions, paste(nm, r, sep='|'))
    }
}
message('Conditions: ', length(conditions))

Bhat = matrix(0, nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))
Shat = matrix(1000, nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))
dglm.qval = matrix(1, nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))
human.symbols.list = list(); dglm.results.list = list()

for (nm in names(objs)) {
    dglm.results = objs[[nm]]$array
    human.sym = objs[[nm]]$human_symbols
    beta.all = dglm.results[, 'beta', , drop=FALSE]
    extreme = apply(beta.all, 1, function(x) any(abs(x) > 100, na.rm=TRUE))
    if (sum(extreme) > 0) {
        dglm.results = dglm.results[!extreme, , , drop=FALSE]
        human.sym = human.sym[!extreme]
    }
    human.symbols.list[[nm]] = human.sym
    dglm.results.list[[nm]] = dglm.results
    genes.ct = dimnames(dglm.results)[[1]]
    for (r in dimnames(dglm.results)[[3]]) {
        cond = paste(nm, r, sep='|')
        if (!cond %in% conditions) next
        beta = dglm.results[, 'beta', r]; bvar = dglm.results[, 'bvar', r]; qval = dglm.results[, 'qval', r]
        idx = genes.ct[genes.ct %in% all.genes]
        no.info = is.na(bvar[idx]) | bvar[idx] <= 0
        Bhat[idx, cond] = ifelse(is.na(beta[idx]) | no.info, 0, beta[idx])
        Shat[idx, cond] = ifelse(no.info, 1000, if (opt$shat_mode=='sqrt') sqrt(bvar[idx]) else bvar[idx])
        dglm.qval[idx, cond] = ifelse(is.na(qval[idx]) | no.info, 1, qval[idx])
    }
}
human.symbols = setNames(all.genes, all.genes)
for (nm in names(human.symbols.list)) {
    sym = human.symbols.list[[nm]]
    for (g in names(sym)) {
        if (g %in% names(human.symbols) && human.symbols[g] == g && !is.na(sym[g]) && sym[g] != g && sym[g] != '') {
            human.symbols[g] = sym[g]
        }
    }
}
message('Final dims: ', nrow(Bhat), ' genes x ', ncol(Bhat), ' conditions')

data = mash_set_data(Bhat, Shat)

# ── determine strong subset per --strong_method ────────────────────────────
if (opt$strong_method == 'qval_recurrence') {
    min.cond = if (!is.null(opt$qval_min_conditions)) opt$qval_min_conditions else ceiling(ncol(Bhat) / 3)
    message('Strong method: qval_recurrence (q<0.05 in >=', min.cond, ' conditions)')
    strong = which(apply(dglm.qval, 1, function(x) sum(x < 0.05, na.rm=TRUE)) >= min.cond)
} else if (opt$strong_method == 'mash1by1') {
    message('Strong method: mash1by1 (lfsr<', opt$strong_thresh, ' in >=1 condition)')
    m.1by1 = mash_1by1(data)
    strong = get_significant_results(m.1by1, thresh=opt$strong_thresh)
} else if (opt$strong_method == 'canonical_lfdr') {
    message('Strong method: canonical_lfdr -- fitting preliminary canonical-only joint mash (expensive)...')
    U.c.prelim = cov_canonical(data)
    start = Sys.time()
    m.c = mash(data, Ulist=U.c.prelim, outputlevel=2)
    message('  preliminary canonical mash elapsed: ', format(Sys.time() - start))
    saveRDS(m.c, file=file.path(opt$out_checkpoints, 'prelim_canonical_joint_mash.rds'))
    strong = tryCatch({
        get_significant_results(m.c, thresh=opt$strong_thresh, sig_fn=ashr::get_lfdr)
    }, error = function(e) {
        message('lfdr-based extraction failed (', conditionMessage(e), '), falling back to default lfsr-based')
        get_significant_results(m.c, thresh=opt$strong_thresh)
    })
} else {
    stop('Unknown --strong_method: ', opt$strong_method)
}

message('Strong subset: ', length(strong), ' / ', nrow(Bhat), ' genes (',
        round(100*length(strong)/nrow(Bhat), 3), '%)')
if (length(strong) < 2) stop('Strong subset too small (<2 genes) to proceed with cov_pca/cov_ed')

# ── data-driven + canonical covariances, COMBINED ───────────────────────────
message('cov_pca + cov_ed on strong subset...')
start = Sys.time()
U.pca = cov_pca(data, min(5, length(strong)-1), subset=strong)
U.ed = cov_ed(data, U.pca, subset=strong)
message('  elapsed: ', format(Sys.time() - start))

U.c = cov_canonical(data)
Ulist = c(U.ed, U.c)
message('Total mixture components (data-driven + canonical): ', length(Ulist))

message('Running full mash() on all genes...')
start = Sys.time()
m = mash(data, Ulist=Ulist)
message('  elapsed: ', format(Sys.time() - start))
message('log likelihood: ', format(get_loglik(m), digits=10))

out.file = file.path(opt$out_checkpoints, 'combined_dglm_mashr_results_strong0.05_lfsr0.2.rds')
saveRDS(list(
    mash=m, Bhat=Bhat, Shat=Shat, dglm_qval=dglm.qval, human_symbols=human.symbols,
    conditions=conditions, dglm_results=dglm.results.list, term=''
), file=out.file)
message('Saved: ', out.file)
message('done.')
