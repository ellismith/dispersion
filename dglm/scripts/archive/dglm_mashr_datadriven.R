#!/usr/bin/env Rscript
#
# dglm_mashr_datadriven.R
#
# Full, proper data-driven mashr fit (cov_pca + cov_ed), following the
# colleague-provided script's structure -- no --fast/canonical-only
# shortcut, no two-stage random-subset speed trick. Strong subset for
# learning covariance patterns is defined via mash_1by1 + get_significant_
# results (the "Option 2" method already validated on this data: 47/15568
# genes, vs. the original q-value-recurrence rule's 7/15568).
#
# Reads per-subcluster *_dglm_results.rds from --checkpoints (same as
# dglm_mashr.R --mode combined), applies the same autosome+X and extreme-
# beta filters as the rest of the pipeline, then:
#   1. mash_1by1(data) -> strong subset
#   2. cov_pca(data, 5, subset=strong); cov_ed(data, U.pca, subset=strong)
#   3. mash(data, U.ed) -- run directly on ALL genes, once
#
# Saves output in the same structure dglm_mashr.R's combined mode uses, so
# dglm_fdr_combined.R can read it unmodified.
#
# Usage:
#   Rscript dglm_mashr_datadriven.R --checkpoints <filtered per-subcluster dir> --out_checkpoints <new dir>

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

suppressMessages({
  library(optparse)
  library(mashr)
})

option_list = list(
    make_option('--checkpoints', type='character', required=TRUE,
                help='directory of per-subcluster *_dglm_results.rds (e.g. the 300-cell filtered dir)'),
    make_option('--out_checkpoints', type='character', required=TRUE,
                help='where to save the result (kept separate from any existing mashr output)'),
    make_option('--shat_mode', type='character', default='sqrt'),
    make_option('--gene_keep_list', type='character',
                default='/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv'),
    make_option('--strong_thresh', type='double', default=0.05,
                help='mash_1by1 lfsr threshold for the strong subset (default 0.05; use 0.1 for a larger, less thin subset)')
)
opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$out_checkpoints, showWarnings=FALSE, recursive=TRUE)

autosome.x.genes = NULL
if (nchar(opt$gene_keep_list) > 0) {
    if (!file.exists(opt$gene_keep_list)) stop('--gene_keep_list not found: ', opt$gene_keep_list)
    autosome.x.genes = read.csv(opt$gene_keep_list, stringsAsFactors=FALSE)$ensembl_gene_id
    message('Loaded gene keep-list: ', length(autosome.x.genes), ' autosome+X genes')
}

files = list.files(opt$checkpoints, pattern='_dglm_results\\.rds$', full.names=TRUE)
if (length(files) == 0) stop('No *_dglm_results.rds files found in ', opt$checkpoints)
message('Found ', length(files), ' subcluster RDS file(s)')

objs = lapply(files, readRDS)
names(objs) = sub('_dglm_results\\.rds$', '', basename(files))

all.genes = Reduce(union, lapply(objs, function(x) dimnames(x$array)[[1]]))
message('Total unique genes (pre-filter): ', length(all.genes))
if (!is.null(autosome.x.genes)) {
    n.before = length(all.genes)
    all.genes = intersect(all.genes, autosome.x.genes)
    message('Autosome+X filter: dropping ', n.before - length(all.genes), ' of ', n.before, ' genes')
}

conditions = character(0)
for (nm in names(objs)) {
    a = objs[[nm]]$array
    for (r in dimnames(a)[[3]]) {
        b = a[, 'beta', r]
        if (sum(!is.na(b) & b != 0) > 0) conditions = c(conditions, paste(nm, r, sep='|'))
    }
}
message('Total conditions (unit x region): ', length(conditions))

Bhat = matrix(0, nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))
Shat = matrix(1000, nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))
dglm.qval = matrix(1, nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))
human.symbols.list = list()
dglm.results.list = list()

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
        beta = dglm.results[, 'beta', r]
        bvar = dglm.results[, 'bvar', r]
        qval = dglm.results[, 'qval', r]
        idx = genes.ct[genes.ct %in% all.genes]
        no.info = is.na(bvar[idx]) | bvar[idx] <= 0
        Bhat[idx, cond] = ifelse(is.na(beta[idx]) | no.info, 0, beta[idx])
        Shat[idx, cond] = ifelse(no.info, 1000, if (opt$shat_mode == 'sqrt') sqrt(bvar[idx]) else bvar[idx])
        dglm.qval[idx, cond] = ifelse(is.na(qval[idx]) | no.info, 1, qval[idx])
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

message('Final dims: ', nrow(Bhat), ' genes x ', ncol(Bhat), ' conditions')

data = mash_set_data(Bhat, Shat)

message('Computing strong subset via mash_1by1...')
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1, thresh=opt$strong_thresh)
message('Strong subset: ', length(strong), ' / ', nrow(Bhat), ' genes (',
        round(100 * length(strong) / nrow(Bhat), 3), '%)')

message('cov_pca + cov_ed on strong subset...')
start = Sys.time()
U.pca = cov_pca(data, 5, subset=strong)
U.ed = cov_ed(data, U.pca, subset=strong)
message('  elapsed: ', format(Sys.time() - start))

message('Running full mash() on all genes (data-driven covariances, no --fast, no random-subset shortcut)...')
start = Sys.time()
m = mash(data, Ulist=U.ed)
message('  elapsed: ', format(Sys.time() - start))
message('log likelihood: ', format(get_loglik(m), digits=10))

out.file = file.path(opt$out_checkpoints, 'combined_dglm_mashr_results_strong0.05_lfsr0.2.rds')
saveRDS(list(
    mash          = m,
    Bhat          = Bhat,
    Shat          = Shat,
    dglm_qval     = dglm.qval,
    human_symbols = human.symbols,
    conditions    = conditions,
    dglm_results  = dglm.results.list,
    term          = ''
), file=out.file)
message('Saved: ', out.file)
message('done.')
