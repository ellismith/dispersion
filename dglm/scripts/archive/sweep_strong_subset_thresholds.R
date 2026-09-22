#!/usr/bin/env Rscript
#
# sweep_strong_subset_thresholds.R
#
# Checks a range of thresholds for two strong-subset definitions, WITHOUT
# running the expensive full data-driven mash fit:
#   1. mash_1by1 lfsr threshold (one mash_1by1 fit, cheap to re-query at
#      different thresholds afterward)
#   2. raw DGLM q-value threshold (free -- just counting the qval matrix
#      we already build, no mashr fit needed at all)
# "Strong" = significant in at least one condition, at each threshold.
#
# Usage:
#   Rscript sweep_strong_subset_thresholds.R --checkpoints <dir>

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

suppressMessages({
  library(optparse)
  library(mashr)
})

option_list = list(
    make_option('--checkpoints', type='character', required=TRUE),
    make_option('--shat_mode', type='character', default='sqrt'),
    make_option('--gene_keep_list', type='character',
                default='/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv')
)
opt = parse_args(OptionParser(option_list=option_list))

autosome.x.genes = NULL
if (nchar(opt$gene_keep_list) > 0) {
    autosome.x.genes = read.csv(opt$gene_keep_list, stringsAsFactors=FALSE)$ensembl_gene_id
}

files = list.files(opt$checkpoints, pattern='_dglm_results\\.rds$', full.names=TRUE)
objs = lapply(files, readRDS)
names(objs) = sub('_dglm_results\\.rds$', '', basename(files))

all.genes = Reduce(union, lapply(objs, function(x) dimnames(x$array)[[1]]))
if (!is.null(autosome.x.genes)) all.genes = intersect(all.genes, autosome.x.genes)

conditions = character(0)
for (nm in names(objs)) {
    a = objs[[nm]]$array
    for (r in dimnames(a)[[3]]) {
        b = a[, 'beta', r]
        if (sum(!is.na(b) & b != 0) > 0) conditions = c(conditions, paste(nm, r, sep='|'))
    }
}

Bhat = matrix(0, nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))
Shat = matrix(1000, nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))
dglm.qval = matrix(1, nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))

for (nm in names(objs)) {
    dglm.results = objs[[nm]]$array
    beta.all = dglm.results[, 'beta', , drop=FALSE]
    extreme = apply(beta.all, 1, function(x) any(abs(x) > 100, na.rm=TRUE))
    if (sum(extreme) > 0) dglm.results = dglm.results[!extreme, , , drop=FALSE]
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

message('Genes: ', nrow(Bhat), ' | Conditions: ', ncol(Bhat))

# ── option 2: raw DGLM q-value sweep (free, no mashr fit needed) ──────────
message('')
message('=== Raw DGLM q-value threshold sweep (significant in >=1 condition) ===')
q.thresholds = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
for (t in q.thresholds) {
    n = sum(apply(dglm.qval, 1, function(x) any(x < t, na.rm=TRUE)))
    message('  q < ', t, ': ', n, ' / ', nrow(Bhat), ' genes (', round(100*n/nrow(Bhat), 3), '%)')
}

# ── option 1: mash_1by1 lfsr sweep (one fit, cheap to re-query) ───────────
message('')
message('Fitting mash_1by1 (one-time cost for this sweep)...')
data = mash_set_data(Bhat, Shat)
start = Sys.time()
m.1by1 = mash_1by1(data)
message('  elapsed: ', format(Sys.time() - start))

message('')
message('=== mash_1by1 LFSR threshold sweep (significant in >=1 condition) ===')
lfsr.thresholds = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
for (t in lfsr.thresholds) {
    strong = get_significant_results(m.1by1, thresh=t)
    message('  lfsr < ', t, ': ', length(strong), ' / ', nrow(Bhat), ' genes (', round(100*length(strong)/nrow(Bhat), 3), '%)')
}

message('done.')
