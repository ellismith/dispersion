#!/usr/bin/env Rscript
#
# sweep_recurrence_thresholds.R
#
# One mash_1by1 fit, then sweeps combinations of (per-condition lfsr
# threshold) x (minimum number of conditions a gene must clear it in) to
# find a strong-subset rule that actually captures cross-condition
# recurring signal, not just single-condition idiosyncratic hits.

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')
suppressMessages({library(optparse); library(mashr)})

option_list = list(
    make_option('--checkpoints', type='character', required=TRUE),
    make_option('--shat_mode', type='character', default='sqrt'),
    make_option('--gene_keep_list', type='character',
                default='/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv')
)
opt = parse_args(OptionParser(option_list=option_list))

autosome.x.genes = read.csv(opt$gene_keep_list, stringsAsFactors=FALSE)$ensembl_gene_id
files = list.files(opt$checkpoints, pattern='_dglm_results\\.rds$', full.names=TRUE)
objs = lapply(files, readRDS)
names(objs) = sub('_dglm_results\\.rds$', '', basename(files))
all.genes = intersect(Reduce(union, lapply(objs, function(x) dimnames(x$array)[[1]])), autosome.x.genes)

conditions = character(0)
for (nm in names(objs)) {
    a = objs[[nm]]$array
    for (r in dimnames(a)[[3]]) {
        if (sum(!is.na(a[, 'beta', r]) & a[, 'beta', r] != 0) > 0) conditions = c(conditions, paste(nm, r, sep='|'))
    }
}

Bhat = matrix(0, nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))
Shat = matrix(1000, nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))
for (nm in names(objs)) {
    dglm.results = objs[[nm]]$array
    beta.all = dglm.results[, 'beta', , drop=FALSE]
    extreme = apply(beta.all, 1, function(x) any(abs(x) > 100, na.rm=TRUE))
    if (sum(extreme) > 0) dglm.results = dglm.results[!extreme, , , drop=FALSE]
    genes.ct = dimnames(dglm.results)[[1]]
    for (r in dimnames(dglm.results)[[3]]) {
        cond = paste(nm, r, sep='|')
        if (!cond %in% conditions) next
        beta = dglm.results[, 'beta', r]; bvar = dglm.results[, 'bvar', r]
        idx = genes.ct[genes.ct %in% all.genes]
        no.info = is.na(bvar[idx]) | bvar[idx] <= 0
        Bhat[idx, cond] = ifelse(is.na(beta[idx]) | no.info, 0, beta[idx])
        Shat[idx, cond] = ifelse(no.info, 1000, if (opt$shat_mode=='sqrt') sqrt(bvar[idx]) else bvar[idx])
    }
}

message('Genes: ', nrow(Bhat), ' | Conditions: ', ncol(Bhat))
data = mash_set_data(Bhat, Shat)
message('Fitting mash_1by1 (one-time cost)...')
start = Sys.time()
m.1by1 = mash_1by1(data)
message('  elapsed: ', format(Sys.time() - start))
lfsr = get_lfsr(m.1by1)

lfsr.thresholds = c(0.05, 0.1, 0.2, 0.3)
min.conditions = c(1, 2, 3, 5, 10)

message('')
message('=== genes clearing (lfsr < T) in >= N conditions ===')
message(sprintf('%-12s %8s %8s %8s %8s %8s', 'lfsr <', '>=1', '>=2', '>=3', '>=5', '>=10'))
for (t in lfsr.thresholds) {
    n.sig.per.gene = apply(lfsr, 1, function(x) sum(x < t, na.rm=TRUE))
    counts = sapply(min.conditions, function(m) sum(n.sig.per.gene >= m))
    message(sprintf('%-12s %8d %8d %8d %8d %8d', t, counts[1], counts[2], counts[3], counts[4], counts[5]))
}
message('done.')
