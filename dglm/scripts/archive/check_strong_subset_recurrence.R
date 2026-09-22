#!/usr/bin/env Rscript
#
# check_strong_subset_recurrence.R
#
# For the mash_1by1-selected strong subset, checks whether those genes show
# RECURRING significance across multiple conditions (real shared signal --
# what cov_pca/cov_ed actually need) or just idiosyncratic single-condition
# hits (which wouldn't help learn a meaningful covariance/sharing pattern).

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')
suppressMessages({library(optparse); library(mashr)})

option_list = list(
    make_option('--checkpoints', type='character', required=TRUE),
    make_option('--shat_mode', type='character', default='sqrt'),
    make_option('--gene_keep_list', type='character',
                default='/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv'),
    make_option('--strong_thresh', type='double', default=0.1)
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

data = mash_set_data(Bhat, Shat)
message('Fitting mash_1by1...')
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1, thresh=opt$strong_thresh)
message('Strong subset: ', length(strong), ' genes')

lfsr.strong = get_lfsr(m.1by1)[strong, , drop=FALSE]
n.sig.conditions = apply(lfsr.strong, 1, function(x) sum(x < opt$strong_thresh, na.rm=TRUE))

message('')
message('=== how many conditions is each "strong" gene actually significant in? ===')
print(table(n.sig.conditions))
message('')
message('Genes significant in only 1 condition (idiosyncratic): ', sum(n.sig.conditions == 1),
        ' / ', length(strong), ' (', round(100*sum(n.sig.conditions==1)/length(strong), 1), '%)')
message('Genes significant in >=5 conditions (real recurring signal): ', sum(n.sig.conditions >= 5),
        ' / ', length(strong), ' (', round(100*sum(n.sig.conditions>=5)/length(strong), 1), '%)')
message('done.')
