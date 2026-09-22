#!/usr/bin/env Rscript
source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(mashr)

option_list = list(
    make_option('--checkpoints', type='character', required=TRUE),
    make_option('--shat_mode', type='character', default='sqrt')
)
opt = parse_args(OptionParser(option_list=option_list))

files = list.files(
    opt$checkpoints,
    pattern='_dglm_results\\.rds$',
    full.names=TRUE
)
if (length(files) == 0) stop('No *_dglm_results.rds files found.')

objs = lapply(files, readRDS)
names(objs) = sub('_dglm_results\\.rds$', '', basename(files))
all.genes = Reduce(union, lapply(objs, function(x) dimnames(x$array)[[1]]))

conditions = character(0)
for (nm in names(objs)) {
    a = objs[[nm]]$array
    for (r in dimnames(a)[[3]]) {
        b = a[, 'beta', r]
        if (sum(!is.na(b) & b != 0) > 0) {
            conditions = c(conditions, paste(nm, r, sep='|'))
        }
    }
}

Bhat = matrix(0, nrow=length(all.genes), ncol=length(conditions),
              dimnames=list(all.genes, conditions))
Shat = matrix(1000, nrow=length(all.genes), ncol=length(conditions),
              dimnames=list(all.genes, conditions))
qval = matrix(1, nrow=length(all.genes), ncol=length(conditions),
              dimnames=list(all.genes, conditions))

for (nm in names(objs)) {
    a = objs[[nm]]$array
    genes = dimnames(a)[[1]]

    for (r in dimnames(a)[[3]]) {
        cond = paste(nm, r, sep='|')
        if (!cond %in% conditions) next

        beta = a[, 'beta', r]
        bvar = a[, 'bvar', r]
        q = a[, 'qval', r]
        no.info = is.na(bvar) | bvar <= 0

        Bhat[genes, cond] = ifelse(is.na(beta) | no.info, 0, beta)
        Shat[genes, cond] = ifelse(
            no.info,
            1000,
            if (opt$shat_mode == 'sqrt') sqrt(bvar) else bvar
        )
        qval[genes, cond] = ifelse(is.na(q) | no.info, 1, q)
    }
}

message('Genes: ', nrow(Bhat))
message('Conditions: ', ncol(Bhat))

current.strong = which(
    apply(qval, 1, function(x) sum(x < 0.05, na.rm=TRUE)) >=
    (ncol(Bhat) / 3)
)

message(
    'Current q-value rule: q < 0.05 in >= ',
    ceiling(ncol(Bhat) / 3),
    ' conditions'
)
message(
    'Current strong subset: ', length(current.strong),
    ' / ', nrow(Bhat),
    ' genes (', round(100 * length(current.strong) / nrow(Bhat), 3), '%)'
)

data = mash_set_data(Bhat, Shat)

message('Running mash_1by1 for strong-set definition 2...')
m.1by1 = mash_1by1(data)
onebyone.strong = get_significant_results(m.1by1, thresh=0.05)

message(
    'One-by-one strong subset (lfsr < 0.05 in >=1 condition): ',
    length(onebyone.strong),
    ' / ', nrow(Bhat),
    ' genes (', round(100 * length(onebyone.strong) / nrow(Bhat), 3), '%)'
)

set.seed(seed)
random.subset = sample(seq_len(nrow(Bhat)), ceiling(nrow(Bhat) / 2))

report_overlap = function(label, strong, random) {
    n.overlap = length(intersect(strong, random))
    message(label)
    message('  Strong genes: ', length(strong))
    message('  Random genes: ', length(random))
    message(
        '  Overlap: ', n.overlap,
        ' genes; ', round(100 * n.overlap / max(1, length(strong)), 2),
        '% of strong subset; ',
        round(100 * n.overlap / length(random), 2),
        '% of random subset'
    )
}

report_overlap('Overlap: q-value rule vs random half', current.strong, random.subset)
report_overlap('Overlap: one-by-one rule vs random half', onebyone.strong, random.subset)

message('Running canonical-only joint mashr for strong-set definition 3...')
U.c = cov_canonical(data)

start = Sys.time()
m.c = mash(data, Ulist=U.c, outputlevel=2)
message('Canonical mash elapsed: ', format(Sys.time() - start))

# save the fit immediately, before attempting extraction -- so a failure in
# get_significant_results below doesn't cost another full mash() run to recover
m.c.file = file.path(opt$checkpoints, 'diag_canonical_joint_mash_fit.rds')
saveRDS(m.c, file=m.c.file)
message('Saved canonical mash fit (checkpoint, in case extraction below fails): ', m.c.file)

canonical.strong = tryCatch({
    get_significant_results(m.c, thresh=0.05, sig_fn=ashr::get_lfdr)
}, error = function(e) {
    message('get_significant_results with sig_fn=get_lfdr failed: ', conditionMessage(e))
    message('Falling back to default sig_fn (lfsr-based) instead')
    get_significant_results(m.c, thresh=0.05)
})

message(
    'Canonical joint-mash strong subset (lfdr < 0.05 in >=1 condition): ',
    length(canonical.strong),
    ' / ', nrow(Bhat),
    ' genes (', round(100 * length(canonical.strong) / nrow(Bhat), 3), '%)'
)

report_overlap('Overlap: canonical joint-mash rule vs random half',
               canonical.strong, random.subset)

message(
    'Overlap: one-by-one vs canonical joint-mash strong subsets: ',
    length(intersect(onebyone.strong, canonical.strong)),
    ' genes'
)

write.csv(
    data.frame(
        gene=rownames(Bhat),
        current_qvalue_rule=seq_len(nrow(Bhat)) %in% current.strong,
        onebyone_lfsr_rule=seq_len(nrow(Bhat)) %in% onebyone.strong,
        canonical_joint_lfdr_rule=seq_len(nrow(Bhat)) %in% canonical.strong,
        in_random_half=seq_len(nrow(Bhat)) %in% random.subset
    ),
    file=file.path(opt$checkpoints, 'mashr_strong_subset_diagnostic.csv'),
    row.names=FALSE
)

message('Saved: mashr_strong_subset_diagnostic.csv')
