#!/usr/bin/env Rscript
# dglm_fdr_combined.R
# Extracts mashr results from combined run, applies global BH-FDR, saves master TSV.

source('/scratch/easmit31/variability/dglm/scripts/_include_options.R')
library(optparse)
library(mashr)

option_list = list(
    make_option('--checkpoints', type='character', default='/scratch/easmit31/variability/dglm/checkpoints')
)
opt = parse_args(OptionParser(option_list=option_list))

message('Loading combined mashr results')
obj           = readRDS(file.path(opt$checkpoints, 'combined_dglm_mashr_results.rds'))
m             = obj$mash
Bhat          = obj$Bhat
human.symbols = obj$human_symbols
conditions    = obj$conditions

if (is.null(m)) stop('mashr object is NULL')

mash.beta = get_pm(m)
mash.lfsr = get_lfsr(m)
message('lfsr range: ', paste(round(range(mash.lfsr, na.rm=TRUE), 4), collapse=' '))
message('Sig at lfsr<0.2: ', sum(mash.lfsr < 0.2, na.rm=TRUE))

# also extract raw DGLM pvals for global FDR
# build master from Bhat + per-condition raw pvals from dglm_results list
dglm.list = obj$dglm_results
all.rows   = list()

for (cond in conditions) {
    parts  = strsplit(cond, '\\|')[[1]]
    ct     = parts[1]
    region = parts[2]

    if (!ct %in% names(dglm.list)) next
    dglm.results = dglm.list[[ct]]
    if (!region %in% dimnames(dglm.results)[[3]]) next

    genes   = dimnames(dglm.results)[[1]]
    pval    = dglm.results[,'pval', region]
    beta    = dglm.results[,'beta', region]
    bvar    = dglm.results[,'bvar', region]

    # mashr posteriors for this condition
    genes.in.mashr = rownames(mash.beta)
    genes.common   = intersect(genes, genes.in.mashr)

    all.rows[[length(all.rows)+1]] = data.frame(
        ensembl_id  = genes,
        symbol      = human.symbols[genes],
        cell_type   = ct,
        region      = region,
        beta        = beta,
        bvar        = bvar,
        pvalue      = pval,
        mash_beta   = ifelse(genes %in% genes.in.mashr,
                             mash.beta[genes, cond], NA),
        mash_lfsr   = ifelse(genes %in% genes.in.mashr,
                             mash.lfsr[genes, cond], NA),
        stringsAsFactors = FALSE
    )
}

master = do.call(rbind, all.rows)
master = master[!is.na(master$pvalue),]
master$qvalue = p.adjust(master$pvalue, method='fdr')

message('Total tests: ', nrow(master))
message('Sig DGLM q<0.05: ',  sum(master$qvalue    < 0.05,  na.rm=TRUE))
message('Sig mashr lfsr<0.2: ', sum(master$mash_lfsr < 0.2,   na.rm=TRUE))
message('Sig mashr lfsr<0.05:', sum(master$mash_lfsr < 0.05,  na.rm=TRUE))

out = file.path(opt$checkpoints, 'master_dglm_combined.tsv')
write.table(master, out, sep='\t', row.names=FALSE, quote=FALSE)
message('Saved: ', out)
message('done.')
