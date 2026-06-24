#!/usr/bin/env Rscript
# dglm_fdr.R
#
# Pools all DGLM p-values across all genes x regions x cell types,
# applies a single global BH-FDR correction, saves master table.
# Analogous to fdr_correct.py for gene_variance pipeline.
#
# Usage:
#   Rscript dglm_fdr.R

source('/scratch/easmit31/variability/dglm/scripts/_include_options.R')

library(optparse)

option_list = list(
    make_option('--checkpoints', type='character', default='/scratch/easmit31/variability/dglm/checkpoints'),
    make_option('--outdir',      type='character', default='/scratch/easmit31/variability/dglm/checkpoints')
)
opt = parse_args(OptionParser(option_list=option_list))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ── pool all results ───────────────────────────────────────────────────────
all.results = list()

for (ct in cell.type.levels) {
    rds.file = file.path(opt$checkpoints, paste0(ct, '_dglm_mashr_results.rds'))
    if (!file.exists(rds.file)) next
    message('Loading ', ct)
    obj          = readRDS(rds.file)
    dglm.results = obj$dglm_results
    human.sym    = obj$human_symbols
    regions      = obj$regions

    # filter extreme betas
    beta.mat = dglm.results[,'beta', regions, drop=FALSE][,1,, drop=FALSE]
    extreme  = apply(beta.mat, 1, function(x) any(abs(x) > 100, na.rm=TRUE))
    if (sum(extreme) > 0) {
        dglm.results = dglm.results[!extreme,,, drop=FALSE]
        human.sym    = human.sym[!extreme]
    }

    # load gene names
    gene.names.file = file.path(opt$checkpoints, paste0(ct, '_gene_names.csv'))
    gene.names.df   = if (file.exists(gene.names.file)) {
        read.csv(gene.names.file, row.names=1, stringsAsFactors=FALSE)
    } else NULL

    for (r in regions) {
        ensembl.ids = dimnames(dglm.results)[[1]]
        pval        = dglm.results[,'pval', r]
        beta        = dglm.results[,'beta', r]
        bvar        = dglm.results[,'bvar', r]

        # get gene symbols
        symbols = sapply(ensembl.ids, function(g) {
            sym = human.sym[g]
            if (!is.na(sym) && sym != g && sym != '') return(sym)
            if (!is.null(gene.names.df) && g %in% rownames(gene.names.df)) {
                nm = gene.names.df[g, 'external_gene_name']
                if (!is.na(nm) && nm != g && nm != '') return(nm)
            }
            return(g)
        })

        all.results[[length(all.results)+1]] = data.frame(
            ensembl_id = ensembl.ids,
            symbol     = symbols,
            cell_type  = ct,
            region     = r,
            beta       = beta,
            bvar       = bvar,
            pvalue     = pval,
            stringsAsFactors = FALSE
        )
    }
}

master = do.call(rbind, all.results)
message('Total tests: ', nrow(master))

# drop NA pvalues
n.before = nrow(master)
master   = master[!is.na(master$pvalue),]
message('Dropped ', n.before - nrow(master), ' NA pvalues')

# single global BH-FDR correction
master$qvalue = p.adjust(master$pvalue, method='fdr')
message('Sig at q<0.05: ', sum(master$qvalue < 0.05))
message('Sig at q<0.2:  ', sum(master$qvalue < 0.2))
message('Sig positive:  ', sum(master$qvalue < 0.05 & master$beta > 0))
message('Sig negative:  ', sum(master$qvalue < 0.05 & master$beta < 0))

# save
out.file = file.path(opt$outdir, 'master_dglm_globalfdr.tsv')
write.table(master, out.file, sep='\t', row.names=FALSE, quote=FALSE)
message('Saved: ', out.file)
message('done.')
