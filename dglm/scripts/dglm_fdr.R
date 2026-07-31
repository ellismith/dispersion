#!/usr/bin/env Rscript
# dglm_fdr.R
# Pools all DGLM p-values across all genes x regions x cell types,
# applies a single global BH-FDR correction, saves master table.
# Reads from _dglm_mashr_results.rds if available, else _dglm_results.rds.

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)

option_list = list(
    make_option('--checkpoints', type='character', default='/scratch/easmit31/dispersion/dglm/checkpoints'),
    make_option('--outdir',      type='character', default='/scratch/easmit31/dispersion/dglm/checkpoints')
)
opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

all.results = list()

for (ct in cell.type.levels) {
    # prefer mashr RDS (has dglm_results), fall back to raw DGLM RDS (has array)
    mashr.file = file.path(opt$checkpoints, paste0(ct, '_dglm_mashr_results.rds'))
    raw.file   = file.path(opt$checkpoints, paste0(ct, '_dglm_results.rds'))

    if (file.exists(mashr.file)) {
        obj          = readRDS(mashr.file)
        dglm.results = obj$dglm_results
        human.sym    = obj$human_symbols
        regions      = obj$regions
    } else if (file.exists(raw.file)) {
        obj          = readRDS(raw.file)
        dglm.results = obj$array
        human.sym    = obj$human_symbols
        regions      = region.levels[region.levels %in% dimnames(dglm.results)[[3]]]
        regions      = regions[apply(dglm.results[,'beta', regions, drop=FALSE], 2,
                                     function(x) sum(!is.na(x) & x != 0) > 0)]
    } else {
        message('Skipping ', ct, ' — no RDS file found')
        next
    }
    message('Loading ', ct, ' from ', ifelse(file.exists(mashr.file), 'mashr RDS', 'raw DGLM RDS'))

    # filter extreme betas
    beta.mat = dglm.results[,'beta', regions, drop=FALSE][,1,, drop=FALSE]
    extreme  = apply(beta.mat, 1, function(x) any(abs(x) > 100, na.rm=TRUE))
    if (sum(extreme) > 0) {
        dglm.results = dglm.results[!extreme,,, drop=FALSE]
        human.sym    = human.sym[!extreme]
    }

    gene.names.file = file.path(opt$checkpoints, paste0(ct, '_gene_names.csv'))
    gene.names.df   = if (file.exists(gene.names.file)) {
        read.csv(gene.names.file, row.names=1, stringsAsFactors=FALSE)
    } else NULL

    for (r in regions) {
        ensembl.ids = dimnames(dglm.results)[[1]]
        pval        = dglm.results[,'pval', r]
        beta        = dglm.results[,'beta', r]
        bvar        = dglm.results[,'bvar', r]

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

master        = do.call(rbind, all.results)
master        = master[!is.na(master$pvalue),]
master$qvalue = p.adjust(master$pvalue, method='fdr')

message('Total tests: ', nrow(master))
message('Sig at q<0.05: ', sum(master$qvalue < 0.05))
message('Sig at q<0.2:  ', sum(master$qvalue < 0.2))

out.file = file.path(opt$outdir, 'master_dglm_globalfdr.tsv')
write.table(master, out.file, sep='\t', row.names=FALSE, quote=FALSE)
message('Saved: ', out.file)
message('done.')
