#!/usr/bin/env Rscript
# dglm_fdr_age_sex_disp.R
# Global BH-FDR for age_sex dispersion model.
# Saves two master TSVs: one for age effect, one for sex effect.

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')
library(optparse)

option_list = list(
    make_option('--checkpoints', type='character',
                default='/scratch/easmit31/dispersion/dglm/age_sex_dispersion'),
    make_option('--outdir', type='character',
                default='/scratch/easmit31/dispersion/dglm/age_sex_dispersion')
)
opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

age.results = list()
sex.results = list()

for (ct in cell.type.levels) {
    rds.file = file.path(opt$checkpoints, paste0(ct, '_dglm_results.rds'))
    if (!file.exists(rds.file)) { message('Skipping ', ct); next }
    message('Loading ', ct)
    obj          = readRDS(rds.file)
    arr          = obj$array
    human.sym    = obj$human_symbols
    regions      = region.levels[apply(arr[,'beta_age',, drop=FALSE], 3,
                                       function(x) sum(!is.na(x)) > 0)]

    gene.names.file = file.path(opt$checkpoints, paste0(ct, '_gene_names.csv'))
    gene.names.df   = if (file.exists(gene.names.file))
        read.csv(gene.names.file, row.names=1, stringsAsFactors=FALSE) else NULL

    for (r in regions) {
        ensembl.ids = dimnames(arr)[[1]]
        symbols = sapply(ensembl.ids, function(g) {
            sym = human.sym[g]
            if (!is.na(sym) && sym != g && sym != '') return(sym)
            if (!is.null(gene.names.df) && g %in% rownames(gene.names.df)) {
                nm = gene.names.df[g, 'external_gene_name']
                if (!is.na(nm) && nm != g && nm != '') return(nm)
            }
            return(g)
        })

        # filter extreme betas
        keep = abs(arr[,'beta_age',r]) <= 100 & abs(arr[,'beta_sex',r]) <= 100
        keep[is.na(keep)] = FALSE

        age.results[[length(age.results)+1]] = data.frame(
            ensembl_id = ensembl.ids[keep],
            symbol     = symbols[keep],
            cell_type  = ct, region = r,
            beta       = arr[keep,'beta_age',r],
            bvar       = arr[keep,'bvar_age',r],
            pvalue     = arr[keep,'pval_age',r],
            stringsAsFactors = FALSE)

        sex.results[[length(sex.results)+1]] = data.frame(
            ensembl_id = ensembl.ids[keep],
            symbol     = symbols[keep],
            cell_type  = ct, region = r,
            beta       = arr[keep,'beta_sex',r],
            bvar       = arr[keep,'bvar_sex',r],
            pvalue     = arr[keep,'pval_sex',r],
            stringsAsFactors = FALSE)
    }
}

save_master = function(results, fname) {
    master        = do.call(rbind, results)
    master        = master[!is.na(master$pvalue),]
    master$qvalue = p.adjust(master$pvalue, method='fdr')
    message(fname, ': ', nrow(master), ' tests, sig q<0.05: ', sum(master$qvalue < 0.05))
    write.table(master, file.path(opt$outdir, fname), sep='\t', row.names=FALSE, quote=FALSE)
}

save_master(age.results, 'master_dglm_age_disp.tsv')
save_master(sex.results, 'master_dglm_sex_disp.tsv')
message('done.')
