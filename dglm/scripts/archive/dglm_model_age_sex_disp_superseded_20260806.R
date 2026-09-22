#!/usr/bin/env Rscript
# dglm_model_age_sex_disp.R
# Same as dglm_model.R but dispersion model is ~ age + sex
# Saves both age and sex dispersion coefficients.
# Mean model: age + sex + mean_n_umi + n_cells
# Dispersion model: ~ age + sex

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(parallel)
library(doParallel)
library(dglm)
library(reshape2)
library(abind)

option_list = list(
    make_option('--cell_type', type='character', help='cell type to analyze'),
    make_option('--outdir',    type='character', help='output directory')
)
opt = parse_args(OptionParser(option_list=option_list))

cell.type = opt$cell_type
outdir    = opt$outdir
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

covs      = c('age', 'sex', 'mean_n_umi', 'n_cells')
disp.covs = c('age', 'sex')
message('Cell type: ', cell.type)
message('Mean covariates: ', paste(covs, collapse=' + '))
message('Dispersion covariates: ', paste(disp.covs, collapse=' + '))
message('Output dir: ', outdir)

# ── load ortholog map ──────────────────────────────────────────────────────
orth = read.csv('/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv')
orth = orth[orth$Human.homology.type == 'ortholog_one2one' & orth$Human.gene.name != '',
            c('Gene.stable.ID','Human.gene.name')]
orth = orth[!duplicated(orth$Gene.stable.ID),]
rownames(orth) = orth$Gene.stable.ID

gene.names.file = file.path(outdir, paste0(cell.type, '_gene_names.csv'))
gene.names.df   = if (file.exists(gene.names.file)) read.csv(gene.names.file, row.names=1) else NULL

get_symbol = function(g) {
    if (g %in% rownames(orth)) return(orth[g,'Human.gene.name'])
    if (!is.null(gene.names.df) && g %in% rownames(gene.names.df)) {
        nm = gene.names.df[g,'external_gene_name']
        if (!is.na(nm) && nm != g && nm != '') return(nm)
    }
    return(g)
}

# ── load metadata ──────────────────────────────────────────────────────────
meta.all           = read.csv(file.path(outdir, paste0(cell.type, '_metadata.csv')),
                               stringsAsFactors=FALSE)
meta.all$sex       = as.factor(meta.all$sex)
meta.all$age       = as.numeric(meta.all$age)
meta.all$mean_n_umi = as.numeric(meta.all$mean_n_umi)
meta.all$n_cells    = as.numeric(meta.all$n_cells)

missing = setdiff(covs, colnames(meta.all))
if (length(missing) > 0) stop('Missing columns in metadata: ', paste(missing, collapse=', '))

# ── run DGLM per region ───────────────────────────────────────────────────
out        = vector('list', length(region.levels))
names(out) = region.levels

for (region in region.levels) {
    pb.file = file.path(outdir, paste0(cell.type, '_', region, '_pseudobulk.csv'))
    if (!file.exists(pb.file)) next

    e.this = read.csv(pb.file, row.names=1, check.names=FALSE)
    m      = meta.all[meta.all$region == region,]
    m      = m[match(colnames(e.this), m$animal_id),]
    rownames(m) = m$animal_id
    if (nrow(m) < 5) next

    c.this      = covs[sapply(covs, function(cv) length(unique(m[[cv]])) > 1)]
    c.disp.this = disp.covs[sapply(disp.covs, function(cv) length(unique(m[[cv]])) > 1)]
    message('  ', region, ': ', nrow(m), ' animals | mean: ', paste(c.this, collapse='+'),
            ' | disp: ', paste(c.disp.this, collapse='+'))

    m.this = m[, unique(c(c.this, c.disp.this)), drop=FALSE]
    clus   = makeCluster(n.cores)
    registerDoParallel(cores=n.cores)
    clusterExport(clus, varlist=c('e.this','m.this','c.this','c.disp.this'), envir=environment())

    out[[region]] = t(parApply(clus, as.matrix(e.this), 1, function(y) {
        require(dglm)
        d   = m.this; d$e = y
        res = try(dglm(
            as.formula(paste('e ~', paste(c.this, collapse=' + '))),
            as.formula(paste('~', paste(c.disp.this, collapse=' + '))),
            family=gaussian(), dlink='log', data=d
        ), silent=TRUE)
        if ('try-error' %in% class(res)) {
            c(NA,NA,NA, NA,NA,NA)
        } else {
            disp.sum = coef(summary(res)$dispersion.summary)
            age.row  = if ('age' %in% rownames(disp.sum)) disp.sum['age', c(1,2,4)] else c(NA,NA,NA)
            sex.row  = if (any(grepl('sex', rownames(disp.sum)))) {
                sex.nm = rownames(disp.sum)[grep('sex', rownames(disp.sum))[1]]
                disp.sum[sex.nm, c(1,2,4)]
            } else c(NA,NA,NA)
            c(age.row, sex.row)
        }
    }))
    stopCluster(clus)
    colnames(out[[region]]) = c('beta_age','bvar_age','pval_age',
                                 'beta_sex','bvar_sex','pval_sex')
    message('  ', region, ': done')
}

# ── assemble array ────────────────────────────────────────────────────────
genes  = Reduce(union, lapply(out[!sapply(out, is.null)], rownames))
stats  = c('beta_age','bvar_age','pval_age','qval_age',
           'beta_sex','bvar_sex','pval_sex','qval_sex')
arr    = array(NA, dim=c(length(genes), length(stats), length(region.levels)),
               dimnames=list(genes, stats, region.levels))

for (r in region.levels) {
    if (is.null(out[[r]])) next
    g = rownames(out[[r]])
    for (s in c('beta_age','bvar_age','pval_age','beta_sex','bvar_sex','pval_sex')) {
        arr[g,s,r] = out[[r]][,s]
    }
    arr[g,'qval_age',r] = p.adjust(out[[r]][,'pval_age'], method='fdr')
    arr[g,'qval_sex',r] = p.adjust(out[[r]][,'pval_sex'], method='fdr')
}

human.symbols        = sapply(genes, get_symbol)
names(human.symbols) = genes

saveRDS(list(
    array         = arr,
    human_symbols = human.symbols,
    run_info      = list(
        covariates      = covs,
        disp_covariates = disp.covs,
        date            = Sys.time(),
        cell_type       = cell.type,
        outdir          = outdir
    )
), file=file.path(outdir, paste0(cell.type, '_dglm_results.rds')))

message('Saved to: ', outdir)
message('done.')
