#!/usr/bin/env Rscript
# dglm_model.R
# Fits DGLM per gene per region for a given cell type.
# --covariates age_sex : mean model = age + sex
# --covariates full    : mean model = age + sex + mean_n_umi + n_cells

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(parallel)
library(doParallel)
library(dglm)
library(reshape2)
library(abind)

option_list = list(
    make_option('--cell_type',   type='character', help='cell type to analyze'),
    make_option('--outdir',      type='character', help='output directory'),
    make_option('--covariates',  type='character', default='full',
                help='age_sex or full (age+sex+mean_n_umi+n_cells)')
)
opt = parse_args(OptionParser(option_list=option_list))

cell.type = opt$cell_type
outdir    = opt$outdir
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# set covariates based on flag
if (opt$covariates == 'age_sex') {
    covs = c('age', 'sex')
} else {
    covs = c('age', 'sex', 'mean_n_umi', 'n_cells')
}
message('Cell type: ', cell.type)
message('Covariates: ', paste(covs, collapse=' + '))
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
meta.all     = read.csv(file.path(outdir, paste0(cell.type, '_metadata.csv')),
                         stringsAsFactors=FALSE)
meta.all$sex = as.factor(meta.all$sex)
meta.all$age = as.numeric(meta.all$age)
if ('mean_n_umi' %in% colnames(meta.all)) meta.all$mean_n_umi = as.numeric(meta.all$mean_n_umi)
if ('n_cells'    %in% colnames(meta.all)) meta.all$n_cells    = as.numeric(meta.all$n_cells)

# verify required covariates present
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

    c.this = covs[sapply(covs, function(cv) length(unique(m[[cv]])) > 1)]
    message('  ', region, ': ', nrow(m), ' animals | covariates: ', paste(c.this, collapse=', '))

    m.this = m[, c.this, drop=FALSE]
    clus   = makeCluster(n.cores)
    registerDoParallel(cores=n.cores)
    clusterExport(clus, varlist=c('e.this','m.this','c.this','predictor'), envir=environment())

    out[[region]] = t(parApply(clus, as.matrix(e.this), 1, function(y) {
        require(dglm)
        d   = m.this; d$e = y
        res = try(dglm(
            as.formula(paste('e ~', paste(c.this, collapse=' + '))),
            as.formula(paste('~', predictor)),
            family=gaussian(), dlink='log', data=d
        ), silent=TRUE)
        if ('try-error' %in% class(res)) c(NA,NA,NA) else
            coef(summary(res)$dispersion.summary)[predictor, c(1,2,4)]
    }))
    stopCluster(clus)
    colnames(out[[region]]) = c('beta','bvar','pval')
    message('  ', region, ': done')
}

# ── assemble array ────────────────────────────────────────────────────────
genes = Reduce(union, lapply(out[!sapply(out, is.null)], rownames))
arr   = array(NA, dim=c(length(genes), 4, length(region.levels)),
              dimnames=list(genes, c('beta','bvar','pval','qval'), region.levels))

for (r in region.levels) {
    if (is.null(out[[r]])) next
    g = rownames(out[[r]])
    arr[g,'beta',r] = out[[r]][,'beta']
    arr[g,'bvar',r] = out[[r]][,'bvar']
    arr[g,'pval',r] = out[[r]][,'pval']
    arr[g,'qval',r] = p.adjust(out[[r]][,'pval'], method='fdr')
}

human.symbols      = sapply(genes, get_symbol)
names(human.symbols) = genes

saveRDS(list(
    array         = arr,
    human_symbols = human.symbols,
    run_info      = list(
        covariates = covs,
        date       = Sys.time(),
        cell_type  = cell.type,
        outdir     = outdir
    )
), file=file.path(outdir, paste0(cell.type, '_dglm_results.rds')))

message('Saved to: ', outdir)
message('Covariates used: ', paste(covs, collapse=' + '))
message('done.')
