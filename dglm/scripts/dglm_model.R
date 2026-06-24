#!/usr/bin/env Rscript
# dglm_model.R
#
# Fits DGLM per gene per region for a given cell type.
# Directly adapted from Chiou et al. 2022 Nat Neurosci dglm_model.R.
#
# SAME as Chiou et al.:
#   - DGLM formulation: mean ~ covariates, dispersion ~ age, gaussian/log
#   - Covariate dropping for uniform variables
#   - Parallel execution with parApply
#   - Output: genes x stats x regions array with per-region FDR
#
# DEVIATIONS from Chiou et al.:
#   - Input: pseudobulk CSV (mean across cells) not bulk RNA-seq library
#   - Covariates: sequencing_run_id/sample_reads/n_umi instead of RIN/batch/reads
#   - Analysis per cell type x region, not per region only
#   - Gene symbol mapping uses external_gene_name + human ortholog (one2one)
#
# Usage:
#   Rscript dglm_model.R --cell_type microglia

source('/scratch/easmit31/variability/dglm/scripts/_include_options.R')

library(optparse)
library(parallel)
library(doParallel)
library(dglm)
library(reshape2)
library(abind)

option_list = list(
    make_option('--cell_type',   type='character', help='cell type to analyze'),
    make_option('--checkpoints', type='character', default='/scratch/easmit31/variability/dglm/checkpoints')
)
opt = parse_args(OptionParser(option_list=option_list))

cell.type = opt$cell_type
message('Cell type: ', cell.type)

dir.create(opt$checkpoints, showWarnings=FALSE, recursive=TRUE)

# ── load ortholog map ──────────────────────────────────────────────────────
# Priority: human one2one ortholog > external_gene_name > ensembl_id
orth = read.csv('/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv')
orth = orth[orth$Human.homology.type == 'ortholog_one2one' & orth$Human.gene.name != '',
            c('Gene.stable.ID','Human.gene.name')]
orth = orth[!duplicated(orth$Gene.stable.ID),]
rownames(orth) = orth$Gene.stable.ID

# load external gene names from var table
gene.names.file = file.path(opt$checkpoints, paste0(cell.type, '_gene_names.csv'))
if (file.exists(gene.names.file)) {
    gene.names.df = read.csv(gene.names.file, row.names=1, stringsAsFactors=FALSE)
} else {
    gene.names.df = NULL
    message('Warning: gene names CSV not found, will use ensembl IDs')
}

get_symbol = function(ensembl.id) {
    if (ensembl.id %in% rownames(orth)) return(orth[ensembl.id, 'Human.gene.name'])
    if (!is.null(gene.names.df) && ensembl.id %in% rownames(gene.names.df)) {
        nm = gene.names.df[ensembl.id, 'external_gene_name']
        if (!is.na(nm) && nm != '' && nm != ensembl.id) return(nm)
    }
    return(ensembl.id)
}

# ── load metadata ──────────────────────────────────────────────────────────
meta.all      = read.csv(file.path(opt$checkpoints, paste0(cell.type, '_metadata.csv')),
                          stringsAsFactors=FALSE)
meta.all$sex  = as.factor(meta.all$sex)
meta.all$age  = as.numeric(meta.all$age)
meta.all$sequencing_run_id = as.factor(meta.all$sequencing_run_id)

# ── initialize output — same structure as Chiou et al. ────────────────────
out        = vector('list', length(region.levels))
names(out) = region.levels

for (i in 1:length(region.levels)) {
    region  = region.levels[i]
    message('Now analyzing region ', region)

    pb.file = file.path(opt$checkpoints, paste0(cell.type, '_', region, '_pseudobulk.csv'))
    if (!file.exists(pb.file)) {
        message('  no pseudobulk file, skipping')
        next
    }

    # load pseudobulk matrix (genes x animals)
    # analogous to Chiou et al. e.keep[keep.genes[[region]], libraries]
    e.this = read.csv(pb.file, row.names=1, check.names=FALSE)
    message('  shape: ', nrow(e.this), ' genes x ', ncol(e.this), ' animals')

    # get metadata for this region — match animal order to columns
    m      = meta.all[meta.all$region == region,]
    m      = m[match(colnames(e.this), m$animal_id),]
    rownames(m) = m$animal_id

    if (nrow(m) < 5) {
        message('  too few animals, skipping')
        next
    }

    # drop covariates uniform across dataset — same as Chiou et al.
    c.this = model.covariates[sapply(model.covariates, function(cv) {
        length(unique(m[[cv]])) > 1
    })]
    message('  covariates: ', paste(c.this, collapse=', '))
    message('  n animals: ', nrow(m))

    m.this = m[, c.this, drop=FALSE]

    # ── parallel DGLM per gene — identical to Chiou et al. ────────────────
    clus = makeCluster(n.cores)
    registerDoParallel(cores=n.cores)
    clusterExport(clus, varlist=c('e.this','m.this','c.this','predictor'),
                  envir=environment())

    out[[region]] = t(parApply(clus, as.matrix(e.this), 1, function(y) {
        require(dglm)
        d   = m.this
        d$e = y
        # identical DGLM call to Chiou et al.
        results = try(dglm(
            as.formula(paste('e ~', paste(c.this, collapse=' + '))),
            as.formula(paste('~', predictor)),
            family  = gaussian(),
            dlink   = 'log',
            data    = d
        ), silent=TRUE)
        if ('try-error' %in% class(results)) {
            c(NA, NA, NA)
        } else {
            coef(summary(results)$dispersion.summary)[predictor, c(1, 2, 4)]
        }
    }))
    stopCluster(clus)

    colnames(out[[region]]) = c('beta', 'bvar', 'pval')
    message('  done: ', nrow(out[[region]]), ' genes')
}

# ── assemble 3D array — identical to Chiou et al. ─────────────────────────
regions.dimnames = list(
    genes   = Reduce(union, lapply(out[!sapply(out, is.null)], rownames)),
    outputs = c('beta', 'bvar', 'pval'),
    regions = region.levels
)
regions.dim     = unlist(lapply(regions.dimnames, length))
regions.numeric = numeric(Reduce(`*`, regions.dim))
regions.numeric[!regions.numeric] = NA
regions.array   = array(unname(regions.numeric),
                         dim=unname(regions.dim),
                         dimnames=unname(regions.dimnames))

for (i in 1:length(out)) {
    if (is.null(out[[i]])) next
    foo      = reshape2::melt(out[[i]])
    foo$Var3 = names(out)[i]
    j        = as.matrix(foo[, paste0('Var', 1:3)])
    regions.array[j] = foo$value
    rm(foo)
}

# ── add per-region FDR — identical to Chiou et al. ────────────────────────
regions.array = abind(
    regions.array,
    array(NA, dim=list(dim(regions.array)[1], 1, length(region.levels)),
          dimnames=append(dimnames(regions.array)[c(1,3)], list('qval'), after=1)),
    along=2
)
for (r in region.levels) {
    regions.array[,'qval',r] = p.adjust(regions.array[,'pval',r], method='fdr')
}

# ── gene symbol mapping ────────────────────────────────────────────────────
gene.ids      = dimnames(regions.array)[[1]]
human.symbols = sapply(gene.ids, get_symbol)
names(human.symbols) = gene.ids

# ── save ──────────────────────────────────────────────────────────────────
out.file = file.path(opt$checkpoints, paste0(cell.type, '_dglm_results.rds'))
saveRDS(list(array=regions.array, human_symbols=human.symbols), file=out.file)
message('Saved: ', out.file)
message('done.')
