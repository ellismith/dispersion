#!/usr/bin/env Rscript

# dglm_model_louvain_astrocytes.R
#
# Louvain-level DGLM analysis for Astrocytes only.
#
# This is intentionally separate from dglm_model.R.
# It does not modify or overwrite the original cell-type-level results.
#
# Expected input files in --datadir:
#
#   Astrocytes_metadata.csv
#   Astrocytes_louvain0_<region>_pseudobulk.csv
#   Astrocytes_louvain1_<region>_pseudobulk.csv
#   Astrocytes_louvain2_<region>_pseudobulk.csv
#
# Each pseudobulk file should contain:
#   rows    = genes
#   columns = animal_id values
#
# Metadata should contain:
#   animal_id
#   region
#   sex
#   age
#   louvain
#   mean_n_umi, if using --covariates full
#   n_cells, if using --covariates full
#
# Existing .rds results are skipped and never overwritten.

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(parallel)
library(doParallel)
library(dglm)
library(abind)

option_list = list(
    make_option(
        '--datadir',
        type='character',
        help='Directory containing original metadata and Louvain-level pseudobulk files'
    ),

    make_option(
        '--outdir',
        type='character',
        help='New output directory for Louvain-level results'
    ),

    make_option(
        '--cell_type',
        type='character',
        default='Astrocytes',
        help='Cell type; this script is restricted to Astrocytes'
    ),

    make_option(
        '--louvain_col',
        type='character',
        default='louvain',
        help='Metadata column containing Louvain labels'
    ),

    make_option(
        '--covariates',
        type='character',
        default='full',
        help='Mean-model covariates: age_sex or full'
    ),

    make_option(
        '--dispersion',
        type='character',
        default='age',
        help='Dispersion-model covariates: age or age_sex'
    )
)

opt = parse_args(OptionParser(option_list=option_list))

if (opt$cell_type != 'Astrocytes') {
    stop(
        "This test script is restricted to cell_type='Astrocytes'. ",
        "Use --cell_type Astrocytes."
    )
}

if (!opt$covariates %in% c('age_sex', 'full')) {
    stop("--covariates must be 'age_sex' or 'full'")
}

if (!opt$dispersion %in% c('age', 'age_sex')) {
    stop("--dispersion must be 'age' or 'age_sex'")
}

datadir = normalizePath(opt$datadir, mustWork=TRUE)
outdir  = opt$outdir

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

cell.type   = opt$cell_type
louvain.col = opt$louvain_col

message('Input directory: ', datadir)
message('Output directory: ', outdir)
message('Cell type: ', cell.type)
message('Louvain column: ', louvain.col)

# -------------------------------------------------------------------------
# Mean-model covariates
# -------------------------------------------------------------------------

if (opt$covariates == 'age_sex') {
    covs = c('age', 'sex')
} else {
    covs = c('age', 'sex', 'mean_n_umi', 'n_cells')
}

# Do not use an externally defined "predictor" object here.
# This preserves the intended age-only dispersion model.
disp.covs = if (opt$dispersion == 'age_sex') {
    c('age', 'sex')
} else {
    c('age')
}

message('Mean covariates: ', paste(covs, collapse=' + '))
message('Dispersion covariates: ', paste(disp.covs, collapse=' + '))

# -------------------------------------------------------------------------
# Load ortholog map
# -------------------------------------------------------------------------

orth = read.csv(
    '/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv',
    stringsAsFactors=FALSE
)

orth = orth[
    orth$Human.homology.type == 'ortholog_one2one' &
    orth$Human.gene.name != '',
    c('Gene.stable.ID', 'Human.gene.name')
]

orth = orth[!duplicated(orth$Gene.stable.ID), ]
rownames(orth) = orth$Gene.stable.ID

gene.names.file = file.path(
    datadir,
    paste0(cell.type, '_gene_names.csv')
)

gene.names.df = if (file.exists(gene.names.file)) {
    read.csv(
        gene.names.file,
        row.names=1,
        stringsAsFactors=FALSE
    )
} else {
    NULL
}

get_symbol = function(g) {
    if (g %in% rownames(orth)) {
        return(orth[g, 'Human.gene.name'])
    }

    if (!is.null(gene.names.df) &&
        g %in% rownames(gene.names.df) &&
        'external_gene_name' %in% colnames(gene.names.df)) {

        nm = gene.names.df[g, 'external_gene_name']

        if (!is.na(nm) && nm != g && nm != '') {
            return(nm)
        }
    }

    return(g)
}

# -------------------------------------------------------------------------
# Load metadata
# -------------------------------------------------------------------------

metadata.file = file.path(
    datadir,
    paste0(cell.type, '_metadata.csv')
)

if (!file.exists(metadata.file)) {
    stop('Metadata file not found: ', metadata.file)
}

meta.all = read.csv(
    metadata.file,
    stringsAsFactors=FALSE,
    check.names=FALSE
)

required.base = c(
    'animal_id',
    'region',
    'sex',
    'age',
    louvain.col
)

missing.base = setdiff(required.base, colnames(meta.all))

if (length(missing.base) > 0) {
    stop(
        'Missing required metadata columns: ',
        paste(missing.base, collapse=', ')
    )
}

meta.all$animal_id = as.character(meta.all$animal_id)
meta.all$region    = as.character(meta.all$region)
meta.all$sex       = as.factor(meta.all$sex)
meta.all$age       = as.numeric(meta.all$age)

if ('mean_n_umi' %in% colnames(meta.all)) {
    meta.all$mean_n_umi = as.numeric(meta.all$mean_n_umi)
}

if ('n_cells' %in% colnames(meta.all)) {
    meta.all$n_cells = as.numeric(meta.all$n_cells)
}

meta.all$louvain_value = as.character(meta.all[[louvain.col]])

louvain.levels = sort(
    unique(meta.all$louvain_value[!is.na(meta.all$louvain_value)])
)

if (length(louvain.levels) == 0) {
    stop('No Louvain levels found in metadata column: ', louvain.col)
}

message(
    'Louvain levels found: ',
    paste(louvain.levels, collapse=', ')
)

# -------------------------------------------------------------------------
# Check required covariates before starting
# -------------------------------------------------------------------------

missing.covs = setdiff(
    union(covs, disp.covs),
    colnames(meta.all)
)

if (length(missing.covs) > 0) {
    stop(
        'Missing model covariates in metadata: ',
        paste(missing.covs, collapse=', ')
    )
}

# -------------------------------------------------------------------------
# Analyze one Louvain subtype at a time
# -------------------------------------------------------------------------

for (louvain in louvain.levels) {

    subtype.name = paste0(cell.type, '_louvain', louvain)

    result.file = file.path(
        outdir,
        paste0(subtype.name, '_dglm_results.rds')
    )

    if (file.exists(result.file)) {
        message(
            'Skipping existing result: ',
            result.file
        )
        next
    }

    message('')
    message('==============================================')
    message('Starting: ', subtype.name)
    message('==============================================')

    meta.sub = meta.all[
        meta.all$louvain_value == louvain,
        ,
        drop=FALSE
    ]

    if (nrow(meta.sub) == 0) {
        message('No metadata rows for Louvain level ', louvain)
        next
    }

    out = vector('list', length(region.levels))
    names(out) = region.levels

    for (region in region.levels) {

        pb.file = file.path(
            datadir,
            paste0(
                cell.type,
                '_louvain',
                louvain,
                '_',
                region,
                '_pseudobulk.csv'
            )
        )

        if (!file.exists(pb.file)) {
            message(
                '  Missing pseudobulk file; skipping: ',
                pb.file
            )
            next
        }

        e.this = read.csv(
            pb.file,
            row.names=1,
            check.names=FALSE
        )

        m = meta.sub[
            meta.sub$region == region,
            ,
            drop=FALSE
        ]

        if (ncol(e.this) == 0 || nrow(e.this) == 0) {
            message('  ', region, ': empty pseudobulk matrix; skipping')
            next
        }

        sample.ids = colnames(e.this)

        m = m[match(sample.ids, m$animal_id), , drop=FALSE]

        if (any(is.na(m$animal_id))) {
            missing.samples = sample.ids[is.na(m$animal_id)]

            stop(
                'Metadata missing for pseudobulk samples in ',
                subtype.name,
                ', region ',
                region,
                ': ',
                paste(missing.samples, collapse=', ')
            )
        }

        rownames(m) = m$animal_id

        if (nrow(m) < 5) {
            message(
                '  ',
                region,
                ': fewer than 5 animals; skipping'
            )
            next
        }

        c.this = covs[
            sapply(
                covs,
                function(cv) {
                    length(unique(m[[cv]])) > 1
                }
            )
        ]

        c.disp.this = disp.covs[
            sapply(
                disp.covs,
                function(cv) {
                    length(unique(m[[cv]])) > 1
                }
            )
        ]

        if (length(c.this) == 0) {
            message(
                '  ',
                region,
                ': no variable mean-model covariates; skipping'
            )
            next
        }

        if (length(c.disp.this) == 0) {
            message(
                '  ',
                region,
                ': no variable dispersion covariates; skipping'
            )
            next
        }

        message(
            '  ',
            region,
            ': ',
            nrow(m),
            ' animals | mean: ',
            paste(c.this, collapse='+'),
            ' | dispersion: ',
            paste(c.disp.this, collapse='+')
        )

        m.this = m[
            ,
            unique(c(c.this, c.disp.this)),
            drop=FALSE
        ]

        clus = makeCluster(n.cores)
        registerDoParallel(cores=n.cores)

        clusterExport(
            clus,
            varlist=c(
                'e.this',
                'm.this',
                'c.this',
                'c.disp.this'
            ),
            envir=environment()
        )

        if (opt$dispersion == 'age') {

            result.region = try(
                t(
                    parApply(
                        clus,
                        as.matrix(e.this),
                        1,
                        function(y) {
                            require(dglm)

                            d = m.this
                            d$e = y

                            res = try(
                                dglm(
                                    as.formula(
                                        paste(
                                            'e ~',
                                            paste(
                                                c.this,
                                                collapse=' + '
                                            )
                                        )
                                    ),
                                    as.formula(
                                        paste(
                                            '~',
                                            paste(
                                                c.disp.this,
                                                collapse=' + '
                                            )
                                        )
                                    ),
                                    family=gaussian(),
                                    dlink='log',
                                    data=d
                                ),
                                silent=TRUE
                            )

                            if ('try-error' %in% class(res)) {
                                return(c(NA, NA, NA))
                            }

                            disp.sum = coef(
                                summary(res)$dispersion.summary
                            )

                            if (c.disp.this[1] %in%
                                rownames(disp.sum)) {
                                return(
                                    disp.sum[
                                        c.disp.this[1],
                                        c(1, 2, 4)
                                    ]
                                )
                            }

                            return(c(NA, NA, NA))
                        }
                    )
                ),
                silent=TRUE
            )

            if ('try-error' %in% class(result.region)) {
                result.region = NULL
            }

            if (!is.null(result.region)) {
                colnames(result.region) = c(
                    'beta',
                    'bvar',
                    'pval'
                )
                out[[region]] = result.region
            }

        } else {

            result.region = try(
                t(
                    parApply(
                        clus,
                        as.matrix(e.this),
                        1,
                        function(y) {
                            require(dglm)

                            d = m.this
                            d$e = y

                            res = try(
                                dglm(
                                    as.formula(
                                        paste(
                                            'e ~',
                                            paste(
                                                c.this,
                                                collapse=' + '
                                            )
                                        )
                                    ),
                                    as.formula(
                                        paste(
                                            '~',
                                            paste(
                                                c.disp.this,
                                                collapse=' + '
                                            )
                                        )
                                    ),
                                    family=gaussian(),
                                    dlink='log',
                                    data=d
                                ),
                                silent=TRUE
                            )

                            if ('try-error' %in% class(res)) {
                                return(
                                    c(
                                        NA, NA, NA,
                                        NA, NA, NA
                                    )
                                )
                            }

                            disp.sum = coef(
                                summary(res)$dispersion.summary
                            )

                            age.row = if (
                                'age' %in% rownames(disp.sum)
                            ) {
                                disp.sum[
                                    'age',
                                    c(1, 2, 4)
                                ]
                            } else {
                                c(NA, NA, NA)
                            }

                            sex.row = if (
                                any(grepl(
                                    'sex',
                                    rownames(disp.sum)
                                ))
                            ) {
                                sex.name = rownames(disp.sum)[
                                    grep(
                                        'sex',
                                        rownames(disp.sum)
                                    )[1]
                                ]

                                disp.sum[
                                    sex.name,
                                    c(1, 2, 4)
                                ]
                            } else {
                                c(NA, NA, NA)
                            }

                            return(c(age.row, sex.row))
                        }
                    )
                ),
                silent=TRUE
            )

            if ('try-error' %in% class(result.region)) {
                result.region = NULL
            }

            if (!is.null(result.region)) {
                colnames(result.region) = c(
                    'beta_age',
                    'bvar_age',
                    'pval_age',
                    'beta_sex',
                    'bvar_sex',
                    'pval_sex'
                )
                out[[region]] = result.region
            }
        }

        stopCluster(clus)
        message('  ', region, ': done')
    }

    valid.out = out[!sapply(out, is.null)]

    if (length(valid.out) == 0) {
        message(
            'No valid regions for ',
            subtype.name,
            '; no result written'
        )
        next
    }

    genes = Reduce(
        union,
        lapply(valid.out, rownames)
    )

    if (opt$dispersion == 'age') {

        stats = c(
            'beta',
            'bvar',
            'pval',
            'qval'
        )

        arr = array(
            NA,
            dim=c(
                length(genes),
                length(stats),
                length(region.levels)
            ),
            dimnames=list(
                genes,
                stats,
                region.levels
            )
        )

        for (region in region.levels) {

            if (is.null(out[[region]])) {
                next
            }

            g = rownames(out[[region]])

            arr[g, 'beta', region] =
                out[[region]][, 'beta']

            arr[g, 'bvar', region] =
                out[[region]][, 'bvar']

            arr[g, 'pval', region] =
                out[[region]][, 'pval']

            arr[g, 'qval', region] =
                p.adjust(
                    out[[region]][, 'pval'],
                    method='fdr'
                )
        }

    } else {

        stats = c(
            'beta_age',
            'bvar_age',
            'pval_age',
            'qval_age',
            'beta_sex',
            'bvar_sex',
            'pval_sex',
            'qval_sex'
        )

        arr = array(
            NA,
            dim=c(
                length(genes),
                length(stats),
                length(region.levels)
            ),
            dimnames=list(
                genes,
                stats,
                region.levels
            )
        )

        for (region in region.levels) {

            if (is.null(out[[region]])) {
                next
            }

            g = rownames(out[[region]])

            for (s in c(
                'beta_age',
                'bvar_age',
                'pval_age',
                'beta_sex',
                'bvar_sex',
                'pval_sex'
            )) {
                arr[g, s, region] =
                    out[[region]][, s]
            }

            arr[g, 'qval_age', region] =
                p.adjust(
                    out[[region]][, 'pval_age'],
                    method='fdr'
                )

            arr[g, 'qval_sex', region] =
                p.adjust(
                    out[[region]][, 'pval_sex'],
                    method='fdr'
                )
        }
    }

    human.symbols = sapply(genes, get_symbol)
    names(human.symbols) = genes

    result.object = list(
        array=arr,
        human_symbols=human.symbols,
        run_info=list(
            covariates=covs,
            disp_covariates=disp.covs,
            dispersion_mode=opt$dispersion,
            date=Sys.time(),
            cell_type=cell.type,
            louvain=louvain,
            louvain_column=louvain.col,
            input_dir=datadir,
            output_dir=outdir
        )
    )

    # Never overwrite an existing result. This is checked again immediately
    # before writing in case another job created the file meanwhile.
    if (file.exists(result.file)) {
        message(
            'Result appeared during this run; skipping write: ',
            result.file
        )
        next
    }

    saveRDS(
        result.object,
        file=result.file
    )

    message('Saved: ', result.file)
}

message('')
message('Louvain-level Astrocyte DGLM analysis complete.')
