#!/usr/bin/env Rscript

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

suppressPackageStartupMessages({
    library(optparse)
    library(parallel)
    library(doParallel)
    library(dglm)
    library(edgeR)
})

option_list = list(
    make_option('--cell_type',   type='character', help='cell type to analyze'),
    make_option('--outdir',      type='character', help='output directory'),
    make_option('--covariates',  type='character', default='full',
                help='age_sex or full (age+sex+mean_n_umi+n_cells) -- mean model'),
    make_option('--dispersion',  type='character', default='age',
                help='age or age_sex -- dispersion model')
)
opt = parse_args(OptionParser(option_list=option_list))

if (is.null(opt$cell_type) || is.null(opt$outdir)) {
    stop('--cell_type and --outdir are required')
}

if (!opt$covariates %in% c('age_sex', 'full')) {
    stop("--covariates must be 'age_sex' or 'full'")
}

if (!opt$dispersion %in% c('age', 'age_sex')) {
    stop("--dispersion must be 'age' or 'age_sex'")
}

cell.type = opt$cell_type
outdir    = opt$outdir
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

covs = if (opt$covariates == 'age_sex') {
    c('age', 'sex')
} else {
    c('age', 'sex', 'mean_n_umi', 'n_cells')
}

disp.covs = if (opt$dispersion == 'age_sex') c('age', 'sex') else c('age')

message('Cell type: ', cell.type)
message('Mean covariates: ', paste(covs, collapse=' + '))
message('Dispersion covariates: ', paste(disp.covs, collapse=' + '))
message('Response: TMM-normalized log2-CPM from CPM-filtered pseudobulk counts')
message('Output dir: ', outdir)

orth = read.csv('/scratch/easmit31/data/human-macaque-orthologs/ensembl113_mmul10_macaque_human.csv')
orth = orth[
    orth$Human.homology.type == 'ortholog_one2one' &
    orth$Human.gene.name != '',
    c('Gene.stable.ID', 'Human.gene.name')
]
orth = orth[!duplicated(orth$Gene.stable.ID), ]
rownames(orth) = orth$Gene.stable.ID

gene.names.file = file.path(outdir, paste0(cell.type, '_gene_names.csv'))
gene.names.df = if (file.exists(gene.names.file)) {
    read.csv(gene.names.file, row.names=1)
} else {
    NULL
}

get_symbol = function(g) {
    if (g %in% rownames(orth)) return(orth[g, 'Human.gene.name'])
    if (!is.null(gene.names.df) && g %in% rownames(gene.names.df)) {
        nm = gene.names.df[g, 'external_gene_name']
        if (!is.na(nm) && nm != g && nm != '') return(nm)
    }
    g
}

meta.all = read.csv(
    file.path(outdir, paste0(cell.type, '_metadata.csv')),
    stringsAsFactors=FALSE
)

meta.all$sex = as.factor(meta.all$sex)
meta.all$age = as.numeric(meta.all$age)

if ('mean_n_umi' %in% colnames(meta.all)) {
    meta.all$mean_n_umi = as.numeric(meta.all$mean_n_umi)
}

if ('n_cells' %in% colnames(meta.all)) {
    meta.all$n_cells = as.numeric(meta.all$n_cells)
}

missing = setdiff(union(covs, disp.covs), colnames(meta.all))
if (length(missing) > 0) {
    stop('Missing columns in metadata: ', paste(missing, collapse=', '))
}

out = vector('list', length(region.levels))
names(out) = region.levels

mean.out = vector('list', length(region.levels))
names(mean.out) = region.levels

for (region in region.levels) {
    filtered.file = file.path(
        outdir,
        'filtered',
        paste0(cell.type, '_', region, '_filtered_cutoff0.5.csv')
    )

    filtered.file.alt = file.path(
        outdir,
        paste0(cell.type, '_', region, '_filtered_cutoff0.5.csv')
    )

    raw.file = file.path(
        outdir,
        paste0(cell.type, '_', region, '_pseudobulk.csv')
    )

    pb.file = if (file.exists(filtered.file)) {
        filtered.file
    } else if (file.exists(filtered.file.alt)) {
        filtered.file.alt
    } else if (file.exists(raw.file)) {
        warning('No filtered file found for ', region,
                '; using raw pseudobulk counts: ', raw.file)
        raw.file
    } else {
        next
    }

    counts = read.csv(pb.file, row.names=1, check.names=FALSE)
    counts = as.matrix(counts)
    storage.mode(counts) = 'numeric'

    # Restrict to autosome + X genes (drops Y-linked and mitochondrial genes,
    # which have known confound/technical-noise issues for this analysis).
    autosome_x_file = '/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv'
    if (file.exists(autosome_x_file)) {
        ax = read.csv(autosome_x_file, stringsAsFactors=FALSE)
        id_col = ax[[grep('ENSMMUG|ensembl|gene_id', colnames(ax), ignore.case=TRUE, value=TRUE)[1]]]
        keep_genes = intersect(rownames(counts), id_col)
        n_before = nrow(counts)
        counts = counts[keep_genes, , drop=FALSE]
        message('    autosome+X filter: ', n_before, ' -> ', nrow(counts), ' genes')
    } else {
        warning('autosome_x_genes.csv not found -- proceeding WITHOUT sex-chromosome/MT gene filter')
    }

    m = meta.all[meta.all$region == region, , drop=FALSE]
    m = m[match(colnames(counts), m$animal_id), , drop=FALSE]
    rownames(m) = m$animal_id

    if (!identical(colnames(counts), rownames(m))) {
        stop('Metadata/count alignment failed for ', cell.type, ' ', region)
    }

    if (nrow(m) < 5) next

    y = DGEList(counts=counts)
    y = calcNormFactors(y, method='TMM')
    logCPM = cpm(y, log=TRUE, prior.count=0.5)

    c.this = covs[sapply(covs, function(cv) length(unique(m[[cv]])) > 1)]
    c.disp.this = disp.covs[sapply(disp.covs, function(cv) length(unique(m[[cv]])) > 1)]

    if (!('age' %in% c.disp.this)) {
        warning(region, ': age has no variation; skipping')
        next
    }

    message(
        '  ', region, ': ', nrow(m), ' animals | ', nrow(counts), ' genes',
        ' | mean: ', paste(c.this, collapse=' + '),
        ' | disp: ', paste(c.disp.this, collapse=' + '),
        ' | input: ', basename(pb.file)
    )

    m.this = m[, unique(c(c.this, c.disp.this)), drop=FALSE]

    clus = makeCluster(n.cores)
    registerDoParallel(clus)

    clusterExport(
        clus,
        varlist=c('logCPM', 'm.this', 'c.this', 'c.disp.this'),
        envir=environment()
    )

    results = t(parApply(clus, logCPM, 1, function(expr) {
        suppressPackageStartupMessages(library(dglm))

        d = m.this
        d$e = as.numeric(expr)

        res = try(
            suppressWarnings(
                dglm(
                    as.formula(paste('e ~', paste(c.this, collapse=' + '))),
                    as.formula(paste('~', paste(c.disp.this, collapse=' + '))),
                    family=gaussian(),
                    dlink='log',
                    data=d
                )
            ),
            silent=TRUE
        )

        if (inherits(res, 'try-error')) {
            return(c(
                disp_beta_age=NA, disp_se_age=NA, disp_pval_age=NA,
                disp_beta_sex=NA, disp_se_sex=NA, disp_pval_sex=NA,
                mean_beta_age=NA, mean_se_age=NA, mean_pval_age=NA,
                mean_beta_sex=NA, mean_se_sex=NA, mean_pval_sex=NA
            ))
        }

        disp.tab = try(coef(summary(res)$dispersion.summary), silent=TRUE)
        mean.tab = try(coef(summary(res)), silent=TRUE)

        get_stat = function(tab, term, p.column) {
            if (inherits(tab, 'try-error') ||
                !(term %in% rownames(tab)) ||
                !all(c('Estimate', 'Std. Error', p.column) %in% colnames(tab))) {
                return(c(NA_real_, NA_real_, NA_real_))
            }

            c(
                unname(tab[term, 'Estimate']),
                unname(tab[term, 'Std. Error']),
                unname(tab[term, p.column])
            )
        }

        disp.age = get_stat(disp.tab, 'age', 'Pr(>|z|)')

        sex.term = grep('^sex', rownames(disp.tab), value=TRUE)
        disp.sex = if (length(sex.term) > 0) {
            get_stat(disp.tab, sex.term[1], 'Pr(>|z|)')
        } else {
            c(NA_real_, NA_real_, NA_real_)
        }

        mean.age = get_stat(mean.tab, 'age', 'Pr(>|t|)')

        mean.sex.term = grep('^sex', rownames(mean.tab), value=TRUE)
        mean.sex = if (length(mean.sex.term) > 0) {
            get_stat(mean.tab, mean.sex.term[1], 'Pr(>|t|)')
        } else {
            c(NA_real_, NA_real_, NA_real_)
        }

        c(
            disp_beta_age=disp.age[1],
            disp_se_age=disp.age[2],
            disp_pval_age=disp.age[3],
            disp_beta_sex=disp.sex[1],
            disp_se_sex=disp.sex[2],
            disp_pval_sex=disp.sex[3],
            mean_beta_age=mean.age[1],
            mean_se_age=mean.age[2],
            mean_pval_age=mean.age[3],
            mean_beta_sex=mean.sex[1],
            mean_se_sex=mean.sex[2],
            mean_pval_sex=mean.sex[3]
        )
    }))

    stopCluster(clus)

    rownames(results) = rownames(logCPM)

    if (opt$dispersion == 'age') {
        out[[region]] = results[, c('disp_beta_age', 'disp_se_age', 'disp_pval_age'), drop=FALSE]
        colnames(out[[region]]) = c('beta', 'se', 'pval')
    } else {
        out[[region]] = results[, c(
            'disp_beta_age', 'disp_se_age', 'disp_pval_age',
            'disp_beta_sex', 'disp_se_sex', 'disp_pval_sex'
        ), drop=FALSE]
        colnames(out[[region]]) = c(
            'beta_age', 'se_age', 'pval_age',
            'beta_sex', 'se_sex', 'pval_sex'
        )
    }

    mean.out[[region]] = results[, c(
        'mean_beta_age', 'mean_se_age', 'mean_pval_age',
        'mean_beta_sex', 'mean_se_sex', 'mean_pval_sex'
    ), drop=FALSE]
    colnames(mean.out[[region]]) = c('beta_age', 'se_age', 'pval_age', 'beta_sex', 'se_sex', 'pval_sex')

    message('  ', region, ': done')
}

non.null.out = out[!sapply(out, is.null)]

if (length(non.null.out) == 0) {
    stop('No regions produced results. Check --outdir and input file locations.')
}

genes = Reduce(union, lapply(non.null.out, rownames))

if (opt$dispersion == 'age') {
    stats = c('beta', 'bvar', 'pval', 'qval')

    arr = array(
        NA,
        dim=c(length(genes), length(stats), length(region.levels)),
        dimnames=list(genes, stats, region.levels)
    )

    for (r in region.levels) {
        if (is.null(out[[r]])) next

        g = rownames(out[[r]])

        arr[g, 'beta', r] = out[[r]][, 'beta']
        arr[g, 'bvar', r] = out[[r]][, 'se']
        arr[g, 'pval', r] = out[[r]][, 'pval']
        arr[g, 'qval', r] = p.adjust(out[[r]][, 'pval'], method='BH')
    }
} else {
    stats = c(
        'beta_age', 'bvar_age', 'pval_age', 'qval_age',
        'beta_sex', 'bvar_sex', 'pval_sex', 'qval_sex'
    )

    arr = array(
        NA,
        dim=c(length(genes), length(stats), length(region.levels)),
        dimnames=list(genes, stats, region.levels)
    )

    for (r in region.levels) {
        if (is.null(out[[r]])) next

        g = rownames(out[[r]])

        arr[g, 'beta_age', r] = out[[r]][, 'beta_age']
        arr[g, 'bvar_age', r] = out[[r]][, 'se_age']
        arr[g, 'pval_age', r] = out[[r]][, 'pval_age']

        arr[g, 'beta_sex', r] = out[[r]][, 'beta_sex']
        arr[g, 'bvar_sex', r] = out[[r]][, 'se_sex']
        arr[g, 'pval_sex', r] = out[[r]][, 'pval_sex']

        arr[g, 'qval_age', r] = p.adjust(out[[r]][, 'pval_age'], method='BH')
        arr[g, 'qval_sex', r] = p.adjust(out[[r]][, 'pval_sex'], method='BH')
    }
}

mean.stats = c('beta_age', 'se_age', 'pval_age', 'qval_age', 'beta_sex', 'se_sex', 'pval_sex', 'qval_sex')

mean.arr = array(
    NA,
    dim=c(length(genes), length(mean.stats), length(region.levels)),
    dimnames=list(genes, mean.stats, region.levels)
)

for (r in region.levels) {
    if (is.null(mean.out[[r]])) next

    g = rownames(mean.out[[r]])

    mean.arr[g, 'beta_age', r] = mean.out[[r]][, 'beta_age']
    mean.arr[g, 'se_age', r] = mean.out[[r]][, 'se_age']
    mean.arr[g, 'pval_age', r] = mean.out[[r]][, 'pval_age']
    mean.arr[g, 'qval_age', r] = p.adjust(mean.out[[r]][, 'pval_age'], method='BH')

    mean.arr[g, 'beta_sex', r] = mean.out[[r]][, 'beta_sex']
    mean.arr[g, 'se_sex', r] = mean.out[[r]][, 'se_sex']
    mean.arr[g, 'pval_sex', r] = mean.out[[r]][, 'pval_sex']
    mean.arr[g, 'qval_sex', r] = p.adjust(mean.out[[r]][, 'pval_sex'], method='BH')
}

human.symbols = sapply(genes, get_symbol)
names(human.symbols) = genes

saveRDS(
    list(
        array = arr,
        mean_age_array = mean.arr,
        human_symbols = human.symbols,
        run_info = list(
            covariates = covs,
            disp_covariates = disp.covs,
            dispersion_mode = opt$dispersion,
            response = 'TMM-normalized log2-CPM from CPM-filtered pseudobulk counts',
            dispersion_pval_source = 'summary(res)$dispersion.summary, hardcoded Gamma dispersion=2, Wald z-test -- REVERTED for comparison, NOT the corrected method',
            bvar_contents = 'raw dispersion standard errors for mashr Shat; bvar is a legacy field name',
            mean_age_array_contents = 'mean-model age beta, SE, p-value, and BH q-value',
            date = Sys.time(),
            cell_type = cell.type,
            outdir = outdir
        )
    ),
    file=file.path(outdir, paste0(cell.type, '_dglm_results.rds'))
)

message('Saved to: ', outdir)
message('done.')
