#!/usr/bin/env Rscript
# dglm_mashr.R
#
# Runs mashr on DGLM results.
# Two modes:
#   --mode per_ct (default): one mashr per cell type across regions (original approach)
#   --mode combined: one mashr across all cell types x regions simultaneously
#
# SAME as Chiou et al. (per_ct mode):
#   - Strong subset, random subset, PCA+ED covariances, canonical fallback
#
# DEVIATIONS (both modes):
#   - estimate_null_correlation_simple (API change)
#   - z-score scaling (pseudobulk compression)
#   - tryCatch fallbacks
#   - Single-region cell types skip mashr (per_ct only)
#
# Usage:
#   Rscript dglm_mashr.R --cell_type microglia
#   Rscript dglm_mashr.R --mode combined

source('/scratch/easmit31/variability/dglm/scripts/_include_options.R')

library(optparse)
library(mashr)
library(abind)

option_list = list(
    make_option('--mode',        type='character', default='per_ct',
                help='per_ct or combined'),
    make_option('--cell_type',   type='character', default=NULL,
                help='cell type (required for per_ct mode)'),
    make_option('--checkpoints', type='character',
                default='/scratch/easmit31/variability/dglm/checkpoints')
)
opt = parse_args(OptionParser(option_list=option_list))

if (opt$mode == 'per_ct' && is.null(opt$cell_type)) {
    stop('--cell_type required for per_ct mode')
}

# ── shared mashr fitting function ─────────────────────────────────────────
run_mashr = function(Bhat, Shat, dglm.qval, label) {
    message('Running mashr for: ', label)
    message('  dims: ', nrow(Bhat), ' genes x ', ncol(Bhat), ' conditions')

    # z-score scaling
    Zhat   = Bhat / Shat
    Zhat[Shat == 1000] = 0
    z.sd   = sd(Zhat[Shat < 1000], na.rm=TRUE)
    message('  z-score SD (scaling factor): ', round(z.sd, 4))
    Bhat.sc = Bhat / z.sd
    Shat.sc = Shat / z.sd

    mash.data = mash_set_data(Bhat.sc, Shat.sc)

    # strong subset: sig in >= 1/3 of conditions
    strong.subset = which(apply(dglm.qval, 1, function(x) {
        sum(x < fsr.cutoff, na.rm=TRUE)
    }) >= max(1, ncol(Bhat) / 3))
    message('  strong subset: ', length(strong.subset))

    set.seed(seed)
    random.subset = sample(1:nrow(Bhat.sc), ceiling(nrow(Bhat.sc)/2))

    # null correlation
    temp = mash_set_data(Bhat.sc[random.subset,], Shat.sc[random.subset,])
    Vhat = tryCatch({
        message('  estimating null correlation')
        estimate_null_correlation_simple(temp)
    }, error = function(e) {
        message('  null correlation failed: ', e$message)
        NULL
    })
    rm('temp')

    if (!is.null(Vhat)) {
        mash.random = mash_set_data(Bhat.sc[random.subset,], Shat.sc[random.subset,], V=Vhat)
    } else {
        mash.random = mash_set_data(Bhat.sc[random.subset,], Shat.sc[random.subset,])
    }

    U.c.r = cov_canonical(mash.random)
    Ulist = tryCatch({
        if (length(strong.subset) > 0 && ncol(Bhat) >= 2) {
            if (!is.null(Vhat)) {
                mash.strong = mash_set_data(Bhat.sc[strong.subset,], Shat.sc[strong.subset,], V=Vhat)
            } else {
                mash.strong = mash_set_data(Bhat.sc[strong.subset,], Shat.sc[strong.subset,])
            }
            n.pcs = min(5, ncol(Bhat)-1, length(strong.subset)-1)
            if (n.pcs < 1) stop('not enough PCs')
            U.pca = cov_pca(mash.strong, n.pcs)
            U.ed  = cov_ed(mash.strong, U.pca)
            message('  using data-driven + canonical covariances')
            c(U.ed, U.c.r)
        } else {
            stop('no strong subset or too few conditions')
        }
    }, error = function(e) {
        message('  data-driven covariances failed, using canonical only: ', e$message)
        U.c.r
    })

    message('  fitting mash on random subset...')
    now = Sys.time()
    m.r = tryCatch(
        mash(mash.random, Ulist=Ulist, outputlevel=1),
        error = function(e) {
            message('  mash fit failed, retrying canonical only: ', e$message)
            mash(mash.random, Ulist=U.c.r, outputlevel=1)
        }
    )
    message('  time: ', format(Sys.time()-now))

    message('  applying to all data...')
    now = Sys.time()
    m = tryCatch(
        mash(mash.data, g=get_fitted_g(m.r), fixg=TRUE),
        error = function(e) {
            message('  final mash failed: ', e$message)
            NULL
        }
    )
    message('  time: ', format(Sys.time()-now))
    if (!is.null(m)) message('  log likelihood: ', format(get_loglik(m), digits=10))
    return(list(mash=m, z_scale=z.sd))
}

# ═══════════════════════════════════════════════════════════════════════════
# PER_CT MODE (original)
# ═══════════════════════════════════════════════════════════════════════════

if (opt$mode == 'per_ct') {
    cell.type = opt$cell_type
    message('Mode: per_ct | Cell type: ', cell.type)

    in.file       = file.path(opt$checkpoints, paste0(cell.type, '_dglm_results.rds'))
    dglm.obj      = readRDS(in.file)
    dglm.results  = dglm.obj$array
    human.symbols = dglm.obj$human_symbols

    Bhat = dglm.results[,'beta',]
    Shat = sqrt(dglm.results[,'bvar',])
    Bhat[is.na(Bhat)] = 0
    Shat[is.na(Shat)] = 1000

    regions.with.data = colnames(Bhat)[apply(Bhat, 2, function(x) sum(x != 0) > 0)]
    message('Regions with data: ', paste(regions.with.data, collapse=', '))
    Bhat = Bhat[, regions.with.data, drop=FALSE]
    Shat = Shat[, regions.with.data, drop=FALSE]

    # filter extreme betas
    extreme = apply(Bhat, 1, function(x) any(abs(x) > 100, na.rm=TRUE))
    if (sum(extreme) > 0) {
        message('Filtering ', sum(extreme), ' genes with extreme beta')
        Bhat         = Bhat[!extreme,, drop=FALSE]
        Shat         = Shat[!extreme,, drop=FALSE]
        dglm.results = dglm.results[!extreme,,, drop=FALSE]
        human.symbols = human.symbols[!extreme]
    }

    out.file = file.path(opt$checkpoints, paste0(cell.type, '_dglm_mashr_results.rds'))

    if (length(regions.with.data) < 2) {
        message('Only one region — skipping mashr')
        saveRDS(list(mash=NULL, dglm_results=dglm.results,
                     human_symbols=human.symbols, regions=regions.with.data, z_scale=1),
                file=out.file)
        message('Saved: ', out.file)
        quit(save='no', status=0)
    }

    dglm.qval = dglm.results[,'qval', regions.with.data, drop=FALSE]
    result    = run_mashr(Bhat, Shat, dglm.qval, cell.type)

    saveRDS(list(mash=result$mash, dglm_results=dglm.results,
                 human_symbols=human.symbols, regions=regions.with.data,
                 z_scale=result$z_scale),
            file=out.file)
    message('Saved: ', out.file)
}

# ═══════════════════════════════════════════════════════════════════════════
# COMBINED MODE — all cell types x regions in one mashr
# ═══════════════════════════════════════════════════════════════════════════

if (opt$mode == 'combined') {
    message('Mode: combined — pooling all cell types x regions')

    # collect all genes across all cell types
    all.genes = character(0)
    for (ct in cell.type.levels) {
        rds.file = file.path(opt$checkpoints, paste0(ct, '_dglm_results.rds'))
        if (!file.exists(rds.file)) next
        obj = readRDS(rds.file)
        all.genes = union(all.genes, dimnames(obj$array)[[1]])
    }
    message('Total unique genes: ', length(all.genes))

    # build condition labels: ct|region
    conditions = character(0)
    for (ct in cell.type.levels) {
        rds.file = file.path(opt$checkpoints, paste0(ct, '_dglm_results.rds'))
        if (!file.exists(rds.file)) next
        obj     = readRDS(rds.file)
        regions = dimnames(obj$array)[[3]]
        for (r in regions) {
            beta.col = obj$array[,'beta',r]
            if (sum(!is.na(beta.col) & beta.col != 0) > 0) {
                conditions = c(conditions, paste(ct, r, sep='|'))
            }
        }
    }
    message('Total conditions (ct x region): ', length(conditions))

    # initialize Bhat and Shat
    Bhat      = matrix(0,    nrow=length(all.genes), ncol=length(conditions),
                       dimnames=list(all.genes, conditions))
    Shat      = matrix(1000, nrow=length(all.genes), ncol=length(conditions),
                       dimnames=list(all.genes, conditions))
    dglm.qval = matrix(1,    nrow=length(all.genes), ncol=length(conditions),
                       dimnames=list(all.genes, conditions))

    # collect human symbols and dglm results
    human.symbols.list = list()
    dglm.results.list  = list()

    for (ct in cell.type.levels) {
        rds.file = file.path(opt$checkpoints, paste0(ct, '_dglm_results.rds'))
        if (!file.exists(rds.file)) next
        obj          = readRDS(rds.file)
        dglm.results = obj$array
        human.sym    = obj$human_symbols

        # filter extreme betas
        beta.all = dglm.results[,'beta',, drop=FALSE]
        extreme  = apply(beta.all, 1, function(x) any(abs(x) > 100, na.rm=TRUE))
        if (sum(extreme) > 0) {
            dglm.results = dglm.results[!extreme,,, drop=FALSE]
            human.sym    = human.sym[!extreme]
        }

        human.symbols.list[[ct]] = human.sym
        dglm.results.list[[ct]]  = dglm.results
        genes.ct = dimnames(dglm.results)[[1]]

        for (r in dimnames(dglm.results)[[3]]) {
            cond = paste(ct, r, sep='|')
            if (!cond %in% conditions) next
            beta = dglm.results[,'beta', r]
            bvar = dglm.results[,'bvar', r]
            qval = dglm.results[,'qval', r]

            # fill in genes present in this ct
            idx = genes.ct[genes.ct %in% all.genes]
            Bhat[idx, cond] = ifelse(is.na(beta[idx]), 0,    beta[idx])
            Shat[idx, cond] = ifelse(is.na(bvar[idx]) | bvar[idx] <= 0, 1000,
                                     sqrt(bvar[idx]))
            dglm.qval[idx, cond] = ifelse(is.na(qval[idx]), 1, qval[idx])
        }
    }

    # build unified human symbols (priority: first non-ensembl symbol found)
    human.symbols = setNames(all.genes, all.genes)
    for (ct in names(human.symbols.list)) {
        sym = human.symbols.list[[ct]]
        for (g in names(sym)) {
            if (g %in% names(human.symbols) &&
                human.symbols[g] == g &&
                !is.na(sym[g]) && sym[g] != g && sym[g] != '') {
                human.symbols[g] = sym[g]
            }
        }
    }

    # run mashr
    result = run_mashr(Bhat, Shat, dglm.qval, 'combined')

    out.file = file.path(opt$checkpoints, 'combined_dglm_mashr_results.rds')
    saveRDS(list(
        mash          = result$mash,
        Bhat          = Bhat,
        Shat          = Shat,
        dglm_qval     = dglm.qval,
        human_symbols = human.symbols,
        conditions    = conditions,
        dglm_results  = dglm.results.list,
        z_scale       = result$z_scale
    ), file=out.file)
    message('Saved: ', out.file)
}

message('done.')
