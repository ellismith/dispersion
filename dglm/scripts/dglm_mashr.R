#!/usr/bin/env Rscript
# dglm_mashr.R
#
# Runs mashr on DGLM results.
# Two modes:
#   --mode per_ct (default): one mashr per cell type across regions (original approach)
#   --mode combined: one mashr across all cell types x regions simultaneously
#
# Auto-detects single-term vs two-term dglm_model.R output from the array's
# column names:
#   'beta'/'bvar'/'pval'/'qval'             -> single term (--dispersion age)
#   'beta_age'/'beta_sex'/... (both)        -> two terms (--dispersion age_sex)
#     -> runs mashr TWICE (once per term), saves two output files
# No flag needed here -- whatever dglm_model.R actually produced is what
# this reads.
#
# SAME as Chiou et al. (per_ct mode):
#   - Strong subset, random subset, PCA+ED covariances, canonical fallback
#
# DEVIATIONS (both modes):
#   - tryCatch fallbacks (paper ran this once by hand; we run it unattended
#     across 12 cell types x 3 model versions)
#   - Single-region cell types skip mashr (per_ct only)
#
# --shat_mode: how dglm's raw 'bvar' output becomes mashr's Shat input.
#   bvar is confirmed to be the dispersion submodel's Std. Error directly
#   (from dglm's own summary object, extracted identically in Chiou's
#   original dglm_model.R and this one -- confirmed by direct comparison
#   against Chiou's uploaded source).
#     raw  (default): Shat = bvar, used as-is.
#     sqrt          : Shat = sqrt(bvar), matches Chiou et al.'s own
#                      dglm_mashr.R line 10 exactly (confirmed against
#                      their uploaded source).
#   A raw-p-value calibration check (fraction of |beta/bvar| exceeding the
#   1.96 z-threshold, which should be ~5% under the null) found: raw bvar
#   gives ~25% (too liberal), sqrt(bvar) gives ~0.05% (far too
#   conservative). This script does not attempt to correct either -- it
#   exists so both of Chiou's-method-vs-less-conservative can be run and
#   compared directly on identical DGLM input.
#
# --random_subset_n / --gridmult: added for large condition-count runs
# (e.g. louvain combined mode, 100+ conditions), where the initial mash()
# fit on a random subset becomes the bottleneck -- both cost roughly
# (genes in subset) x (mixture components) x (conditions)^3, and neither
# knob was previously exposed. Defaults exactly reproduce prior behavior.
#
# Usage:
#   Rscript dglm_mashr.R --cell_type microglia
#   Rscript dglm_mashr.R --mode combined
#   Rscript dglm_mashr.R --mode combined --shat_mode sqrt
#   Rscript dglm_mashr.R --mode combined --checkpoints dglm/disp_age_sex__mean_full  # auto two-term
#   Rscript dglm_mashr.R --mode combined --fast --random_subset_n 8000 --gridmult 2.5  # large condition count

source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')

library(optparse)
library(mashr)
library(abind)

option_list = list(
    make_option('--mode',        type='character', default='per_ct',
                help='per_ct or combined'),
    make_option('--cell_type',   type='character', default=NULL,
                help='cell type (required for per_ct mode)'),
    make_option('--checkpoints', type='character',
                default='/scratch/easmit31/dispersion/dglm/checkpoints'),
    make_option('--shat_mode',   type='character', default='raw',
                help="raw (default; Shat=bvar as-is) or sqrt (Shat=sqrt(bvar), matches Chiou et al.'s original script exactly)"),
    make_option('--fast',        action='store_true', default=FALSE,
                help='skip data-driven covariance estimation (cov_pca/cov_ed), use canonical covariances only. Much faster on large condition counts (100+), at the cost of some statistical refinement -- data-driven covariances pick up correlation patterns canonical ones miss. Recommended for combined-mode runs with many conditions where the full fit is taking too long.'),
    make_option('--random_subset_n', type='integer', default=NULL,
                help='cap the size of the random gene subset used for the initial mash() fit. Default (NULL) preserves existing behavior: half of all genes. A representative random sample for estimating mixture weights does not need to scale with total gene count -- a few thousand is standard. Pass a smaller value (e.g. 5000-10000) for large condition-count runs (100+ conditions) where this step is the bottleneck.'),
    make_option('--gridmult',    type='double', default=sqrt(2),
                help='mash() grid multiplier: larger values use a coarser (fewer-point) grid over effect-size scales, trading resolution for speed. Default sqrt(2) matches mashr\'s own default and prior verified runs. Try 2-3 for large condition-count runs where component count is the bottleneck -- this affects canonical-shape components too, unlike --fast which only skips data-driven ones.'),
    make_option('--gene_keep_list', type='character',
                default='/scratch/easmit31/dispersion/dglm/autosome_x_genes.csv',
                help='CSV (one column: ensembl_gene_id) of genes to keep before mashr -- matches Chiou et al.\'s autosome+X filter (chromosome_name %in% c(1:20,\'X\')), built via build_autosome_x_gene_list.py since biomaRt is unavailable here. Applied to Bhat/Shat/dglm.results/human.symbols together, alongside the extreme-beta filter, so the saved RDS stays internally consistent. Pass an empty string to skip (not recommended -- only for reproducing pre-filter results).'),
    make_option('--only_term',   type='character', default=NULL,
                help='for two-term (age_sex) output: restrict to just this term (age or sex) instead of running both. Use to resume after one term succeeded and the other failed/errored, without redoing the already-saved one. Ignored for single-term output.')
)
opt = parse_args(OptionParser(option_list=option_list))

if (opt$mode == 'per_ct' && is.null(opt$cell_type)) {
    stop('--cell_type required for per_ct mode')
}
if (!opt$shat_mode %in% c('raw', 'sqrt')) {
    stop("--shat_mode must be 'raw' or 'sqrt'")
}
message('Shat mode: ', opt$shat_mode,
        if (opt$shat_mode == 'sqrt') ' (Chiou et al. exact method)' else ' (bvar as-is)')

# ── autosome+X gene filter (matches Chiou et al.: chromosome_name %in%
# c(1:20,'X'), i.e. drop Y/MT/unplaced scaffolds) ──────────────────────────
autosome.x.genes = NULL
if (!is.null(opt$gene_keep_list) && nchar(opt$gene_keep_list) > 0) {
    if (!file.exists(opt$gene_keep_list)) {
        stop('--gene_keep_list not found: ', opt$gene_keep_list,
             ' -- run build_autosome_x_gene_list.py first, or pass --gene_keep_list ""',
             ' to explicitly skip this filter.')
    }
    autosome.x.genes = read.csv(opt$gene_keep_list, stringsAsFactors=FALSE)$ensembl_gene_id
    message('Loaded gene keep-list: ', length(autosome.x.genes), ' autosome+X genes from ', opt$gene_keep_list)
} else {
    message('--gene_keep_list not set -- skipping autosome+X filter (Y/MT/scaffold genes will remain in the analysis)')
}

# ── shat transform helper ──────────────────────────────────────────────────
apply_shat_mode = function(bvar.vals) {
    if (opt$shat_mode == 'sqrt') sqrt(bvar.vals) else bvar.vals
}

# ── detect which term(s) this checkpoint dir's dglm output has ────────────
detect_terms = function(stat.names) {
    if ('beta' %in% stat.names) return('')
    terms = unique(gsub('^beta_', '', grep('^beta_', stat.names, value=TRUE)))
    if (length(terms) == 0) stop('Could not find a beta column in array -- unrecognized dglm_model.R output format')
    terms
}

cols_for_term = function(term) {
    if (term == '') {
        list(beta='beta', bvar='bvar', pval='pval', qval='qval', suffix='')
    } else {
        list(beta=paste0('beta_',term), bvar=paste0('bvar_',term),
             pval=paste0('pval_',term), qval=paste0('qval_',term),
             suffix=paste0('_', term))
    }
}

# ── shared mashr fitting function ─────────────────────────────────────────
run_mashr = function(Bhat, Shat, dglm.qval, label, fast=FALSE,
                      random_subset_n=NULL, gridmult=sqrt(2)) {
    message('Running mashr for: ', label)
    message('  dims: ', nrow(Bhat), ' genes x ', ncol(Bhat), ' conditions')

    mash.data = mash_set_data(Bhat, Shat)

    strong.subset = which(apply(dglm.qval, 1, function(x) {
        sum(x < strong.subset.qval.cutoff, na.rm=TRUE)
    }) >= max(1, ncol(Bhat) / 3))
    message('  strong subset (q<', strong.subset.qval.cutoff, '): ', length(strong.subset))

    set.seed(seed)
    n.random = ceiling(nrow(Bhat)/2)
    if (!is.null(random_subset_n)) {
        n.random = min(n.random, random_subset_n)
    }
    random.subset = sample(1:nrow(Bhat), n.random)
    message('  random subset size: ', n.random,
            if (!is.null(random_subset_n)) ' (capped by --random_subset_n)' else ' (half of all genes, default)')

    temp = mash_set_data(Bhat[random.subset,], Shat[random.subset,])
    Vhat = tryCatch({
        message('  estimating null correlation')
        estimate_null_correlation_simple(temp)
    }, error = function(e) {
        message('  null correlation failed: ', e$message)
        NULL
    })
    rm('temp')

    if (!is.null(Vhat)) {
        mash.random = tryCatch(
            mash_set_data(Bhat[random.subset,], Shat[random.subset,], V=Vhat),
            error = function(e) {
                message('  mash_set_data with estimated V failed (', e$message, ') -- falling back to no V')
                Vhat <<- NULL
                mash_set_data(Bhat[random.subset,], Shat[random.subset,])
            }
        )
    } else {
        mash.random = mash_set_data(Bhat[random.subset,], Shat[random.subset,])
    }

    U.c.r = cov_canonical(mash.random)
    if (fast) {
        message('  --fast: skipping data-driven covariance estimation, using canonical only')
        Ulist = U.c.r
    } else {
        Ulist = tryCatch({
            if (length(strong.subset) > 0 && ncol(Bhat) >= 2) {
                if (!is.null(Vhat)) {
                    mash.strong = mash_set_data(Bhat[strong.subset,], Shat[strong.subset,], V=Vhat)
                } else {
                    mash.strong = mash_set_data(Bhat[strong.subset,], Shat[strong.subset,])
                }
                n.pcs = min(5, ncol(Bhat)-1, length(strong.subset)-1)
                if (n.pcs < 1) stop('not enough PCs')
                message('  computing data-driven covariances (cov_pca/cov_ed) -- this is the expensive step on large condition counts')
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
    }
    message('  total mixture components: ', length(Ulist))
    message('  gridmult: ', gridmult,
            if (gridmult != sqrt(2)) ' (non-default -- coarser grid than mashr default sqrt(2))' else ' (default)')

    message('  fitting mash on random subset (', length(random.subset), ' genes x ', length(Ulist), ' components -- this can take a while on large condition counts)...')
    now = Sys.time()
    m.r = tryCatch(
        mash(mash.random, Ulist=Ulist, outputlevel=1, gridmult=gridmult),
        error = function(e) {
            message('  mash fit failed, retrying canonical only: ', e$message)
            mash(mash.random, Ulist=U.c.r, outputlevel=1, gridmult=gridmult)
        }
    )
    message('  time: ', format(Sys.time()-now))

    message('  applying to all data (', nrow(Bhat), ' genes)...')
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
    return(list(mash=m))
}

# ═══════════════════════════════════════════════════════════════════════════
# PER_CT MODE (original)
# ═══════════════════════════════════════════════════════════════════════════

if (opt$mode == 'per_ct') {
    cell.type = opt$cell_type
    message('Mode: per_ct | Cell type: ', cell.type)

    in.file       = file.path(opt$checkpoints, paste0(cell.type, '_dglm_results.rds'))
    dglm.obj      = readRDS(in.file)
    dglm.results.orig = dglm.obj$array
    human.symbols.orig = dglm.obj$human_symbols

    terms = detect_terms(dimnames(dglm.results.orig)[[2]])
    message('Detected term(s): ', paste(if (length(terms)==1 && terms=='') '(single)' else terms, collapse=', '))
    if (!is.null(opt$only_term)) {
        if (!opt$only_term %in% terms) stop('--only_term ', opt$only_term, ' not among detected terms: ', paste(terms, collapse=', '))
        terms = opt$only_term
        message('--only_term: restricting to just "', terms, '"')
    }

    for (term in terms) {
        cols = cols_for_term(term)
        message('--- term: ', if (term=='') '(single)' else term, ' ---')

        dglm.results  = dglm.results.orig
        human.symbols = human.symbols.orig

        Bhat = dglm.results[, cols$beta, ]
        Shat = apply_shat_mode(dglm.results[, cols$bvar, ])
        Bhat[is.na(Bhat)] = 0
        Shat[is.na(Shat)] = 1000

        regions.with.data = colnames(Bhat)[apply(Bhat, 2, function(x) sum(x != 0) > 0)]
        message('Regions with data: ', paste(regions.with.data, collapse=', '))
        Bhat = Bhat[, regions.with.data, drop=FALSE]
        Shat = Shat[, regions.with.data, drop=FALSE]

        extreme = apply(Bhat, 1, function(x) any(abs(x) > 100, na.rm=TRUE))
        if (sum(extreme) > 0) {
            message('Filtering ', sum(extreme), ' genes with extreme beta')
            Bhat          = Bhat[!extreme,, drop=FALSE]
            Shat          = Shat[!extreme,, drop=FALSE]
            dglm.results  = dglm.results[!extreme,,, drop=FALSE]
            human.symbols = human.symbols[!extreme]
        }

        if (!is.null(autosome.x.genes)) {
            keep.mask = rownames(Bhat) %in% autosome.x.genes
            message('Autosome+X filter: dropping ', sum(!keep.mask), ' of ', length(keep.mask), ' genes (Y/MT/scaffold)')
            Bhat          = Bhat[keep.mask,, drop=FALSE]
            Shat          = Shat[keep.mask,, drop=FALSE]
            dglm.results  = dglm.results[keep.mask,,, drop=FALSE]
            human.symbols = human.symbols[keep.mask]
        }

        out.file = file.path(opt$checkpoints,
            paste0(cell.type, '_dglm_mashr_results', cols$suffix, '_strong',
                   strong.subset.qval.cutoff, '_lfsr', fsr.cutoff, '.rds'))

        if (length(regions.with.data) < 2) {
            message('Only one region — skipping mashr')
            saveRDS(list(mash=NULL, dglm_results=dglm.results,
                         human_symbols=human.symbols, regions=regions.with.data),
                    file=out.file)
            message('Saved: ', out.file)
            next
        }

        dglm.qval = dglm.results[, cols$qval, regions.with.data, drop=FALSE]
        result    = run_mashr(Bhat, Shat, dglm.qval, paste0(cell.type,
                    if (term != '') paste0(' (', term, ')') else ''), fast=opt$fast,
                    random_subset_n=opt$random_subset_n, gridmult=opt$gridmult)

        saveRDS(list(mash=result$mash, dglm_results=dglm.results,
                     human_symbols=human.symbols, regions=regions.with.data),
                file=out.file)
        message('Saved: ', out.file)
    }
}

# ═══════════════════════════════════════════════════════════════════════════
# COMBINED MODE — all cell types x regions in one mashr (per term)
# ═══════════════════════════════════════════════════════════════════════════

if (opt$mode == 'combined') {

    rds.files    = list.files(opt$checkpoints, pattern='_dglm_results\\.rds$')
    ct.discovered = sub('_dglm_results\\.rds$', '', rds.files)
    message('Discovered ', length(ct.discovered), ' unit(s) with DGLM output in ',
            opt$checkpoints, ': ', paste(ct.discovered, collapse=', '))
    if (length(ct.discovered) == 0) stop('No *_dglm_results.rds files found in ', opt$checkpoints)

    first.array = readRDS(file.path(opt$checkpoints, rds.files[1]))$array
    terms = detect_terms(dimnames(first.array)[[2]])
    message('Mode: combined — pooling all units x regions')
    message('Detected term(s): ', paste(if (length(terms)==1 && terms=='') '(single)' else terms, collapse=', '))
    if (!is.null(opt$only_term)) {
        if (!opt$only_term %in% terms) stop('--only_term ', opt$only_term, ' not among detected terms: ', paste(terms, collapse=', '))
        terms = opt$only_term
        message('--only_term: restricting to just "', terms, '"')
    }

    build_combined_matrices = function(term) {
        cols = cols_for_term(term)

        all.genes = character(0)
        for (ct in ct.discovered) {
            rds.file = file.path(opt$checkpoints, paste0(ct, '_dglm_results.rds'))
            if (!file.exists(rds.file)) next
            obj = readRDS(rds.file)
            all.genes = union(all.genes, dimnames(obj$array)[[1]])
        }
        message('[', term, '] Total unique genes (pre-filter): ', length(all.genes))

        if (!is.null(autosome.x.genes)) {
            n.before = length(all.genes)
            all.genes = intersect(all.genes, autosome.x.genes)
            message('[', term, '] Autosome+X filter: dropping ', n.before - length(all.genes),
                    ' of ', n.before, ' genes (Y/MT/scaffold)')
        }
        message('[', term, '] Total unique genes (post-filter): ', length(all.genes))

        conditions = character(0)
        for (ct in ct.discovered) {
            rds.file = file.path(opt$checkpoints, paste0(ct, '_dglm_results.rds'))
            if (!file.exists(rds.file)) next
            obj     = readRDS(rds.file)
            regions = dimnames(obj$array)[[3]]
            for (r in regions) {
                b = obj$array[, cols$beta, r]
                if (sum(!is.na(b) & b != 0) > 0) conditions = c(conditions, paste(ct, r, sep='|'))
            }
        }
        message('[', term, '] Total conditions (unit x region): ', length(conditions))

        Bhat      = matrix(0,    nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))
        Shat      = matrix(1000, nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))
        dglm.qval = matrix(1,    nrow=length(all.genes), ncol=length(conditions), dimnames=list(all.genes, conditions))

        human.symbols.list = list()
        dglm.results.list  = list()

        for (ct in ct.discovered) {
            rds.file = file.path(opt$checkpoints, paste0(ct, '_dglm_results.rds'))
            if (!file.exists(rds.file)) next
            obj          = readRDS(rds.file)
            dglm.results = obj$array
            human.sym    = obj$human_symbols

            beta.all = dglm.results[, cols$beta,, drop=FALSE]
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
                beta = dglm.results[, cols$beta, r]
                bvar = dglm.results[, cols$bvar, r]
                qval = dglm.results[, cols$qval, r]

                idx = genes.ct[genes.ct %in% all.genes]
                no.info = is.na(bvar[idx]) | bvar[idx] <= 0
                Bhat[idx, cond]      = ifelse(is.na(beta[idx]) | no.info, 0, beta[idx])
                Shat[idx, cond]      = ifelse(no.info, 1000, apply_shat_mode(bvar[idx]))
                dglm.qval[idx, cond] = ifelse(is.na(qval[idx]) | no.info, 1, qval[idx])
            }
        }

        human.symbols = setNames(all.genes, all.genes)
        for (ct in names(human.symbols.list)) {
            sym = human.symbols.list[[ct]]
            for (g in names(sym)) {
                if (g %in% names(human.symbols) && human.symbols[g] == g &&
                    !is.na(sym[g]) && sym[g] != g && sym[g] != '') {
                    human.symbols[g] = sym[g]
                }
            }
        }

        list(Bhat=Bhat, Shat=Shat, dglm.qval=dglm.qval, conditions=conditions,
             human.symbols=human.symbols, dglm.results.list=dglm.results.list)
    }

    for (term in terms) {
        cols = cols_for_term(term)
        message('══════════════════════════════════════════')
        message('Term: ', if (term=='') '(single)' else term)
        message('══════════════════════════════════════════')

        mats   = build_combined_matrices(term)
        result = run_mashr(mats$Bhat, mats$Shat, mats$dglm.qval,
                            paste0('combined', if (term != '') paste0(' (', term, ')') else ''), fast=opt$fast,
                            random_subset_n=opt$random_subset_n, gridmult=opt$gridmult)

        out.file = file.path(opt$checkpoints,
            paste0('combined_dglm_mashr_results', cols$suffix, '_strong',
                   strong.subset.qval.cutoff, '_lfsr', fsr.cutoff, '.rds'))
        saveRDS(list(
            mash          = result$mash,
            Bhat          = mats$Bhat,
            Shat          = mats$Shat,
            dglm_qval     = mats$dglm.qval,
            human_symbols = mats$human.symbols,
            conditions    = mats$conditions,
            dglm_results  = mats$dglm.results.list,
            term          = term
        ), file=out.file)
        message('Saved: ', out.file)
    }
}

message('done.')
