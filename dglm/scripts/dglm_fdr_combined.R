#!/usr/bin/env Rscript
# dglm_fdr_combined.R
# Extracts mashr results from combined run, applies global BH-FDR, saves
# master TSV(s).
#
# Auto-detects single-term vs two-term output (age-only dispersion vs
# age+sex dispersion) from which mashr result files actually exist in
# --checkpoints, and loops accordingly -- same pattern as dglm_mashr.R's
# own single/two-term auto-detection, so this can't drift out of sync with
# what that script actually produces.
#   single-term: combined_dglm_mashr_results_strong..._lfsr....rds
#     -> master_dglm_combined.tsv (unchanged filename/behavior)
#   two-term:    combined_dglm_mashr_results_age_strong..._lfsr....rds
#                combined_dglm_mashr_results_sex_strong..._lfsr....rds
#     -> master_dglm_age_disp_combined.tsv, master_dglm_sex_disp_combined.tsv
#
# Previously this script only handled the single-term case (hardcoded
# 'beta'/'bvar'/'pval' column names when reading the raw dglm array, and a
# single fixed mashr-result filename to load) -- would have failed on
# two-term (disp_age_sex__mean_full) output twice over: first on the
# mismatched mashr filename (dglm_mashr.R actually saves separate
# _age/_sex-suffixed files for two-term runs), then on the missing plain
# 'beta' column in the two-term dglm array even if that were fixed. A
# dglm_fdr_age_sex_disp_combined.R was referenced in prior notes as the
# intended fix but never actually existed on disk (confirmed via ls) --
# fixed here instead by extending this script, rather than adding a fourth
# FDR script that duplicates logic and can drift out of sync again --
# matches how dglm_mashr.R itself was unified into one script rather than
# kept as per-case duplicates.
source('/scratch/easmit31/dispersion/dglm/scripts/_include_options.R')
library(optparse)
library(mashr)

option_list = list(
    make_option('--checkpoints', type='character', default='/scratch/easmit31/dispersion/dglm/checkpoints')
)
opt = parse_args(OptionParser(option_list=option_list))

# ── same term-column-naming convention as dglm_mashr.R's cols_for_term() ──
cols_for_term = function(term) {
    if (term == '') {
        list(beta='beta', bvar='bvar', pval='pval', suffix='')
    } else {
        list(beta=paste0('beta_',term), bvar=paste0('bvar_',term),
             pval=paste0('pval_',term), suffix=paste0('_', term))
    }
}

mashr_file_for = function(suffix) {
    file.path(opt$checkpoints,
        paste0('combined_dglm_mashr_results', suffix, '_strong',
               strong.subset.qval.cutoff, '_lfsr', fsr.cutoff, '.rds'))
}

# ── detect single-term vs two-term from which mashr output files exist ────
single.file = mashr_file_for('')
age.file    = mashr_file_for('_age')
sex.file    = mashr_file_for('_sex')

if (file.exists(age.file) && file.exists(sex.file)) {
    terms = c('age','sex')
    message('Detected two-term (age+sex) combined mashr output')
} else if (file.exists(single.file)) {
    terms = ''
    message('Detected single-term combined mashr output')
} else {
    stop('No combined mashr results found in ', opt$checkpoints,
         ' matching expected filenames.\n  checked: ', single.file,
         '\n       or: ', age.file, ' + ', sex.file)
}

process_term = function(term) {
    cols       = cols_for_term(term)
    mashr.file = mashr_file_for(cols$suffix)
    message('══════════════════════════════════════')
    message('Term: ', if (term=='') '(single)' else term)
    message('  file: ', mashr.file)

    obj           = readRDS(mashr.file)
    m             = obj$mash
    human.symbols = obj$human_symbols
    conditions    = obj$conditions
    dglm.list     = obj$dglm_results
    if (is.null(m)) stop('mashr object is NULL for term: ', if (term=='') '(single)' else term)

    mash.beta = get_pm(m)
    mash.lfsr = get_lfsr(m)
    message('  lfsr range: ', paste(round(range(mash.lfsr, na.rm=TRUE), 4), collapse=' '))
    message('  Sig at lfsr<0.2: ', sum(mash.lfsr < 0.2, na.rm=TRUE))

    genes.in.mashr = rownames(mash.beta)
    all.rows = list()
    for (cond in conditions) {
        parts  = strsplit(cond, '\\|')[[1]]
        ct     = parts[1]
        region = parts[2]
        if (!ct %in% names(dglm.list)) next
        dglm.results = dglm.list[[ct]]
        if (!region %in% dimnames(dglm.results)[[3]]) next

        # restrict to genes actually present in the mash object -- any gene
        # here but not there was deliberately excluded upstream (autosome+X
        # filter applied to Bhat/Shat in combined mode, but NOT to the raw
        # dglm.results arrays also saved in the RDS), so it should be
        # excluded from the master TSV entirely, not kept with NA mash
        # columns. This also avoids a real crash: ifelse() evaluates BOTH
        # branches for every row regardless of the condition, so
        # mash.beta[genes, cond] was being evaluated even for genes not in
        # mash.beta's rownames -- matrix indexing by name (unlike vector
        # indexing) errors on any missing name instead of returning NA,
        # which is what threw "subscript out of bounds".
        genes.avail = dimnames(dglm.results)[[1]]
        genes       = intersect(genes.avail, genes.in.mashr)
        n.dropped   = length(genes.avail) - length(genes)
        if (n.dropped > 0) {
            message('    [', ct, '|', region, '] dropping ', n.dropped,
                    ' gene(s) in raw DGLM output but not in the mash object (filtered upstream)')
        }

        beta = dglm.results[genes, cols$beta, region]
        bvar = dglm.results[genes, cols$bvar, region]
        pval = dglm.results[genes, cols$pval, region]

        all.rows[[length(all.rows)+1]] = data.frame(
            ensembl_id  = genes,
            symbol      = human.symbols[genes],
            cell_type   = ct,
            region      = region,
            beta        = beta,
            bvar        = bvar,
            pvalue      = pval,
            mash_beta   = mash.beta[genes, cond],
            mash_lfsr   = mash.lfsr[genes, cond],
            stringsAsFactors = FALSE
        )
    }

    master        = do.call(rbind, all.rows)
    master        = master[!is.na(master$pvalue),]
    master$qvalue = p.adjust(master$pvalue, method='fdr')
    message('  Total tests: ', nrow(master))
    message('  Sig DGLM q<0.05: ',    sum(master$qvalue    < 0.05, na.rm=TRUE))
    message('  Sig mashr lfsr<0.2: ', sum(master$mash_lfsr < 0.2,  na.rm=TRUE))
    message('  Sig mashr lfsr<0.05: ',sum(master$mash_lfsr < 0.05, na.rm=TRUE))

    out.name = if (term == '') 'master_dglm_combined.tsv' else paste0('master_dglm_', term, '_disp_combined.tsv')
    out = file.path(opt$checkpoints, out.name)
    write.table(master, out, sep='\t', row.names=FALSE, quote=FALSE)
    message('  Saved: ', out)
}

for (term in terms) process_term(term)

message('done.')
