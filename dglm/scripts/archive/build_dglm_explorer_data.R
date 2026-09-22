#!/usr/bin/env Rscript
#
# build_dglm_explorer_data.R
#
# Builds the DATA JSON blob consumed by dglm_explorer.html. Now produces
# TWO model entries: V2 (dispersion~age) and V3 (dispersion~age+sex, age
# term), both at louvain resolution, sqrt shat mode, per-cell-type mashr.
# Merges every cell type's combined mashr result for DISPLAY into one
# unified unit list per model -- each subcluster's numbers still come from
# its own cell type's independent mash fit.
#
# Usage:
#   Rscript build_dglm_explorer_data.R --out /scratch/easmit31/dispersion/dglm/dglm_explorer_data.json

suppressMessages({
  library(optparse)
  library(mashr)
  library(jsonlite)
})

option_list = list(
    make_option('--base_dir', type='character',
                default='/scratch/easmit31/dispersion/dglm/disp_age__raw_counts'),
    make_option('--cutoff', type='character', default='0.5'),
    make_option('--out', type='character', required=TRUE)
)
opt = parse_args(OptionParser(option_list=option_list))

THRESH_OPTIONS = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5)

build_model = function(ckpt_suffix, mashr_pattern, model_label) {
    cell_type_dirs = list.dirs(opt$base_dir, recursive=FALSE)
    mashr_files = c()
    for (d in cell_type_dirs) {
        ckpt = file.path(d, paste0('dglm_checkpoints_cutoff', opt$cutoff, ckpt_suffix))
        f = list.files(ckpt, pattern=mashr_pattern, full.names=TRUE)
        if (length(f) > 0) mashr_files = c(mashr_files, f[1])
    }
    message('[', model_label, '] Found ', length(mashr_files), ' cell type(s): ',
            paste(basename(dirname(dirname(mashr_files))), collapse=', '))
    if (length(mashr_files) == 0) {
        return(list(label=model_label, resolution='louvain', available=FALSE,
                    units=list(), regions=list(), bar=list(), volcano=list(), heatmap=list()))
    }

    all_mash_beta = list(); all_mash_lfsr = list()
    all_raw_beta  = list(); all_raw_qval  = list()
    human_symbols_global = c()
    all_units = c(); all_regions = c()

    for (f in mashr_files) {
        obj = readRDS(f)
        mash_beta = mashr::get_pm(obj$mash)
        mash_lfsr = mashr::get_lfsr(obj$mash)
        raw_beta  = obj$Bhat
        raw_qval  = obj$dglm_qval
        conds     = obj$conditions
        human_symbols_global[names(obj$human_symbols)] = obj$human_symbols
        for (cond in conds) {
            parts = strsplit(cond, '\\|')[[1]]
            unit = parts[1]; region = parts[2]
            all_units = union(all_units, unit)
            all_regions = union(all_regions, region)
            all_mash_beta[[cond]] = mash_beta[, cond]
            all_mash_lfsr[[cond]] = mash_lfsr[, cond]
            all_raw_beta[[cond]]  = raw_beta[, cond]
            all_raw_qval[[cond]]  = raw_qval[, cond]
        }
    }
    all_units = sort(all_units); all_regions = sort(all_regions)
    message('[', model_label, '] units: ', length(all_units), ' | regions: ', paste(all_regions, collapse=', '))

    sym_or_gene = function(g) {
        s = human_symbols_global[g]
        if (is.na(s) || s == '') g else s
    }

    build_bar = function() {
        mash_by_thresh = list()
        for (thresh in THRESH_OPTIONS) {
            by_ct_region = list()
            for (unit in all_units) {
                by_ct_region[[unit]] = list()
                for (region in all_regions) {
                    key = paste(unit, region, sep='|')
                    if (is.null(all_mash_beta[[key]])) next
                    beta = all_mash_beta[[key]]; lfsr = all_mash_lfsr[[key]]
                    sig = !is.na(lfsr) & lfsr < thresh
                    by_ct_region[[unit]][[region]] = list(
                        increase = sum(sig & beta > 0, na.rm=TRUE),
                        decrease = sum(sig & beta <= 0, na.rm=TRUE))
                }
            }
            by_region = list()
            for (region in all_regions) {
                inc = sum(sapply(all_units, function(u) { v = by_ct_region[[u]][[region]]$increase; if (is.null(v)) 0 else v }))
                dec = sum(sapply(all_units, function(u) { v = by_ct_region[[u]][[region]]$decrease; if (is.null(v)) 0 else v }))
                by_region[[region]] = list(increase=inc, decrease=dec)
            }
            mash_by_thresh[[as.character(thresh)]] = list(by_region=by_region, by_celltype_region=by_ct_region)
        }
        raw_by_ct_region = list()
        for (unit in all_units) {
            raw_by_ct_region[[unit]] = list()
            for (region in all_regions) {
                key = paste(unit, region, sep='|')
                if (is.null(all_raw_beta[[key]])) next
                beta = all_raw_beta[[key]]; qval = all_raw_qval[[key]]
                sig = !is.na(qval) & qval < 0.05
                raw_by_ct_region[[unit]][[region]] = list(
                    increase = sum(sig & beta > 0, na.rm=TRUE),
                    decrease = sum(sig & beta <= 0, na.rm=TRUE))
            }
        }
        raw_by_region = list()
        for (region in all_regions) {
            inc = sum(sapply(all_units, function(u) { v = raw_by_ct_region[[u]][[region]]$increase; if (is.null(v)) 0 else v }))
            dec = sum(sapply(all_units, function(u) { v = raw_by_ct_region[[u]][[region]]$decrease; if (is.null(v)) 0 else v }))
            raw_by_region[[region]] = list(increase=inc, decrease=dec)
        }
        list(mash=mash_by_thresh, raw=list(by_region=raw_by_region, by_celltype_region=raw_by_ct_region))
    }

    build_volcano = function(max_background=2000) {
        volcano = list()
        for (unit in all_units) {
            volcano[[unit]] = list()
            for (region in all_regions) {
                key = paste(unit, region, sep='|')
                if (is.null(all_mash_beta[[key]])) next
                beta = all_mash_beta[[key]]; lfsr = all_mash_lfsr[[key]]
                genes = names(beta)
                sig_idx = which(!is.na(lfsr) & lfsr < 0.5)
                bg_idx = setdiff(seq_along(genes), sig_idx)
                if (length(bg_idx) > max_background) bg_idx = sample(bg_idx, max_background)
                keep_idx = c(sig_idx, bg_idx)
                pts = lapply(keep_idx, function(i) list(g = sym_or_gene(genes[i]), b = unname(beta[i]), l = unname(lfsr[i])))
                volcano[[unit]][[region]] = pts
            }
        }
        volcano
    }

    build_heatmap = function() {
        mk_matrix = function(fill_fn) {
            m = matrix(NA_real_, nrow=length(all_units), ncol=length(all_regions))
            for (i in seq_along(all_units)) for (j in seq_along(all_regions)) {
                key = paste(all_units[i], all_regions[j], sep='|')
                m[i, j] = fill_fn(key)
            }
            m
        }
        genes_tested = mk_matrix(function(key) if (is.null(all_raw_beta[[key]])) NA_real_ else length(all_raw_beta[[key]]))
        raw_beta_all = mk_matrix(function(key) { if (is.null(all_raw_beta[[key]])) return(NA_real_); mean(abs(all_raw_beta[[key]]), na.rm=TRUE) })
        raw_beta_sig = mk_matrix(function(key) {
            if (is.null(all_raw_beta[[key]])) return(NA_real_)
            b = all_raw_beta[[key]]; q = all_raw_qval[[key]]; sig = !is.na(q) & q < 0.05
            if (!any(sig)) return(0); mean(abs(b[sig]), na.rm=TRUE)
        })
        by_threshold = list()
        for (thresh in THRESH_OPTIONS) {
            pct_sig = mk_matrix(function(key) { if (is.null(all_mash_lfsr[[key]])) return(NA_real_); l = all_mash_lfsr[[key]]; 100 * sum(!is.na(l) & l < thresh) / length(l) })
            n_increase = mk_matrix(function(key) { if (is.null(all_mash_lfsr[[key]])) return(NA_real_); b = all_mash_beta[[key]]; l = all_mash_lfsr[[key]]; sum(!is.na(l) & l < thresh & b > 0, na.rm=TRUE) })
            n_decrease = mk_matrix(function(key) { if (is.null(all_mash_lfsr[[key]])) return(NA_real_); b = all_mash_beta[[key]]; l = all_mash_lfsr[[key]]; sum(!is.na(l) & l < thresh & b <= 0, na.rm=TRUE) })
            n_net = n_increase - n_decrease
            mash_beta_sig = mk_matrix(function(key) { if (is.null(all_mash_lfsr[[key]])) return(NA_real_); b = all_mash_beta[[key]]; l = all_mash_lfsr[[key]]; sig = !is.na(l) & l < thresh; if (!any(sig)) return(0); mean(b[sig], na.rm=TRUE) })
            mash_beta_mag_sig = mk_matrix(function(key) { if (is.null(all_mash_lfsr[[key]])) return(NA_real_); b = all_mash_beta[[key]]; l = all_mash_lfsr[[key]]; sig = !is.na(l) & l < thresh; if (!any(sig)) return(0); mean(abs(b[sig]), na.rm=TRUE) })
            by_threshold[[as.character(thresh)]] = list(
                pct_sig=list(z=pct_sig, kind='sequential', label='% significant'),
                n_increase=list(z=n_increase, kind='sequential', label='# sig genes, increase'),
                n_decrease=list(z=n_decrease, kind='sequential', label='# sig genes, decrease'),
                n_net=list(z=n_net, kind='diverging', label='Net (increase - decrease)'),
                mash_beta_sig=list(z=mash_beta_sig, kind='diverging', label='Mean mash beta (sig genes)'),
                mash_beta_mag_sig=list(z=mash_beta_mag_sig, kind='sequential', label='Mean |mash beta| (sig genes)'))
        }
        list(units=all_units, regions=all_regions,
             static=list(genes_tested=list(z=genes_tested, kind='sequential', label='Genes tested'),
                         raw_beta_all=list(z=raw_beta_all, kind='sequential', label='Mean |raw beta| (all genes)'),
                         raw_beta_sig=list(z=raw_beta_sig, kind='sequential', label='Mean |raw beta| (raw q<0.05)')),
             by_threshold=by_threshold)
    }

    list(label=model_label, resolution='louvain', available=TRUE,
         units=all_units, regions=all_regions,
         bar=build_bar(), volcano=build_volcano(), heatmap=build_heatmap())
}

message('Building V2 (age)...')
v2 = build_model('', '^combined_dglm_mashr_results_strong.*_lfsr.*\\.rds$',
                  'V2: sqrt Shat, per-cell-type mashr, louvain resolution (dispersion~age)')

message('Building V3 (age_sex, age term)...')
v3 = build_model('_agesex', '^combined_dglm_mashr_results_age_strong.*_lfsr.*\\.rds$',
                  'V3: sqrt Shat, per-cell-type mashr, louvain resolution (dispersion~age+sex, age term)')

result = list(v2_sqrt_louvain = v2, v3_sqrt_louvain_agesex = v3)

write_json(result, opt$out, auto_unbox=TRUE, na='null', digits=6)
message('Saved: ', opt$out)
