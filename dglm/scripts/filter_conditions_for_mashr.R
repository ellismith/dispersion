#!/usr/bin/env Rscript
#
# filter_conditions_for_mashr.R
#
# For each subcluster's *_dglm_results.rds (array: genes x stats x regions),
# drops any region that fails the >=min_cells AND >=min_animals condition
# filter (per assess_condition_sparsity.py's output), and saves a trimmed
# copy into a NEW checkpoints directory -- original files untouched.
# dglm_mashr.R --mode combined then runs unmodified on the new directory.
#
# Usage:
#   Rscript filter_conditions_for_mashr.R --cell_type GABAergic_neurons --dispersion age

suppressMessages(library(optparse))

option_list = list(
    make_option('--cell_type',    type='character', default=NULL),
    make_option('--dispersion',   type='character', default='age'),
    make_option('--sparsity_csv', type='character', default=NULL),
    make_option('--src_base',     type='character',
                default='/scratch/easmit31/dispersion/dglm/disp_age__raw_counts'),
    make_option('--out_base',     type='character',
                default='/scratch/easmit31/dispersion/dglm/disp_age__min300cells_min30animals')
)
opt = parse_args(OptionParser(option_list=option_list))
if (is.null(opt$cell_type)) stop('--cell_type is required')

ckpt_suffix = if (opt$dispersion == 'age_sex') '_cutoff0.5_agesex' else '_cutoff0.5'
src_dir = file.path(opt$src_base, opt$cell_type, paste0('dglm_checkpoints', ckpt_suffix))
out_dir = file.path(opt$out_base, opt$cell_type, paste0('dglm_checkpoints', ckpt_suffix))
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

sparsity_csv = if (!is.null(opt$sparsity_csv)) opt$sparsity_csv else
    file.path(opt$out_base, paste0(opt$cell_type, '_condition_sparsity.csv'))
if (!file.exists(sparsity_csv)) stop('Sparsity CSV not found: ', sparsity_csv)

sparsity = read.csv(sparsity_csv, stringsAsFactors=FALSE)
sparsity$pass_both = as.logical(sparsity$pass_both)
message('Loaded sparsity table: ', nrow(sparsity), ' rows from ', sparsity_csv)

rds_files = list.files(src_dir, pattern='_dglm_results\\.rds$', full.names=TRUE)
message('Found ', length(rds_files), ' subcluster RDS file(s) in ', src_dir)
if (length(rds_files) == 0) stop('No RDS files found in ', src_dir)

total_before = 0
total_after  = 0

for (f in rds_files) {
    subcluster = sub('_dglm_results\\.rds$', '', basename(f))
    if (grepl('^opc-olig_', subcluster)) {
        subcluster = paste0(opt$cell_type, '_', sub('^opc-olig_', '', subcluster))
    }
    obj = readRDS(f)
    arr = obj$array
    all_slots = dimnames(arr)[[3]]
    stat_names = dimnames(arr)[[2]]
    beta_col = if ('beta' %in% stat_names) 'beta' else if ('beta_age' %in% stat_names) 'beta_age' else
        stop('Could not find a beta or beta_age column in array stats: ', paste(stat_names, collapse=', '))
    has_real_data = sapply(all_slots, function(r) {
        b = arr[, beta_col, r]
        any(!is.na(b) & b != 0)
    })
    regions_present = all_slots[has_real_data]
    total_before = total_before + length(regions_present)

    pass_regions = sparsity$region[sparsity$subcluster == subcluster & sparsity$pass_both]
    keep = intersect(regions_present, pass_regions)
    dropped = setdiff(regions_present, keep)
    total_after = total_after + length(keep)

    message(subcluster, ': ', length(regions_present), ' regions -> ', length(keep),
            ' kept', if (length(dropped) > 0) paste0(' (dropped: ', paste(dropped, collapse=', '), ')') else '')

    if (length(keep) == 0) {
        message('  WARNING: ', subcluster, ' has ZERO passing regions -- skipping, not saved')
        next
    }

    obj$array = arr[, , keep, drop=FALSE]
    out_file = file.path(out_dir, basename(f))
    saveRDS(obj, file=out_file)
}

message('══════════════════════════════════════════')
message('Total conditions before: ', total_before, ' | after: ', total_after, ' | removed: ', total_before - total_after)
message('Saved to: ', out_dir)
message('done.')
