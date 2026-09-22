#!/usr/bin/env Rscript
# For each subcluster x region condition, check whether each age bin
# has a median n_cells >= min_median_cells. If any bin fails, that
# condition is flagged and its array values set to NA (excluded from mashr).
# Applied as a post-DGLM, pre-mashr cleaning step.
#
# Usage: Rscript filter_age_bin_cells.R --checkpoints <dir> --metadata_dir <dir> \
#          --out_checkpoints <dir> --min_median_cells 20

suppressMessages(library(optparse))

option_list = list(
  make_option('--checkpoints', type='character'),
  make_option('--metadata_dir', type='character'),
  make_option('--out_checkpoints', type='character'),
  make_option('--min_median_cells', type='integer', default=20)
)
opt = parse_args(OptionParser(option_list=option_list))
dir.create(opt$out_checkpoints, showWarnings=FALSE, recursive=TRUE)

age_breaks = c(0, 7, 12, 25)
age_labels = c("young","mid","old")

files = list.files(opt$checkpoints, pattern='_dglm_results.rds$', full.names=TRUE)
total_flagged = 0; total_conditions = 0

for (f in files) {
  r = readRDS(f)
  arr = r$array
  subcl = r$run_info$cell_type

  meta_file = file.path(opt$metadata_dir, paste0(subcl, '_metadata.csv'))
  if (!file.exists(meta_file)) {
    saveRDS(r, file.path(opt$out_checkpoints, basename(f)))
    next
  }
  meta = read.csv(meta_file, stringsAsFactors=FALSE)
  meta$age_bin = cut(meta$age, breaks=age_breaks, labels=age_labels)

  for (reg in dimnames(arr)[[3]]) {
    total_conditions = total_conditions + 1
    meta_reg = meta[meta$region == reg, ]
    if (nrow(meta_reg) == 0) next

    bin_medians = tapply(meta_reg$n_cells, meta_reg$age_bin, median, na.rm=TRUE)
    bins_with_data = bin_medians[!is.na(bin_medians)]
    failing_bins = names(bins_with_data)[bins_with_data < opt$min_median_cells]

    if (length(failing_bins) > 0) {
      cat(sprintf('  %s | %s: flagging (bins below %d cells median: %s)\n',
                  subcl, reg, opt$min_median_cells, paste(failing_bins, collapse=', ')))
      arr[, , reg] = NA
      total_flagged = total_flagged + 1
    }
  }
  r$array = arr
  saveRDS(r, file.path(opt$out_checkpoints, basename(f)))
}

cat(sprintf('\nTotal conditions flagged: %d / %d (%.1f%%)\n',
            total_flagged, total_conditions, 100*total_flagged/total_conditions))
