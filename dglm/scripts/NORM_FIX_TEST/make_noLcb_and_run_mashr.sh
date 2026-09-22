#!/bin/bash
set -e
CT=$1
SRC=/scratch/easmit31/dispersion/dglm/disp_age__FINAL_autosomeX_300c/${CT}/dglm_checkpoints_cutoff0.5
DEST=/scratch/easmit31/dispersion/dglm/disp_age__FINAL_autosomeX_300c_noLcb/${CT}/dglm_checkpoints_cutoff0.5
mkdir -p $DEST

Rscript -e "
files = list.files('$SRC', pattern='_dglm_results.rds\$', full.names=TRUE)
for (f in files) {
  r = readRDS(f)
  if ('lCb' %in% dimnames(r\$array)[[3]]) {
    keep_regions = setdiff(dimnames(r\$array)[[3]], 'lCb')
    r\$array = r\$array[, , keep_regions, drop=FALSE]
  }
  saveRDS(r, file.path('$DEST', basename(f)))
}
cat('Rebuilt', length(files), 'files for $CT without lCb\n')
"

/home/easmit31/.conda/envs/mashr_env/bin/Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_mashr.R \
  --checkpoints $DEST --approach qval_single --shat_mode raw

cp $DEST/combined_dglm_mashr_results_qval_single_raw_fill1000.rds \
   $DEST/combined_dglm_mashr_results_strong0.05_lfsr0.2.rds

Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_fdr_combined.R --checkpoints $DEST

echo "DONE: ${CT} (noLcb)"
