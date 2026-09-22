#!/bin/bash
set -e
CT=$1
SRC=/scratch/easmit31/dispersion/dglm/disp_age__FINAL_autosomeX_300c/${CT}/dglm_checkpoints_cutoff0.5
DEST=/scratch/easmit31/dispersion/dglm/disp_age__FINAL_autosomeX_300c_QCfiltered/${CT}/dglm_checkpoints_cutoff0.5
mkdir -p $DEST

Rscript -e "
files = list.files('$SRC', pattern='_dglm_results.rds\$', full.names=TRUE)
n_flagged_total = 0
for (f in files) {
  r = readRDS(f)
  arr = r\$array
  regions = dimnames(arr)[[3]]
  beta_mat = matrix(arr[, 'beta', ], nrow=dim(arr)[1], ncol=length(regions), dimnames=list(dimnames(arr)[[1]], regions))
  bvar_mat = matrix(arr[, 'bvar', ], nrow=dim(arr)[1], ncol=length(regions), dimnames=list(dimnames(arr)[[1]], regions))
  bad = which(abs(beta_mat) > 10 | bvar_mat > 10, arr.ind=TRUE)
  if (!is.null(bad) && nrow(bad) > 0) {
    for (i in seq_len(nrow(bad))) {
      g = bad[i,1]; reg_idx = bad[i,2]
      arr[g,'beta',reg_idx]=NA; arr[g,'bvar',reg_idx]=NA; arr[g,'pval',reg_idx]=NA; arr[g,'qval',reg_idx]=NA
    }
    n_flagged_total = n_flagged_total + nrow(bad)
    cat('  ', basename(f), ': flagged', nrow(bad), 'implausible fit(s)\n')
  }
  r\$array = arr
  saveRDS(r, file.path('$DEST', basename(f)))
}
cat('$CT total flagged:', n_flagged_total, '\n')
"

/home/easmit31/.conda/envs/mashr_env/bin/Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_mashr.R \
  --checkpoints $DEST --approach qval_single --shat_mode raw
cp $DEST/combined_dglm_mashr_results_qval_single_raw_fill1000.rds \
   $DEST/combined_dglm_mashr_results_strong0.05_lfsr0.2.rds
Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_fdr_combined.R --checkpoints $DEST
echo "DONE: ${CT}"
