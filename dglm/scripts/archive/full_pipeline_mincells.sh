#!/bin/bash
set -e
CT=$1
MINCELLS=$2
BASE=/scratch/easmit31/dispersion/dglm/disp_age__mincells${MINCELLS}
RAW_DIR=${BASE}/${CT}
FILTERED_DIR=${RAW_DIR}/filtered
STAGE_BASE=${BASE}/${CT}_staged/dglm_run_cutoff0.5
CHECKPOINTS=${BASE}/${CT}_staged/dglm_checkpoints_cutoff0.5
SPARSITY_DIR=${BASE}/sparsity_diagnostics
FILTERED300_BASE=${BASE}_300c
QC_BASE=${BASE}_300c_QCfiltered

mkdir -p $RAW_DIR $STAGE_BASE $CHECKPOINTS $SPARSITY_DIR

echo "=== 1. Pseudobulk (min_cells_per_animal=${MINCELLS}) ==="
source /packages/apps/mamba/2.0.8/etc/profile.d/conda.sh
conda activate mixed_models
python3 /scratch/easmit31/dispersion/dglm/scripts/pseudobulk.py \
  --cell_type ${CT} --outdir ${RAW_DIR} --subcluster_col ct_louvain \
  --min_cells_per_animal ${MINCELLS}

echo "=== 2. CPM filter ==="
conda activate mashr_env
bash /scratch/easmit31/dispersion/dglm/scripts/filter_all_pseudobulk.sh ${RAW_DIR} ${FILTERED_DIR}

echo "=== 3. DGLM ==="
/scratch/easmit31/dispersion/dglm/scripts/run_dglm_for_celltype.sh \
  ${RAW_DIR} ${FILTERED_DIR} 0.5 ${STAGE_BASE} ${CHECKPOINTS} age

echo "=== 4. Condition filter (300c/30 animals) ==="
python3 /scratch/easmit31/dispersion/dglm/scripts/assess_condition_sparsity.py \
  --cell_type ${CT} --base_dir ${RAW_DIR} --out_dir ${SPARSITY_DIR} \
  --min_cells 300 --min_animals 30

Rscript /scratch/easmit31/dispersion/dglm/scripts/filter_conditions_for_mashr.R \
  --cell_type ${CT} --dispersion age \
  --sparsity_csv ${SPARSITY_DIR}/${CT}_condition_sparsity.csv \
  --src_base ${BASE}/${CT}_staged --out_base ${FILTERED300_BASE}

echo "=== 5. mashr ==="
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
FINAL_CHECKPOINTS=${FILTERED300_BASE}/${CT}/dglm_checkpoints_cutoff0.5
Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_mashr.R \
  --checkpoints ${FINAL_CHECKPOINTS} --approach qval_single --shat_mode raw

echo "=== 6. QC filter (|beta|>10 guardrail) ==="
QC_CHECKPOINTS=${QC_BASE}/${CT}/dglm_checkpoints_cutoff0.5
mkdir -p ${QC_CHECKPOINTS}
Rscript -e "
files = list.files('${FINAL_CHECKPOINTS}', pattern='_dglm_results.rds\$', full.names=TRUE)
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
  }
  r\$array = arr
  saveRDS(r, file.path('${QC_CHECKPOINTS}', basename(f)))
}
cat('Flagged:', n_flagged_total, '\n')
"

cp ${FINAL_CHECKPOINTS}/combined_dglm_mashr_results_qval_single_raw_fill1000.rds ${QC_CHECKPOINTS}/ 2>/dev/null || true
Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_mashr.R \
  --checkpoints ${QC_CHECKPOINTS} --approach qval_single --shat_mode raw
cp ${QC_CHECKPOINTS}/combined_dglm_mashr_results_qval_single_raw_fill1000.rds \
   ${QC_CHECKPOINTS}/combined_dglm_mashr_results_strong0.05_lfsr0.2.rds
Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_fdr_combined.R --checkpoints ${QC_CHECKPOINTS}

echo "DONE: ${CT} min_cells_per_animal=${MINCELLS}"
