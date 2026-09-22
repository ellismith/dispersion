#!/bin/bash
# audit_pipeline.sh -- for each cell type, checks data size at every stage:
# raw pseudobulk (files, genes, animals) -> filtered (genes kept, min/max/mean
# across that cell type's subcluster x region files) -> master TSV per model
# (total rows, unique genes, unique subclusters, unique regions).

BASE=/scratch/easmit31/dispersion/dglm/disp_age__raw_counts

printf "%-24s %8s %8s %10s %10s %10s | %8s %8s %8s %8s | %8s %8s %8s %8s\n" \
  "cell_type" "n_raw_f" "raw_gen" "filt_min" "filt_max" "filt_mean" \
  "age_n" "age_sub" "age_reg" "age_gen" "asx_n" "asx_sub" "asx_reg" "asx_gen"

for d in "$BASE"/*/; do
  ct=$(basename "$d")

  n_raw_f=$(ls "$d"/*_pseudobulk.csv 2>/dev/null | wc -l)
  raw_f1=$(ls "$d"/*_pseudobulk.csv 2>/dev/null | head -1)
  raw_genes=$( [ -n "$raw_f1" ] && wc -l < "$raw_f1" | awk '{print $1-1}' || echo "NA" )

  filt_files=$(ls "$d"/filtered/*_filtered_cutoff0.5.csv 2>/dev/null)
  if [ -n "$filt_files" ]; then
    filt_counts=$(for f in $filt_files; do wc -l < "$f" | awk '{print $1-1}'; done)
    filt_min=$(echo "$filt_counts" | sort -n | head -1)
    filt_max=$(echo "$filt_counts" | sort -n | tail -1)
    filt_mean=$(echo "$filt_counts" | awk '{s+=$1; n++} END{if(n>0) printf "%.0f", s/n; else print "NA"}')
  else
    filt_min="NA"; filt_max="NA"; filt_mean="NA"
  fi

  age_tsv="$d/dglm_checkpoints_cutoff0.5/master_dglm_combined.tsv"
  if [ -f "$age_tsv" ]; then
    age_n=$(($(wc -l < "$age_tsv") - 1))
    age_sub=$(tail -n +2 "$age_tsv" | cut -f3 | sort -u | wc -l)
    age_reg=$(tail -n +2 "$age_tsv" | cut -f4 | sort -u | wc -l)
    age_gen=$(tail -n +2 "$age_tsv" | cut -f1 | sort -u | wc -l)
  else
    age_n="NA"; age_sub="NA"; age_reg="NA"; age_gen="NA"
  fi

  asx_tsv="$d/dglm_checkpoints_cutoff0.5_agesex/master_dglm_age_disp_combined.tsv"
  if [ -f "$asx_tsv" ]; then
    asx_n=$(($(wc -l < "$asx_tsv") - 1))
    asx_sub=$(tail -n +2 "$asx_tsv" | cut -f3 | sort -u | wc -l)
    asx_reg=$(tail -n +2 "$asx_tsv" | cut -f4 | sort -u | wc -l)
    asx_gen=$(tail -n +2 "$asx_tsv" | cut -f1 | sort -u | wc -l)
  else
    asx_n="NA"; asx_sub="NA"; asx_reg="NA"; asx_gen="NA"
  fi

  printf "%-24s %8s %8s %10s %10s %10s | %8s %8s %8s %8s | %8s %8s %8s %8s\n" \
    "$ct" "$n_raw_f" "$raw_genes" "$filt_min" "$filt_max" "$filt_mean" \
    "$age_n" "$age_sub" "$age_reg" "$age_gen" "$asx_n" "$asx_sub" "$asx_reg" "$asx_gen"
done
