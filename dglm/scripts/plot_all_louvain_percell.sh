#!/usr/bin/env bash
set -uo pipefail

RSCRIPT="/home/easmit31/.conda/envs/mashr_env/bin/Rscript"
PLOT_SCRIPT="/scratch/easmit31/dispersion/dglm/scripts/dglm_plot_fig4.R"
QTHRESH="0.2"

RAW_ROOT="/scratch/easmit31/dispersion/dglm/disp_age__mean_full_louvain_percell_raw"
SQRT_ROOT="/scratch/easmit31/dispersion/dglm/disp_age__mean_full_louvain_percell_sqrt"

RAW_FIG_ROOT="/scratch/easmit31/dispersion/dglm/figures/louvain_percell_raw"
SQRT_FIG_ROOT="/scratch/easmit31/dispersion/dglm/figures/louvain_percell_sqrt"

run_root() {
    local mode="$1"
    local input_root="$2"
    local output_root="$3"
    local checkpoint_dir
    local cell_type
    local master_tsv
    local figdir

    echo "============================================================"
    echo "Starting ${mode}: ${input_root}"
    echo "============================================================"

    for checkpoint_dir in "$input_root"/*; do
        [[ -d "$checkpoint_dir" ]] || continue

        cell_type="$(basename "$checkpoint_dir")"
        master_tsv="${checkpoint_dir}/master_dglm_combined.tsv"
        figdir="${output_root}/${cell_type}"

        if [[ ! -s "$master_tsv" ]]; then
            echo "[skip] ${mode}/${cell_type}: missing or empty master_dglm_combined.tsv"
            continue
        fi

        mkdir -p "$figdir"

        echo "[run] ${mode}/${cell_type}: region-rollup Fig 4a PNG"
        "$RSCRIPT" "$PLOT_SCRIPT" \
            --checkpoints "$checkpoint_dir" \
            --master_tsv "$master_tsv" \
            --figdir "$figdir" \
            --outfmt png \
            --qthresh "$QTHRESH" \
            --panel_a_only || echo "[warn] ${mode}/${cell_type}: region-rollup plot failed; continuing"

        echo "[run] ${mode}/${cell_type}: faceted subcluster Fig 4a PNG"
        "$RSCRIPT" "$PLOT_SCRIPT" \
            --checkpoints "$checkpoint_dir" \
            --master_tsv "$master_tsv" \
            --figdir "$figdir" \
            --outfmt png \
            --qthresh "$QTHRESH" \
            --by_celltype \
            --panel_a_only || echo "[warn] ${mode}/${cell_type}: faceted plot failed; continuing"
    done

    echo "Finished ${mode}"
}

run_root "raw" "$RAW_ROOT" "$RAW_FIG_ROOT"
run_root "sqrt" "$SQRT_ROOT" "$SQRT_FIG_ROOT"

echo "All available bar-chart PNGs complete."
