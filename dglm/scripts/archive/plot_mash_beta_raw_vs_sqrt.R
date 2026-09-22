#!/usr/bin/env Rscript
# plot_mash_beta_raw_vs_sqrt.R
#
# Direct comparison of mashr's shrunk mash_beta between --shat_mode raw and
# --shat_mode sqrt, for the same genes/conditions -- shows what mashr
# actually did differently with identical raw beta input, which is the
# real question (not whether raw beta itself differs -- it doesn't).

library(ggplot2)

raw = read.table('/scratch/easmit31/dispersion/dglm/disp_age__mean_full/master_dglm_combined.tsv',
                  sep='\t', header=TRUE, stringsAsFactors=FALSE)
sqrt_run = read.table('/scratch/easmit31/dispersion/dglm/disp_age__mean_full_SQRT/master_dglm_combined.tsv',
                       sep='\t', header=TRUE, stringsAsFactors=FALSE)

raw$key = paste(raw$ensembl_id, raw$cell_type, raw$region, sep='|')
sqrt_run$key = paste(sqrt_run$ensembl_id, sqrt_run$cell_type, sqrt_run$region, sep='|')

merged = merge(
    raw[, c('key','cell_type','region','mash_beta','mash_lfsr')],
    sqrt_run[, c('key','mash_beta','mash_lfsr')],
    by='key', suffixes=c('_raw','_sqrt')
)
message('Matched ', nrow(merged), ' gene x condition pairs between the two runs')

p1 = ggplot(merged, aes(mash_beta_raw, mash_beta_sqrt)) +
    geom_point(alpha=0.1, size=0.5) +
    geom_abline(slope=1, intercept=0, color='red', linetype=2) +
    theme_classic(base_size=14) +
    xlab('mash_beta (raw Shat)') + ylab('mash_beta (sqrt Shat)') +
    ggtitle('mash_beta: raw vs sqrt Shat mode, same genes/conditions')
ggsave(p1, file='/scratch/easmit31/dispersion/dglm/figures/mash_beta_raw_vs_sqrt.png', width=8, height=8, dpi=150)
message('Saved: mash_beta_raw_vs_sqrt.png')

# Focus on the two known outlier conditions specifically
for (combo in list(c('ependymal_cells','MB'), c('glutamatergic_neurons','NAc'))) {
    sub = merged[merged$cell_type == combo[1] & merged$region == combo[2], ]
    if (nrow(sub) == 0) next
    p = ggplot(sub, aes(mash_beta_raw, mash_beta_sqrt)) +
        geom_point(alpha=0.3) +
        geom_abline(slope=1, intercept=0, color='red', linetype=2) +
        theme_classic(base_size=14) +
        xlab('mash_beta (raw Shat)') + ylab('mash_beta (sqrt Shat)') +
        ggtitle(paste0(combo[1], ' | ', combo[2], ' — mash_beta raw vs sqrt (n=', nrow(sub), ')'))
    out = paste0('/scratch/easmit31/dispersion/dglm/figures/mash_beta_raw_vs_sqrt_', combo[1], '_', combo[2], '.png')
    ggsave(p, file=out, width=7, height=7, dpi=150)
    message('Saved: ', out)
}
message('done.')
