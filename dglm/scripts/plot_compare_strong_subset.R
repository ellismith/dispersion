#!/usr/bin/env Rscript
# plot_compare_strong_subset.R — old (q<0.05) vs new (q<0.2) strong-subset runs.
# Reads both master tables per cell type; skips types missing either one.
# Significance = mash_lfsr < 0.05; direction = sign(mash_beta); genes keyed by ensembl_id.
# Outputs: figures_915q02/compare_q005_vs_q02_celltype.png, ..._region.png
suppressMessages(library(ggplot2))
ab <- c(GABAergic_neurons="INH", astrocytes="AST", glutamatergic_neurons="EXC", oligodendrocytes="OLIG",
        ependymal_cells="EPEN", medium_spiny_neurons="MSN", midbrain_neurons="MBN", microglia="MGL",
        opc="OPC", vascular_cells="VASC")
runs <- c("q < 0.05 (old)"="disp_age__mask1_915_300c", "q < 0.2 (new)"="disp_age__mask1_915_300c_q0.2")
ct_rows <- list(); rg_rows <- list()
for (ct in names(ab)) {
  fs <- file.path(runs, ct, "dglm_checkpoints_cutoff0.5", "master_dglm_combined.tsv")
  if (!all(file.exists(fs))) { message("skip (not finished): ", ct); next }
  for (i in seq_along(runs)) {
    m <- read.delim(fs[i]); s <- m[!is.na(m$mash_lfsr) & m$mash_lfsr < 0.05, ]
    ct_rows[[length(ct_rows)+1]] <- data.frame(ct=ab[ct], run=names(runs)[i],
      inc=length(unique(s$ensembl_id[s$mash_beta > 0])), dec=length(unique(s$ensembl_id[s$mash_beta < 0])))
    for (r in unique(s$region)) { x <- s[s$region == r, ]
      rg_rows[[length(rg_rows)+1]] <- data.frame(ct=ab[ct], region=r, run=names(runs)[i],
        inc=length(unique(x$ensembl_id[x$mash_beta > 0])), dec=length(unique(x$ensembl_id[x$mash_beta < 0]))) }
  }
}
d <- do.call(rbind, ct_rows); g <- do.call(rbind, rg_rows)
d$run <- factor(d$run, names(runs)); g$run <- factor(g$run, names(runs))
d$ct <- factor(d$ct, sort(unique(as.character(d$ct)))); g$ct <- factor(g$ct, levels(d$ct))
cols <- c("q < 0.05 (old)"="grey65", "q < 0.2 (new)"="#1f3a5f")

p1 <- ggplot(d, aes(x=ct, fill=run)) +
  geom_col(aes(y=inc), position=position_dodge(0.8), width=0.75) +
  geom_col(aes(y=-dec), position=position_dodge(0.8), width=0.75, alpha=0.55) +
  geom_text(aes(y=inc, label=inc), position=position_dodge(0.8), vjust=-0.3, size=3.2) +
  geom_text(aes(y=-dec, label=dec), position=position_dodge(0.8), vjust=1.3, size=3.2) +
  geom_hline(yintercept=0) + scale_fill_manual(values=cols, name=NULL) +
  scale_y_continuous(labels=abs, expand=expansion(mult=0.12)) +
  labs(x=NULL, y="Unique significant genes (lfsr < 0.05)\n\u2191 increasing dispersion   \u2193 decreasing (transparent)",
       title="Strong-subset rule: q < 0.05 vs q < 0.2 (qval_single, \u2265 1 condition)") +
  theme_bw(base_size=14) + theme(legend.position="top", panel.grid.major.x=element_blank())
ggsave("figures_915q02/compare_q005_vs_q02_celltype.png", p1, width=10, height=6.5, dpi=200)

p2 <- ggplot(g, aes(x=region, fill=run)) +
  geom_col(aes(y=inc), position=position_dodge(0.8), width=0.75) +
  geom_col(aes(y=-dec), position=position_dodge(0.8), width=0.75, alpha=0.55) +
  geom_hline(yintercept=0) + facet_wrap(~ct, scales="free", ncol=2) +
  scale_fill_manual(values=cols, name=NULL) + scale_y_continuous(labels=abs) +
  labs(x=NULL, y="Unique significant genes (\u2191 inc / \u2193 dec)", title="By region: q < 0.05 vs q < 0.2") +
  theme_bw(base_size=12) + theme(legend.position="top", axis.text.x=element_text(angle=45, hjust=1))
ggsave("figures_915q02/compare_q005_vs_q02_region.png", p2, width=11, height=2.6*ceiling(nlevels(g$ct)/2)+1.5, dpi=200)
cat("Saved: figures_915q02/compare_q005_vs_q02_celltype.png\nSaved: figures_915q02/compare_q005_vs_q02_region.png\n")
print(d, row.names=FALSE)
