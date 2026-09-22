#!/usr/bin/env Rscript
# 4-quadrant scatter: dispersion vs CCC overlap
# Two panels: left=as sender, right=as receiver
# x-axis: net dispersion (n_sig_inc - n_sig_dec)
# y-axis: net CCC (n_strengthening - n_weakening)
# Quadrants: top-right = more variable + more connected
#            bottom-left = less variable + less connected
#            top-left = less variable + more connected
#            bottom-right = more variable + less connected

suppressMessages({library(ggplot2); library(patchwork)})

ct_colors = c(
  "GABAergic_neurons"="#E07B39", "glutamatergic_neurons"="#E8A857",
  "cerebellar_neurons"="#F2C97E", "medium_spiny_neurons"="#C4A862",
  "midbrain_neurons"="#A07840",  "basket_cells"="#8B5E3C",
  "microglia"="#3B7DD8",         "astrocytes"="#5BA4CF",
  "opc"="#2E5FA3",               "oligodendrocytes"="#1A3A6B",
  "vascular_cells"="#7FBEDB",    "ependymal_cells"="#B8D9E8"
)

df = read.csv("/scratch/easmit31/dispersion/dglm/dispersion_ccc_overlap.csv")
df$net_disp = df$n_sig_inc - df$n_sig_dec
df$net_ccc_sender = df$n_str_sender - df$n_wk_sender
df$net_ccc_receiver = df$n_str_receiver - df$n_wk_receiver

# label quadrants
label_quad = function(x, y) {
  ifelse(x >= 0 & y >= 0, "Inc disp + Strengthening",
  ifelse(x <  0 & y >= 0, "Dec disp + Strengthening",
  ifelse(x >= 0 & y <  0, "Inc disp + Weakening",
                           "Dec disp + Weakening")))
}
df$quad_sender   = label_quad(df$net_disp, df$net_ccc_sender)
df$quad_receiver = label_quad(df$net_disp, df$net_ccc_receiver)

make_panel = function(data, y_col, y_lab, title) {
  data$y = data[[y_col]]
  ggplot(data, aes(x=net_disp, y=y, color=cell_type_broad)) +
    geom_hline(yintercept=0, color="grey50", linetype="dashed") +
    geom_vline(xintercept=0, color="grey50", linetype="dashed") +
    geom_point(alpha=0.7, size=1.8) +
    scale_color_manual(values=ct_colors, name=NULL) +
    annotate("text", x=Inf,  y=Inf,  label="Inc disp\nStrengthening", hjust=1.1, vjust=1.5, size=3, color="grey40") +
    annotate("text", x=-Inf, y=Inf,  label="Dec disp\nStrengthening", hjust=-0.1, vjust=1.5, size=3, color="grey40") +
    annotate("text", x=Inf,  y=-Inf, label="Inc disp\nWeakening",     hjust=1.1, vjust=-0.5, size=3, color="grey40") +
    annotate("text", x=-Inf, y=-Inf, label="Dec disp\nWeakening",     hjust=-0.1, vjust=-0.5, size=3, color="grey40") +
    theme_minimal(base_size=11) +
    theme(legend.position="none") +
    labs(x="Net dispersion change\n(n increasing genes - n decreasing genes)",
         y=y_lab, title=title)
}

p_sender   = make_panel(df, "net_ccc_sender",   "Net CCC as sender\n(strengthening - weakening)",   "As sender")
p_receiver = make_panel(df, "net_ccc_receiver", "Net CCC as receiver\n(strengthening - weakening)", "As receiver")

# shared legend
combined = (p_sender | p_receiver) +
  plot_layout(guides='collect') &
  theme(legend.position='bottom') &
  plot_annotation(
    title="Dispersion-age vs. cell-cell communication overlap (baseline, lfsr<0.05)",
    subtitle="Each dot = one subcluster x region. x=0: equal inc/dec genes. y=0: equal strengthening/weakening interactions.")

ggsave("/scratch/easmit31/dispersion/dglm/figures/dispersion_ccc_scatter.png",
       combined, width=14, height=8, dpi=150)
cat("Saved.\n")
