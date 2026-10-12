# =================================================
# Multipanel cumulative incidence figure
# =================================================

install.packages("patchwork")
install.packages("cowplot")
library(patchwork)
library(cowplot)
library(grid)

# -------------------------------------------------
# Make sure legends are shown for collection
# -------------------------------------------------

cuminc_plot_aml <- cuminc_plot_aml +
  ylab(NULL) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    colour = guide_legend(
      nrow = 2,
      ncol = 3,
      byrow = TRUE,
      override.aes = list(linewidth = 1.6)
    )
  )

cuminc_plot_mds <- cuminc_plot_mds +
  ylab(NULL) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    colour = guide_legend(
      nrow = 2,
      ncol = 3,
      byrow = TRUE,
      override.aes = list(linewidth = 1.6)
    )
  )

cuminc_plot_mpn <- cuminc_plot_mpn +
  ylab(NULL) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    colour = guide_legend(
      nrow = 2,
      ncol = 3,
      byrow = TRUE,
      override.aes = list(linewidth = 1.6)
    )
  )

# -------------------------------------------------
# Combine plots vertically with shared legend
# -------------------------------------------------

combined_plot <- (
  (cuminc_plot_aml + labs(tag = "a")) /
  (cuminc_plot_mds + labs(tag = "b")) /
  (cuminc_plot_mpn + labs(tag = "c"))
) +
  plot_layout(
    ncol = 1,
    guides = "collect",
    axis_titles = "collect",
    heights = c(1, 1, 1)
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(family = "sans", size = 16, margin = margin(l = 4, r = 8), colour = "black"),
    legend.key.width = unit(10, "mm"),
    legend.key.height = unit(4, "mm"),
    legend.spacing.x = unit(3, "mm"),
    legend.spacing.y = unit(1, "mm"),
    legend.margin = margin(0, 0, 0, 0),

    plot.tag = element_text(
      family = "sans",
      face = "bold",
      size = 18,
      colour = "black"
    ),
    plot.tag.position = c(-0.05, 0.99),

    plot.margin = margin(2, 4, 2, 2)
  ) &
  guides(
    colour = guide_legend(
      nrow = 2,
      ncol = 3,
      byrow = TRUE,
      override.aes = list(linewidth = 1.6)
    )
  )

# -------------------------------------------------
# Add one shared centred y-axis title
# -------------------------------------------------

final_multi_panel_cuminc <- ggdraw() +
  draw_label(
    "Cumulative incidence",
    x = 0.05,
    y = 0.5,
    angle = 90,
    hjust = 0.5,
    vjust = 0.5,
    size = 17,
    fontfamily = "sans",
    colour = "black"
  ) +
  draw_plot(
    combined_plot,
    x = 0.07,
    y = 0,
    width = 0.93,
    height = 1
  )

# -------------------------------------------------
# Display final figure
# -------------------------------------------------

final_multi_panel_cuminc

# -------------------------------------------------
# Save outputs
# -------------------------------------------------

ggsave(
  filename = "fig6.pdf",
  plot = final_multi_panel_cuminc,
  width = 89,
  height = 180,
  units = "mm"
)

ggsave(
  filename = "fig6.png",
  plot = final_multi_panel_cuminc,
  width = 89,
  height = 180,
  units = "mm",
  dpi = 600
)