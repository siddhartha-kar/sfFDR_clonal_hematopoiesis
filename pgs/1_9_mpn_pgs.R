# =========================
# Cumulative incidence curves
# =========================

pc_means_mpn <- analysis_df %>%
  summarise(across(starts_with("PC"), ~ mean(.x, na.rm = TRUE)))

pgs_percentile_values_mpn <- quantile(
  analysis_df$pgs_std,
  probs = target_percentiles / 100,
  na.rm = TRUE
)

curve_newdata_mpn <- tibble(
  percentile = target_percentiles,
  pgs_std = as.numeric(pgs_percentile_values_mpn),
  sex = factor(sex_for_curves, levels = levels(analysis_df$sex))
) %>%
  bind_cols(pc_means_mpn[rep(1, length(target_percentiles)), ]) %>%
  mutate(
    curve_label = factor(
      paste0("PGS percentile ", percentile),
      levels = paste0("PGS percentile ", target_percentiles)
    )
  )

sf_curves_mpn <- survfit(cox_model, newdata = curve_newdata_mpn)

curve_df_mpn <- purrr::map_dfr(seq_len(nrow(curve_newdata_mpn)), function(j) {
  tibble(
    attained_age = sf_curves_mpn$time,
    cuminc = 1 - sf_curves_mpn$surv[, j],
    curve_label = curve_newdata_mpn$curve_label[j]
  )
}) %>%
  dplyr::filter(is.finite(attained_age), is.finite(cuminc)) %>%
  dplyr::mutate(cuminc = pmax(pmin(cuminc, 1), 0))

# =========================
# Format HR + p-value label
# =========================

p_val_mpn <- hr_pgs$p_value

if (p_val_mpn == 0) {
  p_label <- "p < 1 %*% 10^-300"
} else {
  exp_val <- floor(log10(p_val_mpn))
  mantissa <- p_val_mpn / (10^exp_val)
  p_label <- sprintf("p == %.2f %%*%% 10^{%d}", mantissa, exp_val)
}

hr_label_expr_mpn <- paste0(
  "'HR for MPN per 1 SD CH PGS' == ", sprintf("%.2f", hr_pgs$HR),
  "~'('*95*'% CI: '*",
  sprintf("%.2f", hr_pgs$lower_95), "*'-'*", sprintf("%.2f", hr_pgs$upper_95),
  "*')'~",
  p_label
)

# =========================
# Plot
# =========================

palette_mpn <- c(
  "#000000",
  "#D55E00",
  "#0072B2",
  "#009E73",
  "#CC79A7",
  "#E69F00"
)
palette_mpn <- palette_mpn[seq_along(levels(curve_df_mpn$curve_label))]

x_rng_mpn <- range(curve_df_mpn$attained_age, na.rm = TRUE)
y_rng_mpn <- range(curve_df_mpn$cuminc, na.rm = TRUE)

x_annot_mpn <- mean(x_rng_mpn)
y_annot_mpn <- y_rng_mpn[1] + 0.95 * diff(y_rng_mpn)

cuminc_plot_mpn <- ggplot(
  curve_df_mpn,
  aes(x = attained_age, y = cuminc, colour = curve_label)
) +
  geom_line(linewidth = 1.2, lineend = "round") +
  scale_colour_manual(values = palette_mpn, name = NULL) +
  scale_x_continuous(
    name = "Age (years)",
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  scale_y_continuous(
    name = NULL,
    labels = scales::percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, NA), clip = "off") +
  annotate(
    "text",
    x = x_annot_mpn,
    y = y_annot_mpn,
    label = hr_label_expr_mpn,
    parse = TRUE,
    hjust = 0.5,
    vjust = 1,
    size = 6.2,
    family = "sans",
    colour = "black",
    lineheight = 1.0
  ) +
  guides(
    colour = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(linewidth = 1.6)
    )
  ) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title = element_text(size = 17, colour = "black"),
    axis.text = element_text(size = 15, colour = "black"),
    axis.line = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks.length = unit(1.2, "mm"),
    legend.position = "none",
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = 16, colour = "black"),
    legend.key.width = unit(10, "mm"),
    legend.key.height = unit(4, "mm"),
    legend.spacing.y = unit(0.5, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 4, 2, 2)
  )