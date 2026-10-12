library(ggplot2)
library(patchwork)

#-----------------------------#
# 1. Read and prepare the data
#-----------------------------#

prepare_data <- function(file, panel_title) {
  dat <- read.delim(
    file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(
    "beta_META", "se_META",
    "beta_UKB_EUR", "se_UKB_EUR"
  )

  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "%s is missing required column(s): %s",
        file,
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  dat[required_cols] <- lapply(dat[required_cols], as.numeric)

  keep <- complete.cases(dat[required_cols]) &
    apply(dat[required_cols], 1, function(x) all(is.finite(x)))

  if (!all(keep)) {
    warning(sum(!keep), " incomplete/non-finite row(s) removed from ", file)
    dat <- dat[keep, , drop = FALSE]
  }

  # Discovery: UKB EUR
  dat$discovery_beta <- dat$beta_UKB_EUR
  dat$discovery_low  <- dat$beta_UKB_EUR - 1.96 * dat$se_UKB_EUR
  dat$discovery_high <- dat$beta_UKB_EUR + 1.96 * dat$se_UKB_EUR

  # Replication: META
  dat$replication_beta <- dat$beta_META
  dat$replication_low  <- dat$beta_META - 1.96 * dat$se_META
  dat$replication_high <- dat$beta_META + 1.96 * dat$se_META

  dat$panel_title <- panel_title
  dat
}

overall <- prepare_data("fig4a.txt", "Overall CH")
dnmt3a <- prepare_data("fig4b.txt", "DNMT3A-CH")

all_limits <- c(
  overall$discovery_low, overall$discovery_high,
  overall$replication_low, overall$replication_high,
  dnmt3a$discovery_low, dnmt3a$discovery_high,
  dnmt3a$replication_low, dnmt3a$replication_high
)

axis_limit <- max(abs(all_limits), na.rm = TRUE)
axis_limit <- ceiling(axis_limit / 0.02) * 0.02
common_limits <- c(-axis_limit, axis_limit)

#-----------------------------#
# 2. Plotting function
#-----------------------------#

make_panel <- function(dat) {
  # Pearson correlation between discovery and replication estimates.
  cor_result <- cor.test(
    dat$discovery_beta,
    dat$replication_beta,
    method = "pearson"
  )

  r_label <- sprintf(
    "italic(r) == %.2f",
    unname(cor_result$estimate)
  )

  if (cor_result$p.value < 0.001) {
    p_scientific <- sprintf("%.1e", cor_result$p.value)
    p_parts <- strsplit(p_scientific, "e", fixed = TRUE)[[1]]
    p_mantissa <- as.numeric(p_parts[1])
    p_exponent <- as.integer(p_parts[2])

    p_label <- sprintf(
      "italic(P) == %.1f %%*%% 10^{%d}",
      p_mantissa,
      p_exponent
    )
  } else {
    p_label <- sprintf(
      "italic(P) == %.2f",
      cor_result$p.value
    )
  }

  n_label <- sprintf(
    "italic(n) == %d",
    nrow(dat)
  )

  ggplot(
    dat,
    aes(x = discovery_beta, y = replication_beta)
  ) +
    # Reference lines
    geom_hline(
      yintercept = 0,
      linewidth = 0.35,
      colour = "grey75"
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.35,
      colour = "grey75"
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.45,
      colour = "grey45"
    ) +

    # 95% confidence intervals for discovery and replication
    geom_errorbar(
      aes(xmin = discovery_low, xmax = discovery_high),
      orientation = "y",
      width = 0,
      linewidth = 0.45
    ) +
    geom_errorbar(
      aes(ymin = replication_low, ymax = replication_high),
      width = 0,
      linewidth = 0.45
    ) +

    # Point estimates
    geom_point(
      shape = 21,
      size = 2.6,
      stroke = 0.55,
      fill = "white"
    ) +

    # Correlation annotations in the upper-left corner.
    annotate(
      "text",
      x = common_limits[1] + 0.05 * diff(common_limits),
      y = common_limits[2] - 0.05 * diff(common_limits),
      label = r_label,
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 2.6,
      family = "Arial"
    ) +
    annotate(
      "text",
      x = common_limits[1] + 0.05 * diff(common_limits),
      y = common_limits[2] - 0.12 * diff(common_limits),
      label = p_label,
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 2.6,
      family = "Arial"
    ) +
    annotate(
      "text",
      x = common_limits[1] + 0.05 * diff(common_limits),
      y = common_limits[2] - 0.19 * diff(common_limits),
      label = n_label,
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = 2.6,
      family = "Arial"
    ) +

    coord_equal(
      xlim = common_limits,
      ylim = common_limits,
      expand = FALSE
    ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 5)
    ) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 5)
    ) +
    labs(
      title = unique(dat$panel_title),
      x = "Discovery beta coefficient",
      y = "Replication beta coefficient"
    ) +
    theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      plot.title = element_text(
        size = 9,
        face = "bold",
        hjust = 0
      ),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7, colour = "black"),
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks.length = grid::unit(1.5, "mm"),
      plot.margin = margin(5, 6, 5, 5)
    )
}

panel_a <- make_panel(overall) +
  labs(
    title = expression(bold("New overall CH lead SNPs")),
    tag = "a"
  ) +
  theme(
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 15
    ),
    plot.tag.position = c(0, 1)
  )

panel_b <- make_panel(dnmt3a) +
  labs(
    title = expression(
      bold("New ") * bolditalic("DNMT3A") * bold("CH lead SNPs")
    ),
    tag = "b"
  ) +
  theme(
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 15
    ),
    plot.tag.position = c(0, 1)
  )

figure <- panel_a + panel_b +
  plot_layout(ncol = 2)

figure

#-----------------------------#
# 3. Export files
#-----------------------------#

ggsave(
  filename = "fig4.pdf",
  plot = figure,
  width = 180,
  height = 88,
  units = "mm",
  device = cairo_pdf
)

ggsave(
  filename = "fig4.tiff",
  plot = figure,
  width = 180,
  height = 88,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)
