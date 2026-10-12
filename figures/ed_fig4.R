
library(tidyverse)
library(patchwork)
library(grid)

read_clean <- function(file) {
  x <- read.delim(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "")
  )
  x[, !is.na(names(x)) & names(x) != "", drop = FALSE]
}

base_family <- "sans"

phenotype_levels <- c("AML", "MDS", "MPN", "AML+MDS+MPN")

phenotype_cols <- c(
  "AML" = "#0072B2",
  "MDS" = "#D55E00",
  "MPN" = "#4D4D4D",
  "AML+MDS+MPN" = "#009E73"
)

fmt_p_label <- function(p) {
  case_when(
    is.na(p) ~ "plain('-')",
    p < 0.001 ~ {
      exponent <- floor(log10(p))
      mantissa <- p / 10^exponent
      sprintf("%.1f %%*%% 10^{%d}", mantissa, exponent)
    },
    TRUE ~ sprintf("%.2f", p)
  )
}

theme_forest <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(
        family = base_family,
        colour = "black"
      ),
      axis.text = element_text(
        colour = "black",
        size = base_size + 2
      ),
      axis.title = element_text(
        size = base_size + 4
      ),
      strip.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.5
      ),
      strip.text = element_text(
        family = base_family,
        colour = "black",
        size = base_size + 3,
        face = "plain",
        margin = margin(2, 2, 2, 2)
      ),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      ),
      panel.grid.major.x = element_line(
        colour = "grey88",
        linewidth = 0.25
      ),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_blank(),

      plot.margin = margin(3, 12, 3, 8)
    )
}

# --------------------------------------------------
# Read and prepare data
# --------------------------------------------------

infile <- "ed_fig4.txt"

dat <- read_clean(infile) |>
  transmute(
    phenotype = factor(
      phenotype,
      levels = phenotype_levels
    ),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p)
  ) |>
  arrange(phenotype) |>
  mutate(
    or = exp(beta),
    l95cl = exp(beta - 1.96 * se),
    u95cl = exp(beta + 1.96 * se),

    y_pos = rev(seq_len(n())),

    or_ci_label = sprintf(
      "%.2f (%.2f, %.2f)",
      or, l95cl, u95cl
    ),

    p_label = fmt_p_label(p),

    panel_title =
      "CH PGS-outcome associations in OFH"
  )

stopifnot(
  all(!is.na(dat$phenotype)),
  all(is.finite(dat$or)),
  all(is.finite(dat$l95cl)),
  all(is.finite(dat$u95cl)),
  all(
    dat$l95cl <= dat$or &
    dat$or <= dat$u95cl
  )
)

y_limits <- c(0.5, nrow(dat) + 0.85)
y_header <- nrow(dat) + 0.52

# --------------------------------------------------
# LEFT PANEL: Forest plot
# --------------------------------------------------

p_forest <- ggplot(
  dat,
  aes(x = or, y = y_pos, colour = phenotype)
) +

  geom_vline(
    xintercept = 1,
    linewidth = 0.3,
    linetype = "dashed",
    colour = "black"
  ) +

  geom_segment(
    aes(
      x = l95cl,
      xend = u95cl,
      yend = y_pos
    ),
    linewidth = 0.55,
    lineend = "butt"
  ) +

  geom_segment(
    aes(
      x = l95cl,
      xend = l95cl,
      y = y_pos - 0.08,
      yend = y_pos + 0.08
    ),
    linewidth = 0.4,
    lineend = "butt"
  ) +

  geom_segment(
    aes(
      x = u95cl,
      xend = u95cl,
      y = y_pos - 0.08,
      yend = y_pos + 0.08
    ),
    linewidth = 0.4,
    lineend = "butt"
  ) +

  geom_point(size = 2.1) +

  facet_wrap(
    ~ panel_title,
    ncol = 1,
    strip.position = "top"
  ) +

  scale_colour_manual(
    values = phenotype_cols,
    guide = "none"
  ) +

  scale_y_continuous(
    limits = y_limits,
    breaks = dat$y_pos,
    labels = as.character(dat$phenotype),
    expand = expansion(mult = c(0, 0))
  ) +

  scale_x_log10(
    breaks = c(0.8, 1, 1.2, 1.5, 1.8),
    labels = c("0.8", "1", "1.2", "1.5", "1.8"),
    limits = c(0.7, 1.95),
    expand = expansion(mult = c(0, 0))
  ) +

  labs(
    x = "Odds ratio (95% CI)",
    y = NULL
  ) +

  coord_cartesian(clip = "off") +

  theme_forest() +

  theme(
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

# --------------------------------------------------
# RIGHT PANEL: OR (95% CI) and P-value table
# --------------------------------------------------

p_table <- ggplot(dat, aes(y = y_pos)) +

  geom_text(
    aes(
      x = 0,
      label = or_ci_label
    ),
    hjust = 0,
    colour = "black",
    size = 3.1,
    family = base_family
  ) +

  geom_text(
    aes(
      x = 0.85,
      label = p_label
    ),
    hjust = 0,
    parse = TRUE,
    colour = "black",
    size = 3.1,
    family = base_family
  ) +

  annotate(
    "text",
    x = 0,
    y = y_header,
    label = "OR (95% CI)",
    hjust = 0,
    fontface = "bold",
    size = 3.2,
    family = base_family,
    colour = "black"
  ) +

  annotate(
    "text",
    x = 0.85,
    y = y_header,
    label = "P-value",
    hjust = 0,
    fontface = "bold",
    size = 3.2,
    family = base_family,
    colour = "black"
  ) +

  facet_wrap(
    ~ panel_title,
    ncol = 1,
    strip.position = "top"
  ) +

  scale_y_continuous(
    limits = y_limits,
    expand = expansion(mult = c(0, 0))
  ) +

  scale_x_continuous(
    limits = c(-0.02, 1.25),
    expand = expansion(mult = c(0, 0))
  ) +

  coord_cartesian(clip = "off") +

  theme_void(base_size = 8) +

  theme(
    text = element_text(
      family = base_family,
      colour = "black"
    ),

    strip.text = element_text(
      family = base_family,
      colour = "white",
      size = 11,
      margin = margin(2, 2, 2, 2)
    ),

    strip.background = element_rect(
      fill = "white",
      colour = NA
    ),

    plot.margin = margin(3, 2, 3, 1)
  )

# --------------------------------------------------
# Combine panels
# --------------------------------------------------

p_combined <- p_forest + p_table +
  plot_layout(
    ncol = 2,
    widths = c(1.5, 0.70)
  )

# --------------------------------------------------
# Export
# --------------------------------------------------

ggsave(
  "ed_fig4.pdf",
  p_combined,
  width = 7.6,
  height = 3.3,
  units = "in",
  device = cairo_pdf
)

ggsave(
  "ed_fig4.png",
  p_combined,
  width = 7.6,
  height = 3.3,
  units = "in",
  dpi = 600
)

p_combined
