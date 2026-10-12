
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

model_levels <- c(
  "CH PGS, age, age2, sex, 10 PCs",
  "CH PGS, LTL PGS, age, age2, sex, 10 PCs",
  "CH PGS, LTL PGS, baseline CH, age, age2, sex, 10 PCs"
)

model_cols <- c(
  "CH PGS, age, age2, sex, 10 PCs" = "#0072B2",
  "CH PGS, LTL PGS, age, age2, sex, 10 PCs" = "#D55E00",
  "CH PGS, LTL PGS, baseline CH, age, age2, sex, 10 PCs" = "#4D4D4D"
)

model_labels <- c(
  "CH PGS, age, age2, sex, 10 PCs" =
    "CH PGS + age + age² + sex + 10 PCs",
  "CH PGS, LTL PGS, age, age2, sex, 10 PCs" =
    "CH PGS + LTL PGS + age + age² + sex + 10 PCs",
  "CH PGS, LTL PGS, baseline CH, age, age2, sex, 10 PCs" =
    "CH PGS + LTL PGS + baseline CH + age + age² + sex + 10 PCs"
)

fmt_p_plotmath <- function(p) {
  case_when(
    is.na(p) ~ "italic(P) == plain('-')",
    p < 0.001 ~ {
      exponent <- floor(log10(p))
      mantissa <- p / 10^exponent
      sprintf("italic(P) == %.1f %%*%% 10^{%d}", mantissa, exponent)
    },
    TRUE ~ sprintf("italic(P) == %.3f", p)
  )
}

theme_forest <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = base_family, colour = "black"),
      axis.text = element_text(colour = "black", size = base_size + 2),
      axis.title = element_text(size = base_size + 4),

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

      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size + 2),
      legend.box.spacing = unit(0.5, "mm"),
      legend.margin = margin(0, 0, 1, 0),
      legend.spacing.y = unit(0.5, "mm"),

      plot.title = element_blank(),
      plot.margin = margin(3, 8, 3, 8)
    )
}

infile <- "ed_fig5.txt"

dat <- read_clean(infile) |>
  transmute(
    cohort = as.character(cohort),
    model = as.character(model),
    beta = as.numeric(beta),
    l95cl = as.numeric(l95cl),
    u95cl = as.numeric(u95cl),
    p = as.numeric(p)
  ) |>
  mutate(
    cohort_label = recode(
      cohort,
      "Full" = "Full cohort",
      "No CH" = "No detectable CH at recruitment",
      .default = cohort
    ),
    cohort_label = factor(
      cohort_label,
      levels = c(
        "Full cohort",
        "No detectable CH at recruitment"
      )
    ),
    model = factor(model, levels = model_levels),
    p_label = fmt_p_plotmath(p)
  ) |>
  arrange(cohort_label, model)

stopifnot(
  all(is.finite(dat$beta)),
  all(is.finite(dat$l95cl)),
  all(is.finite(dat$u95cl)),
  all(dat$l95cl <= dat$beta & dat$beta <= dat$u95cl)
)

ci_min <- min(dat$l95cl, na.rm = TRUE)
ci_max <- max(dat$u95cl, na.rm = TRUE)
ci_span <- ci_max - ci_min

x_left  <- min(0, ci_min - 0.08 * ci_span)
x_p     <- ci_max + 0.22 * ci_span
x_right <- ci_max + 1.05 * ci_span



make_panel <- function(panel_dat, panel_title, show_x_axis = FALSE, show_legend = TRUE) {

  panel_dat <- panel_dat |>
    arrange(model) |>
    mutate(
      y_pos = rev(seq_len(n()))
    )

  p <- ggplot(panel_dat) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.3,
      linetype = "dashed",
      colour = "black"
    ) +

    geom_segment(
      aes(
        x = l95cl,
        xend = u95cl,
        y = y_pos,
        yend = y_pos,
        colour = model
      ),
      linewidth = 0.55,
      lineend = "butt"
    ) +

    geom_segment(
      aes(
        x = l95cl,
        xend = l95cl,
        y = y_pos - 0.08,
        yend = y_pos + 0.08,
        colour = model
      ),
      linewidth = 0.4,
      lineend = "butt"
    ) +

    geom_segment(
      aes(
        x = u95cl,
        xend = u95cl,
        y = y_pos - 0.08,
        yend = y_pos + 0.08,
        colour = model
      ),
      linewidth = 0.4,
      lineend = "butt"
    ) +

    geom_point(
      aes(
        x = beta,
        y = y_pos,
        colour = model
      ),
      size = 2.1
    ) +

    geom_text(
      aes(
        x = x_p,
        y = y_pos,
        label = p_label
      ),
      inherit.aes = FALSE,
      parse = TRUE,
      hjust = 0,
      family = base_family,
      size = 3.1,
      colour = "black"
    ) +

    facet_wrap(
      ~ factor(panel_title, levels = panel_title),
      ncol = 1,
      strip.position = "top"
    ) +

    scale_colour_manual(
      values = model_cols,
      breaks = model_levels,
      labels = model_labels,
      limits = model_levels,
      drop = FALSE
    ) +

    scale_y_continuous(
      breaks = NULL,
      expand = expansion(add = c(0.12, 0.12))
    ) +

    scale_x_continuous(
      limits = c(x_left, x_right),
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0, 0))
    ) +

    guides(
      colour = guide_legend(
        ncol = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = 0.7,
          size = 2.0
        ),
        keywidth = unit(5, "mm"),
        keyheight = unit(3.0, "mm")
      )
    ) +

    labs(
      x = if (show_x_axis) {
        "beta coefficient from linear regression (95% CI)"
      } else {
        NULL
      },
      y = NULL
    ) +

    coord_cartesian(clip = "off") +
    theme_forest() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = if (show_legend) "top" else "none"
    )

  if (!show_x_axis) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      )
  }

  p
}



dat_full <- dat |>
  filter(cohort_label == "Full cohort")

dat_no_ch <- dat |>
  filter(cohort_label == "No detectable CH at recruitment")

p_full <- make_panel(
  panel_dat = dat_full,
  panel_title = "Full cohort",
  show_x_axis = FALSE,
  show_legend = TRUE
)

p_no_ch <- make_panel(
  panel_dat = dat_no_ch,
  panel_title = "No detectable CH at recruitment",
  show_x_axis = TRUE,
  show_legend = FALSE
)


p_forest <- p_full / p_no_ch +
  plot_layout(
    heights = c(nrow(dat_full), nrow(dat_no_ch))
  )



ggsave(
  "ed_fig5.pdf",
  p_forest,
  width = 7.6,
  height = 4.8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  "ed_fig5.png",
  p_forest,
  width = 7.6,
  height = 4.8,
  units = "in",
  dpi = 600
)

p_forest
