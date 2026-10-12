
library(tidyverse)
library(patchwork)
library(gridExtra)
library(grid)
library(gtable)

# =========================
# Helpers
# =========================

label_a_size <- 16
label_b_size <- 16
label_c_size <- 16

label_face <- "bold"

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
body_fontsize <- 9
header_fontsize <- 9

ng_cols <- c(
  "Age, sex" = "#4D4D4D",
  "Age, sex, CH PGS, 10 PCs" = "#0072B2",
  "Age, sex, baseline CH" = "#999999",
  "Age, sex, baseline CH, CH PGS, 10 PCs" = "#D55E00"
)

theme_ng <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = base_family, colour = "black"),
      axis.text = element_text(colour = "black", size = base_size + 3),
      axis.title = element_text(size = base_size + 5),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
      strip.text = element_text(family = base_family, colour = "black", size = base_size + 3, face = "plain"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size + 3),
      plot.title = element_blank(),
      plot.margin = margin(8, 10, 8, 10)
    )
}

fmt_p_plain <- function(p) {
  case_when(
    is.na(p) ~ "–",
    p < 1e-3 ~ {
      exp <- floor(log10(p))
      mant <- round(p / 10^exp, 1)
      paste0(mant, " x 10^", exp)
    },
    TRUE ~ sprintf("%.3f", p)
  )
}

find_cell <- function(tbl, r, c, nm) {
  which(tbl$layout$t == r & tbl$layout$l == c & tbl$layout$name == nm)
}

# =========================
# Read data
# =========================

fig7a <- read_clean("fig7a.txt") |>
  select(cohort, model, outcome, cindex, l95cl, u95cl)

fig7b <- read_clean("fig7b.txt") |>
  mutate(.row_id = row_number()) |>
  select(.row_id, cohort, model, outcome, deltac, l95cl, u95cl, p)

fig7c_raw <- read_clean("fig7c.txt")

if ("HR" %in% names(fig7c_raw) && !"hr" %in% names(fig7c_raw)) {
  fig7c_raw <- fig7c_raw |> rename(hr = HR)
}

fig7c <- fig7c_raw |>
  select(cohort, model, outcome, hr, l95cl, u95cl, p)

# =========================
# Ordering
# =========================

outcome_levels <- c("AML", "MDS", "MPN")
outcome_levels_c <- c("MPN", "MDS", "AML")  

cohort_levels <- c("Full", "No CH", "No CH at recruitment")
cohort_label_levels <- c("Full cohort", "No detectable CH at recruitment")

model_levels <- c(
  "Age, sex",
  "Age, sex, CH PGS, 10 PCs",
  "Age, sex, baseline CH",
  "Age, sex, baseline CH, CH PGS, 10 PCs"
)

fig7a <- fig7a |>
  mutate(
    outcome = factor(outcome, outcome_levels),
    cohort = factor(cohort, cohort_levels),
    cohort_label = case_when(
      cohort == "Full" ~ "Full cohort",
      cohort %in% c("No CH", "No CH at recruitment") ~ "No detectable CH at recruitment",
      TRUE ~ as.character(cohort)
    ),
    cohort_label = factor(cohort_label, cohort_label_levels),
    model = factor(model, model_levels)
  )

fig7b <- fig7b |>
  mutate(
    outcome = factor(outcome, outcome_levels),
    cohort = factor(cohort, cohort_levels),
    delta_ci = ifelse(
      is.na(deltac),
      "–",
      sprintf("%.3f (%.3f–%.3f)", deltac, l95cl, u95cl)
    ),
    p_fmt = fmt_p_plain(p)
  )

fig7c <- fig7c |>
  mutate(
    outcome = factor(outcome, outcome_levels_c),
    cohort = factor(cohort, cohort_levels),
    cohort_label = case_when(
      cohort == "Full" ~ "Full cohort",
      cohort %in% c("No CH", "No CH at recruitment") ~ "No detectable CH at recruitment",
      TRUE ~ as.character(cohort)
    ),
    cohort_label = factor(cohort_label, cohort_label_levels),
    model = factor(model, model_levels),
    p_fmt = fmt_p_plain(p)
  )

# =========================
# Panel a
# =========================

p_a <- ggplot(
  fig7a,
  aes(x = outcome, y = cindex, ymin = l95cl, ymax = u95cl, colour = model)
) +
  geom_errorbar(
    width = 0.16,
    linewidth = 0.35,
    position = position_dodge(width = 0.55)
  ) +
  geom_point(
    size = 1.6,
    position = position_dodge(width = 0.55)
  ) +
  facet_wrap(~ cohort_label, ncol = 1, strip.position = "top") +
  scale_colour_manual(
    values = ng_cols,
    limits = model_levels,
    drop = FALSE
  ) +
  guides(
    colour = guide_legend(
      ncol = 2,
      byrow = TRUE,
      keywidth = unit(4, "mm"),
      keyheight = unit(3, "mm")
    )
  ) +
  scale_y_continuous(
    limits = c(0.50, 0.75),
    breaks = seq(0.50, 0.75, 0.05),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  labs(x = NULL, y = "C-index (95% CI)") +
  theme_ng() +
  theme(
    panel.spacing.y = unit(0.75, "lines"),
    legend.box.spacing = unit(1, "mm"),
    legend.margin = margin(0, 0, 2, 0),
    legend.spacing.y = unit(1.2, "mm")
  )

p_a <- p_a +
  labs(tag = "a") +
  theme(
    plot.tag = element_text(
      family = base_family,
      face = label_face,
      size = label_a_size,
      colour = "black"
    ),
    plot.tag.position = c(0.005, 0.985)
  )

# =========================
# Panel b table
# =========================

tab_b <- fig7b |>
  arrange(.row_id) |>
  transmute(
    Cohort = as.character(cohort),
    Outcome = as.character(outcome),
    Model = model,
    `ΔC-index (95% CI)` = delta_ci,
    P = p_fmt
  )

tbl_theme_b <- ttheme_minimal(
  base_size = body_fontsize,
  base_family = base_family,
  core = list(
    fg_params = list(
      hjust = 0,
      x = 0.02,
      fontsize = body_fontsize,
      col = "black"
    ),
    bg_params = list(
      fill = rep(c("white", "#EEECE2"), length.out = nrow(tab_b)),
      col = NA
    ),
    padding = unit(c(1.6, 2.2), "mm")
  ),
  colhead = list(
    fg_params = list(
      fontface = "bold",
      hjust = 0,
      x = 0.02,
      fontsize = header_fontsize,
      col = "black"
    ),
    bg_params = list(
      fill = "#D8D5C5",
      col = NA
    ),
    padding = unit(c(1.8, 2.2), "mm")
  )
)

g_b <- tableGrob(tab_b, rows = NULL, theme = tbl_theme_b)

g_b$widths <- unit(
  c(0.52, 0.55, 2.55, 1.35, 0.70),
  "in"
)

p_col_index <- which(colnames(tab_b) == "P")

header_idx <- find_cell(g_b, 1, p_col_index, "colhead-fg")

if (length(header_idx) == 1) {
  g_b$grobs[[header_idx]] <- textGrob(
    expression(bolditalic(P)),
    x = 0.02,
    hjust = 0,
    gp = gpar(
      fontsize = header_fontsize,
      fontfamily = base_family,
      col = "black"
    )
  )
}

for (i in seq_len(nrow(tab_b))) {
  val <- tab_b[[p_col_index]][i]

  if (!is.na(val) && grepl(" x 10\\^", val)) {
    parts <- strsplit(val, " x 10\\^")[[1]]
    mant <- as.numeric(parts[1])
    exponent <- as.numeric(parts[2])

    body_idx <- find_cell(g_b, i + 1, p_col_index, "core-fg")

    if (length(body_idx) == 1) {
      g_b$grobs[[body_idx]] <- textGrob(
        bquote(.(mant) %*% 10^.(exponent)),
        x = 0.02,
        hjust = 0,
        gp = gpar(
          fontsize = body_fontsize,
          fontfamily = base_family,
          col = "black"
        )
      )
    }
  }
}

n_rows_b <- nrow(tab_b)

g_b <- gtable_add_grob(
  g_b,
  segmentsGrob(
    x0 = 0, x1 = 1, y0 = 1, y1 = 1,
    gp = gpar(lwd = 1.0, col = "#3A3A3A")
  ),
  t = 1, b = 1, l = 1, r = ncol(g_b)
)

for (i in 1:(n_rows_b + 1)) {
  g_b <- gtable_add_grob(
    g_b,
    segmentsGrob(
      x0 = 0, x1 = 1, y0 = 0, y1 = 0,
      gp = gpar(lwd = 0.7, col = "#3A3A3A")
    ),
    t = i, b = i, l = 1, r = ncol(g_b)
  )
}

p_b <- wrap_elements(g_b) +
  theme(plot.margin = margin(8, 10, 8, 10))

p_b <- p_b +
  labs(tag = "b") +
  theme(
    plot.tag = element_text(
      family = base_family,
      face = label_face,
      size = label_b_size,
      colour = "black"
    ),
    plot.tag.position = c(0.005, 0.985)
  )

# =========================
# Panel c
# =========================

p_c <- ggplot(
  fig7c,
  aes(x = hr, y = outcome, xmin = l95cl, xmax = u95cl, colour = model)
) +
  geom_vline(
    xintercept = 1,
    linewidth = 0.25,
    linetype = "dashed"
  ) +
  geom_errorbarh(
    height = 0.18,
    linewidth = 0.35,
    position = position_dodge(width = 0.5)
  ) +
  geom_point(
    size = 1.6,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~ cohort_label, ncol = 1, strip.position = "top") +
  scale_colour_manual(
    values = ng_cols,
    limits = model_levels,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(1.00, 1.35),
    breaks = seq(1.00, 1.35, 0.05),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  labs(
    x = "Hazard ratio per s.d. CH PGS (95% CI)",
    y = NULL
  ) +
  theme_ng() +
  theme(
    panel.spacing.y = unit(0.75, "lines"),
    legend.position = "none"
  )

p_c <- p_c +
  labs(tag = "c") +
  theme(
    plot.tag = element_text(
      family = base_family,
      face = label_face,
      size = label_c_size,
      colour = "black"
    ),
    plot.tag.position = c(0.005, 0.985)
  )

# =========================
# Assemble and save
# =========================

fig7 <- p_a / plot_spacer() / p_b / plot_spacer() / p_c +
  plot_layout(
    heights = c(1.25, 0.10, 1.30, 0.10, 1.15)
  ) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(6, 6, 6, 6)
    )
  )

ggsave(
  "fig7.pdf",
  fig7,
  width = 7.6,
  height = 10.6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  "fig7.png",
  fig7,
  width = 7.6,
  height = 10.6,
  units = "in",
  dpi = 600
)

fig7