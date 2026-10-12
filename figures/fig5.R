# ============================================
# Inputs:
#   beta_forest.txt (or fig5a.txt)
#   or_forest.txt (or fig5b.txt)
#
# Outputs:
#   fig5a.png / .pdf
#   fig5b.png / .pdf
#   fig5.png  / .pdf
# ============================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(ggtext)
  library(grid)
})

# -----------------------------
# Global theme
# -----------------------------
base_family <- "Helvetica"

theme_set(
  theme_minimal(base_family = base_family, base_size = 9.5)
)

# -----------------------------
# Formatting helpers
# -----------------------------
fmt_num <- function(x, digits = 3) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

fmt_ci <- function(lcl, ucl, digits = 3) {
  paste0("(", fmt_num(lcl, digits), ", ", fmt_num(ucl, digits), ")")
}

fmt_p <- function(p) {
  sapply(p, function(x) {
    if (is.na(x)) return("")
    if (x >= 0.01) return(sprintf("%.3f", round(x, 3)))

    e <- floor(log10(x))
    m <- x / (10^e)
    m_txt <- sprintf("%.1f", round(m, 1))
    paste0(m_txt, " × 10<sup>", e, "</sup>")
  })
}

# -----------------------------
# Read data
# -----------------------------
beta_df <- read_tsv("fig5a.txt", show_col_types = FALSE, trim_ws = TRUE)
or_df   <- read_tsv("fig5b.txt",   show_col_types = FALSE, trim_ws = TRUE)

# -----------------------------
# Data prep
# -----------------------------
prep_panel_a <- function(df) {
  df %>%
    mutate(
      block = case_when(
        outpheno == "Clonal growth rate" ~ "Clonal growth rate",
        grepl("^Plasma", outpheno) ~ "Plasma protein levels",
        outpheno == "Leukocyte telomere length" ~ "Leukocyte telomere length",
        TRUE ~ outpheno
      ),
      block_detail = case_when(
        outpheno == "Plasma FLT3 level"   ~ "FLT3",
        outpheno == "Plasma FLT3LG level" ~ "FLT3LG",
        outpheno == "Plasma KIT level"    ~ "KIT",
        TRUE ~ ""
      )
    )
}

prep_panel_b <- function(df) {
  df %>%
    mutate(
      block = case_when(
        outpheno %in% c("AML", "MDS", "MPN") ~ "Myeloid cancers",
        outpheno == "Mosaic loss of Y" ~ "Mosaic loss of Y",
        TRUE ~ outpheno
      ),
      block_detail = case_when(
        outpheno %in% c("AML", "MDS", "MPN") ~ outpheno,
        TRUE ~ ""
      )
    )
}

add_headers_and_y <- function(df, est_col, block_gap = 0.42) {
  blocks <- unique(df$block)
  out <- list()
  idx <- 1

  for (b in blocks) {
    out[[idx]] <- tibble(
      row_type = "header",
      block = b,
      detail = "",
      rsid = NA_character_,
      ra = NA_character_,
      gene = NA_character_,
      estimate = NA_real_,
      l95cl = NA_real_,
      u95cl = NA_real_,
      p = NA_real_
    )
    idx <- idx + 1

    tmp <- df %>%
      filter(block == b) %>%
      mutate(
        row_type = "data",
        detail   = block_detail,
        estimate = .data[[est_col]]
      ) %>%
      select(row_type, block, detail, rsid, ra, gene, estimate, l95cl, u95cl, p)

    out[[idx]] <- tmp
    idx <- idx + 1
  }

  dat <- bind_rows(out) %>%
    mutate(y = rev(seq_len(n()))) %>%
    group_by(block) %>%
    mutate(row_in_block = cumsum(row_type == "data")) %>%
    ungroup() %>%
    mutate(
      stripe = ifelse(row_type == "data" & row_in_block %% 2 == 0, TRUE, FALSE),
      block_id = match(block, unique(block))
    ) %>%
    mutate(
      # smaller block gap than original -> more compact
      y = y - (block_id - 1) * block_gap,
      gene_lab = ifelse(is.na(gene), NA, paste0("<i>", gene, "</i>")),
      est_lab  = ifelse(row_type == "data", fmt_num(estimate), ""),
      ci_lab   = ifelse(row_type == "data", fmt_ci(l95cl, u95cl), ""),
      p_lab    = ifelse(row_type == "data", fmt_p(p), "")
    )

  dat
}

a_df <- beta_df %>% prep_panel_a() %>% add_headers_and_y(est_col = "beta", block_gap = 0.42)
b_df <- or_df   %>% prep_panel_b() %>% add_headers_and_y(est_col = "or",   block_gap = 0.42)

# -----------------------------
# Background helpers
# -----------------------------
make_background_df <- function(df) {
  row_bg <- df %>%
    filter(row_type == "data", stripe) %>%
    transmute(ymin = y - 0.48, ymax = y + 0.48)

  block_bg <- df %>%
    group_by(block) %>%
    summarise(
      ymin = min(y) - 0.48,
      ymax = max(y) + 0.48,
      .groups = "drop"
    )

  list(row_bg = row_bg, block_bg = block_bg)
}

# -----------------------------
# Column x positions
# -----------------------------
left_cols <- list(
  subgroup = 0.015,
  snp      = 0.155,
  allele   = 0.555,
  gene     = 0.80
)

right_cols <- list(
  est = 0.12,
  ci  = 0.57,
  p   = 0.99
)

# -----------------------------
# Block colours
# -----------------------------
block_cols <- c(
  "Clonal growth rate"         = "#4C78A8",  # muted blue
  "Plasma protein levels"      = "#72B7B2",  # muted teal
  "Leukocyte telomere length"  = "#B279A2",  # muted mauve
  "Myeloid cancers"            = "#E08D79",  # muted salmon
  "Mosaic loss of Y"           = "#8C8C3A"   # muted olive
)

# -----------------------------
# Plot: left table
# -----------------------------
plot_left_table <- function(df) {
  bg <- make_background_df(df)
  top_y <- max(df$y) + 0.78

  ggplot(df, aes(y = y)) +
    geom_rect(
      data = bg$block_bg,
      aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = NA,
      colour = "grey72",
      linewidth = 0.25
    ) +
    geom_rect(
      data = bg$row_bg,
      aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "grey97",
      colour = NA
    ) +
    annotate("text", x = left_cols$snp,    y = top_y, label = "SNP",
             fontface = "bold", hjust = 0, size = 3.25) +
    annotate("text", x = left_cols$allele, y = top_y, label = "CH risk allele",
             fontface = "bold", hjust = 0, size = 3.25) +
    annotate("text", x = left_cols$gene,   y = top_y, label = "Gene",
             fontface = "bold", hjust = 0, size = 3.25) +
    annotate("segment",
             x = 0, xend = 1,
             y = top_y - 0.26, yend = top_y - 0.26,
             linewidth = 0.35, colour = "black") +
    geom_text(
      data = filter(df, row_type == "header"),
      aes(x = left_cols$subgroup, label = block, colour = block),
      hjust = 0, fontface = "bold", size = 3.2,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = block_cols) +
    geom_text(
      data = filter(df, row_type == "data", detail != ""),
      aes(x = left_cols$subgroup, label = detail),
      hjust = 0, size = 2.8, colour = "grey25"
    ) +
    geom_text(
      data = filter(df, row_type == "data"),
      aes(x = left_cols$snp, label = rsid),
      hjust = 0, size = 3.0
    ) +
    geom_text(
      data = filter(df, row_type == "data"),
      aes(x = left_cols$allele, label = ra),
      hjust = 0, size = 3.0
    ) +
    ggtext::geom_richtext(
      data = filter(df, row_type == "data"),
      aes(x = left_cols$gene, label = gene_lab),
      hjust = 0, vjust = 0.5,
      fill = NA, label.color = NA, size = 3.0
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(min(df$y) - 0.22, top_y + 0.08),
      clip = "off"
    ) +
    theme_void() +
    theme(
      text = element_text(family = base_family),
      plot.margin = margin(0.8, 1.0, 0.8, 0.8)
    )
}

# -----------------------------
# Plot: right table
# -----------------------------
plot_right_table <- function(df, est_name = "Beta") {
  bg <- make_background_df(df)
  top_y <- max(df$y) + 0.78

  ggplot(df, aes(y = y)) +
    geom_rect(
      data = bg$block_bg,
      aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = NA,
      colour = "grey72",
      linewidth = 0.25
    ) +
    geom_rect(
      data = bg$row_bg,
      aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "grey97",
      colour = NA
    ) +
    annotate("text", x = right_cols$est, y = top_y, label = est_name,
             fontface = "bold", hjust = 1, size = 3.25) +
    annotate("text", x = right_cols$ci,  y = top_y, label = "95% CI",
             fontface = "bold", hjust = 1, size = 3.25) +
    annotate("text", x = right_cols$p,   y = top_y, label = "P-value",
             fontface = "bold", hjust = 1, size = 3.25) +
    annotate("segment",
             x = 0, xend = 1,
             y = top_y - 0.26, yend = top_y - 0.26,
             linewidth = 0.35, colour = "black") +
    geom_text(
      data = filter(df, row_type == "data"),
      aes(x = right_cols$est, label = est_lab),
      hjust = 1, size = 3.0
    ) +
    geom_text(
      data = filter(df, row_type == "data"),
      aes(x = right_cols$ci, label = ci_lab),
      hjust = 1, size = 3.0
    ) +
    ggtext::geom_richtext(
      data = filter(df, row_type == "data"),
      aes(x = right_cols$p, label = p_lab),
      hjust = 1, vjust = 0.5,
      fill = NA, label.color = NA, size = 3.0
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(min(df$y) - 0.22, top_y + 0.08),
      clip = "off"
    ) +
    theme_void() +
    theme(
      text = element_text(family = base_family),
      plot.margin = margin(0.8, 0.8, 0.8, 1.0)
    )
}

# -----------------------------
# Plot: central forest
# -----------------------------
plot_forest_clean <- function(df, null_value, xlim, breaks, xlab_txt) {
  bg <- make_background_df(df)
  top_y <- max(df$y) + 0.78

  ggplot(filter(df, row_type == "data"), aes(x = estimate, y = y)) +
    geom_rect(
      data = bg$block_bg,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = NA,
      colour = "grey72",
      linewidth = 0.25
    ) +
    geom_rect(
      data = bg$row_bg,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "grey97",
      colour = NA
    ) +
    annotate(
      "segment",
      x = xlim[1], xend = xlim[2],
      y = top_y - 0.26, yend = top_y - 0.26,
      linewidth = 0.35, colour = "black"
    ) +
    geom_vline(
      xintercept = null_value,
      linewidth = 0.38,
      colour = "grey42"
    ) +
    geom_errorbar(
    aes(xmin = l95cl, xmax = u95cl, colour = block),
      orientation = "y",
      width = 0,
      linewidth = 0.50,
      show.legend = FALSE
    ) +
    geom_point(
       aes(fill = block, colour = block),
        size = 2.6,
        shape = 21,
        stroke = 0.50,
        show.legend = FALSE
    ) +
scale_colour_manual(values = block_cols) +
scale_fill_manual(values = block_cols) +
    scale_x_continuous(
      limits = xlim,
      breaks = breaks,
      expand = expansion(mult = c(0.015, 0.015))
    ) +
    coord_cartesian(
      ylim = c(min(df$y) - 0.22, top_y + 0.08),
      clip = "off"
    ) +
    labs(x = xlab_txt, y = NULL) +
    theme_minimal(base_size = 9.5, base_family = base_family) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.30),
      axis.line.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 9.8, margin = margin(t = 2)),
      axis.text.x  = element_text(size = 8.8, colour = "black"),
      axis.ticks.x = element_line(linewidth = 0.25, colour = "black"),
      plot.margin  = margin(0.9, 0.8, 0.6, 0.8)
    )
}

# -----------------------------
# Panel assembly
# -----------------------------
assemble_panel <- function(df, label, null_value, xlim, breaks, xlab_txt, est_name) {
  left  <- plot_left_table(df)
  mid   <- plot_forest_clean(df, null_value, xlim, breaks, xlab_txt)
  right <- plot_right_table(df, est_name)

  core <- plot_grid(
    left, mid, right,
    nrow = 1,
    rel_widths = c(0.445, 0.275, 0.280),
    align = "h",
    axis = "tb"
  )

  ggdraw(core) +
    draw_label(
      label,
      x = 0.0015, y = 0.998,
      hjust = 0, vjust = 1,
      fontface = "bold",
      size = 14,
      fontfamily = base_family
    )
}

# -----------------------------
# Build panels
# -----------------------------
panel_a <- assemble_panel(
  a_df,
  label      = "a",
  null_value = 0,
  xlim       = c(-0.27, 0.32),
  breaks     = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3),
  xlab_txt   = "Beta coefficient",
  est_name   = "Beta"
)

panel_b <- assemble_panel(
  b_df,
  label      = "b",
  null_value = 1,
  xlim       = c(0.90, 2.35),
  breaks     = c(1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2),
  xlab_txt   = "Odds ratio",
  est_name   = "OR"
)

multi_panel <- plot_grid(
  panel_a, panel_b,
  ncol = 1,
  rel_heights = c(1.00, 1.04),
  align = "v"
)

# -----------------------------
# Save outputs
# -----------------------------
ggsave("fig5a.png", panel_a, width = 10.4, height = 4.10, dpi = 600, bg = "white")
ggsave("fig5a.pdf", panel_a, width = 10.4, height = 4.10, device = cairo_pdf, bg = "white")

ggsave("fig5b.png", panel_b, width = 10.4, height = 4.45, dpi = 600, bg = "white")
ggsave("fig5b.pdf", panel_b, width = 10.4, height = 4.45, device = cairo_pdf, bg = "white")

ggsave("fig5.png", multi_panel, width = 10.4, height = 8.05, dpi = 600, bg = "white")
ggsave("fig5.pdf", multi_panel, width = 10.4, height = 8.05, device = cairo_pdf, bg = "white")