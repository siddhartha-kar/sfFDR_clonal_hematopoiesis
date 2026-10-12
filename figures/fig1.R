library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggtext)
library(grid)

# ----------------------------
# Editable paths / filenames
# ----------------------------
infile_1 <- "fig1.txt"
out_dir   <- "."

file_1_png <- file.path(out_dir, "fig1.png")
file_1_pdf <- file.path(out_dir, "fig1.pdf")

# Export size
single_width  <- 16
single_height <- 7.2
single_dpi    <- 600

# ----------------------------
# Shared settings
# ----------------------------
trait_order <- c(
  "WBC", "NEU", "LYM", "MON", "BAS", "EOS",
  "RBC", "HGB", "HCT", "MCH", "MCV", "MCHC", "RDW", "PLT", "MPV"
)

trait_label_map <- c(
  WBC  = "WBC",
  NEU  = "Neutrophils",
  LYM  = "Lymphocytes",
  MON  = "Monocytes",
  BAS  = "Basophils",
  EOS  = "Eosinophils",
  RBC  = "RBC",
  HGB  = "Hemoglobin",
  HCT  = "Hematocrit",
  MCH  = "MCH",
  MCV  = "MCV",
  MCHC = "MCHC",
  RDW  = "RDW",
  PLT  = "Platelets",
  MPV  = "MPV"
)

fill_breaks <- c(
  "p < 5e-8 / beta > 0",
  "p < 5e-8 / beta < 0",
  "p < 0.05 / beta > 0",
  "p < 0.05 / beta < 0"
)

fill_labels <- c(
  "p &lt; 5 × 10<sup>-8</sup>, beta &gt; 0",
  "p &lt; 5 × 10<sup>-8</sup>, beta &lt; 0",
  "p &lt; 0.05, beta &gt; 0",
  "p &lt; 0.05, beta &lt; 0"
)

# ----------------------------
# Helpers
# ----------------------------
build_plot_data <- function(dat, gene_col) {
  dat <- dat %>%
    mutate(
      chr_num = suppressWarnings(as.numeric(chr)),
      chr_label = paste0("chr", chr_num),
      snp_id_plain = paste(chr_label, rsID, CH_risk_increasing_allele, .data[[gene_col]], sep = "_"),
      snp_label = paste(rsID, CH_risk_increasing_allele, .data[[gene_col]])
    )

  pval_cols <- names(dat)[str_detect(names(dat), "^pval_")]
  trait_names <- str_remove(pval_cols, "^pval_")
  beta_cols <- paste0("beta_", trait_names)

  plot_df <- dat %>%
    select(chr_num, chr_label, pos, snp_id_plain, snp_label, all_of(beta_cols), all_of(pval_cols)) %>%
    pivot_longer(
      cols = c(all_of(beta_cols), all_of(pval_cols)),
      names_to = c(".value", "trait"),
      names_pattern = "(beta|pval)_(.*)"
    ) %>%
    mutate(
      tile_class = case_when(
        pval < 5e-8 & beta > 0 ~ "p < 5e-8 / beta > 0",
        pval < 5e-8 & beta < 0 ~ "p < 5e-8 / beta < 0",
        pval < 0.05 & beta > 0 ~ "p < 0.05 / beta > 0",
        pval < 0.05 & beta < 0 ~ "p < 0.05 / beta < 0",
        TRUE ~ "NS"
      ),
      trait = factor(
        trait,
        levels = rev(trait_order),
        labels = rev(trait_label_map[trait_order])
      )
    )

  snp_order_df <- dat %>%
    distinct(chr_num, chr_label, pos, snp_id_plain, snp_label) %>%
    arrange(chr_num, pos) %>%
    mutate(x_pos = row_number())

  plot_df <- plot_df %>%
    mutate(snp_id_plain = factor(snp_id_plain, levels = snp_order_df$snp_id_plain))

  x_label_map <- snp_order_df$snp_label
  names(x_label_map) <- snp_order_df$snp_id_plain

  chr_top_df <- snp_order_df %>%
    group_by(chr_num, chr_label) %>%
    summarise(
      xmin = min(x_pos),
      xmax = max(x_pos),
      xmid = mean(x_pos),
      .groups = "drop"
    )

  chr_sep_df <- chr_top_df %>%
    arrange(chr_num)

  if (nrow(chr_sep_df) > 1) {
    chr_sep_df <- chr_sep_df[1:(nrow(chr_sep_df) - 1), , drop = FALSE] %>%
      mutate(xintercept = xmax + 0.5)
  } else {
    chr_sep_df <- chr_sep_df[0, , drop = FALSE] %>%
      mutate(xintercept = numeric(0))
  }

  list(
    plot_df = plot_df,
    chr_top_df = chr_top_df,
    chr_sep_df = chr_sep_df,
    x_label_map = x_label_map,
    n_traits = length(trait_order)
  )
}

make_tile_plot <- function(data_file, gene_col, pos_color, show_legend = TRUE) {
  dat <- read_tsv(data_file, show_col_types = FALSE)
  obj <- build_plot_data(dat, gene_col)

  top_seg_y <- obj$n_traits + 0.40
  top_lab_y <- obj$n_traits + 1.00

  p <- ggplot(obj$plot_df, aes(x = snp_id_plain, y = trait, fill = tile_class)) +
    geom_vline(
      data = obj$chr_sep_df,
      aes(xintercept = xintercept),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "#d9d9d9"
    ) +
    geom_tile(
      colour = "white",
      linewidth = 0.75,
      width = 0.96,
      height = 0.96
    ) +
    geom_segment(
      data = obj$chr_top_df,
      aes(x = xmin - 0.48, xend = xmax + 0.48, y = top_seg_y, yend = top_seg_y),
      inherit.aes = FALSE,
      linewidth = 0.45,
      colour = "black"
    ) +
    geom_text(
      data = obj$chr_top_df,
      aes(x = xmid, y = top_lab_y, label = chr_label),
      inherit.aes = FALSE,
      size = 5,
      fontface = "plain",
      colour = "black",
      vjust = 0.5
    ) +
    scale_fill_manual(
      values = c(
        "p < 5e-8 / beta > 0" = pos_color,
        "p < 5e-8 / beta < 0" = "#2166ac",
        "p < 0.05 / beta > 0" = "#f4a3a8",
        "p < 0.05 / beta < 0" = "#92c5de",
        "NS" = "#f0f0f0"
      ),
      breaks = fill_breaks,
      labels = fill_labels
    ) +
    scale_x_discrete(
      labels = obj$x_label_map,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_discrete(
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = NULL
    ) +
    coord_cartesian(
      clip = "off",
      ylim = c(0.5, obj$n_traits + 1.35)
    ) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.text.x = element_text(
        angle = 65,
        hjust = 1,
        vjust = 1,
        size = 13,
        colour = "black",
        margin = margin(t = 9)
      ),
      axis.text.y = element_text(
        size = 15,
        face = "plain",
        colour = "black",
        margin = margin(r = 10)
      ),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = ggtext::element_markdown(size = 13.5, colour = "black"),
      legend.key.size = unit(1.1, "cm"),
      legend.margin = margin(t = 2, b = 2),
      legend.spacing.x = unit(0.4, "cm"),
      plot.margin = margin(30, 18, 15, 18)
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(colour = NA)
      )
    )

  p
}

save_single_plot <- function(plot_obj, png_file, pdf_file, width, height, dpi = 600) {
  ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  grDevices::cairo_pdf(
    filename = pdf_file,
    width = width,
    height = height,
    bg = "white"
  )
  print(plot_obj)
  dev.off()
}

# ----------------------------
# Build Figure 1
# ----------------------------
fig1_plot <- make_tile_plot(
  data_file = infile_1,
  gene_col = "nearestGene",
  pos_color = "#b2181a",
  show_legend = TRUE
)

# ----------------------------
# Save Figure 1 outputs
# ----------------------------
save_single_plot(
  fig1_plot,
  file_1_png,
  file_1_pdf,
  single_width,
  single_height,
  single_dpi
)

# Optional preview in interactive session
print(fig1_plot)
