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

ancestry_levels <- c("EUR", "AFR", "EAS", "SAS")
ancestry_labels <- c(
  "EUR" = "European",
  "AFR" = "African",
  "EAS" = "East Asian",
  "SAS" = "South Asian"
)
ancestry_cols <- c(
  "EUR" = "#0072B2",
  "AFR" = "#D55E00",
  "EAS" = "#4D4D4D",
  "SAS" = "#009E73"
)

theme_forest <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = base_family, colour = "black"),
      axis.text = element_text(colour = "black", size = base_size + 2),
      axis.title = element_text(size = base_size + 4),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
      strip.text = element_text(
        family = base_family,
        colour = "black",
        size = base_size + 3,
        face = "plain",
        margin = margin(2, 2, 2, 2)
      ),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_blank(),
      plot.margin = margin(3, 8, 3, 8)
    )
}

infile <- "ed_fig2.txt"

dat_wide <- read_clean(infile) |>
  mutate(
    snp_label = sprintf(
      "'%s - %s - ' * italic('%s')",
      rsID, ea, nearestGene
    )
  )

dat <- bind_rows(
  dat_wide |>
    transmute(snp_label, ancestry = "EUR", or = as.numeric(or_EUR), l95cl = as.numeric(l95cl_EUR), u95cl = as.numeric(u95cl_EUR)),
  dat_wide |>
    transmute(snp_label, ancestry = "AFR", or = as.numeric(or_AFR), l95cl = as.numeric(l95cl_AFR), u95cl = as.numeric(u95cl_AFR)),
  dat_wide |>
    transmute(snp_label, ancestry = "EAS", or = as.numeric(or_EAS), l95cl = as.numeric(l95cl_EAS), u95cl = as.numeric(u95cl_EAS)),
  dat_wide |>
    transmute(snp_label, ancestry = "SAS", or = as.numeric(or_SAS), l95cl = as.numeric(l95cl_SAS), u95cl = as.numeric(u95cl_SAS))
) |>
  mutate(
    snp_label = factor(snp_label, levels = dat_wide$snp_label),
    ancestry = factor(ancestry, levels = ancestry_levels)
  )

p_forest <- ggplot(dat, aes(x = or, y = ancestry, colour = ancestry)) +
  geom_vline(xintercept = 1, linewidth = 0.3, linetype = "dashed", colour = "black") +
  geom_segment(aes(x = l95cl, xend = u95cl, yend = ancestry), linewidth = 0.55, lineend = "butt") +
  geom_segment(aes(x = l95cl, xend = l95cl, y = as.numeric(ancestry) - 0.08, yend = as.numeric(ancestry) + 0.08), linewidth = 0.4, lineend = "butt") +
  geom_segment(aes(x = u95cl, xend = u95cl, y = as.numeric(ancestry) - 0.08, yend = as.numeric(ancestry) + 0.08), linewidth = 0.4, lineend = "butt") +
  geom_point(size = 2.1) +
  facet_wrap(
  ~ snp_label,
  ncol = 1,
  strip.position = "top",
  labeller = label_parsed
  ) +
  scale_colour_manual(values = ancestry_cols, breaks = ancestry_levels, labels = ancestry_labels, limits = ancestry_levels, drop = FALSE) +
  scale_y_discrete(labels = ancestry_labels) +
  scale_x_log10(
    breaks = c(0.5, 1, 2, 5),
    labels = c("0.5", "1", "2", "5"),
    limits = c(0.5, 5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Odds ratio (95% CI)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_forest() +
  theme(
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

ggsave(
  "ed_fig2.pdf",
  p_forest,
  width = 7.2,
  height = 7.3,
  units = "in",
  device = cairo_pdf
)

ggsave(
  "ed_fig2.png",
  p_forest,
  width = 7.2,
  height = 7.3,
  units = "in",
  dpi = 600
)

p_forest
